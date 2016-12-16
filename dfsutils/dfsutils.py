#! /usr/bin/env python
# dfs.py - implements the Double Force Scanning method for elastic network models
# Copyright (C) 2015 Matteo Tiberti <matteo.tiberti@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/li.censes/>.

import prody
from prody import SelectionError
from libdfs import _get_drmsd
import logging as log
import numpy as np
import ConfigParser as cp
import re
import itertools
import collections
import multiprocessing as mp
import time
import socket
import getpass
import os
import traceback
import h5py

atomic_weights = {  'C' : 12.0,
                    'H' : 1.0,
                    'O' : 16.0,
                    'N' : 14.0,
                    'S' : 32.0,
                    'P' : 31.0,
                    'CA': 40.0,
                    'CL': 35.45  }


ALL = "all"
DF = "force_vectors"
DR = "displacements"
F_DR = "displacements_after_fit"
P_XYZ = "perturbed_coordinates"
FP_XYZ = "fitted_perturbed_coordinates"
SCORES = "score_matrix"
RAW_SCORES = "raw_score_matrix"
CD = "conformational_distance_matrixes"
ref_CD = "reference_cdm"
dfs_CD = "dfs_cdm"
SF="scaling_factors"


writable_data = [   ALL,
                    DF,
                    DR,
                    F_DR,
                    P_XYZ,
                    FP_XYZ,
                    SCORES,
                    RAW_SCORES,
                    CD,
                    SF
                ]

class AtomSet:
    def __init__(self, selection):
        self.atom_list = list(selection)
        self.selection = selection

    def __repr__(self):
        return self.selection.__repr__()

    def __len__(self):
        return len(self.atom_list)

    def __iter__(self):
        return self.atom_list

    def next(self):
        self.atom_list.next()

    def set(self):
        return set(self.atom_list)

    def list(self):
        return self.atom_list

class ForceSet:
    def __init__(self, origin, name, *args, **kwargs):
        self.forces = []
        self.origin = origin
        self.details = " ".join(args)
        self.name = name
        self.metadata = {}
        self.options = kwargs

        if args[0] == "fibonacci_lattice":
            times = int(args[1])
            final_r = float(args[2])
            self.metadata["times"] = times
            self.metadata["ordering"] = None

            golden_angle = np.pi * (3 - np.sqrt(5))
            theta = golden_angle * np.arange(times)
            z = np.linspace(1.0 - 1.0/times, 1.0/times-1.0, times)
            r = np.sqrt(1 - z * z)
            points = np.array(( r * np.cos(theta),
                                r * np.sin(theta),
                                z)).T * final_r

            if len(args) == 4:
                if args[3] == 'dfs':
                    points = np.repeat(points, points.shape[0], axis=0)
                    self.metadata["ordering"] = "c"
                elif args[3] == 'ref':
                    points = np.concatenate([points] * points.shape[0])
                    self.metadata["ordering"] = "r"
                else:
                    raise
            self.forces = [ Force(origin=origin, cartesian=point) for point in points ]
            self.random = False
        else:
            raise

    def __repr__(self):
        return "<ForceSet at %s, %s>" % (self.origin, self.details)

    def __str__(self):
        return "%s-%s-%s%d" % (self.name, self.origin.getChid(), self.origin.getResname(), self.origin.getResnum())

    def __getitem__(self, index):
        return self.forces[index]

    def get_vectors(self):
        return np.vstack([f.vector for f in self.forces])

class Force:
    def __init__(self, origin, cartesian=None, spherical=None):
        if (cartesian is None and spherical is None) or (cartesian is not None and spherical is not None):
            raise
        if cartesian is not None:
            self.vector = np.array(cartesian)
        if spherical is not None:
            self.vector = np.array( (spherical[0]*np.sin(spherical[1])*np.cos(spherical[2]), # r sin(theta) cos(phi)
                                    spherical[0]*np.sin(spherical[1])*np.sin(spherical[2]), # r sin(theta) sin(phi)
                                    spherical[0]*np.cos(spherical[1])) ) # r cos(theta)
        self.origin = origin.getCoords()

    def __repr__(self):
        return "<Force on %s; vector is %s>" % (self.origin, self.vector)

class ForceRun:
    def __init__(self):
        self.covariance = None
        self.force_set = None
        self.name = None
        self.F = None
        self.ready = False
        self.results = None

    def prepare(self):
        pass

    def apply_forces(self, results=None):
        if results is not None:
            self.results = results
            return

        if not self.ready or self.F is None:
            raise Exception

        self.results = np.zeros(self.F.shape)
        for i,this_F in enumerate(self.F):
            self.results[i,:] = np.dot(self.covariance, this_F)

    run = apply_forces

class PerturbationResponseRun(ForceRun):
    def __init__(self, force_sets, atoms, covariance, name="", scaling_factors=None):

        ForceRun.__init__(self)

        self.atoms = atoms
        self.force_sets = force_sets
        self.covariance = covariance
        self.name = name
        self.scaling_factors = scaling_factors

    def __repr__(self):
        return "<PerturbationResponseRun %s>" % self.name

    def prepare(self, F=None):
        if F is not None:
            self.F = F
        else:
            # assert no double atoms
            assert(len(set(self.atoms.getResindices())) == self.atoms.getResindices().shape[0])
            # assert consistence between covariance and atoms
            assert(self.covariance.shape[0] == len(self.atoms)*3)

            force_set_sizes = []
            #print self.force_sets, 'AAA'
            for force_set in self.force_sets:
                force_set_sizes.append(len(force_set.forces))

            #print "FSS", force_set_sizes
            # ensure at least one force for force set
            assert(0 not in force_set_sizes)
            # remove all the cases in which there's only one force
            force_set_sizes = filter(lambda x: x != 1, force_set_sizes)

            # times the forces to applied
            if len(force_set_sizes) < 1:
                times_forces_are_applied = 1
                force_set_sizes = [1]
            else:
                times_forces_are_applied = list(set(force_set_sizes))[0]

            # assert that all the force sets have the same size
            if len(set(force_set_sizes)) != 1:
                log.error("The provided force sets are not compatible!")
                exit(1)

            # expand to times_forces_are_applied times the force sets that only have 1 force
            for force_set in self.force_sets:
                if len(force_set.forces) == 1:
                    force_set.forces = [force_set.forces[0]] * times_forces_are_applied

            # inizialize forces arrays
            self.F = np.zeros((times_forces_are_applied, self.covariance.shape[0]))

            # get the array indexes corresponding to the atoms the force will be applied to
            force_atom_idxs = []

            for force_set in self.force_sets:
                force_atom_idxs.append(np.where(self.atoms.getResindices() == force_set.origin.getResindex())[0])

            # fill force arrays as required
            if self.scaling_factors is None:
                self.scaling_factors = np.ones(( len(self.force_sets), len(self.force_sets[0].forces)))
            else:
                #print self.scaling_factors.shape
                assert(self.scaling_factors.shape[1]) == len(self.force_sets[0].forces)
            #print "SF", self.scaling_factors
            #print "sf", self.scaling_factors
            for i in range(times_forces_are_applied):
                #print self.scaling_factors.shape
                #print len(self.force_sets[0].forces)
                    #print self.scaling_factors.shape
                    #print self.force_sets, "FS"
                    #for u,fs in enumerate(self.force_sets[0].forces):
                        #print u, fs
                    #print np.ones((len(self.force_sets[0].forces), len(self.force_sets)))

                for j,force_set in enumerate(self.force_sets):
                    #print self.scale_factor[i,j]
                    #print np.reshape(force_set.forces[i].vector, (3,1)).T
                    #print "force_v", np.reshape(force_set.forces[i].vector, (3,1)).T                    
                    scaled_force = np.reshape(force_set.forces[i].vector, (3,1)).T / self.scaling_factors[j,i]
                    #print "scaled_fv", scaled_force
                    self.F[i, force_atom_idxs[j]*3:(force_atom_idxs[j]+1)*3] = scaled_force
                    #print "F_v", self.F

        self.ready = True

class PerturbationResponseJob:
    def __init__(   self,
                    refs,
                    sets,
                    atoms,
                    original,
                    covariance,
                    fitting_operation=None,
                    fitting_string="name CA",
                    scoring_function=None,
                    name="",
                    do_write=None,
                    write_queue=None,
                    score_kwargs=None,
                    scaling_factors=None):
        # refs, variables, atoms, covariance, name
        self.refs = refs
        self.fitting_operation = fitting_operation
        self.fitting_string = fitting_string
        self.scoring_function = scoring_function
        self.sets = sets
        self.atoms = atoms
        self.covariance = covariance
        self.name = name
        self.runs = []
        self.ready = False
        self.complete = False
        self.results = None
        self.scores = None
        self.do_write = do_write
        self.write_queue = write_queue
        self.details_fh = None
        if score_kwargs is None:
            self.score_kwargs = {}
        else:
            self.score_kwargs = score_kwargs
        self.scaling_factors = scaling_factors

        self.indices = atoms.getResindices()
    def __repr__(self):
        return "<PerturbationResponseRun %s; done: %s>" % (self.name, str(self.scores is None))

    def is_in_details(self, name):
        def f(x):
            if name in x:
                return True
            return None

        return self.details_fh.visit(f)

        #self.save_displacements = save_displacements
        #self.save_displacement_magnitudes = save_displacement_magnitudes
        #self.save_force_vectors = save_force_vectors

    def prepare(self):
        log.debug("Preparing Perturbation Job")

        if self.scaling_factors is not None:
            indices = self.atoms.getResindices()

            ref_idxs = [ref.origin.getResindex() for ref in self.refs]
            ref_scale_idx = [np.where(x == indices)[0][0] for x in ref_idxs]
            ref_scaling_factors = np.array([self.scaling_factors[r] for r in ref_scale_idx])
            sf_shape = int(np.sqrt(self.scaling_factors.shape[1]))
        else:
            ref_scaling_factors = None

        self.ref_run = PerturbationResponseRun(self.refs,
                                               self.atoms,
                                               self.covariance,
                                               name="%s" % (",".join(map(str, self.refs))),
                                               scaling_factors=ref_scaling_factors)

        ref_combinations = [[j] for j in self.refs]

        combinations = itertools.product(*(ref_combinations + self.sets))

        for combination in combinations:

            if self.scaling_factors is not None:
                this_scaling_factors = np.zeros((len(combination), self.scaling_factors.shape[1]))
                for r_i,r in enumerate(ref_scaling_factors):
                    this_scaling_factors[r_i,:] = r

                sets = combination[len(ref_combinations):]
                sets_idxs = [s.origin.getResindex() for s in sets]
                sets_scale_idx = [np.where(x == indices) for x in sets_idxs]
                for s_i,s in enumerate(sets_scale_idx):
                    #print r_i+s_i+1, "RISI"
                    this_scaling_factors[r_i+1+s_i,:] = self.scaling_factors[s,:].reshape((sf_shape,sf_shape)).T.flatten()

                #sets_scaling_factors = np.array([ref_scaling_factors] + [self.scaling_factors[r][0] for r in sets_scale_idx])
            else:
                this_scaling_factors = None

            #print "RF", ref_scaling_factors.shape
            #print "SSF", this_scaling_factors.shape

            #print sets_scaling_factors
            self.runs.append(PerturbationResponseRun(combination,
                                                     self.atoms,
                                                     self.covariance,
                                                     name="%s" % (",".join(map(str, combination))),
                                                     scaling_factors = this_scaling_factors))
                                                     
        self.ready = True
        log.info("Preparation done.")

    def run(self):

        if self.do_write is not None:
            write_DF = DF in self.do_write
            write_DR = DR in self.do_write
            write_F_DR = F_DR in self.do_write
            write_P_XYZ = P_XYZ in self.do_write
            write_FP_XYZ = FP_XYZ in self.do_write
            write_RAW_SCORES = RAW_SCORES in self.do_write
            write_CD = CD in self.do_write
        else:
            write_DF = False
            write_DR = False
            write_F_DR = False
            write_P_XYZ = False
            write_FP_XYZ = False
            write_RAW_SCORES = False
            write_CD = False

        ref_idxs = [ref.origin.getResindex() for ref in self.refs]
        ref_combinations = [[j] for j in self.refs]
        combinations = itertools.product(*(ref_combinations + self.sets))

        self.scores = np.zeros([len(self.indices) for i in self.sets])
        self.max_scores = np.zeros([len(self.indices) for i in self.sets])
        self.min_scores = np.zeros([len(self.indices) for i in self.sets])
        displaced_atoms = self.atoms.copy()

        self.ref_run.prepare()
        self.ref_run.run()

        fitted_ref_structures = prody.Ensemble(title="displaced_fitted_structures")
        fitted_ref_structures.setAtoms(self.atoms)

        if write_F_DR:
            F_DR_values = np.zeros(self.ref_run.results.shape)
        if write_FP_XYZ:
            FP_XYZ_values = np.zeros((self.ref_run.results.shape[0],self.ref_run.results.shape[1]/3,3))
        if write_P_XYZ:
            P_XYZ_values = np.zeros((self.ref_run.results.shape[0],self.ref_run.results.shape[1]/3,3))

        for d,displacements in enumerate(self.ref_run.results):
            displaced_atoms.setCoords(self.atoms.getCoords() + displacements.reshape((displacements.shape[0]/3,3)))
            if self.fitting_operation is not None:
                fitting_original_atoms = self.atoms.select(self.fitting_string)
                fitting_displaced_atoms = displaced_atoms.select(self.fitting_string)
                transformation = self.fitting_operation(fitting_displaced_atoms, fitting_original_atoms).transformation
                this_ref_fit = self.fitting_operation(displaced_atoms, self.atoms, transformation=transformation)
                fitted_ref_structures.addCoordset(this_ref_fit.transformed())
                if write_F_DR:
                    F_DR_values[d,:] = (this_ref_fit.transformed().getCoords()- self.atoms.getCoords()).flatten()
                if write_FP_XYZ:
                    FP_XYZ_values[d,:] = this_ref_fit.transformed().getCoords()
            else:
                fitted_ref_structures.addCoordset(displaced_atoms)

            if write_P_XYZ:
                P_XYZ_values[d,:] = displaced_atoms.getCoords()

        if self.do_write is not None:
            if write_DR:
                self.write_queue.put((DR, ref_idxs, self.ref_run.results))
            if write_DF:
                self.write_queue.put((DF, ref_idxs, self.ref_run.F))
            if write_FP_XYZ:
                self.write_queue.put((FP_XYZ, ref_idxs, FP_XYZ_values))
            if write_F_DR:
                self.write_queue.put((F_DR, ref_idxs, F_DR_values))
            if write_P_XYZ:
                self.write_queue.put((P_XYZ, ref_idxs, P_XYZ_values))

            del self.ref_run.results
            del self.ref_run.F

        for r,run in enumerate(self.runs):

            this_comb = combinations.next()
            this_ref = this_comb[:len(self.refs)]
            this_perturbation = this_comb[-len(self.sets):]
            fitted_structures = prody.Ensemble(title="displaced_fitted_structures")
            fitted_structures.setAtoms(self.atoms)

            this_idxs = [this.origin.getResindex() for this in this_comb]        

            run.prepare()
            run.run()

            for d,displacements in enumerate(run.results):

                displaced_atoms.setCoords(self.atoms.getCoords() + displacements.reshape((displacements.shape[0]/3,3)))

                if self.fitting_operation is not None:
                    fitting_original_atoms = self.atoms.select(self.fitting_string)
                    fitting_displaced_atoms = displaced_atoms.select(self.fitting_string)
                    transformation = self.fitting_operation(fitting_displaced_atoms, fitting_original_atoms).transformation
                    this_fit = self.fitting_operation(displaced_atoms, self.atoms, transformation=transformation)
                    fitted_structures.addCoordset(this_fit.transformed())
                    if write_F_DR:
                        F_DR_values[d,:] = (this_ref_fit.transformed().getCoords() - self.atoms.getCoords()).flatten()
                    if write_FP_XYZ:
                        FP_XYZ_values[d,:] = this_fit.transformed().getCoords()
                else:
                    fitted_structures.addCoordset(displaced_atoms)

                if write_P_XYZ:
                    P_XYZ_values[d,:] = displaced_atoms.getCoords()

            ref_score_idxs = [ np.where(c.origin.getResindex() == self.indices) for c in this_ref ]
            this_score_idxs = [ np.where(c.origin.getResindex() == self.indices) for c in this_comb ]

            self.score_kwargs["ref_idxs"] = ref_idxs
            self.score_kwargs["this_idxs"] = this_idxs
            self.score_kwargs["ref_score_idxs"] = ref_score_idxs
            self.score_kwargs["this_score_idxs"] = this_score_idxs
            
            score_idxs = [ np.where(c.origin.getResindex() == self.indices) for c in this_perturbation ]
            score = self.scoring_function(  self.atoms,
                                            fitted_ref_structures,
                                            fitted_structures,
                                            metadata=[fs.metadata for fs in run.force_sets],
                                          )

            if r == 0:
                ref_cds = score.ref_cds
            else:
                score.ref_cds = ref_cds

            avg_s, min_s, max_s, details = score.get(**self.score_kwargs)

            self.scores[np.array(score_idxs)] = avg_s
            self.min_scores[np.array(score_idxs)] = min_s
            self.max_scores[np.array(score_idxs)] = max_s

            if self.do_write is not None:
                if write_DR:
                    self.write_queue.put((DR, this_idxs, run.results))
                if write_DF:
                    self.write_queue.put((DF, this_idxs, run.F))
                if write_FP_XYZ:
                    self.write_queue.put((FP_XYZ, this_idxs, FP_XYZ_values))
                if write_F_DR:
                    self.write_queue.put((F_DR, this_idxs, F_DR_values))
                if write_P_XYZ:
                    self.write_queue.put((P_XYZ, this_idxs, P_XYZ_values))
                if write_RAW_SCORES:
                    self.write_queue.put((RAW_SCORES, this_idxs, details[RAW_SCORES]))
                if write_CD:
                        self.write_queue.put((ref_CD, this_idxs, details[CD][0]))
                        self.write_queue.put((dfs_CD, this_idxs, details[CD][1]))

            del run.F
            del run.results

class Fitting:
    def __init__(self, mobile, target, transformation=None, *args, **kwargs):
        self.mobile = mobile
        self.target = target
        self._transformed = None

        if transformation is not None:
            self.transformation = transformation
        else:
            self.compute_transformation()

    def compute_transformation(self):
        pass

    def apply_transformation(self):
        self._transformed = self.transformation.apply(self.mobile.copy())

    def transformed(self):
        if self._transformed is None:
            self.apply_transformation()
        return self._transformed


class LsqFitting(Fitting):
    def compute_transformation(self):
        self.transformation = prody.calcTransformation(self.mobile, self.target)

def parse_force_string(s, type=None):
    args = s.strip().split()
    if len(args) < 2:
        raise
    if type:
        args.append(type)
    return args

def _worker(work_package):
    try:
        idxs,job = work_package
        log.info("Now doing %s - starting PerturbationResponse" % str(idxs))
        job.prepare()
        job.run()
        log.info("Job %s done" % str(idxs))
        if job.do_write is not None:
            job.write_queue.put('DONE')
            log.info("Waiting for write queue to flush")
            while job.write_queue.qsize() > 100:
                time.sleep(1)
                continue
        return (idxs, job.scores, job.min_scores, job.max_scores)
    except:
        traceback.print_exc()
        raise

def _scribe(write_queue, details_fname, runsn, dtype='float32', precision=-1, metadata=None):
    DONE = 'DONE'
    dead_procs = 0

    fh = h5py.File(details_fname,'w')
    log.info("Scribe starting on queue %d" % id(write_queue))

    if metadata is not None:
        for k,v in metadata.iteritems():
            fh.attrs[k] = v

    for data_piece in iter(write_queue.get, 'NEVER'):
        if data_piece == DONE:
            dead_procs += 1
            log.debug("DONE signal received! We're at %d/%d" % (dead_procs,runsn))
            if dead_procs == runsn and write_queue.empty():
                log.debug("All processes done! Scribe will finish and die")
                fh.close()
                return
            else:
                continue
        if type(data_piece[1]) is str:
            this_group = data_piece[1]
        else:
            this_group = "%s" % ("/".join(map(str,data_piece[1])))
        fh.require_group(this_group)

        if precision < 0:
            fh[this_group].create_dataset(data_piece[0], data=data_piece[2], dtype=dtype)
        else:
            fh[this_group].create_dataset(data_piece[0], data=data_piece[2], dtype=dtype, scaleoffset=precision)

class DFSJob:
    def __init__(   self,
                    structure,
                    atoms,
                    original,
                    covariance,
                    config_file,
                    fitting_operation=None,
                    fitting_string="",
                    scoring_function=None,
                    scoring_function_name=None,
                    name="",
                    do_write=None,
                    output_fname="",
                    score_kwargs=None,
                    options={}):

        self.name = name
        self.atoms = atoms
        self.structure = structure
        self.jobs = []
        self.job_results = None
        self.do_write = do_write
        self.options = options

        if self.do_write is not None:
            if fitting_operation is None and FP_XYZ in do_write:
                if P_XYZ in do_write:
                    log.warning("no fitting option specified; fitted coordinates will not be saved.")
                else:
                    log.warning("saving fitted coordinates was requested, but no fitting option specified. Will save perturbed non-fitted coordinates instead.")
                    self.do_write.append(P_XYZ)

        reference_force_sets, force_sets, config_options = self.parse_config_file(config_file, self.structure, self.atoms)

        #print "valssss", reference_force_sets.values()
        ref_combinations = itertools.product(*reference_force_sets.values())
        #print "refsss", reference_force_sets
        #print "listone",list(ref_combinations)

        #combinations =     itertools.product(*force_sets.values())

        dimensions = len(force_sets.keys()) + len(reference_force_sets.keys())

        indices = atoms.getResindices()

        self.score_matrix =     np.zeros([indices.shape[0] for i in range(dimensions)])
        self.min_score_matrix = np.zeros([indices.shape[0] for i in range(dimensions)])
        self.max_score_matrix = np.zeros([indices.shape[0] for i in range(dimensions)])

        #print "combsss", list(ref_combinations)

        #log.info("Ok, we will write details to file")
        manager = mp.Manager()
        self.write_queue = manager.Queue()

        if self.options["scale_forces_by_reference_cd"]:
            log.info("Computing scale factors ...")
            scale_runs = [PerturbationResponseRun(comb,
                                                atoms,
                                                covariance,
                                                name="scale run",
                                                scaling_factors=None) for comb in ref_combinations]

            for s in scale_runs:
                s.prepare()
                s.run()
            log.info("Done!")

            if scoring_function_name == "drmsd":
                cd_function = _get_drmsd
            else:
                raise

            self.scaling_factors, self.scaling_cds, self.scaling_coords = get_scaling_factors(original,
                                                 scale_runs,
                                                 cd_function,
                                                 exclude_application_site=True,
                                                 fitting_operation=None,
                                                 fitting_string=None)
        else:
            self.scaling_factors = None

        #print "SCALA", scaling_factors,
        #print scaling_factors.shape

        ref_combinations = itertools.product(*reference_force_sets.values())

        for comb in ref_combinations:
            matrix_ref_idxs = [np.where(c.origin.getResindex() == indices)[0][0] for c in comb]

            self.jobs.append([ matrix_ref_idxs, PerturbationResponseJob(comb,
                                                                        force_sets.values(),
                                                                        atoms,
                                                                        atoms,
                                                                        covariance,
                                                                        fitting_operation,
                                                                        fitting_string,
                                                                        scoring_function,
                                                                        "",
                                                                        do_write,
                                                                        self.write_queue,
                                                                        scaling_factors = self.scaling_factors,
                                    score_kwargs = score_kwargs)])

    def parse_config_file(self, config_file, structure, atoms):
        # XXX check uniqueness of names

        reference_atom_sets_section = "first sites atom sets"
        reference_force_sets_section = "first sites force sets"

        dfs_atom_sets_section = "second sites atom sets"
        dfs_force_sets_section = "second sites force sets"

        options_section = "options"

        recognized_option_types = {}
        recognized_option_defaults = {}
        options = recognized_option_defaults.copy()

        required_sections = set([reference_atom_sets_section,
                                 reference_force_sets_section,
                                 dfs_atom_sets_section,
                                 dfs_force_sets_section])
        
        reference_atom_sets = {}
        reference_force_sets = {}
        dfs_atom_sets = {}
        dfs_force_sets = {}

        parser = cp.ConfigParser()
        parser.SECTCRE = re.compile(r"\[ *(?P<header>[^]]+?) *\]") # To remove trailing white space error thing
        
        if type(config_file) is str:   
            config_file_stream = open(config_file, 'r')
        else:
            config_file_stream = config_file
            
        parser.readfp(config_file_stream)

        if parser.has_section(options_section):
            for name, s in parser.items(options_section):
                if name not in recognized_option_types.keys():
                    log.warning("Option %s in file %s not recognized! It will be ignored." % (config_file, name))
                    continue
                if recognized_option_types[name] == "bool":
                    options[name] = parser.getboolean(options_section, name)
                else:
                    options[name] = parser.get(options_section, name)

        if len(required_sections & set(parser.sections())) == 4: # combinations mode # XXX all of none these sections must be present
            for name, s in parser.items(reference_atom_sets_section):
                try:
                    this_set = AtomSet(structure.select(s))
                except:
                    log.error("Selection string \"%s\" for set %s was invalid." % (s, name))
                    exit(1)

                if not this_set:
                    log.error("Your selection string %s resulted in no atoms being selected in your structure." % s )
                    exit(1)

                if not this_set.selection in atoms:
                    log.error("Your selection string \"%s\" resulted in a selection which is not included in the ANM model." % s )
                    exit(1)

                reference_atom_sets[name] = this_set
                reference_force_sets[name] = []
                
                print reference_atom_sets, "AAZ"

            for name, s in parser.items(reference_force_sets_section):
                if name not in reference_atom_sets.keys():
                    log.error("Force set %s was defined, but no correspondent atom sets found." % name)
                    exit(1)
                try:
                    args = parse_force_string(s, type="ref")
                except:
                    log.error("Force set %s wasn't defined properly." % name)
                    exit(1)

                for atom in reference_atom_sets[name].atom_list:
                    try:
                        reference_force_sets[name].append(ForceSet(atom, name, *args))
                    except:
                        log.error("Couldn't calculate force set for atom %s in set %s" % (atom, name))
                        exit(1)

            for name, s in parser.items(dfs_atom_sets_section):
                try:
                    this_set = AtomSet(structure.select(s))
                except:
                    log.error("Selection string \"%s\" for set %s was invalid." % (s, name))
                    exit(1)

                if not this_set:
                    log.error("Your selection string %s resulted in no atoms being selected in your structure." % s )
                    exit(1)

                if not this_set.selection in atoms:
                    log.error("Your selection string \"%s\" resulted in a selection which is not included in the ANM model." % s )
                    exit(1)

                dfs_atom_sets[name] = this_set
                dfs_force_sets[name] = []

            for name, s in parser.items(dfs_force_sets_section):
                if name not in dfs_atom_sets.keys():
                    log.error("Force set %s was defined, but no correspondent atom sets found." % name)
                    exit(1)
                try:
                    args = parse_force_string(s, type="dfs")
                except:
                    log.error("Force set %s wasn't defined properly." % name)
                    exit(1)

                for atom in dfs_atom_sets[name].atom_list:
                    try:
                        dfs_force_sets[name].append(ForceSet(atom, name, *args))
                    except:
                        log.error("Couldn't calculate force set for atom %s in set %s" % (atom, name))
                        exit(1)

            return reference_force_sets, dfs_force_sets, options

        else:
            # XXX implement custom cases
            return None
        #print self.score_matrix

    def get_score_matrix(self):
        for idxs, avg_results, min_results, max_results in self.job_results:
            self.score_matrix[idxs] = avg_results
            self.min_score_matrix[idxs] = min_results
            self.max_score_matrix[idxs] = max_results

        np.fill_diagonal(self.score_matrix, 0.)
        np.fill_diagonal(self.min_score_matrix, 0.)
        np.fill_diagonal(self.max_score_matrix, 0.)

        return (self.score_matrix, self.min_score_matrix, self.max_score_matrix)

    def run(self, nprocs=1, details_output_fname=None, metadata=None, precision=5, read_details=False):

        minproc = np.min((len(self.jobs),nprocs))

        log.info("Running parallel DFS. %d procs will be used for %d jobs" % (minproc,len(self.jobs)))
        if nprocs > 1:
            pool = mp.Pool( processes=minproc )
        else:
            pool = None

        if details_output_fname is not None:

            scribe_p = mp.Process(target=_scribe, name="scribe_process", kwargs={   "write_queue" : self.write_queue,
                                                                                    "details_fname" : details_output_fname,
                                                                                    "runsn" : len(self.jobs),
                                                                                    "precision" : precision,
                                                                                    "metadata" : metadata })

            scribe_p.start()
        if self.do_write is not None and self.scaling_factors is not None:    
            if SF in self.do_write:
                log.info("Writing scaling factors to file ...")
                self.write_queue.put((SF, "/", self.scaling_factors))
            #self.write_queue.put(("scaling_cds", "/", self.scaling_cds))
            #for i,c in enumerate(self.scaling_coords):
                #self.write_queue.put(("%d" %i, "/scaling_coords/", c))
            
            

        if pool:
            results = pool.map_async(_worker, self.jobs)

            pool.close()

            if details_output_fname is not None:
                scribe_p.join()
            else:
                pool.join()

            self.job_results = results.get()

        else:
            results = [ _worker(j) for j in self.jobs ]

            self.job_results = results
class Score:
    def __init__(self, original, ref_structures, structures, metadata=None, **kwargs):
        self.ref_structures = ref_structures
        self.structures = structures
        self.original = original
        self.metadata = metadata
        self.ref_cds = None
        self.ref_original = self.original.copy()
        self.this_original = self.original.copy()

    def score(self, multi=False, details=None, **kwargs):

        self.ref_structures.setCoords(self.ref_original)
        self.structures.setCoords(self.this_original)

        if kwargs["exclude_sites"]:
            idxs = self.ref_original.getResindices()
            #print list(self.ref_original)
            #print "XXX", kwargs["this_score_idxs"]

            mask = idxs[np.logical_and.reduce([ idxs != s[0] for s in kwargs["this_score_idxs"]])]
            
            #print "mask", mask            
            
            self.ref_original = self.ref_original[mask]
            self.this_original = self.this_original[mask]

            ref_structures = prody.Ensemble()
            structures = prody.Ensemble()
            ref_structures.setAtoms(self.ref_original)
            structures.setAtoms(self.this_original)

            #print self.ref_structures, "GGG"
            #print "idxs", idxs

            mask = idxs[np.logical_and.reduce([ idxs != s for s in np.array(kwargs["this_score_idxs"]).T[0,0,:]])]        
            #print "mask2", mask    

            
            #print '----'            
            
            #print "this_idxs", kwargs["this_idxs"]
            #print "this_score_idxs", kwargs["this_score_idxs"]            
            #print self.ref_structures.getCoordsets().shape
            #print len(ref_structures.getAtoms())
            
            ref_structures.addCoordset(self.ref_structures.getCoordsets()[:,mask,:])
            structures.addCoordset(self.structures.getCoordsets()[:,mask,:])
 
            self.ref_structures = ref_structures
            self.structures = structures

            self.ref_cds = None

        if multi:
            return self._score_single(kwargs=kwargs)
        else:
            return self._score_multi(kwargs=kwargs)

    def _score_single(self, details=None, **kwargs):
        return None

    def _score_multi(self, details=None, **kwargs):
        return None

    get = score

class ScoreRMSD(Score):
    
    name = "rmsd"

    def _score(self, **kwargs):
        self.ref_structures.setCoords(self.ref_original)
        self.structures.setCoords(self.this_original)

        metadata_available = np.sum([map(bool,self.metadata)]) != 0
        if self.metadata is not None and metadata_available:
            if len(self.metadata) > 2: #XXX support more than 2
                raise
            if np.sum([map(bool,self.metadata)]) == len(self.metadata):

                if not ((self.metadata[0]["ordering"] == 'c' and self.metadata[1]["ordering"] == 'r') or (self.metadata[0]["ordering"] == 'r' and self.metadata[1]["ordering"] == 'c')):
                    raise
                structures_per_ref = int(np.sqrt(len(self.structures)))
            else:
                structures_per_ref = len(self.structures)
        else:
            structures_per_ref = len(self.structures)

        assert( len(self.structures) % len(self.ref_structures) == 0 )
        #reference_rmsds = self.ref_structures.getRMSDs().repeat(len(self.structures)/len(self.ref_structures)) # 1 1 1 2 2 2 3 3 3 ...
        #print "referenzio"
        if not self.ref_cds:
            ref_original_c    = self.ref_original.getCoordsets()
            ref_structures_c  = self.ref_structures.getCoordsets()
            if ref_original_c.shape[0] == 1:
                ref_original_c = np.repeat(ref_original_c, ref_structures_c.shape[0], axis=0)


            reference_rmsds = self.ref_structures.getRMSDs()

            reference_rmsds = reference_rmsds.reshape(( reference_rmsds.shape[0]/structures_per_ref,
                                                        structures_per_ref )).T
            self.ref_cds = reference_rmsds

        dfs_rmsds = self.structures.getRMSDs().reshape(self.ref_cds.shape)

        scores = (self.ref_cds - dfs_rmsds)/self.ref_cds

        details =   {
                        RAW_SCORES : scores,
                        CD : [self.ref_cds, dfs_rmsds]
                    }

        return (np.average(np.max(scores, axis=1), ),
                np.min(    np.max(scores, axis=1), ),
                np.max(    np.max(scores, axis=1), ), details)

    _score_multi = _score
    _score_single = _score


class ScoreDRMSD(Score):

    name = "drmsd"

    def __init__(self, *args, **kwargs):
        Score.__init__(self, *args,**kwargs)

        try:
            from libdfs import _get_drmsd
            self._get_drmsd = _get_drmsd
            #log.info("Fast Cython implementation will be used!")
        except ImportError:
            log.warning("Fast compiled version of DRMSD calculation not found; standard python will be used. Expect high memory usage")
            self._get_drmsd = self._get_drmsd_slow

    def _score(self, **kwargs):

        metadata_available = np.sum([map(bool,self.metadata)]) != 0

        if self.metadata is not None and metadata_available:
            if len(self.metadata) > 2: # support more than 2
                raise
            if np.sum([map(bool,self.metadata)]) == len(self.metadata):

                if not ((self.metadata[0]["ordering"] == 'c' and self.metadata[1]["ordering"] == 'r') or (self.metadata[0]["ordering"] == 'r' and self.metadata[1]["ordering"] == 'c')):
                    raise
                structures_per_ref = int(np.sqrt(len(self.structures)))
            else:
                structures_per_ref = len(self.structures)
        else:
            structures_per_ref = len(self.structures)


        if not self.ref_cds:
            ref_original_c    = self.ref_original.getCoordsets()
            ref_structures_c  = self.ref_structures.getCoordsets()
            if ref_original_c.shape[0] == 1:
                ref_original_c = np.repeat(ref_original_c, ref_structures_c.shape[0], axis=0)

            reference_rmsds = self._get_drmsd(ref_original_c, ref_structures_c)

            reference_rmsds = reference_rmsds.reshape(( reference_rmsds.shape[0]/structures_per_ref,
                                                        structures_per_ref, ))
            self.ref_cds = reference_rmsds

        this_original_c = self.this_original.getCoordsets()
        structures_c = self.structures.getCoordsets()
        if this_original_c.shape[0] == 1:
            this_original_c = np.repeat(this_original_c, ref_structures_c.shape[0], axis=0)

        dfs_rmsds = self._get_drmsd(this_original_c, structures_c).reshape(self.ref_cds.shape)

        scores = (self.ref_cds - dfs_rmsds)/self.ref_cds

        details =   {
                        RAW_SCORES : scores,
                        CD : [self.ref_cds, dfs_rmsds]
                    }

        return (np.average(np.max(scores, axis=0), ),
                np.min(    np.max(scores, axis=0), ),
                np.max(    np.max(scores, axis=0), ), details)


    _score_multi = _score
    _score_single = _score

    def _get_drmsd_slow(self, ref, ens, **kwargs):

        eln = float(ens.shape[1])*(float(ens.shape[1])-1)/2.0

        assert(len(ref.shape) == 3)
        assert(ref.shape[0] == 1 or ref.shape[0] == ens.shape[0])

        idx_a, idx_b = np.triu_indices(ens.shape[1], k=1)

        ref_a = ref[:,idx_a,:]
        ref_b = ref[:,idx_b,:]

        delta_ref = np.linalg.norm(ref_a-ref_b, axis=2)

        del ref_a, ref_b

        ens_a = ens[:,idx_a,:]
        ens_b = ens[:,idx_b,:]

        del idx_a, idx_b

        delta_ens = np.linalg.norm(ens_a-ens_b, axis=2)

        del ens_a, ens_b,

        return np.sqrt(np.sum((delta_ref - delta_ens)**2, axis=1)/eln)

def get_scaling_factors(original,
                        runs,
                        cd_function,
                        exclude_application_site=True,
                        fitting_operation=None,
                        fitting_string=None,
                        target_cd=1.0):

    out = np.zeros((len(original), runs[0].results.shape[0]))
    coords_out = []

    atoms = original

    for r in runs:
        this_original = original.copy()
        if exclude_application_site:
            this_idxs = [ np.where(atoms.getResindices() == fs.origin.getResindex())[0][0] for fs in r.force_sets ]
            ref_original = this_original.select(" and ".join(["(not resindex %s)"%s for s in this_idxs ]))
            ref_original_c = ref_original.getCoordsets()
        else:
            ref_original_c = this_original.getCoordsets()

        ref_original_c = np.repeat(ref_original_c, r.results.shape[0], axis=0)

        fitted_ref_structures = prody.Ensemble(title="displaced_fitted_structures")
        fitted_ref_structures.setAtoms(atoms)

        for d,displacements in enumerate(r.results):
            displaced_atoms = original.copy()

            displaced_atoms.setCoords(atoms.getCoords() + displacements.reshape((displacements.shape[0]/3,3)))

            if fitting_operation is not None:
                fitting_original_atoms = atoms.select(fitting_string)
                fitting_displaced_atoms = displaced_atoms.select(fitting_string)
                transformation = fitting_operation(fitting_displaced_atoms, fitting_original_atoms).transformation
                this_ref_fit = fitting_operation(displaced_atoms, atoms, transformation=transformation)
                fitted_ref_structures.addCoordset(this_ref_fit.transformed())
            else:
                fitted_ref_structures.addCoordset(displaced_atoms)

        if exclude_application_site:
            filtered_fitted_ref_structures = prody.Ensemble(title="displaced_fitted_structures")
            filtered_fitted_ref_structures.setAtoms(ref_original)
            filtered_fitted_ref_structures.addCoordset(np.delete(fitted_ref_structures.getCoordsets(), this_idxs, axis=1))
        else:
            filtered_fitted_ref_structures = fitted_ref_structures

        coords_out.append(fitted_ref_structures.getCoordsets())
        out[this_idxs,:] = cd_function(ref_original_c, filtered_fitted_ref_structures.getCoordsets())

    out = np.array(out)
    #return target_cd / out, out, np.array(coords_out)
    return out, out, np.array(coords_out)


def assign_masses(atoms):
    masses = []
    for elm in atoms.getElements():
        try:
            masses.append(atomic_weights[elm])
        except:
            log.warning("Element %s was not recognized; Standard mass of 12 will be used" % elm)
            masses.append(12)
    atoms.setMasses(masses)

def inertia_tensor(atoms):
    return np.sum( [ m*(np.dot(c,c)*np.identity(3)-np.outer(c,c)) for c,m in zip(atoms.getCoords(),atoms.getMasses()) ], axis=0)

def inertia_rotation_matrix(inertia_tensor):
    return np.linalg.eig(inertia_tensor)[1].T

def safe_savetxt(fname, data, fmt="%1.5f"):
    try:
        np.savetxt(fname, data, fmt)
    except:
        log.error("Couldn't save file %s; exiting" % fname)
        exit(1)

def safe_savenpy(fname, data):
    try:
        np.save(fname, data)
    except:
        log.error("Couldn't save file %s; exiting" % fname)
        exit(1)

def safe_savenmd(fname, eigvs, atoms):
    try:
        prody.dynamics.writeNMD(fname, eigvs, atoms)
    except:
        log.error("Couldn't write file %s; exiting" % fname)
        exit(1)
