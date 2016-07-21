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
from dfsutils import *
from dfsutils.libdfs import _get_drmsd
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
import argparse

print ScoreRMSD

scores =    {   "rmsd": ScoreRMSD,
                "drmsd": ScoreDRMSD,
            }

fitting_operations =    {   "lsq_fit" : LsqFitting,
                            "none" : None
                        }

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("PDB", type=str, help="input PDB file")
    parser.add_argument("-s", "--selection", dest="selection_string", type=str, default="name CA", help="CHARMM-style selection of atoms to be used in the model")
    parser.add_argument("-v", "--verbose", dest="verbose", help="toggle verbose mode")
    parser.add_argument("-m", "--write-hessian", dest="hessian_output", nargs='?', const="hessian.dat", type=str, default=None, help="Save Hessian matrix (default: matrix.dat)")
    parser.add_argument("-M", "--load-hessian", dest="hessian_input", type=str, default=None, help="Load hessian matrix instead of calculating it")
    parser.add_argument("-i", "--write-covariance", dest="covariance_output", nargs='?', const="covariance.dat", type=str, default=None, help="Save covariance matrix (default: covariance.dat)")
    parser.add_argument("-I", "--load-covariance", dest="covariance_input", type=str, default=None, help="Load covariance matrix instead of calculating it")
    parser.add_argument("-b", "--eigenvectors", dest="nm_eigvec_output", nargs='?', const=1, type=str, default="nmodes_eigenvectors.dat", help="Save normal modes eigenvectors as a text file")
    #parser.add_argument("-l", "--eigenvalues", dest="nm_eigval_output", nargs='?', const=1, type=str, default="nmodes_eigenvalues.dat", help="Save normal modes eigenvalues as a text file")
    parser.add_argument("-D", "--nmd", dest="nmd_output", nargs='?', const=1, type=str, default="modes.nmd", help="Save normal mode eigenvectors as a NMD file")
    parser.add_argument("-f", "--dfs", dest="dfs_config", default=None, help="Configuration file for DFS runs")
    parser.add_argument("-r", "--ref-dfs", dest="ref_dfs_config", default=None, help="Configuration file for DFS runs (reference)")
    parser.add_argument("-c", "--cutoff", dest="anm_cutoff", default=15.0, type=float, help="Cut-off for generating the topology of ANM (in A; default is 15.0)")
    parser.add_argument("-g", "--gamma", dest="anm_gamma", default=1.0, type=float, help="Spring force constant for ANM; default is 1.0)")
    parser.add_argument("-F", "--fitting", dest="fitting_operation", default="none", action="store", type=str, choices=fitting_operations.keys(), help="perform least square fitting between the displaced and original geometries")
    parser.add_argument("-a", "--fit-selection", dest="fit_selection", default="name CA", help="Fit the selection in --selection on this group of atoms")
    parser.add_argument("-d", "--conformational-difference", dest="conformational_difference", choices=scores.keys(), help="type of conformational score to use (default: drmsd)", default="drmsd")
    parser.add_argument("-n", "--np", dest="np", default=1, type=int, action="store", help="Number of cores to be used at the same time for calculation")
    parser.add_argument("-w", "--write", dest="write", nargs='*', type=str, choices=writable_data, help="Choose which data should be saved in the output file. 'all' (default behaviour) just saves everything.")
    parser.add_argument("-o", "--output", dest="output_fname", default="details.h5", type=str, action="store", help="Filename to be used as output in case --write is specified")
    parser.add_argument("-t", "--force_mode", dest="force_mode", nargs='?', type=str, default="both", choices=["fixed_force, fixed_displacement, both"], help="help text I will write")
    parser.add_argument("-S", "--output-scores", dest="output_scores", default="scores", type=str,  action="store", help="Filename to be used for scores")
    parser.add_argument("-q", "--precision", dest="precision", default=-1, type=int, help="Number of decimal places to be used for compression (-1 for lossless)")
    parser.add_argument("-x", "--include-application-sites", dest="exclude_sites", action='store_false', default=True, help="Include force application sites in the calculation of scores")
    #XXX add option: number of normal modes

    log.basicConfig(level=log.DEBUG)

    args = parser.parse_args()

    if args.write is None:
        do_write = None
    elif args.write == [] or (ALL in args.write):
        if len(args.write) > 1:
            log.warning("all was specified in option --write, together with others; all will be used.")
        do_write = writable_data[1:]
    else:
        do_write = args.write
    if do_write is not None:
        log.info("this information will be written to file: %s" % (", ".join(list(do_write))))
        import h5py

    if args.covariance_input and args.hessian_input:
        log.error("You cannot specify covariance matrix and hessian at the same time. Exiting...")
        exit(1)

    # 1. Parse PDB
    try:
        structure = prody.parsePDB(args.PDB)
    except:
        log.error("Couldn't parse your model file (%s); Exiting..." % args.PDB)
        exit(1)

    assign_masses(structure)

    # 2. select atoms for building the model
    log.info("Selection string is %s" % args.selection_string)
    try:
        atoms = structure.select(args.selection_string)
    except SelectionError:
        log.error("Your selection string was invalid. Exiting...")
        exit(1)
    if not atoms:
        log.error("Your selection string resulted in no atoms being selected in your structure. Exiting...")
        exit(1)
    if len(atoms) < 2:
        log.error("Your selection string resulted in less than two atoms being selected in your structure. Exiting...")
        exit(1)

    # XXX: Check 1 atom per residue

    if args.covariance_input is not None:
        try:
            covariance = np.load(args.covariance_input)
        except:
            log.error("File %s is not readable or in wrong format. Exiting...")
            exit(1)
        assert(covariance.shape[0] == len(atoms)*3)
        log.warning("The covariance matrix has been loaded from file; therefore gamma and cut-off options will be ignored.")
    else:
        hessian = None
        if args.hessian_input is not None:
            for load_function in [np.load, np.loadtxt]:
                try:
                    hessian = load_function(args.hessian_input)
                    break
                except:
                    continue
            else:
                log.error("File %s is not readable or in wrong format. Exiting..." % args.hessian_input)
                exit(0)

            assert(hessian.shape[0] == 3*len(atoms))
            log.warning("Hessian has been loaded from file; therefore gamma and cut-off options will be ignored.")

        anm = prody.dynamics.ANM(atoms)

        # Build Hessian

        if hessian is None:
            anm.buildHessian(atoms, cutoff=args.anm_cutoff, gamma=args.anm_gamma)
        else:
            anm.setHessian(hessian)

        # Write Hessian to text file
        if args.hessian_output:
            safe_savenpy(args.hessian_output, anm.getHessian())

        # Calculate NMs
        anm.calcModes(n_modes = anm.numDOF())

        # Write eigenvectors to NMD file
        safe_savenmd(args.nmd_output, anm, atoms)

        covariance = anm.getCovariance()
        if args.covariance_output:
            safe_savenpy(args.covariance_output, covariance)

    if do_write is not None:

        attrs = {}
        attrs['hostname'] = socket.gethostname()
        attrs['user'] = getpass.getuser()
        attrs['date'] = time.strftime("%c")
        attrs['directory'] = os.getcwd()

        attrs['input_PDB'] = args.PDB
        attrs['atom_selection'] = args.selection_string
        attrs['gamma'] = args.anm_gamma
        attrs['cutoff'] = args.anm_cutoff

        attrs['cdm'] = args.conformational_difference
        attrs['fitting'] = args.fitting_operation
        attrs['fitting_selection'] = args.fit_selection
        attrs['available_data'] = ",".join(do_write)


    else:
        fh = None
        attrs = None

    score_kwargs = { "exclude_sites" : args.exclude_sites,
                    }

    if args.force_mode == "fixed_force" or args.force_mode == "both":
        if do_write is not None:
            output_fname = "%s_%s" % ("fixed_force", args.output_fname)
        else:
            output_fname = None
        output_scores = "%s_%s" % ("fixed_force", args.output_scores)

        dfs_job = DFSJob(
                            structure,
                            atoms,
                            atoms,
                            np.matrix(covariance),
                            args.ref_dfs_config,
                            args.dfs_config,
                            fitting_operation   = fitting_operations[args.fitting_operation],
                            fitting_string      = args.fit_selection,
                            scoring_function    = scores[args.conformational_difference],
                            scoring_function_name = args.conformational_difference,
                            do_write            = do_write,
                            score_kwargs        = score_kwargs,
                            options             = {"scale_forces_by_reference_cd" : False}
                        )

        dfs_job.run(details_output_fname=output_fname, metadata=attrs, nprocs=args.np)        

        score_matrix, score_matrix_min, score_matrix_max = dfs_job.get_score_matrix()

        if do_write and SCORES in do_write:
            fh = h5py.File(output_fname, 'r+')
            fh.create_dataset(SCORES, data=score_matrix)
            fh.close()

        np.savetxt("%s.dat" % output_scores, score_matrix, fmt="%.5f",)

    if args.force_mode == "fixed_displacements" or args.force_mode == "both":
        if do_write is not None:
            output_fname = "%s_%s" % ("fixed_displacements", args.output_fname)
        else:
            output_fname = None
        output_scores = "%s_%s" % ("fixed_displacements", args.output_scores)        

        dfs_job = DFSJob(
                        structure,
                        atoms,
                        atoms,
                        np.matrix(covariance),
                        args.ref_dfs_config,
                        args.dfs_config,
                        fitting_operation   = fitting_operations[args.fitting_operation],
                        fitting_string      = args.fit_selection,
                        scoring_function    = scores[args.conformational_difference],
                        scoring_function_name = args.conformational_difference,
                        do_write            = do_write,
                        score_kwargs        = score_kwargs,
                        options             = {"scale_forces_by_reference_cd" : True}
                    )

        dfs_job.run(details_output_fname=output_fname, metadata=attrs, nprocs=args.np)

        score_matrix, score_matrix_min, score_matrix_max = dfs_job.get_score_matrix()
       
	if do_write and SCORES in do_write:
            fh = h5py.File(output_fname, 'r+')
            fh.create_dataset(SCORES, data=score_matrix)
            fh.close()

        np.savetxt("%s.dat" % output_scores, score_matrix, fmt="%.5f",)

    log.info("All done! Exiting...")
