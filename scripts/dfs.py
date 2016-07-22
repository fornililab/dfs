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
import argparse
from argparse import RawTextHelpFormatter
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
import StringIO


intro = """dfs implements the Dual Force Scanning method for the identification 
of compensatory mutations in proteins as detailed in:

<paper>

The method uses Linear Response Theory on Anisotropic Network Models (ANM)
to calculate a rescuability score R between pairs of residues, one that is 
considered a Pathogenic Mutation (PM) to be compensated and another that is a 
Secondary  Mutation site, which represents a potential compensatory site. 
Positive values of R indicate a compensatory effect. By default, every residue 
is tested against each other both as a PM and as a SM; the outcome is a NxN 
matrix of rescuability score, where N is the number of residues. By default 
all pairs are calculated.

Running with the default options, dfs just requires a single PDB input file. 
It is strongly advised to edit the PDB file so that it contains only
the atoms that will be used in the calculation (i.e. the CA atoms of each
residue), in order to remove possible spurious "CA" atoms (as for instance 
calcium ions or atoms from cofactors) that would be otherwise considered as 
protein residues by the software.

After this, running dfs is just a matter of:

./dfs input.pdb

the DFS method is run in two variants: one in which all the forces are kept
with the same magnitude (fixed force) and the other in which magnitues are 
changed as to obtain constant displacements (fixed displacements). Both of them
are run by default and thus the software return two ASCII .dat matrix files
containing one or the other, and named accordingly.

While the default options should be adequate for most cases, it is possible
to change them as follows. Further description of the available options is
presented after them.
"""

outro = """In the DFS method, the ANM is first calculated from the protein 
structure. In ANM, each residue is represented by the CA (alpha carbon) atom
alone. You can select which atoms should be included in the model by using
option -s, which uses CHARMM-style selections (as in ProDy) to select which
atoms/residues need to be included. By default, all atoms named "CA" take part
in the model. The parameters for the model (distance cut-off to connect atoms
in the model and the spring constant of the springs between them) can be 
supplied using -c and -g. Notice that, while -c can somehow change the outcome,
varying -g only has the effect of scaling the displacements linearly and thus
does not substantially affect the rescuability score. It is possible to provide
an external hessian or covariance matrix instead of calculating it in the
program (options -M and -I) as well as save the hessian or covariance matrices
for further reference (options -m and -i).

For DFS, forces need to be defined in terms of i) number of forces ii) their 
direction iii) their magnitude. These parameters can be set through 
configuration files as detailed in the examples directory. It should be noted 
that, according to our tests, the score values converge already with 12 force 
application when their direction is generated with the Fibonacci lattice 
approach. The user can supply two configuration files: one for the forces 
applied to the pathogenic mutations (option -f) and one for the forces applied 
to secondary sites (option -r). If one or both these files are not supplied,
default options will be used as detailed at the bottom of this text.

Once the forces have been applied, the calculation of the rescuability score
depends on the conformational distance between the original, unperturbed
structure and the perturbed structures. To this extent, two conformation
measures are provided through option -d: root mean square deviation (RMSD) 
and distance root mean square deviation (dRMSD). RMSD is usually considered
more intuitive and it's simpler to calculate, however it has the disadvantage
of requiring best-fit superimposition between the two structures, and sometimes
isolating the protein groups or domains that give rise to meaningful 
superimpositions can be tricky (especially when one part of the protein changes 
much more than the other e.g. in a multi-domain protein). Using option -F
it is possible to require for least-square superimposition to be performed
prior to RMSD calculation, and options -a allows to specify which selection
(in CHARMM style selection string) needs to be used for fitting. No 
superimpisition is required for dRMSD, which is used by default.
It should noted that the atoms on which the force is applied are by default not
considerend in the calculation as their displacement is not really meaningful 
for the calculation of the rescuability score. It is however possible to 
include them in the calculation by using option -x.
 
Finally, the rescuability score can be calculated according to two different
stragies, which can be selected using option -t. One is to just keep the force 
magnitude constant for all atoms and directions (-t fixed_forces), and the 
other is to normalize force magnitudes so that they give rise to constant
conformational distance values (close to 1A). Both modes can also be performed
in the same run, one after the other (-t both), and the runtime will be 
roughly twice as long than in the single case. This is the default.

As far as output formats are concerned, dfs writes two types of output files.
One is the ASCII numeric NxN matrix, as previously described. The basename
for this file can be changed with option -S. The other
is a binary file in the HDF5 portable format (https://www.hdfgroup.org/) that 
contains relevant details about the calculation. By default, this file is not
written; the user needs to specify the output filename with option -o for 
that to happen. The content of the file can be tuned with option -w, which
determines what will be actually written in the file. The following type of 
data are supported:

* force_vectors: the force vectors used in Linear Response Theory
displacements: atomic displacements induced by the force vectors
displacements_after_fit: atomic displacements between the perturbed and non-
* perturbed structure after least-square fitting before them has been 
performed. Only meaningful if fitting is performed.
* perturbed_coordinates: atomic coordinates of the original structure to which
atomic displacements have been applied. This is the default option for -w.
* fitted_perturbed_coordinates: atomic coordinates of the original structure to 
which atomic displacements and least square fitting have been applied.
* score_matrix: Rescuability score matrix
* raw_score_matrix: score matrix per every residue pair, before the final R 
value is calculated (see the paper for details)
* conformational_distance_matrixes: conformational distance values on which
the raw_score_matrix is calculated.
* scaling_factors: scaling factors used for force magnitudes
* all: all of the above

More than one option can be specified with -w (e.g. "-w force_vectors 
displacements" is valid). Option -q helps setting the accuracy/space trade-off
of the compression algorithm that is used to store data, in particular by 
by setting the number of decimal places that will be reliably stored in the 
file (3 is usually enough). Setting this option to -1 means lossless 
compression. It should be noted that, as calculations are usually 
relatively fast respect to file writing disk speeds, it often happens that file 
writing often lags behind computation and thus time has to be spent on waiting 
for the write operations to finish before proceeding on the calculation, 
in order to avoid unnecessary build-up of system memory usage.

Calculations can be sped up by using more than one core at the time by means of
option --np, which sets the number of cores to be used. If -o and -w are enabled
an additional process will be spawned to write the compressed .h5 file, and 
this should be taken into account when setting --np. 

By default, DFS will be run with the following options:

ANM model
    atom selection (-s): 'name CA'
    distance cut-off (-c): 15.0 A
    force constant (-g): 1.0
Force application
    number of forces: 12
    directions of forces: Fibonacci lattice
    Force magnitude: 10.0 UNIT
    selection of atoms on which forces will be applied: 'name CA'
Calculation of the rescuability score
    conformational distance measure (-d): dRMSD (distance RMSD)
    fitting operation (-F): none
    fitting operation selection string (-a): none
    include force application site for the calculation of dRMSD (-x): no
    force mode (-t): constant force and constant displacement
 
"""


scores =    {   "rmsd": ScoreRMSD,
                "drmsd": ScoreDRMSD,
            }

fitting_operations =    {   "lsq_fit" : LsqFitting,
                            "none" : None
                        }

default_ref_config_string="""[ atom sets ]

refset1:   name CA

[ force sets ]

refset1:   fibonacci_lattice 12 10.818 r"""

default_config_string="""[ atom sets ]

set2:   name CA

[ force sets ]

set2:   fibonacci_lattice 12 10.818 c"""

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=intro,
                                     epilog=outro,
                                     formatter_class=RawTextHelpFormatter)

    parser.add_argument("PDB", type=str, help="input PDB file")
    parser.add_argument("-s", "--selection", dest="selection_string", type=str, 
                        default="name CA", help="CHARMM-style selection of \
                        atoms to be used to build the ANM model. The final \
                        selection should contain at most one CA atom per \
                        residue.")
    parser.add_argument("-v", "--verbose", dest="verbose", 
                        help="toggle verbose mode")
    parser.add_argument("-c", "--cutoff", dest="anm_cutoff", default=15.0, \
                        type=float, help="Distance cut-off for generating the \
                        topology of ANM (in A; default is 15.0)")
    parser.add_argument("-g", "--gamma", dest="anm_gamma", default=0.1, 
                        type=float, help="Spring force constant for ANM; \
                        default is 1.0)")
    parser.add_argument("-m", "--write-hessian", dest="hessian_output", 
                        nargs='?', const="hessian.dat", type=str, default=None, 
                        help="Save Hessian matrix of the ANM model (default: \
                        matrix.dat)")
    parser.add_argument("-M", "--load-hessian", dest="hessian_input", type=str, 
                        default=None, help="Load hessian matrix from file \
                        instead of calculating it")
    parser.add_argument("-i", "--write-covariance", dest="covariance_output", 
                        nargs='?', const="covariance.dat", type=str, 
                        default=None, help="Save the covariance matrix from \
                        the ANM used to calculate the effect of force \
                        application (default: covariance.dat)")
    parser.add_argument("-I", "--load-covariance", dest="covariance_input", 
                        type=str, default=None, help="Load covariance \
                        matrix instead of calculating it. The hessian matrix \
                        will not be calculated if this flag is used.")
    #parser.add_argument("-b", "--eigenvectors", dest="nm_eigvec_output", nargs='?', const=1, type=str, default="nmodes_eigenvectors.dat", help="Save normal modes eigenvectors as a text file")
    #parser.add_argument("-l", "--eigenvalues", dest="nm_eigval_output", nargs='?', const=1, type=str, default="nmodes_eigenvalues.dat", help="Save normal modes eigenvalues as a text file")
    #parser.add_argument("-D", "--nmd", dest="nmd_output", nargs='?', const=1, 
    #                    type=str, default="modes.nmd", help="Save normal mode 
    #                    eigenvectors as a NMD file")
    parser.add_argument("-f", "--dfs", dest="dfs_config", default=None, 
                        help="Configuration file for DFS runs (Secondary \
                        sites. Defaults will be used if such a file is not \
                        specified.)")
    parser.add_argument("-r", "--ref-dfs", dest="ref_dfs_config", default=None,
                        help="Configuration file for DFS runs (Pathogenic \
                        Mutation sites). Defaults will be used if such a file \
                        is not specified.")
    parser.add_argument("-F", "--fitting", dest="fitting_operation", 
                        default="none", action="store", type=str, 
                        choices=fitting_operations.keys(), help="perform least\
                         square fitting between the displaced and original \
                         geometries")
    parser.add_argument("-a", "--fit-selection", dest="fit_selection", 
                        default="name CA", help="selection string that selects\
                        the group that will be used to superimpose the \
                        modified structures to the original one.")
    parser.add_argument("-d", "--conformational-distance", 
                        dest="conformational_distance", choices=scores.keys(), 
                        help="which conformational distance measure will be \
                        used (default: drmsd)", default="drmsd")
    parser.add_argument("-n", "--np", dest="np", default=1, type=int, 
                        action="store", help="Number of cores to be used \
                        for calculation")
    parser.add_argument("-w", "--write", dest="write", nargs='*', type=str, 
                        choices=writable_data, help="Choose which data should \
                        be saved in the output details file. 'all' just saves \
                        everything. Default is %s. This option will be ignored\
                         if -o is not supplied" % P_XYZ, default=P_XYZ)
    parser.add_argument("-o", "--output", dest="output_fname", 
                        default=None, type=str, action="store", 
                        help="Filename to be used as output in case -w/--write\
                         is specified")
    parser.add_argument("-t", "--force_mode", dest="force_mode", nargs='?', 
                        type=str, default="both", choices=["fixed_forces, \
                        fixed_displacementsfo, both"], help="Force mode to be \
                        used (default: both)")
    parser.add_argument("-S", "--output-scores", dest="output_scores", 
                        default="scores", type=str,  action="store", 
                        help="Basename for the output rescuability score matrix")
    parser.add_argument("-q", "--precision", dest="precision", default=-1, 
                        type=int, help="Number of decimal places to be used \
                        in .h5 detail files (-1 for lossless)")
    parser.add_argument("-x", "--include-application-sites", 
                        dest="exclude_sites", action='store_false', 
                        default=True, help="Include force application sites \
                        in the calculation of scores (by default it doesn't)")
    #XXX add option: number of normal modes

    log.basicConfig(level=log.DEBUG)

    args = parser.parse_args()

    if args.output_fname is None:
        do_write = None
    elif args.write == [] or (ALL in args.write):
        if len(args.write) > 1:
            log.warning("all was specified in option --write, together with others; all will be used.")
        do_write = writable_data[1:]
    else:
        do_write = args.write
    if do_write is not None:
        log.info("this information will be written to file: %s" % (", ".join(list(do_write))))

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
        #safe_savenmd(args.nmd_output, anm, atoms)

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

        attrs['cdm'] = args.conformational_distance
        attrs['fitting'] = args.fitting_operation
        attrs['fitting_selection'] = args.fit_selection
        attrs['available_data'] = ",".join(do_write)


    else:
        fh = None
        attrs = None



    score_kwargs = { "exclude_sites" : args.exclude_sites,
                    }

    if args.force_mode == "fixed_force" or args.force_mode == "both":

        if args.ref_dfs_config is None:
            log.info("Default settings will be used for forces on pathogenic \
    mutation sites as no input file was provided.")
            ref_config = StringIO.StringIO(default_ref_config_string)
        else:
            ref_config = args.ref_dfs_config            
        if args.dfs_config is None:
            log.info("Default settings will be used for forces on secondary \
    mutation sites as no input file was provided.")
            config =     StringIO.StringIO(default_config_string)
        else:
            config = args.dfs_config
            
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
                            ref_config,
                            config,
                            fitting_operation   = fitting_operations[args.fitting_operation],
                            fitting_string      = args.fit_selection,
                            scoring_function    = scores[args.conformational_distance],
                            scoring_function_name = args.conformational_distance,
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

        if args.ref_dfs_config is None:
            log.info("Default settings will be used for forces on pathogenic \
    mutation sites as no input file was provided.")
            ref_config = StringIO.StringIO(default_ref_config_string)
        else:
            ref_config = args.ref_dfs_config            
        if args.dfs_config is None:
            log.info("Default settings will be used for forces on secondary \
    mutation sites as no input file was provided.")
            config =     StringIO.StringIO(default_config_string)
        else:
            config = args.dfs_config

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
                        ref_config,
                        config,
                        fitting_operation   = fitting_operations[args.fitting_operation],
                        fitting_string      = args.fit_selection,
                        scoring_function    = scores[args.conformational_distance],
                        scoring_function_name = args.conformational_distance,
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
