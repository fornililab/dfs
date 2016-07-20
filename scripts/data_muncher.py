    # data_muncher.py - Write data in .h5 files from DFS runs in easily readable formats
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

import h5py
import numpy as np
import logging as log
import prody as pd
from perturbation_magic import *

def traverse_tree(fh, current="/", data_types=[]):
    if not data_types:
        yield None
        return

    this_keys = filter(unicode.isdigit, fh[current].keys())

    for key in this_keys:
        for t in traverse_tree(fh, current=os.path.join(current,key), data_types=data_types):
            yield t

    for d in data_types:
        if d in fh[current].keys():
            yield (current, d, np.array(fh[current][d]))

def traverse_tree_best(fh, data_type, current="/"):
    if data_type == MAX_SCORE_FP_XYZ:
        coord_type = FP_XYZ
    elif data_type == MAX_SCORE_P_XYZ:
        coord_type = P_XYZ
    else:
        raise TypeError

    this_keys = filter(unicode.isdigit, fh[current].keys())


    for key in this_keys:
        for t in traverse_tree_best(fh,
                                    data_type = data_type,
                                    current=os.path.join(current,key)):
            yield t

    try:
        this_rs = np.array(fh[current][RAW_SCORES])
    except:
        yield None
        return

    current_ref = current.split("/")[1]

    this_coords_ref = np.array(fh[current_ref][coord_type])
    this_coords = np.array(fh[current][coord_type])

    n_elms = int(np.sqrt(this_coords.shape[0]))

    where = np.array(np.where(this_rs == np.max(this_rs)))
    where_idx = where[0][0]*n_elms + where[1][0]

    max_coords = this_coords[where_idx,:,:]
    ref_max_coords = this_coords_ref[where_idx,:,:]

    all_coords = np.array([ref_max_coords, max_coords])

    yield (current, data_type, all_coords)

def write_data(structure, atoms, atom_indices, d, fmt="%.5f", highest=False, lowest=False):
    current, data_type, data = d

    print data_type, 
    
    if current == '/':
        return

    this_idxs = map(int, current.split("/")[1:])
    #print np.where(atom_indices==this_idxs[0])[0]
    #print atoms.getIndices()[np.where(atom_indices==idx)[0]
    
    idx_matching = []
    for idx in this_idxs:
        idx_matching.append( atoms.getIndices()[np.where(atom_indices==idx)[0]][0] )

    this_atoms = [ structure[s] for s in idx_matching ]

    label = ",".join(["%s-%s%d" % (a.getChid(), a.getResname(), a.getResnum()) for a in this_atoms])

    if data_type in [DF, DR, F_DR, RAW_SCORES, ref_CD, dfs_CD]:
        if data_type == DF:
            fname = "force_vectors_%s.npz" % label
    
        elif data_type == DR:
            fname = "displacements_%s.npz" % label
    
        elif data_type == F_DR:
            fname = "fitted_displacements_%s.npz" % label

        elif data_type == RAW_SCORES:
            fname = "raw_scores_%s.npz" % label

        elif data_type == ref_CD:
            fname = "reference_cd_%s.npz" % label

        elif data_type == dfs_CD:
            fname = "dfs_cd_%s.npz" % label
            
        np.savetxt(fname, data, fmt=fmt)

    elif data_type in [P_XYZ, FP_XYZ]:
        if data_type == P_XYZ:
            fname = "displaced_structures_%s.dcd" % label

        elif data_type == FP_XYZ:
            fname = "fitted_displaced_structures_%s.dcd" % label

        ensemble = prody.Ensemble()

        lendata = int(np.sqrt(len(data)))    
        for i in range(lendata):
            for j in range(lendata):
                ensemble.addCoordset(data[j*lendata+i])

        prody.writeDCD(fname, ensemble)

    elif data_type in [MAX_SCORE_P_XYZ, MAX_SCORE_FP_XYZ]:

        if data_type == MAX_SCORE_P_XYZ:
            fname = "displaced_structures_%s_max.pdb" % label
            fname_ref = "displaced_structures_%s_max_ref.pdb" % label

        elif data_type == MAX_SCORE_FP_XYZ:
            fname = "fitted_displaced_structures_%s_max.pdb" % label
            fname_ref = "fitted_displaced_structures_%s_max_ref.pdb" % label

        this_structure = structure.copy()
        this_structure_ref = structure.copy()

        this_structure_ref.setCoords(data[0])
        this_structure.setCoords(data[1])

        prody.writePDB(fname_ref, this_structure_ref)
        prody.writePDB(fname, this_structure)


if __name__ == "__main__":

    import argparse

    ALL = "all"
    DF = "force_vectors"
    DR = "displacements"
    F_DR = "displacements_after_fit"
    P_XYZ = "perturbed_coordinates"
    FP_XYZ = "fitted_perturbed_coordinates"
    MAX_SCORE_P_XYZ = "max_score_perturbed_coordinates"
    MAX_SCORE_FP_XYZ = "max_score_fitted_perturbed_coordinates"
    SCORES = "score_matrix"
    RAW_SCORES = "raw_score_matrix"
    CD = "conformational_distance_matrixes"
    ref_CD = "reference_cdm"
    dfs_CD = "dfs_cdm"
    SF = "scaling_factors"
    
    data_types =    [
                        ALL,
                        DF,
                        DR,
                        F_DR,
                        P_XYZ,
                        FP_XYZ,
                        SCORES,
                        RAW_SCORES,
                        CD,
                        MAX_SCORE_P_XYZ,
                        MAX_SCORE_FP_XYZ,
                        SF
                    ]

    argparser = argparse.ArgumentParser()

    argparser.add_argument("HDF5", type=str, help="input HDF5 file")
    argparser.add_argument("-p", "--pdb", type=str, dest="PDB", help="initial PDB structure", default=None)
    argparser.add_argument("-v", "--verbose", dest="verbose", help="toggle verbose mode")
    argparser.add_argument("-w", "--write", dest="data_types", nargs='*', type=str, choices=data_types, help="Choose which data should be saved in the output file. 'all' (default behaviour) just saves everything.", default=[])
    #argparser.add_argument("-u", "--highest", dest="highest", default=None, action='store', help="Per run, save a PDB file with the structure having the highest score")
    #argparser.add_argument("-l", "--lowest", dest="lowest", default=None, action='store', help="Per run, save a PDB file with the structure having the lowest score")

    args = argparser.parse_args()

    if args.data_types == [] or (ALL in args.data_types):
        if len(args.data_types) > 1:
            log.warning("all was specified in option --write, together with others; all will be used.")
        this_data_types = data_types[1:]
    else:
        this_data_types = args.data_types
        
    if CD in this_data_types:
        this_data_types.extend([ref_CD, dfs_CD])

    log.basicConfig(level=log.DEBUG)

    hf = h5py.File(args.HDF5, 'r')

    print "File %s metadata:" % args.HDF5
    for k,v in hf.attrs.iteritems():
        print "\t%s\t\t:%s" % (k,v)

    if not args.PDB and not args.data_types:
        exit(0)
    sf_path = "/%s" % SF
    if sf_path in hf and SF in this_data_types:
        sf = 1.0/np.array(hf[sf_path])
	sf_shape = np.int(np.sqrt(sf.shape[1]))
        np.savetxt("scaling_factors.dat", sf[:,0:sf_shape], fmt="%.5f")
        this_data_types.remove(SF)
        
    available_data = hf.attrs['available_data'].split(",")

    if MAX_SCORE_P_XYZ in this_data_types and (RAW_SCORES not in available_data or P_XYZ not in available_data):
        log.error("%s requires %s and %s to be present in details file. Exiting ..." % (MAX_SCORE_P_XYZ, RAW_SCORES, P_XYZ))
    if MAX_SCORE_FP_XYZ in this_data_types and (RAW_SCORES not in available_data or FP_XYZ not in available_data):
        log.error("%s requires %s and %s to be present in details file. Exiting ..." % (MAX_SCORE_FP_XYZ, RAW_SCORES, FP_XYZ))

    # 1. Parse PDB
    try:
        structure = prody.parsePDB(args.PDB)
    except:
        log.error("Couldn't parse your model file (%s); Exiting..." % args.PDB)
        exit(1)

    # 2. select atoms for building the model
    try:
        atoms = structure.select(hf.attrs['atom_selection'])

    except prody.SelectionError:
        log.error("Your selection string was invalid. Exiting...")
        exit(1)
    if not atoms:
        log.error("Your selection string resulted in no atoms being selected in your structure. Exiting...")
        exit(1)
    if len(atoms) < 2:
        log.error("Your selection string resulted in less than two atoms being selected in your structure. Exiting...")
        exit(1)

    atom_indices = atoms.getResindices()

    if P_XYZ in this_data_types:
        p_structures = prody.Ensemble()

    if FP_XYZ in this_data_types:
        fp_structures = prody.Ensemble()

    if MAX_SCORE_P_XYZ in this_data_types:
        for d in traverse_tree_best(hf, data_type = MAX_SCORE_P_XYZ):
            if d:
                write_data(structure, atoms, atom_indices, d, highest=True, lowest=False)
        this_data_types.remove(MAX_SCORE_P_XYZ)

    if MAX_SCORE_FP_XYZ in this_data_types:
        for d in traverse_tree_best(hf, data_type = MAX_SCORE_FP_XYZ):
            if d:
                write_data(structure, atoms, atom_indices, d, highest=True, lowest=False)
        this_data_types.remove(MAX_SCORE_FP_XYZ)                

    if this_data_types:
        for d in traverse_tree(hf, data_types = this_data_types):
            write_data(structure, atoms, atom_indices, d)


