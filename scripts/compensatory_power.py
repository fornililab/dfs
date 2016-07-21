import numpy as np
import prody as pd
from prody import LOGGER as log
import logging as log

def get_contacts(s, dist_co):
	return np.sum(pd.buildDistMatrix(s) <= dist_co, axis=0, dtype=np.float)-1

def get_raw_index(m, score_co, axis=0):
	return np.sum(m >= score_co, axis=axis, dtype=np.float)

if __name__ == '__main__':

	import argparse

	parser = argparse.ArgumentParser(description="Calculate \"Rescue ability power\" or \"rescuability power\" index")
	parser.add_argument(dest="score_matrix", type=str, action='store', help="Score matrix filename")
	parser.add_argument(dest="pdb", type=str, action='store', help="PDB file filename")
	parser.add_argument('-t', '--score-type', type=str, dest="type", action='store', choices=["rescue_ability", "rescuability"], default="rescue_ability", help="Use rescue_ability or rescuability (default: rescue_ability)")	
	parser.add_argument('-s', '--selection-string', dest="sel_string", action='store', type=str, default="name CA", help=("Selection string for the PDB file (default: name CA)"))
	parser.add_argument('-c','--score-cutoff', dest="score_cutoff", type=float, action='store', default=0.0, help="Cut-off to detect rescuing positions in the matrix (default: 0.0)")
	parser.add_argument('-d','--distance-cutoff', dest="distance_cutoff", type=float, action='store', default=15.0, help="Distance cut-off for the contact matrix (A; detault: 15.0)")	
	parser.add_argument('-o', '--output', type=str, dest="outfile", action='store', default="normalized_score.dat", help="Name for the output file (default: normalized_score.dat)")
	parser.add_argument('-p', '--pdb-output', type=str, dest="pdb_outfile", action='store', default=None, help="Name for the output PDB file (default: don't save)")
	parser.add_argument('-x', '--pdb-percentile-output', type=float, dest="percentile", action='store', default=None, help="Percentile values to be written in occupancy column")


	args = parser.parse_args()

	if args.type == 'rescue_ability':
		axis = 0
	else:
		axis = 1

	try: 
		score_matrix = np.loadtxt(args.score_matrix)
	except IOError:
		log.error("Matrix file not found or not readable. Exiting...")
		exit(1)
	except ValueError:
		log.error("Matrix file is in the wrong format. Exiting...")
		exit(1)
	
	structure = pd.parsePDB(args.pdb)

	if structure is None:
		log.error("Matrix file not found or not readable. Exiting...")

	try:
		selection = structure.select(args.sel_string)
	except pd.SelectionError:
		log.error("The selection string is not valid. Exiting...")

	contacts = get_contacts(selection, args.distance_cutoff)

	raw_index = get_raw_index(score_matrix, args.score_cutoff, axis=axis)

	ra_power = raw_index/contacts

	try:
		np.savetxt(args.outfile, ra_power, fmt="%.3f")
	except:
		log.error("Couldn't write file %s, Exiting..." % args.outfile)

	res_indices = [ a.getResindex() for a in structure]
	
	betas = []
	this_i = -1
	k = -1
	for i,idx in enumerate(res_indices):
		if idx != this_i:
			k += 1
			this_i = idx
		betas.append(ra_power[k])

        if args.percentile:
	        this_i = -1
        	k = -1
                occupancies = []
                cutoff = np.percentile(ra_power, args.percentile, interpolation='midpoint')
                cutoff_bool = ra_power >= cutoff
                for i,idx in enumerate(res_indices):
                        if idx != this_i:
                                k += 1
                                this_i = idx
                        occupancies.append(np.float(cutoff_bool[k]))
        
                

	if args.pdb_outfile is not None:
                if args.percentile:
		        try:
			        pd.writePDB(args.pdb_outfile, structure, beta=betas, occupancy=occupancies)
		        except:
			        log.error("Couldn't write file %s, Exiting..." % args.pdb_outfile)
                else:
                        try:
                                pd.writePDB(args.pdb_outfile, structure, beta=betas)
                        except:
                                log.error("Couldn't write file %s, Exiting..." % args.pdb_outfile)



