import numpy as np
import prody as pd
from prody import LOGGER as log

def parse_mutfile(fname, structure):

    try:
        muts = np.loadtxt(fname, dtype=np.int)
        return muts
    except ValueError:
        pass

    try:
        muts = np.loadtxt(fname, dtype=str)
    except ValueError:
        log.error("Couldn't load mutation list file. Exiting...")
        exit(1)

    out = []
    for m in muts:
        selection = structure.select("chain %s and resnum %s and name CA" % (m[0], m[1:]))
        if selection is None:
            log.warning("Mutation %s will be skipped!" % str(m))
            continue
        this_out = selection.getResindices()
        assert(this_out.shape[0] == 1)
        out.append(this_out[0])

    return np.array(out, dtype=np.int)

def parse_range_option(range_string, validity_range=(-np.inf, np.inf)):
    tmp = range_string.strip().split(",")
    tmp = map(float, tmp)

    if tmp[0] < validity_range[0] or tmp[0] > validity_range[1]:
        raise ValueError

    if len(tmp) == 1:
        return np.array(map(float, tmp), dtype=np.float)
    elif len(tmp) == 3:
        if tmp[1] < validity_range[0] or tmp[1] > validity_range[1]:
            raise ValueError
        if tmp[1] >= tmp[0]:
            return np.arange(tmp[0], tmp[1], tmp[2])
        else:
            return np.arange(tmp[0], tmp[1], tmp[2])[::-1]
    else:
        raise IndexError

def get_contacts(s, dist_co):
    return np.sum(pd.buildDistMatrix(s) <= dist_co, axis=0, dtype=np.float)-1

def get_raw_index(m, score_co, axis=0):
    return np.sum(m >= score_co, axis=axis, dtype=np.float)

if __name__ == '__main__':

    import argparse
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    parser = argparse.ArgumentParser(description="Calculate enrichment for dfs and verified compensatory mutations")
    parser.add_argument(dest="score_matrix", type=str, action='store', help="Score matrix file")
    parser.add_argument(dest="pdb", type=str, action='store', help="PDB file filename")
    parser.add_argument(dest="mutlist", action='store', help="Mutation list file")
    parser.add_argument('-c','--score-cutoff', dest="score_cutoff", type=str, action='store', default="0.0", help="Cut-off to detect rescuing positions in the matrix (default: 0.0)")
    parser.add_argument('-d','--distance-cutoff', dest="distance_cutoff", type=float, action='store', default=15.0, help="Distance cut-off for the contact matrix (A; detault: 15.0)")
    parser.add_argument('-p','--percentile-cutoff', dest="percentile_cutoff", type=str, action='store', default="5.0", help="topX value for enrichment score (default: 5.0)")
    parser.add_argument('-t', '--score-type', type=str, dest="type", action='store', choices=["rescue_ability", "rescuability"], default="rescue_ability", help="Use rescue_ability or rescuability (default: rescue_ability)")
    parser.add_argument('-s', '--selection-string', dest="sel_string", action='store', type=str, default="name CA", help=("Selection string for the PDB file (default: name CA)"))
    parser.add_argument('-o', '--output', type=str, dest="outfile", action='store', default="enrichment_matrix.dat", help="Name for the output file (default: normalized_score.dat)")
    parser.add_argument('-f', '--plot-matrix', type=str, dest="plot_outfile", action='store', default="enrichment_matrix.pdf", help="Name for the plot file (default: enrichment_matrix.pdf)")
    parser.add_argument('-P', '--percentile-score',  dest="score_percentile", action='store_true', default=False, help="Use percentile to calculate score cut-offs as well")
    parser.add_argument('-q', '--positive',  dest="positive", action='store_true', default=False, help="Use only positive values from the matrix. Works only in combination with -P")
    parser.add_argument('-n', '--title',  dest="title", action='store', default=None, type=str, help="Name of the system - to be used as title on plots")    
    
    args = parser.parse_args() 

    if args.type == 'rescue_ability':
        axis = 0
    else:
        axis = 1

    try:
        score_matrix = np.loadtxt(args.score_matrix)
    except IOError:
        log.error("Matrix file not found or not readable). Exiting...")
    except ValueError:
        log.error("Matrix file is in the wrong format. Exiting...")

    structure = pd.parsePDB(args.pdb)

    if structure is None:
        log.error("Matrix file not found or not readable. Exiting...")

    try:
        selection = structure.select(args.sel_string)
    except pd.SelectionError:
        log.error("The selection string is not valid. Exiting...")

    muts = parse_mutfile(args.mutlist, structure)
    print muts

    contacts = get_contacts(selection, args.distance_cutoff)

    ra_powers = []                
        
    try:
        percentile_range = parse_range_option(args.percentile_cutoff, validity_range=(0,100))
        #percentile_range = 100.0 - percentile_range
        #percentile_range = percentile_range[::-1]
    except ValueError:
        log.error("topX range must be within 0 and 100. Exiting...")
        exit(1)
        
    if not args.score_percentile:
        score_cutoff_range = parse_range_option(args.score_cutoff)
    else:
        score_cutoff_range = percentile_range
        score_cutoffs = []
        

    for sco in score_cutoff_range:
        if args.score_percentile:
            if args.positive:
                this_score_matrix = score_matrix[score_matrix >= 0.0]
            else:
                this_score_matrix = score_matrix
            sco = np.percentile(this_score_matrix.flatten(), 100.0-sco, interpolation='midpoint')
            score_cutoffs.append(sco)

        raw_index = get_raw_index(score_matrix, sco, axis=axis)

        ra_powers.append(raw_index/contacts)

    if not args.score_percentile:
        log.info("Will use score cut-offs: %s" % (", ".join(map(str, score_cutoff_range))))
    else:
        #print ", ".join(["%.1f (%.3f)" % (i,j) for i,j in zip(percentile_range, score_cutoffs)])
        log.info("Will use topX scores and related score cut-offs: %s " % ", ".join(["%.1f (%.3f)" % (i,j) for i,j in zip(percentile_range, score_cutoffs)]))
        
    log.info("Will use topX values: %s" % (", ".join(["%.1f" % i for i in percentile_range])))
    
    enrichment_matrix = np.zeros((percentile_range.shape[0], score_cutoff_range.shape[0]))

    for i,p in enumerate(percentile_range):
        for j,rp in enumerate(ra_powers):
            perc = np.percentile(rp, 100-p, interpolation='midpoint')
            top = np.where(rp >= perc)[0]
            #print p, score_cutoff_range[j]
            #print "mut in top", np.float(np.intersect1d(top, muts).shape[0])
            #print "size top", np.float(top.shape[0])
            #print "size mut", np.float(muts.shape[0])
            #print "size all", np.float(rp.shape[0])
            #print "--"
        
            enrichment_matrix[i,j] = (np.float(np.intersect1d(top, muts).shape[0]) /\
                                        np.float(top.shape[0]))/\
                                        (np.float(muts.shape[0])\
                                        / np.float(rp.shape[0]))

    try:
        np.savetxt(args.outfile, enrichment_matrix, fmt="%.5f")
    except:
        log.error("Couldn't write file %s, Exiting..." % args.outfile)

    if not args.score_percentile:
        score_cutoffs = np.array(score_cutoff_range)
        zero_line = score_cutoffs[np.where(score_cutoffs == 0.0)][0]

    else:
        score_cutoffs = np.array(score_cutoffs)
        zero_line = 100.0 - score_cutoff_range[np.sum(score_cutoffs < 0)+1]

    imshow_extent = (score_cutoff_range[0], score_cutoff_range[-1],
                     percentile_range[-1], percentile_range[0])
    fig = plt.imshow(enrichment_matrix, interpolation='none', extent=imshow_extent,
               aspect="auto", vmin=0.0, vmax=2.0, origin='upper')
    plt.title(args.title)        

    plt.axvline(x=zero_line, ymin=-100, ymax=100, color='black')
    plt.axhline(y=5.0, xmin=-100, xmax=100, color='black')    

    plt.gca().invert_yaxis()

    if args.score_percentile:
        plt.gca().invert_xaxis()

    plt.colorbar(fig)
    plt.ylabel("top X%")
    plt.xlabel("Score cut-off")
    plt.savefig(args.plot_outfile)
    
    pos_idx = np.where(score_cutoffs >= 0)[0]
    neg_idx = np.where(score_cutoffs <  0)[0]
        
    if pos_idx.shape[0] > 0:

        pos_em = enrichment_matrix[:,pos_idx]
        
        top_pos_value = np.where(pos_em == np.max(pos_em))
        
        for tpv in zip(*top_pos_value):
            if not args.score_percentile:
                log.info("Highest enrichment value for positive score cut-offs: %.5f for score \
cut-off %.5f and topX score %.5f" % (pos_em[tpv[0],tpv[1]],
                                         score_cutoffs[pos_idx][tpv[1]],
                                         percentile_range[tpv[0]]))
            else:    
                log.info("Highest enrichment value for positive score cut-offs: %.5f for \
percentile %.3f (score is %.5f) and topX score %.5f" % (pos_em[tpv[0],tpv[1]],
                                            score_cutoff_range[pos_idx][tpv[1]],
                                            score_cutoffs[pos_idx][tpv[1]],
                                            percentile_range[tpv[0]]))
    else:
        log.info("No positive score cut-offs.")
        
    if neg_idx.shape[0] > 0:
    
        neg_em = enrichment_matrix[:,neg_idx]        
        
        top_neg_value = np.where(neg_em == np.max(neg_em))
        #print top_neg_value

        for tnv in zip(*top_neg_value):
        
            if not args.score_percentile:
                log.info("Highest enrichment value for positive score cut-offs: %.5f for score \
cut-off %.5f and topX% score %.5f" % (neg_em[tnv[0],tnv[1]],
                                                              score_cutoffs[neg_idx][tnv[1]],
                                                              percentile_range[tnv[0]]))
            else:
               log.info("Highest enrichment value for negative score cut-offs: %.5f for \
percentile %.3f (score is %.5f) and topX score %.5f" % (neg_em[tnv[0],tnv[1]],
                                                         score_cutoff_range[neg_idx][tnv[1]],
                                                         score_cutoffs[neg_idx][tnv[1]],
                                                         percentile_range[tnv[0]]))
    else:
        log.info("No negative score cut-offs.")
