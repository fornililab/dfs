######################################################################
# protein: third immunoglobulin binding domain of Protein G (GB3) 
# PDB ID: 2LUM

######################################################################
# path and configuration settings
compensatory_power=`which compensatory_power`

# input matrix
matrix="../T01/fixed_force_scores_T01.dat"

# structure used in the dfs step
pdb="../T00/original.pdb"

# input file settings
input_settings="-m $matrix -p $pdb"

######################################################################
# output settings
# percentile for top compensatory mutations
percentile="-x 95"

# output prefix
outprefix="T02"

# output filename
outputscore="-o $outprefix.fixed_force_normalized_score.dat"

# output pdb filename
outputpdb="-P $outputprefix.fixed_force_normalized_score.pdb" 

# output_settings
output_settings="$outputscore $outputpdb $percentile"

######################################################################
# check availability of input file
if [ ! -e $matrix ]
then
        echo "TEST T02: NOT EXECUTED - REQUIRES T01"
        exit
fi

######################################################################
# command line
$compensatory_power $input_settings $output_settings >& log.T02.log

######################################################################
# test check
cp T02.fixed_force_normalized_score.dat ../check/T02.out
diff ../check/T02.out ../check/T02.chk > ../check/T02.diff
if [ -s ../check/T02.diff ]
then
        echo "TEST T02: FAILED"
else
        echo "TEST T02: PASSED"
fi
