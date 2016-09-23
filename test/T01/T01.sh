######################################################################
# protein: third immunoglobulin binding domain of Protein G (GB3) 
# PDB ID: 2LUM

######################################################################
# path and configuration settings
dfs=`which dfs`

# forces configuration files
config=../T00/forces.ini

# configuration file settings
cfg_settings="-f $config"

# structure to be used to generate the Anisotropic Network Model (ANM)
pdb="../T00/original.pdb"

# input file settings
input_settings="-p $pdb"

######################################################################
# ANM settings
# distance cut-off for ANM
c=15.00000

# force constant for ANM harmonic potential (gamma)
g=0.10000

# number of cores to be used 
np=2

# anm settings
anm_settings="-c $c -g $g --np $np"

######################################################################
# output settings
# write options - what to write (-w) and in which file
write="-w perturbed_coordinates raw_score_matrix -o details.h5"

# presision for details.h5 file
precision="-q 3"

# output suffix
outsuffix="-S scores_T01"

# output_settings
output_settings="$write $precision $outsuffix"
######################################################################
# command line
$dfs $cfg_settings $anm_settings -s "name CA" $output_settings $input_settings >& log.T01.log

######################################################################
# test check
cksum fixed_*_scores_T01.dat > ../check/T01.out
diff ../check/T01.out ../check/T01.chk > ../check/T01.diff
if [ -s ../check/T01.diff ]
then
        echo "TEST T01: FAILED"
else
        echo "TEST T01: PASSED"
fi
