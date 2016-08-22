######################################################################
# protein: third immunoglobulin binding domain of Protein G (GB3) 
# PDB ID: 2LUM

######################################################################
# path and configuration settings
dfs=`which dfs`

# pathogenic mutation site (PM) configuration files
pm_config=../T00/config_file_pathogenic_mut.ini
# secondary mutation site (PM) configuration files
sm_config=../T00/config_file_secondary_mut.ini

# configuration file settings
cfg_settings="-r $pm_config -f $sm_config"

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
