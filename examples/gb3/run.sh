# this is an example for a DFS run, which uses the third immunoglobulin 
# binding domain of Protein G (GB3) as a technical test (PDB ID 2LUM).
# We will define several parameters here # and then add them to the 
# command line. Each parameter is commented.

# dfs script - to be changed according to its position in your system
dfs=$HOME/.local/bin/dfs

# structure to be used to generate the Anisotropic Network Model (ANM)
pdb="original.pdb"

# distance cut-off for ANM
c=15.00000

# force constant for ANM harmonic potential (gamma)
g=0.10000

# pathogenic mutation site (PM) and Secondary mutation site (SM) 
# configuration files. Read the files themselves for information
# on their format and syntax. 
pm_config=config_file_pathogenic_mut.ini
sm_config=config_file_secondary_mut.ini

# number of cores to be used 
np=8

# write options - what to write (-w) and in which file
write="-w perturbed_coordinates raw_score_matrix -o details.h5"

# presision (number of decimal places) to be used for the data
# stored in the details.h5 file
precision="-q 3"

# command line
$dfs -r $pm_config -f $sm_config -c $c -g $g -s "name CA" $write $precision --np $np $pdb



