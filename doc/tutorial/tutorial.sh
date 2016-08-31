# this is an example for a complete application of DFS, and as in "tests" 
# uses the third immunoglobulin binding domain of Protein G (GB3). This includes
# running the dfs script to get the rescuability matrices and then the 
# compensatory_power script to calculate the rescaubility scores.


#################################################################################
# 1 - use the Dual Force Scanning method to generate rescuability scores

# First we run the dfs script that will generate the rescuability
# score matrices. We will not save any additional details file
# for the moment as it's not strictly required for the successive
# analyses.

# dfs script - to be changed according to its position in your system
scripts_dir=../../scripts_src

dfs=$scripts_dir"/dfs"
compp=$scripts_dir"/compensatory_power"

# structure to be used to generate the Anisotropic Network Model (ANM)
pdb="1TSR.pdb"

# selection of atoms to be used for the ANM model. We're going to use
# CA atoms from chain A from p53's multi-chain PDB 1TSR

selection="chain A and protein and name CA"

# distance cut-off for ANM
c=15.0

# force constant for ANM harmonic potential (gamma)
g=1.0

# pathogenic mutation site (PM) and Secondary mutation site (SM) 
# configuration files. Read the files themselves for information
# on their format and syntax. 
pm_config=config_file_pathogenic_mut.ini
sm_config=config_file_secondary_mut.ini

# number of cores to be used 
np=8

# command line
echo now running: $dfs -p $pdb -r $pm_config -f $sm_config -c $c -g $g -s $selection --np $np -v
$dfs -p $pdb -r $pm_config -f $sm_config -c $c -g $g -s $selection --np $np

# Two matrix files have been saved: one that contains outcome of the calculation
# in constant force mode and the other in constant displacement mode.


#################################################################################
# 2 - calculate rescuability power index

# We feed both matrices calculated in the previous step to the compensatory_power
# script in order to calculate the CP index on the average matrix:

echo "now running: $compp -m fixed_force_scores.dat fixed_displacements_scores.dat -p $pdb -s $selection"
$compp -m fixed_force_scores.dat fixed_displacements_scores.dat -p $pdb -s $selection

# The output file normalized_score.dat contains the per-residue Compensatory
# power score.
