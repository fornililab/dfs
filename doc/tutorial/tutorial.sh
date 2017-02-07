# this is an example for a complete application of DFS, and as in "tests" 
# uses the p53 DNA Binding domain. This includes running the dfs script 
# to get the rescuability matrices and then the compensatory_power script 
# to calculate the rescaubility scores.


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

pymol_binary=$(which pymol)
pymol_runfile="view_compensatory_power.pml"

# structure to be used to generate the Anisotropic Network Model (ANM)
pdb="1TSR.pdb"

# selection of atoms to be used for the ANM model. We're going to use
# CA atoms from chain A from p53's multi-chain PDB 1TSR

selection="chain A and protein and name CA"

# distance cut-off for ANM
c=15.0

# force constant for ANM harmonic potential (gamma)
g=0.1

# Forces configuration file. Read the file itself for information
# on its format and syntax. 
forces_config=forces.ini

# number of cores to be used 
np=2

# command line
echo now running: $dfs -p $pdb -f $forces_config -c $c -g $g -s $selection --np $np -v
$dfs -p $pdb -f $forces_config -c $c -g $g -s $selection --np $np

# Two matrix files have been saved: one that contains outcome of the calculation
# in constant force mode and the other in constant displacement mode.


#################################################################################
# 2 - calculate rescuability power index

# We feed both matrices calculated in the previous step to the compensatory_power
# script in order to calculate the CP index on the average matrix:

echo "now running: $compp -m fixed_force_scores.dat fixed_displacements_scores.dat -p $pdb -s $selection" -P compensatory_power.pdb
$compp -m fixed_force_scores.dat fixed_displacements_scores.dat -p $pdb -s $selection -P compensatory_power.pdb

# The output file normalized_score.dat contains the per-residue compensatory 
# power score. The same information is also stored in the B-factor field of 
# the output pdb file "compensatory_power.pdb", only for those residues for
# which the DFS score was calculated (i.e. those of chain A). 
# We can visualize the compensatory power scores on the structure using the
# compensatory_power.pdb file and PyMOL, simply by running the following 
# command in the current directory (PyMOL must be installed and available
# on command line):

#################################################################################
# 3 - visualize rescuability power using PyMOL

echo "now running: $pymol_binary -d run $pymol_runfile"
$pymol_binary -d run $pymol_runfile

# the color scale goes from blue to red, with blue representing lower values
# and red representing higher ones. 
