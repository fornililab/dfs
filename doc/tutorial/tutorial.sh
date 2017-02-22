# This is an example of a complete application of DFS on p53 (DNA binding 
# domain). It includes three steps: 
#
# 1. dfs run to generate the rescuability score matrices
# 2. compensatory_power run to calculate the compensatory power 
# 3. visualisation of the compensatory power mapped onto the 3D 
#    structure with pymol

#################################################################################
# 1 - use the Double Force Scanning method to generate rescuability score matrices

# First we run the dfs script that will generate the rescuability
# score matrices. We will not save any  HDF5 details file
# for the moment as it is not strictly required for the subsequent
# analyses. This run should take around 30 minutes on 4 Xeon cores @3.70GHz.

# dfs scripts directory - to be changed according to the location of the 
# dfs scripts in your filesystem
scripts_dir=../../scripts_src

dfs=$scripts_dir"/dfs"
compp=$scripts_dir"/compensatory_power"

pymol_binary=$(which pymol)
pymol_runfile="view_compensatory_power.pml"

# input structure
pdb="1TSR.pdb"

# selection of atoms to be used for the ANM model 
# (CA atoms from chain A in this case). 
selection="chain A and protein and name CA"

# distance cut-off for ANM (A)
c=15.0

# force constant for ANM harmonic potential (gamma, kcal/mol/A^2)
g=0.1

# Forces configuration file. The forces.ini file included in this 
# tutorial also contains a description of its format. 
forces_config=forces.ini

# number of cores to be used 
np=4

# command line
echo now running: $dfs -p $pdb -f $forces_config -c $c -g $g -s $selection --np $np -v
$dfs -p $pdb -f $forces_config -c $c -g $g -s $selection --np $np

# Two matrix files are saved, containing the rescuability score matrix 
# calculated in Fixed Force (fixed_force_scores.dat) and Fixed RMSD 
# (fixed_displacements_scores.dat) mode.

#################################################################################
# 2 - calculate rescuability power index
# We feed both matrices calculated in the previous step to the 
# compensatory_power script in order to calculate the compensatory 
# power for each residue:

echo "now running: $compp -m fixed_force_scores.dat fixed_displacements_scores.dat -p $pdb -s $selection -P compensatory_power.pdb"
$compp -m fixed_force_scores.dat fixed_displacements_scores.dat -p $pdb -s $selection -P compensatory_power.pdb

# The output file normalized_score.dat contains the per-residue compensatory 
# power values. The same information is also stored in the B-factor field of 
# the output pdb file "compensatory_power.pdb", only for those residues for
# which the DFS score was calculated (chain A residues in our example).

#################################################################################
# 3 - visualize rescuability power using PyMOL

# We can visualize the compensatory power on the structure using the
# compensatory_power.pdb file generated in the previous setp and PyMOL. 
# We just run the following command in the current directory (PyMOL must 
#be installed and available on command line):

echo "now running: $pymol_binary -d run $pymol_runfile"
$pymol_binary -d run $pymol_runfile
