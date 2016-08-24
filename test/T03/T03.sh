######################################################################
# protein: third immunoglobulin binding domain of Protein G (GB3) 
# PDB ID: 2LUM

######################################################################
# path and configuration settings
dotM=`which dotM`

# input h5 file
filename="../T01/fixed_force_details.h5"

# structure used in the dfs step
pdb="../T00/original.pdb"

# input file settings
input_settings="-f $filename -p $pdb"

######################################################################
# residue settings

# reference pathogenic mutation site
pathogenic_site="A5"

# secondary site
secondary_site="A10"

# residue settings
residue_settings="-R $pathogenic_site -D $secondary_site"

######################################################################
# output settings

# output prefix
outprefix="T03"

# output filename
outputRMSIP="$outprefix.fixed_force_rmsip.txt"

# output_settings
output_settings="-o $outputRMSIP"

######################################################################
# check availability of input file
if [ ! -e $filename ]
then
        echo "TEST T03: NOT EXECUTED - REQUIRES T01"
        exit
fi

######################################################################
# command line
$dotM $input_settings $residue_settings $output_settings >& log.T03.log

######################################################################
# test check
cp $outputRMSIP ../check/T03.out
diff ../check/T03.out ../check/T03.chk > ../check/T03.diff
if [ -s ../check/T03.diff ]
then
        echo "TEST T03: FAILED"
else
        echo "TEST T03: PASSED"
fi
