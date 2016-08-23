The DFS package implements the Dual Force Scanning method for the identification 
of compensatory mutations in proteins as detailed in:

<paper>

The method uses Linear Response Theory on Anisotropic Network Models (ANM) to 
calculate a rescuability score R between pairs of residues, one that is 
considered a Pathogenic Mutation (PM) to be compensated and another that is a 
Secondary  Mutation site, which represents a potential compensatory site. 
Positive values of R indicate a compensatory effect. 

The installation package includes a Python module (dfsutils) that comprises the 
main classes and functions, as well as a number of user-accessible scripts that
can be run by the user in order to run calculations using the DFS method as well
as analyzing the results. 

The following scripts are included in the package:

    * dfs - the main program, used to run calculations using the DFS method
    * compensatory_power - calculate per-residue compensatory power index
    * data_muncher - write data from details HDF5 file to more common formats
    * dotM - calculate RMSIP between vector spaces from DFS

Please check the files in the "doc" directory for relevant 
documentation and instructions for each script. A tutorial is also included.
Each script also hes a summary help section which is accessible by running:

	script --help 

Happy simulating!

