# DFS

The **D**ouble **F**orce **S**canning (**DFS**) is a computational method to
predict compensatory mutations in proteins. The effect of mutations are
represented by forces calculated using Linear Response Theory on Anisotropic
Network Models (ANM).

Given a Pathogenic Mutation (PM) all the residues in the protein are scanned to
identify secondary sites for candidate Compensatory Mutations (CM). For each
pair of residues (PM and secondary site) a rescuability score R is calculated.
Secondary sites with positive values of R indicate a putative compensatory
effect. Finally for each secondary site a comprehensive value of Compensatory
Power (CP) is calculated. Residues can be ranked according to their CP.

## Availability

The program is made available under the GNU Public License for academic
scientific purposes and under the condition that proper acknowledgement is made
to the authors of the program in publications resulting from the use of the
program. Please see [LICENSE](LICENSE) for details.

DFS is available on [GitHub](https://github.research.its.qmul.ac.uk/eex058/dfs
"DFS").  The distribution includes a Python module (dfsutils) that comprises the
main classes and functions, as well as the following scripts:

 * **dfs**: runs the Double Force Scanning calculations
 * **compensatory_power**: calculates per-residue Compensatory Power index
 * **data_muncher**: extracts data from the compressed HDF5 file
 * **dotM**: calculates the RMSIP between force and motion vectors

A summary help can be obtained with the --help option, e.g.:

`dfs --help`

## Installation

Please refer to the [INSTALL](INSTALL) file in the distribution.

A [test](test) directory is provided to check the functionallity of your
installation.

## Documentation

Please refer to the [doc](doc) directory in the distribution.
A step-by-step [tutorial](doc/tutorial) is also included.
 
## Citation

If you publish results produced with the DFS or develop methods based on the DFS
code, please cite the following paper:

<reference>
