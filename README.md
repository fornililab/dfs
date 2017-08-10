# DFS

**D**ouble **F**orce **S**canning (**DFS**) is a computational method to
predict rescue sites in proteins.

DFS uses external forces applied to (first-site,second-site) pairs of residues
to mimic the effect of a pathogenic mutation (first site) and a candidate 
rescue mutation (second site). The effect of forces on the protein structure
is calculated using the Linear Response Theory combined with the Anisotropic
Network Model. The compensatory effect is quantified through a rescuability 
index, which is > 0 when the structural perturbation induced by the forces at
the two sites is smaller than the one produced by a single force at the first
site. The rescuability indices of a given second site are combined to generate
an overall estimate of its compensatory power, which can be used to predict if
the residue is a rescue site.

## Availability

The program is made available under the GNU Public License for academic
scientific purposes and under the condition that proper acknowledgement is made
to the authors of the program in publications resulting from the use of the
program. Please see [LICENSE](LICENSE) for details.

DFS is available on [GitHub](tobedefined "DFS").  The distribution includes a 
Python module (dfsutils) that comprises the main classes and functions, as well
as the following scripts:

 * **dfs**: runs the Double Force Scanning calculations
 * **compensatory_power**: calculates per-residue Compensatory Power
 * **data_muncher**: extracts data from the compressed HDF5 file
 * **dotM**: calculates the RMSIP between force and displacement vectors

A summary help can be obtained with the --help option, e.g.:

`dfs --help`

## Installation

Please refer to the [INSTALL](INSTALL) file in the distribution.

A [test](test) directory is provided to check your
installation.

## Documentation

Please refer to the [doc](doc) directory in the distribution.
A step-by-step [tutorial](doc/tutorial) is also included.
 
## Citation

If you publish results produced with the DFS or develop methods based on the DFS
code, please cite the following paper:

	M. Tiberti, A. Pandini, F. Fraternali, A. Fornili,
        "In silico identification of rescue sites by double force scanning",
        Bioinformatics, accepted

