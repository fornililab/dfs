<p align="center">
  <img src="https://afornililab.files.wordpress.com/2017/02/dfs_logo_600ppi.png" width="256" />
</p>

## Double Force Scanning

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
program. Please see the LICENSE file for details.

DFS is available on [GitHub](https://github.com/fornililab/dfs/). The distribution includes a Python module (dfsutils) that comprises the main classes and functions, as well as the following scripts:

 * **dfs**: runs the Double Force Scanning calculations
 * **compensatory_power**: calculates per-residue Compensatory Power
 * **data_muncher**: extracts data from the compressed output HDF5 file
 * **dotM**: calculates the RMSIP between force and displacement vectors

A summary help can be obtained with the --help option, e.g.:

    dfs --help

The latest release of the dfs package is available [here](https://github.com/fornililab/dfs/releases/latest)

## Installation

Please refer to the INSTALL file in the distribution.

A test directory is also provided to check your installation.

## Documentation

Please refer to the doc directory in the distribution.

A step-by-step tutorial is included, together with a user manual.
 
## Citation

If you publish results produced with the DFS or develop methods based on the DFS
code, please cite the following paper:

    M. Tiberti, A. Pandini, F. Fraternali, A. Fornili, "In silico identification of rescue sites by double force scanning", Bioinformatics, accepted
