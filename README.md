# PorinPredict

*In silico* prediction of OprD porin-loss in *Pseudomonas aeruginosa* isolates from assembled genomes

## Overview

PorinPredict is a command-line tool for the standardized detection of OprD-inactivating mutations in *Pseudomonas aeruginosa* genomes.
OprD loss typically leads to carbapenem resistance.

The program and associated databases can be downloaded from https://github.com/MBiggel/PorinPredict.

For the complementary detection of carbapenemases consider using
[AMRFinderPlus](https://github.com/ncbi/amr#ncbi-antimicrobial-resistance-gene-finder-amrfinderplus),
[RGI](https://card.mcmaster.ca/analyze/rgi),
[Resfinder](https://cge.cbs.dtu.dk/services/ResFinder/)
or similar tools. 


## Installation

### Requirements

* Linux or macOS
* Python3 and pandas 
* DIAMOND >=v2.0.0
* BLAST+
* R and packages dplyr and Biostrings

PorinPredict expects DIAMOND, BLASTn, and Rscript to be available in $PATH.

### Installation of dependencies via conda

To install required dependencies in a separate [Conda](https://conda.io/docs/install/quick.html) environment, run the following commands:
```
git clone https://github.com/MBiggel/PorinPredict
cd PorinPredict
conda env create -n porinpredict --file porinpredict.yaml

conda activate porinpredict
/path/to/porinpredict.py --version
```


## Run PorinPredict

```
conda activate porinpredict
/path/to/porinpredict.py -i /path/to/genome.fasta -o /path/to/output_directory/ --threads 8
```

To analyze multiple assemblies and create a summary table, you can use a for loop in combination with the "--summarize" flag:
```
for i in /path/to/input_directory/*.fasta; do /path/to/porinpredict.py -i $i -o /path/to/output_directory/ --summarize ; done
```

Assemblies for a test run are included in the repository folder PorinPredict/test_assemblies:
```
for i in ./test_assemblies/*.fasta; do ./porinpredict.py -i $i -o test_run --summarize -t 8; done
```

## Output Files

Filename | Description
---------|------------
`PREFIX_PA_PorinPredict.tsv` | PorinPredict results and classification
`PorinPredict_results_table.tsv` | PorinPredict collated results (when run with the "--summarize" flag)

## Considerations

* PorinPredict relies on high-quality genome assemblies. Before running PorinPredict, we recommend to confirm *Pseudomonas aeruginosa* species affiliation using e.g. [rMLST](https://pubmlst.org/species-id) and to assess the assembly quality using [CheckM](https://github.com/Ecogenomics/CheckM/wiki), [QUAST](http://quast.sourceforge.net/), or similar tools. In low-quality assemblies, absence of *oprD* may be caused by technical reasons such as an insufficient coverage.

* In addition to inactivating mutations in OprD and promoter disruptions, PorinPredict reports missense mutations. Specific amino acid substitutions associated with carbapenem resistance are described in our [publication](https://doi.org/10.1128/spectrum.03588-22)

* 15 intact OprD variants are currently included in the database. With the increasing number of available genomes, the database will be supplemented with additional variants of carbapenem-susceptible isolates.

## Citation

Biggel, M., Johler, S., Roloff, T., Tschudin-Sutter, S., Bassetti, S., Siegemund, M., Egli, A., Stephan, R., & Seth-Smith, H. M. B. (2023). PorinPredict: In Silico Identification of OprD Loss from WGS Data for Improved Genotype-Phenotype Predictions of *P. aeruginosa* Carbapenem Resistance. Microbiology spectrum, e0358822. https://doi.org/10.1128/spectrum.03588-22

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
