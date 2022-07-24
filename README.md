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
conda env create -n PorinPredict --file PorinPredict.yaml

conda activate PorinPredict
/path/to/PorinPredict.py --version
```
For a faster installation, use [mamba](https://github.com/mamba-org/mamba) instead of conda.

## Run PorinPredict

```
conda activate PorinPredict
/path/to/PorinPredict.py -i /path/to/genome.fasta -o /path/to/output_directory/ [other options]
```

To analyze multiple assemblies and create a summary table, you can use a for loop in combination with the "--summarize" flag:
```
for i in /path/to/input_directory/*.fasta; do /path/to/PorinPredict.py -i $i -o /path/to/output_directory/ --summarize ; done
```

Assemblies for a test run are included in PorinPredict/test_assemblies:
```
for i in ./test_assemblies/*.fasta; do ./PorinPredict.py -i $i -o test_run --summarize ; done
```

## Output Files

Filename | Description
---------|------------
`PREFIX_PA_PorinPredict.tsv` | PorinPredict results and classification
`PorinPredict_results_table.tsv` | PorinPredict collated results (when run with the "--summarize" flag)

## Considerations

* PorinPredict relies on high-quality genome assemblies. Before running PorinPredict, we recommend to confirm *Pseudomonas aeruginosa* species affiliation using e.g. [rMLST](https://pubmlst.org/species-id) and to assess the assembly quality using [CheckM](https://github.com/Ecogenomics/CheckM/wiki), [QUAST](http://quast.sourceforge.net/), or similar tools. In low-quality assemblies, an absent *oprD* gene may be due to technical reasons such as an insufficient coverage In our study, the following criteria were used to define low-quality assemblies: N50 < 15 kb, assembly length <5.9 Mb, CheckM completeness <97 %, or CheckM contamination >3 %.

* In addition to inactivating mutations in OprD and promoter disruptions, PorinPredict reports missense mutations. Specific amino acid substitutions associated with carbapenem resistance are described in our publication: LINK

* 15 intact OprD variants are currently included in the database. With the increasing number of available genomes, the database will be supplemented with additional variants of carbapenem-susceptible isolates.

## Citation

LINK

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
