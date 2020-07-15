# GenomicModelCreator

This repo has scripts designed to train an XGBoost (and other) models based on input genomic data.  By default it accepts both one-hot-encoded alignments and fasta files as genomic features for training.  The model builder was originally designed to do antimicrobial resistance, but extended train on various metadata as well.  

## Requirements

There are two major prerequisites to run these scripts:
1. KMC must be installed and within the paths.
2. Python and the required libraries should also be installed.

### KMC Setup

The KMC Directory contains a kmc.sh script that is designed as a shortcut to run both kmc and kmc_dump in tandem.  The KMC executible package (kmc, kmc_dump, and kmc_tools) will still need to be downloaded from [here](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about).  Once that is downloaded and/or compiled, it's easiest to just copy the executibles into the KMC folder and run:

```bash
PATH="$PATH:<path_to_this_directory>/KMC"
```

### Python Packages

This repo leans on several python packages to work properly and must be installed.  It's recommended to use Anaconda to do some of this work as it creates a "container" for your Python setup.  Additionally these can also be installed on top of your existing Python if you so choose to.  The following packages need to be installed:
- numpy ([website](https://numpy.org), [anaconda](https://anaconda.org/anaconda/numpy))
- sklearn ([website](https://scikit-learn.org/stable/), [anaconda](https://anaconda.org/anaconda/scikit-learn))
- scipy ([website](https://www.scipy.org), [anaconda](https://anaconda.org/anaconda/scipy))
- xgboost ([website](https://xgboost.readthedocs.io/en/latest/), [anaconda](https://anaconda.org/conda-forge/xgboost))

## Typical Input Data Formatting

All input and options are sent to the script through options.  Running the script with the *-h* option will show all available options and descriptions.  

Regardless of whether or not alignments/binary or fasta (k-mer) inputs are to be used with the script, the source of the predictive labels is to be sent to the script as a tabular *.tab* file that is formatted as follows:

```
genome_id	<test_cond1>	<test_cond2>	...	<test_cond3>	label
```

The genome ID is the genome ID that corresponds to the alginment or fasta file.  The test conditions could be an antibiotic and testing standards for example.  The label would be what is to be predicted and trained on.  Anything in angled brackets (<>) is optional

There are two input formats for genome features, fasta and alignment.  Nucleotide fasta files should be placed into a directory with their files named genome_id.fasta and the directory name passed to the script.  The type of fasta file can be anything (full contig, genes, etc.) but should be consistant accross all files.  The genome_ids must match those found in the supplied tabular file.  At runtime the script will read these files and run KMC.

Alignments are sent as a tab delimited file name whose file is formatted as follows:

```
genome_id	alignment
```

The genome_id needs to match those found in the supplied tabular file.  The alignment is a string of 0's and 1's representing the one-hot encoding of the genome's alignment to one another.  The 0's and 1's should not be separated by any characters.  Note that in theory it is possible to just train a model one one-hot encoded data using this method.  

