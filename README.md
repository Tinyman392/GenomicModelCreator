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

## Model Building Script

The *buildModel.py* script is used to build a model.  These options can be seen by running *buildModel.py -h*.  The options are described below:
- -f --fasta_dir : Specify a directory containing fasta files to train with.  There is no default for this option.
- -t --tabular_file : Specity the file containing the genomic metadata and testing conditions to use.  There is no default for this option.
- -H --header : Specify whether the tabular file (-t) contains a header.  This is specified as either "True" or "False".  The default for this is "False".
- -T --temp_dir : Specify the temporary directory to use while training.  This director will be filled during training and can be used for debugging in the event of a crash.  It doesn't completely empty after the script finishes.  The default value for this is "temp".
- -o --out_dir : Specify the output directory to store the model after training.  If stats are computed (-S), they will be stored in here as well.  The default value for this is "model".
- -n --threads : Specify the number of threads to train with.  Note that if you're running a SciKit Learn model (like Random Forest), this may use a large amount of RAM (total RAM use < single-thread use * number of threads).  Some SciKit Learn models don't support multithreading.  For XGBoost you can specify as many threads as your machine has.  The default value for this is "1".
- -d --depth : Specify the maximum tree depth for XGBoost.  The default value for this is "4".
- -k --kmer_size : Specify the kmer size to use if using fasta files for input.  Larger values for k will result in longer train times and RAM use.  Default value for this is "10".  Values as low as 8 should have low bearing on overall accuracy, values as high as 15 can be used if feature importance is needed.
- -K --kmc_dir : Specify the location of previously run KMC output.  If this is specified, then the KMC will not be rerun for anything.  Note if none is supplied, the default location for KMC output will be "/dev/shm/kmc\<pid\>" where the pid is the pid of the process.
- -P --presence_absence : Specify whether or not to use presence vs absence of a kmer or kmer counts.  Defaults to "False".  Setting this option to true for k > 12 may be benifical for RAM conservation.
- -i --individual : Specify whether or not to run individual models for \<test_cond1\> which is typically the antibiotic.  The default value for this is "False".
- -e --enumerate_classes : Specify whether the classes in the tabular (-t) file are already enumerated or if they require enumeration.  For SIR/SR models, do not enumerate classes.  Do not enumerate classes for regressions.  Default value is "False".
- -a --folds_to_run : Specify the number of folds to run in the CV.  This must be <= the total folds to run (-A).  Defualts to "5".
- -A --total_folds : Specify the total folds to run.  The dataset will be split up into this many parts and one chosen for each fold ran will be the test set, one as the validation set, and the rest used for training.  
- -c --classification : Specify whether or not this model should be a classification model.  Defaults to "False".
- -j --SvNS : Specify whether or not to run an S vs NS model for SIR.  Defaults to "False".
- -J --noI : Specify whether or not to run an S vs R model for SIR.  Defaults to "False".
- -m --model_params : Specify any additional model parameters to run with.  An eta or silent marker (or other XGB option) can be specified here.  Must be passed in as a Python hash.  Defaults to "{'eta':0.0625, 'silent':1}".  
- -N --num_rounds : Specify the number of rounds to boost/number of trees to make.  Defaults to "1000".
- -C --cleanup : Specify whether or not to cleanup the temp directory.  Useful to debugging.  Defaults to "True".
- -S --compute_stats : Specify what stats to compute on the model's output.  We currently support rawAcc, w1Acc, r2, VMEME, confMatrix, and classReport.  Combinations for AMRcls, AMRreg, and cls are supported as well.  Future support for heatmap will be added at a later date.  Separate multiple with commas.  There is no default for this option.





