# entirety

entirety is a tool to assess the completeness of your epigenetic dataset. Can be used with ChIP-seq, ATAC-seq or WGBS data.

For questions on installation or usage, please open an issue or contact Jessica Shields (j.m.shields@exeter.ac.uk).

## 1. Install

You can install entirety via ```git clone``` and then by running ```pip install -e ./[dev]```. The program requires python>=3.7 and pip>=23.

## 2. Basic usage
entirety is a single command line tool which takes as input a set of BED files, and outputs a metric and plot indicating how complete the dataset is.

Basic usage is shown below. See [detailed usage](https://github.com/aspides-js/entirety/edit/main/README.md#3-detailed-usage) below for more info.

```
entirety \
  -f <reads.bed> \
  [-i <increment>] \
  -o <outdir> \
  -g <genome> \
  --gsize \
  -s <seed> \
  --paired \
  --force-overwrite
```

## 3. Detailed usage

#### entirety

Required parameters:

* `--files/-f <reads.bed>`: BED files of the dataset. Must be aligned files. Should include full path to file.

Optional parameters:

* `-o <outdir>`: Path to output directory where structure <outdir>/1_subsample <outdir>/2_binarise and output files will be created. Directory will be created if it does not exist. Default: current working directory.  
* `--increment/-i <int>`: Amount of reads/read pairs by which to incrementally subsample. e.g. If left to default the whole dataset will be subsampled to 50000000, 100000000, 150000000 etc. Default: 50000000. 
* `--genome/-g`: A two column tab delimited file with the first column being the chromosome and the second being the chromosome length. Genome assemblies hg18, hg19, hg38, mm9, mm10, rn5, rn6, danRer7, danRer10, dm3, dm6, ce6, and ce10 can be accessed by their genome assembly name (e.g. hg19). For other assemblies these can be
obtained with the fetchChromSizes script available from the UCSC browser http://hgdownload.cse.ucsc.edu/admin/exe/
specifying the desired assembly and redirecting the output to a text file. Default: hg38.
* `--gsize <int>`: Size of genome using. Required if specifying own genome chromosome length file. Default is FALSE
* `--seed <int>`: The seed used for initiating randomization operations. Default: randomly generated.
* `--paired`: If specified data will be treated as paired end. Default: false.
* `--force-overwrite`: If specified, files and directories in outdir will be overwritten. Default: false.

## 4. Usage example
In this example, we use the files test_s1.bed test_s2.bed as input.

`entirety -f test_s1.bed test_s2.bed --outdir tmp/ --increment 25 --genome hg37 --paired `

This command first concatenates the files into downsampled.0.bed. This file is then downsampled, binarisation run on each of these and proportion of marks found calculated. The command outputs a text file with the metrics, and plots the results.
