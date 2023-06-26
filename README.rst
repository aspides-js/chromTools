chromTools
==========

chromTools complete is a tool to assess the completeness of your epigenetic dataset. Can be used with ChIP-seq, ATAC-seq or WGBS data, required files are aligned BED files.

For detailed documentation, please visit the `documentation pages <https://aspides-js.github.io/chromTools/>`_.

For questions on installation or usage, please open an issue or contact Jessica Shields (j.m.shields@exeter.ac.uk).

Install
-------

You can install chromTools via :code:`git clone` and then by running :code:`pip install -e ./[dev]`. The program requires python>=3.7 and pip>=23.

Basic usage
-----------

**chromTools complete** is a single command line tool which takes as input a set of BED files, and outputs a metric and plot indicating how complete the dataset is.

.. code-block:: console

   $ chromTools complete \
      --files/-f <reads.bed> \
      --control/-c <controlreads.bed> \
      --increment/-i <int> \
      --outdir/-o <outdir> \
      --genome/-g <genome>  \
      --gsize <int> \
      --seed/-s <int>  \
      --paired \
      --force-overwrite


Required parameters
~~~~~~~~~~~~~~~~~~~

:code:`--files/-f <reads.bed>`: BED files of the dataset. Must be aligned files. Should include full path to file.

Optional parameters
~~~~~~~~~~~~~~~~~~~

:code:`--control/-c <controlreads.bed>`: Optional flag containing BED-formatted control files. Must include full path. If unspecified, flag is set to FALSE.

:code:`--outdir/-o <outdir>`: Path to output directory where structure :code:`<outdir>/1_subsample/ <outdir>/2_binarise/` and output files will be created. Directory will be created if it does not exist. Default: current working directory.

:code:`--increment/-i <int>``: Amount of reads/read pairs by which to incrementally subsample. e.g. If left to default the whole dataset will be subsampled to 50000000, 100000000, 150000000 etc. Default: 50000000.

:code:`--genome/-g``: A two column tab delimited file with the first column being the chromosome and the second being the chromosome length. Genome assemblies hg18, hg19, hg38, mm9, mm10, rn5, rn6, danRer7, danRer10, dm3, dm6, ce6, and ce10 can be accessed by their genome assembly name (e.g. hg19). For other assemblies these can be obtained with the fetchChromSizes script available from the UCSC browser http://hgdownload.cse.ucsc.edu/admin/exe/ specifying the desired assembly and redirecting the output to a text file. Default: hg38.

:code:`--gsize <int>`: Size of genome using. Required if specifying own genome chromosome length file. Default: FALSE

:code:`--seed <int>`: The seed used for initiating randomization operations. Default: randomly generated.

:code:`--paired`: If specified data will be treated as paired end. Default: FALSE.

:code:`--force-overwrite`: If specified, files and directories in outdir will be overwritten. Default: FALSE.


Example
-------

In this example, we use the files test_s1.bed test_s2.bed as input.

:code:`chromTools complete -f test_s1.bed test_s2.bed --outdir tmp/ --increment 25 --genome hg37 --paired`

This command first concatenates the files into downsampled.0.bed. This file is then downsampled, binarisation run on each of these and proportion of marks found calculated. The command outputs a text file with the metrics, and plots the results.
