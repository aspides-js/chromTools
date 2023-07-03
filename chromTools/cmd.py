#!/usr/bin/python


# ------------------------------------
# modules
# ------------------------------------

import gzip as gz
import multiprocessing as mp
import pathlib
import shutil
import sys
import tempfile
import time
from subprocess import check_output

import lmfit
import matplotlib.pyplot as plt
import mmh3
import numpy as np
import pandas as pd

from chromTools.chmm_cmd import make_binary_data_from_bed
from chromTools.validate import assert_compressed, chmm_validator, macs_validator


# ------------------------------------
# Main function
# ------------------------------------
def run(options):
    """The Main function pipeline for chromTools

    Args:
            options (Namespace object): Command line options
    """
    ## Concatenating
    cat_bed(options.files, options.control, options.subdir, options.info)
    total, nfile = wc(
        options.increment, options.subdir, options.info, options.warn, options.paired
    )

    ## Downsampling
    start_time = time.time()
    options.start_time = start_time
    options.info("Downsampling...")
    options.info(f"CPU number: {str(mp.cpu_count())}")

    pool = mp.Pool()
    args = [
        (n, options, total) for n in range(1, nfile)
    ]  # nfile should be number calculated by wc()

    r = {}  # initiate empty dictionary
    for res in pool.starmap(downsample, args):
        r.setdefault(res[0], [])
        r[res[0]].append(res[1])
        options.info(f"--- {(time.time() - start_time)} seconds ---")

    ## Binarising
    options.info("CHMM binarising...")
    args = [
        (n, options) for n in range(0, nfile)
    ]  # nfile should be number calculated by wc()
    r["0"] = [total]
    for res in pool.starmap(use_chmm, args):
        print("activated")
        r.setdefault(res[0], [])
        r[res[0]].append(res[1])
    options.info(f"--- {(time.time() - start_time)} seconds ---")

    print(r)
    param_write(r, options.outdir)
    param_plot(r, options.outdir)

    print("Complete")


# --------------------------------------------------------------------------------#


def cat_bed(files, control, subdir, info):
    """
    Concatenate compressed or uncompressed files into downsampled.0.bed.

    :param files: Files to concatenate.
    :type files: str
    :param control: Control files to concatenate.
    :type control: str or None
    :param subdir: Path to the subsample directory.
    :type subdir: str
    :param info: Logging function for informational messages.

    """
    info("Concatenating files...")
    start_time = time.time()
    with open(pathlib.Path(subdir / "downsampled.0.bed"), "wb") as wfd:
        for f in files:
            if assert_compressed(f):
                print("File is compressed")
                with gz.open(f, "rb") as fd:
                    shutil.copyfileobj(fd, wfd)
            else:
                with open(f, "rb") as fd:
                    shutil.copyfileobj(fd, wfd)
    if control != False:
        with open(pathlib.Path(subdir / "downsampled.ctrl.bed"), "wb") as wfd:
            for f in control:
                if assert_compressed(f):
                    print("File is compressed")
                    with gz.open(f, "rb") as fd:
                        shutil.copyfileobj(fd, wfd)
                else:
                    with open(f, "rb") as fd:
                        shutil.copyfileobj(fd, wfd)

    info(f"--- {(time.time() - start_time)} seconds ---")


def wc(increment, subdir, info, warn, paired):
    """Count total read number

    :param increment: The amount to increase each subsampled file by.
    :type increment: int
    :param subdir: Path to the subsample directory.
    :type subdir: str
    :param info: Logging function for informational messages.
    :param warn: Logging function for warning messages.
    :param paired: Indicates if the reads are paired end.
    :type paired: bool

    :return: A tuple containing the total number of reads and the number of files that will be generated.
    :rtype: tuple[int, int]

    """
    start_time = time.time()
    info("Calculating total read number...")
    total = int(check_output(["wc", "-l", f"{subdir}/downsampled.0.bed"]).split()[0])

    if paired:
        total = total / 2

    info(f"--- {(time.time() - start_time)} seconds ---")

    nfile = int(total / increment)

    # give warning if nfile is very high
    if nfile > 100:
        warn(f"Number of downsampled files will be {nfile}")

    if total == 0:
        warn(
            "Total number of lines is equal to 0. Are your input files empty? Terminating."
        )
        sys.exit(1)

    return total, nfile


# --------------------------------------------------------------------------------#


def params(proportion):
    """Calculate a value between a range to reflect proportion of reads kept

    Args:
            proportion (int): Proportion of reads to subsample from whole dataset

    Returns:
            maxHashValue (int): Threshold above which reads are discarded
    """
    max_size = sys.maxsize
    min_size = -sys.maxsize - 1
    maxRange = max_size - min_size
    maxHashValue = min_size + round(maxRange * proportion)
    return maxHashValue


def discard(maxHashValue, seed, line):
    """Generate a random hash from readname using seed. If number above proportional cut-off (True), discard read.

    Args:
            maxHashValue (int): Threshold above which reads are discarded
            seed (int): Random seed
            line (str): Read/line in file

    Returns:
            bool: Boolean specifying if readname is below or above discard threshold
    """
    readname = line.split("\t")[3].rsplit("/")[
        0
    ]  # extract readname, remove everything after '/' (read pair if paired)
    hashInt = mmh3.hash64(readname, seed)
    return hashInt[0] > maxHashValue


def downsample(n, options, total):
    """Downsample a bed file. For each read pair, a random value is assigned between the range.
    The proportion is used to calculate a maximum acceptable value within the range. Records
    whose value is below the limit are written to outfile, records whose hash value is above
    the limit are discarded.

    :param n: Numerical descriptor of file
    :type n: int
    :param total: Total number of reads/read pairs in downsampled.0.bed
    :type total: int
    :param options: Command line arguments
    :type options: Namespace object

    :return: A tuple containing the numerical descriptor of file and the number of reads.
    :rtype: tuple[str, float]

    """
    proportion = (options.increment * n) / total
    outfile = pathlib.Path(options.subdir / f"downsampled.{n}.bed")
    reads = 0
    a = params(proportion)
    with open(outfile, "w") as outf:
        with open(pathlib.Path(options.subdir / "downsampled.0.bed"), "r") as f:
            for line in f:
                if discard(a, options.seed, line):
                    continue
                else:
                    outf.write(line)
                    reads += 1
    if options.paired:
        reads = reads / 2
    return str(n), reads


# --------------------------------------------------------------------------------#


def use_macs(n, options):
    """Uses macs3 peak calling functionality to identify areas in the genome that
    have been enriched with aligned reads.

    Args:
            n (int): Numerical descriptor of file
            options (Namespace object): Command line arguments

    Returns:
            n (int): Numerical descriptor of file
            (int): Sum of counted peaks / total chromosome length
    """
    options = macs_validator(n, options)

    # if set by user, use alternative tempdir
    if options.tempdir:
        tempfile.tempdir = options.tempdir

    if options.paired:
        (treat, control) = load_frag_files_options(options)
    else:
        (treat, control) = load_tag_files_options(options)

    t0 = treat.total
    t1 = t0
    options.d = options.tsize

    peakdetect = PeakDetect(treat=treat, control=control, opt=options)

    peakdetect.call_peaks()

    # filter out low fe peaks
    # peakdetect.peaks.filter_fc( fc_low = options.fecutoff )

    ofhd_bed = open(options.peakBed, "w")
    peakdetect.peaks.write_to_narrowPeak(
        ofhd_bed,
        name_prefix=b"%s_peak_",
        name=options.name.encode(),
        score_column="qscore",
        trackline=False,
    )
    ofhd_bed.close()

    return str(n), count_peakmark(options)


def count_peakmark(options):
    """Count the number of base pairs overlapping peaks per chromosome

    Args:
            options (Namespace object): Command line arguments

    Returns:
            (int): Sum of counted peaks / total chromosome length
    """
    g = chr_len(options.genome)
    with open(options.peakBed, "r") as f:
        for line in f:
            chrom = line.split("\t")[0]
            num = int(line.split("\t")[2]) - int(line.split("\t")[1])
            g[chrom].append(num)

    if 0 in [sum(x[1:]) for x in g.values()]:
        options.warn(
            "0 values in chromosome(s) peak count may indicate incorrect chromosome length file."
        )

    return (sum([sum(x[1:]) for x in g.values()])) / (
        sum([x[0] for x in g.values()])
    )  # sum of counted peaks / total chromosome length


def chr_len(genome):
    """Generate dictionary from genome chromosome length files

    Args:
            genome (str): Path to genome chromosome length file

    Returns:
            g (dict): {CHR number : length}
    """
    g = {}
    with open(genome) as f:
        for line in f:
            (key, val) = line.split()
            g.setdefault(key, [])
            g[key].append(int(val))
    return g


# --------------------------------------------------------------------------------#


def use_chmm(n, options):
    options = chmm_validator(options)
    count, total = make_binary_data_from_bed(n, options)
    options.info("this is happening")
    return str(n), count / total


# --------------------------------------------------------------------------------#


def param_write(r, outdir):
    """Write the output to text file (tsv)

    Args:
            r (dict): {Subsampled file : [<Number of reads>, <Proportion of genome bound>] }
            outdir (str): Output directory
    """
    with open(pathlib.Path(outdir / "completeness.txt"), "w") as f:
        for key, value in r.items():
            f.write(f"{key}\t{value[0]}\t{value[1]}\n")


def param_plot(r, outdir):
    """Plot the output to graph and call the Micheal-Menten function

    Args:
            r (dict): {Subsampled file : [<Number of reads>, <Proportion of genome bound>] }
            outdir (str): Output directory
    """
    df = pd.DataFrame(r)
    plt.figure(figsize=(10, 6), tight_layout=True)
    plt.plot(df.loc[0].tolist(), df.loc[1].tolist(), "s-", color="#06846a")
    plt.xlabel("Number of Reads")
    plt.ylabel("Proportion of marks")
    plt.savefig(pathlib.Path(outdir / "completeplot.jpg"))

    mm(df, outdir)


# -------------------------------------------------------------------------------#


def v(s, Vm, Km):
    return (Vm * s) / (Km + s)


def residuals(p, x, y):
    Vm = p["Vm"]
    Km = p["Km"]
    fi = v(x, Vm, Km)
    return y - fi


def mm(df, outdir):
    """Calculate the Michealis-Menten kinetics and plot to graph

    Args:
            df (pandas dataframe): [<Subsampled file>] [<Number of reads>] [<Proportion of genome bound>]
            outdir (str): Output directory
    """
    data = np.array(df)

    params = lmfit.Parameters()
    params.add("Vm", value=1, min=0, max=5_000_000_000)
    params.add("Km", value=1, min=0, max=5_000_000_000)

    result = lmfit.minimize(residuals, params, args=(data[0], data[1]))

    fm = np.linspace(0, max(data[0]), 100)
    plt.figure(figsize=(10, 6), tight_layout=True)
    plt.scatter(df.loc[0].tolist(), df.loc[1].tolist(), color="k")
    plt.plot(fm, v(fm, result.params["Vm"].value, result.params["Km"].value), "k")
    plt.xlabel("[S] (reads)")
    plt.ylabel("v (proportion)")
    plt.axhline(y=result.params["Vm"].value, linestyle="-", color="#06846a")
    if result.params["Km"].value < max(data[0]):
        plt.axvline(x=result.params["Km"].value, linestyle="-", color="#06846a")

    plt.title(label=f'Vm: {result.params["Vm"].value}')
    plt.savefig(pathlib.Path(outdir / "mmplot.jpg"))

    with open(pathlib.Path(outdir / "mm.txt"), "w") as f:
        f.write(f'{result.params["Vm"].value}\t{result.params["Km"].value}')


# -------------------------------------------------------------------------------#
