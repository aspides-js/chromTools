#!/usr/bin/python

# ------------------------------------
# Modules
# ------------------------------------

import errno
import logging
import os
import random
import sys
from pathlib import Path

from chromTools.constants import GENOME as gnm

# ------------------------------------
# Misc function
# ------------------------------------


def args_validator(options):
    """
    Validate command line arguments and set additional necessary parameters.

    :param options: Command line arguments.
    :type options: Namespace object.
    :raises FileNotFoundError: If paths to files or directories are inaccessible.
    :return: Validated command line arguments.
    :rtype: Namespace object.
    """
    ## logging
    logging.basicConfig(
        level=20,
        format="%(levelname)-5s @ %(asctime)s: %(message)s ",
        datefmt="%a, %d %b %Y %H:%M:%S",
        stream=sys.stderr,
        filemode="w",
    )

    options.error = logging.critical
    options.warn = logging.warning
    options.debug = logging.debug
    options.info = logging.info

    ## outdir
    if options.outdir == "":
        options.outdir = Path.cwd()

    ## genome
    if not options.gsize:
        try:
            options.gsize = gnm[options.genome][1]
            options.genome = (
                Path(__file__).parent / "chromsize" / gnm[options.genome][0]
            )
            options.info(f"Path to genome chromosome sizes is: {options.genome}")
        except:
            options.warn(
                f'If not using available genome shortcuts, effective genome size should also be specified with the --gsize flag. Available shortcuts for genome chromosome sizes are: {", ".join(list(gnm.keys()))}'
            )
            sys.exit(1)

    else:
        options.genome = Path(__file__).parent / "chromsize" / options.genome
        options.gsize = int(options.gsize)
        options.info(f"Path to genome chromosome sizes is: {options.genome}")
        options.info(f"Effective genome size is: {options.gsize}")

    ## files (path)
    for f in options.files:
        if not Path(f).exists():
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f)
        # incorrect data format


    ## ctrl files (path)
    if options.control != False:
        for f in options.control:
            if not Path(f).exists():
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), f)

    ## outdir
    options.outdir = Path(options.outdir)
    options.subdir = Path(options.outdir / "1_subsample/")
    options.bindir = Path(options.outdir / "2_binarised/")

    if options.force:
        try:
            options.subdir.mkdir(parents=True, exist_ok=True)
            options.bindir.mkdir(parents=True, exist_ok=True)
        except:
            options.warn(
                f"Output directories in {options.outdir} could not be accessed. Terminating program."
            )
            sys.exit(1)
    else:
        try:
            options.subdir.mkdir(parents=True, exist_ok=False)
            options.bindir.mkdir(parents=True, exist_ok=False)
        except:
            options.warn(
                f"Output directories in {options.outdir} could not be accessed or already exist. Please use --force-overwrite if you wish to overwrite output files. Terminating program."
            )
            sys.exit(1)

    ## seed
    if options.seed == None:
        options.seed = random.randint(0, 5_000_000_000)
        options.info(f"RANDOM SEED: {options.seed}")
    else:
        options.info(f"RANDOM SEED: {options.seed}")

    ## paired
    options.info(f"PE_MODE: {options.paired}")

    return options


def assert_compressed(f):
    """
    Test whether the file is gz compressed (first two bytes are 1f 8b).

    :param f: Path to the file.
    :type f: str
    :return: Boolean specifying whether the file can be opened or not.
    :rtype: bool
    """
    with open(f, "rb") as f:
        return f.read(2) == b'\x1f\x8b'


def benchmark(outdir, step, interval_time, timestr):
    bench_f = Path(outdir, f".bench_{timestr}.log")
    with open(bench_f, "a") as bench_f:
        bench_f.write(f"{step}\t{interval_time}\n")


# -----------------------------------------------------------------------------------#


def chmm_validator(options):
    """
    Validate and configure the options for the ChromHMM binarization.

    :param options: An object containing the pipeline arguments.
    :param szchromlengthfile: The file containing the chromosome lengths (genome file).
    :param szoutputbinarydir: The directory to store the output binary files.
    :param szmarkdir: The subdirectory for mark-related files.
    :param szcontroldir: The subdirectory for control-related files.
    :param bpairend: A flag indicating whether the input data is paired-end (obtained from options.paired).
    :param nshift: The shift value (default: 100).
    :param noffsetleft: The amount to subtract from the left coordinate to make it 0-based inclusive (default: 0).
    :param noffsetright: The value to add to the right coordinate for 0-based, non-inclusive bed files (default: 1).
    :param npseudocountcontrol: An integer pseudocount added uniformly to each bin in the control data (default: 1).
    :param dpoissonthresh: The Poisson tail probability threshold (default: 0.0001).
    :param dfoldthresh: The fold threshold (default: 0).
    :param bcontainsthresh: A flag indicating the behavior of the Poisson cutoff (default: True).
    :param dcountthresh: The absolute signal threshold (default: 0).
    :param nflankwidthcontrol: The flank width for control data (default: 5).
    :type options: Namespace
    :return: The validated and configured options.
    :rtype: Namespace
    """
    options.szchromlengthfile = options.genome
    options.szoutputbinarydir = options.bindir
    options.szmarkdir = options.subdir
    options.szcontroldir = options.subdir

    options.bpairend = options.paired
    options.nshift = 100  # chmm default
    options.noffsetleft = 0  # the amount that should be subtracted from the left coordinate so it is 0-based inclusive
    options.noffsetright = 1  # based on bed files being 0-based but not inclusive
    options.npseudocountcontrol = 1  # An integer pseudocount that is uniformly added to every bin in the control data in order to smooth the control data from 0. The default value is 1.

    options.dpoissonthresh = 0.0001  # chmm default
    options.dfoldthresh = 0  # chmm default
    options.bcontainsthresh = True  # i think chmm default - if true poisson cut off should be highest that still contains dpoissonthresh probability and if false requires strictly greater
    options.dcountthresh = 0  # chmm default
    options.nflankwidthcontrol = 5  # chmm default

    return options
