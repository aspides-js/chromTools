#!/usr/bin/python

# ------------------------------------
# Modules
# ------------------------------------

import random
import logging
import pytest
from pathlib import Path
from subprocess import check_output

from chromTools.complete_cmd import cat_bed, discard, wc

# ------------------------------------
# Test happy path
# ------------------------------------

def test_concatenate():
    """
    Test that the function cat_bed creates a new file which is the concat of input
    """
    subdir = Path("test_tmp")
    subdir.mkdir(parents=True, exist_ok=True)

    files = ["data/1_alignments/ENCFF111ZBO.bed", "data/1_alignments/ENCFF378LKV.bed", "data/1_alignments/ENCFF497RWO.bed"]
    total_ifs = [int(check_output(["wc", "-l", f"{x}"]).split()[0]) for x in files]

    # setup
    logging.basicConfig(level=20)
    warn = logging.warning
    info = logging.info

    cat_bed(files = files, control = False, subdir = subdir, info = info, warn = warn)

    outfile = Path(subdir, "subsampled.0.bed")
    total_of = int(check_output(["wc", "-l", f"{outfile}"]).split()[0])

    assert outfile.is_file() 
    assert total_of == sum(total_ifs)


def test_no_discard_paired():
    """
    Test that the function discard() keeps paired samples
    """
    ## set up
    maxHashValue = -9038904596117680288  # params(0.01)
    seed = 10

    R1 = "chr1\t10060\t10160\tSIM:chr1:10060:141:0:303:312/1\t31\t+\n"
    R2 = "chr1\t10102\t10202\tSIM:chr1:10060:141:0:303:312/2\t31\t-\n"

    a = discard(maxHashValue, seed, R1)
    b = discard(maxHashValue, seed, R2)

    assert a == b


def test_discard_diff_single():
    """
    This may fail sometimes as there is an element of probability but on average should pass.
    """
    ## set up
    maxHashValue = -9038904596117680288  # params(0.01)

    R1 = "chr1\t10060\t10160\tSIM:chr1:10060:141:0:303:312\t31\t+\n"
    R2 = "chr5\t10113\t10213\tSIM:chr5:10089761:151:0:127:706\t30\t-\n"

    x = 0
    boolean = []
    while x < 100:
        seed = random.randint(0, 5_000_000_000)
        x += 1
        a = discard(maxHashValue, seed, R1)
        b = discard(maxHashValue, seed, R2)
        boolean.append(a != b)

    assert True in boolean


# def test_mm_calc():


# ------------------------------------
# Test edge cases
# ------------------------------------

def test_empty_input():
    """
    Test that the function wc() raises a SystemExit error when given an empty file
    """
    subdir = Path("test_")
    subdir.mkdir(parents=True, exist_ok=True)
    files = ["test_empty.bed"]
    increment = 50

    # setup
    logging.basicConfig(level=20)
    warn = logging.warning
    info = logging.info

    cat_bed(files = files, control = False, subdir = subdir, info = info, warn = warn)

    # assert that 
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        wc(increment, subdir, info, warn, False)

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1


def test_incorrect_input():
    """
    Test that the function cat_bed() raises a SystemExit error when given a
    file with fewer than four tab-delimited columns
    """
    subdir = Path("test_")
    subdir.mkdir(parents=True, exist_ok=True)
    files = ["test_wrong.bam"]

    # setup
    logging.basicConfig(level=20)
    warn = logging.warning
    info = logging.info

    # assert that 
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cat_bed(files = files, control = False, subdir = subdir, info = info, warn = warn)

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1
    
