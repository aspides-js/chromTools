#!/usr/bin/python


# ------------------------------------
# modules
# ------------------------------------

import random
import logging

import pytest

from chromTools.complete_cmd import cat_bed, discard, wc

# --------------------------------------------------------------------------------#
## happy path


def test_cat_bed():
    files = ["test/test_s1.bed", "test/test_s2.bed"]
    subdir = "test/test_"
    increment = 50

    # setup
    logging.basicConfig(level=20)
    warn = logging.warning
    info = logging.info

    cat_bed(files, subdir, info)

    total, nfile = wc(increment, subdir, info, warn, False)
    assert total == 150


def test_discard_paired():
    ## set up
    maxHashValue = -9038904596117680288  # params(0.01)
    seed = 10

    R1 = "chr1\t10060\t10160\tSIM:chr1:10060:141:0:303:312/1\t31\t+\n"
    R2 = "chr1\t10102\t10202\tSIM:chr1:10060:141:0:303:312/2\t31\t-\n"

    a = discard(maxHashValue, seed, R1)
    b = discard(maxHashValue, seed, R2)

    assert a == b


def test_discard_single():
    ## set up
    maxHashValue = -9038904596117680288  # params(0.01)
    seed = 10

    R1 = "chr1\t10060\t10160\tSIM:chr1:10060:141:0:303:312\t31\t+\n"
    R2 = "chr1\t10102\t10202\tSIM:chr1:10060:141:0:303:312\t31\t-\n"

    a = discard(maxHashValue, seed, R1)
    b = discard(maxHashValue, seed, R2)

    assert a == b


def test_discard_singlediff():
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


# def test_


# --------------------------------------------------------------------------------#
## edge cases
def test_genome_no_size():
    pass


def test_wrong_genome():
    pass


def test_empty_data():
    files = ["test/test_empty.bed"]
    subdir = "test/test_"
    increment = 50

    # setup
    logging.basicConfig(level=20)
    warn = logging.warning
    info = logging.info

    cat_bed(files, subdir, info)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        wc(increment, subdir, info, warn, False)

    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 1


# ---------------------------------------------------------------------------#

## chmm binarise


def test_load_grid():
    pass
