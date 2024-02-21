#!/usr/bin/python

# ------------------------------------
# Modules
# ------------------------------------

from libc.stdlib cimport malloc, free, atoi, calloc
from libc.string cimport strtok, strncpy, strlen
from libc.stdio cimport FILE, fopen, fclose, fgets, fwrite, printf
from cpython.string cimport PyString_FromStringAndSize
import cython

import numpy as np
import mmh3
import sys
cimport numpy as np

from chromTools.mmh cimport murmurhash_32

#import chromTools.complete_cmd

# ------------------------------------
# Misc function(s)
# ------------------------------------

cdef int c_strand_check(bytes szstrand, int szbeginline, int szendline, int noffsetleft, int noffsetright, int nbinsize, int nshift):
    """
    Determine the bin index based on strand information.

    This function calculates the bin index for a given genomic region, taking into account the strand information.

    :param bytes szstrand: The strand information as bytes, either b"+\\n" for positive strand or b"-\\n" for negative strand.
    :param int szbeginline: The starting position of the genomic region.
    :param int szendline: The ending position of the genomic region.
    :param int noffsetleft: Offset for the left coordinate so it is 0-based inclusive
    :param int noffsetright: Offset for the right coordinate based on bed files being 0-based but not inclusive
    :param int nbinsize: The size of the bins. default 200
    :param int nshift: Shift value for calculating the bin index.

    :return: The calculated bin index based on the provided parameters.
    :rtype: int
    """
    cdef int nbin 
    if szstrand == b"+\n":
        nbin = (szbeginline - noffsetleft + nshift) // nbinsize
    elif szstrand == b"-\n":
        nbin = (szendline - noffsetright - nshift) // nbinsize
    return nbin


cpdef read_to_grid(file_path, 
        dict hmchrom,
        int nchromcol, 
        int nstrandcol, 
        int nbegincol, 
        int nendcol, 
        int noffsetleft, 
        int noffsetright, 
        int nshift, 
        int nmark,
        int nbinsize, 
        np.ndarray[np.int64_t, ndim=3] grid, 
        list bpresent):
    """
    Reads data from a file and populates a grid with the data.

    This function reads the specified file line by line, processes each line, and adds the data to the grid.
    The grid is a three-dimensional array that represents each chromosome, the lengths of each chr divided by binsize, 
    and the number of marks (always 1 in chromTools_complete).

    :param str file_path: The path to the file to be read.
    :param dict hmchrom: A dictionary containing mappings of chromosome names to chromosome indices.
    :param int nchromcol: The column index for the chromosome information in the file.
    :param int nstrandcol: The column index for the strand information in the file.
    :param int nbegincol: The column index for the start position of the data in the file.
    :param int nendcol: The column index for the end position of the data in the file.
    :param int noffsetleft: The amount to subtract from the left coordinate to make it 0-based inclusive.
    :param int noffsetright: The value to add to the right coordinate for 0-based, non-inclusive bed files.
    :param int nshift: The shift value to apply to the read data.
    :param int nmark: The mark index for the data in the grid.
    :param int nbinsize: The number of base pairs in a bin.
    :param numpy.ndarray grid: A three-dimensional numpy array for storing the data.
    :param list bpresent: A list that keeps track of the presence of data for each chromosome.

    :raises FileNotFoundError: If the specified file_path does not exist or cannot be opened.
    :raises ValueError: If an invalid strand is encountered in the file.

    :returns: The updated grid with the data read from the file and a list indicating the presence of data for each chromosome.
    :rtype: numpy.ndarray, list
    """
    cdef bytes szchrom, szstrand
    cdef int nbin, nchrom, nbeginline, nendline
    cdef int nmaxindex
    cdef list szLine
    cdef int shape = grid.shape[1]

    # Open the file
    cdef FILE *file = fopen(file_path.encode(), "r")
    
    # Error handling if file cannot be opened
    if file is NULL:
        raise FileNotFoundError(f"File '{file_path}' not found.")
    
    # Read the file line by line
    cdef char line[1024]
    while fgets(line, sizeof(line), file) is not NULL:
        # Process the line as needed
        szLine = line.split(b"\t")
        szchrom = szLine[nchromcol]
        objInt = hmchrom.get(szchrom)

        if objInt is not None:
            nchrom = objInt
            szstrand = szLine[nstrandcol]
            szbeginline = atoi(szLine[nbegincol])
            szendline = atoi(szLine[nendcol])

            # run check for positive/negative strand to output nbin
            nbin = c_strand_check(szstrand, szbeginline, szendline, noffsetleft, noffsetright, nbinsize, nshift)

            if nbin >= 0 and nbin < shape:
                grid[nchrom, nbin, nmark] += 1
                bpresent[nchrom] = 1

    # Close the file
    fclose(file)
    return grid, bpresent



# --------------------------------------------------------------------------------#                         

cpdef int c_subsample(str file_path, str outf_path, int a, int seed):
    """
    Subsample reads from a file based on a hash threshold.

    This function reads a file line by line, extracts the readname, hashes it using the MurmurHash3_32 algorithm, 
    and writes the original line to a new file if the hash value is below a specified threshold, calculated 
    based on the proportion of lines to be retained. If number above proportional cut-off (True), discard read.


    :param str file_path: Path to the input file.
    :param str outf_path: Path to the output file.
    :param int a: Hash threshold; reads with hash values above this threshold will be skipped.
    :param int seed: Seed for the MurmurHash3_32 algorithm.

    :return: The number of reads written to the output file.
    :rtype: int
    """
    # Open the file
    cdef FILE *file = fopen(file_path.encode(), "r")
    cdef FILE *outf = fopen(outf_path.encode(), "w")

    cdef char* c_readname 
    cdef bytes readname
    cdef int hashInt, c_readname_len

    # Read the file line by line    
    cdef char line[1024]
    cdef int reads = 0
    while fgets(line, sizeof(line), file) is not NULL:
        # extract readname from line (4th element, strip indication of readnames 1 and 2)
        readname = line.split(b"\t")[3].rsplit(b"/")[0]
        
        # convert to c-type string
        c_readname = readname
        c_readname_len = strlen(c_readname)

        # use c function murmurhash to hash string
        hashInt = murmurhash_32(c_readname, c_readname_len, seed)
        
        # if int is above threshold then skip, otherwise write
        if hashInt > a: 
            continue
        else:
            fwrite(line, sizeof(char), len(line), outf)
            reads += 1

    fclose(file)
    fclose(outf)
    return reads    

