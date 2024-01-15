#!/usr/bin/python

# ------------------------------------
# Modules
# ------------------------------------

from libc.stdlib cimport malloc, free
from libc.stdio cimport FILE, fopen, fclose, fgets, fwrite
from cpython.string cimport PyString_FromStringAndSize


import numpy as np
import mmh3
cimport numpy as np

import chromTools.complete_cmd

# ------------------------------------
# Misc function(s)
# ------------------------------------

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
    The grid is a three-dimensional array that represents the data in a structured manner based on the provided
    parameters.

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
    cdef char* szchrom
    cdef char* szstrand
    cdef int nbin, nchrom
    cdef int nmaxindex
    cdef list szLine

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
            if szstrand == b"+\n":
                nbin = (
                    int(szLine[nbegincol]) - noffsetleft + nshift
                ) // nbinsize
            elif szstrand == b"-\n":
                nbin = (
                    int(szLine[nendcol]) - noffsetright - nshift
                ) // nbinsize
            else:
                raise ValueError(f"{szstrand} is an invalid strand!")
            if nbin >= 0 and nbin < grid.shape[1]:
                grid[nchrom, nbin, nmark] += 1
                bpresent[nchrom] = 1


    # Close the file
    fclose(file)
    return grid, bpresent


cpdef int c_subsample(str file_path, str outf_path, long a, long seed):
    """Generate a random hash from readname using seed. If number above proportional cut-off (True), discard read.

    Args:
            maxHashValue (int): Threshold above which reads are discarded
            seed (int): Random seed
            line (str): Read/line in file

    Returns:
            bool: Boolean specifying if readname is below or above discard threshold
    """
    # Open the file
    cdef FILE *file = fopen(file_path.encode(), "r")
    cdef FILE *outf = fopen(outf_path.encode(), "w")

    cdef bytes readname
    cdef long hashInt

    # Error handling if file cannot be opened
    if file is NULL:
        raise FileNotFoundError(f"File '{file_path}' not found.")
        
    # Read the file line by line
    cdef char line[1024]
    cdef int reads = 0

    while fgets(line, sizeof(line), file) is not NULL:
        readname = line.split(b"\t")[3].rsplit(b"/")[0]
        hashInt = mmh3.hash64(readname, seed)[0]
        if hashInt > a: #c_discard(a, seed line):
            continue
        else:
            fwrite(line, sizeof(char), len(line), outf)
            reads += 1

    fclose(file)
    fclose(outf)
    return reads    


# cdef c_discard(long maxHashValue, long seed, char line):

#     cdef bytes readname = line.split(b"\t")[3].rsplit(b"/")[
#         0
#     ]  # extract readname, remove everything after '/' (read pair if paired)

#     if len(readname) < 2:
#         raise ValueError(
#             f"Readname length is {len(readname)}! Check fourth column of input files contains a valid readname."
#         )

#     cdef long hashInt = mmh3.hash64(readname, 10)[0]

#     return hashInt > maxHashValue


