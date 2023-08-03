#!/usr/bin/python

# ------------------------------------
# Modules
# ------------------------------------

from libc.stdlib cimport malloc, free
from libc.stdio cimport FILE, fopen, fclose, fgets
from cpython.string cimport PyString_FromStringAndSize


import numpy as np
cimport numpy as np

# ------------------------------------
# Misc function(s)
# ------------------------------------

def read_to_grid(file_path, 
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
