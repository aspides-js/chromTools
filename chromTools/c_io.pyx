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

from chromTools.mmh cimport murmurhash

#import chromTools.complete_cmd

# ------------------------------------
# Misc function(s)
# ------------------------------------

cdef int c_strand_check(bytes szstrand, int szbeginline, int szendline):
    cdef int nbin 
    if szstrand == b"+\n":
        nbin = (szbeginline - 0 + 100) // 200
    elif szstrand == b"-\n":
        nbin = (szendline - 1 - 100) // 200
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
            nbin = c_strand_check(szstrand, szbeginline, szendline)
            if nbin >= 0 and nbin < shape:
                grid[nchrom, nbin, nmark] += 1
                bpresent[nchrom] = 1

    # Close the file
    fclose(file)
    return grid, bpresent



# --------------------------------------------------------------------------------#                         

cdef int convert_fourth_ele_to_hash(char *remaining_line):
    # memory address is remaining_line
    # *remaining_line is used to access the value (character) that the pointer is pointing to.

    #cython_remaining_line[0] accesses the first character of the Cython char array or memoryview.
    #&cython_remaining_line[0] retrieves the memory address of the first character in the Cython char array or memoryview.
    #cdef char *start_pointer = remaining_line
    cdef int chars_in_ele = 0

    while (remaining_line[chars_in_ele] != b'\0'):
        chars_in_ele+=1
        if (remaining_line[chars_in_ele] == b'/' or remaining_line[chars_in_ele] == b'\t'):
            break
    #cdef char[:chars_in_ele] start_pointer = &remaining_line[0] doesnt worj 
    cdef char* start_pointer = <char *> malloc(chars_in_ele + 1)  # Allocate memory, +1 for null-termination

    # Copy characters into the buffer
    for i in range(chars_in_ele):
        start_pointer[i] = remaining_line[i]

    # Null-terminate the buffer
    start_pointer[chars_in_ele] = b'\0'

    #while (remaining_line != b'\0'):
    #    chars_in_ele+=1
    #    if (remaining_line == b'/' or remaining_line == b'+'):
    #        break
    #    else:
    #        remaining_line+=1

    #// Here you would call the murmur hash with key being 'start_pointer' and len being chars_in_ele
    #// Hopefully should work
    #// Else, here it is copied into a seperate char array
    cdef int hashInt = murmurhash(start_pointer, chars_in_ele, 10)
    #free(start_pointer)
    #free(remaining_line)
    #cdef char *fourth_ele = <char *>calloc(chars_in_ele + 1, sizeof(char))
    #fourth_ele = calloc(chars_in_ele + 1, sizeof(char))
    #strncpy(fourth_ele, start_pointer, chars_in_ele)
    #fourth_ele[chars_in_ele - 1] = b'\0'
    return hashInt





cdef int convert_line_to_hash(char *line):
    cdef int str_len = strlen(line)
    cdef char chr
    cdef int delim_counter = 0
    cdef int hashInt

    for i in range(str_len):
        chr = line[i]
        if (chr == b"\t"):
            delim_counter+=1
        if (delim_counter == 3):
            # Assuming 'line' is a Cython memoryview of characters
            #next_char_ptr = &line[i + 1]

            # Access the character at the next position using the pointer
            #next_char = next_char_ptr[0]
            #fourth_ele = convert_fourth_ele_to_hash(next_char)
            #fourth_ele = convert_fourth_ele_to_hash(&line[i + 1])
            hashInt = convert_fourth_ele_to_hash(&line[i + 1])
            break
    return hashInt 


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

    cdef char* fourth_ele 
    cdef bytes readname
    cdef long hashInt

    # Read the file line by line    
    cdef char line[1024]
    cdef int reads = 0
    print("test print", outf_path)
    while fgets(line, sizeof(line), file) is not NULL:
        #readname = line.split(b"\t")[3].rsplit(b"/")[0]
        hashInt = convert_line_to_hash(line)
        #hashInt = mmh3.hash(readname, seed)
        if hashInt > a: #c_discard(a, seed line):
            continue
        else:
            fwrite(line, sizeof(char), len(line), outf)
            reads += 1

    fclose(file)
    fclose(outf)
    print("done", outf_path)
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
