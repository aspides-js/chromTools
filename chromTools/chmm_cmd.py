#!/usr/bin/python


# ------------------------------------
# Modules
# ------------------------------------

import math
import os
import pathlib
import time
import sys
from itertools import compress

import numba as nb
import numpy as np

import chromTools.c_io

# ------------------------------------
# Main function
# ------------------------------------


def make_binary_data_from_bed(
    n, options
):  
    """Binarize BED data, both directly and with control.

    :param n: File number.
    :type n: int
    :param options: Namespace object.
    :type options: object
    :param control: If not False (default), control contains a list of control files.
    :type control: bool or list
    :param szcontroldir: The directory containing the bed files with the control data. If set to None, no control data is used.
    :type szcontroldir: str
    :param szchromlengthfile: The name of the two-column file containing chromosome and strand information.
    :type szchromlengthfile: str
    :param szmarkdir: The directory containing the bed files with the regular mark data.
    :type szmarkdir: str
    :param nflankwidthcontrol: Specifies the number of bins used in both directions to estimate the background. This attribute is only relevant if control data is being used.
    :type nflankwidthcontrol: int
    :param nshift: The number of bases a read should be shifted in the 5' to 3' direction of a read.
    :type nshift: int
    :param noffsetleft: The amount that should be subtracted from the left coordinate to make it 0-based inclusive.
    :type noffsetleft: int
    :param noffsetright: The amount that should be subtracted from the right coordinate to make it 0-based inclusive.
    :type noffsetright: int
    :param dpoissonthresh: The tail probability threshold on the Poisson distribution.
    :type dpoissonthresh: float
    :param dfoldthresh: The fold threshold required for a present call.
    :type dfoldthresh: float
    :param bcontainsthresh: If True, the Poisson cutoff should be the highest value that still contains the dpoissonthresh probability. If False, it requires strictly greater.
    :type bcontainsthresh: bool
    :param npseudocountcontrol: An integer pseudocount that is uniformly added to every interval to smooth the control data.
    :type npseudocountcontrol: int
    :param nbinsize: The number of base pairs in a bin
    :type nbinsize: int
    :param dcountthresh: The absolute signal threshold for a present call.
    :type dcountthresh: float

    :raises ValueError: Invalid line found in the chromosome length file if fewer than two cols found

    :return: A tuple containing the count of present marks across included chromosomes and the total number of bins
             across included chromosomes.
    :rtype: tuple[int, int]

    """
    start_time = time.time()
    nbinsize = 200
    szcell = f"{n}"
    szmark = "mark"
    hscells = set([szcell])
    hsmarks = set([szmark])
    hmfiles = {f"{szcell}\t{szmark}": [f"subsampled.{n}.bed"]}

    ## reads in the chromosome length information file
    # the first column of this file is the chromosome and the second is the chromsome length
    with open(options.szchromlengthfile) as brchrom:
        allines = brchrom.readlines()
        chroms = []  # stores the chromosome name
        lengths = []  # stores the chromosome length
        hmchrom = {}  # stores a mapping of chromosome to chromosome index
        for szLine in allines:
            szLine = szLine.strip().split()
            if len(szLine) < 2:
                raise ValueError(f"Invalid line found in {options.szchromlengthfile}")
            chrom = szLine[0]
            chroms.append(chrom)
            hmchrom[chrom] = len(chroms) - 1
            lengths.append(int(szLine[1]))

    ## encode to bytes for cython
    chmchrom = {y.encode("utf-8"): hmchrom.get(y) for y in hmchrom.keys()}
    nchroms = len(chroms) # number of chromosomes

    # loads all the marks in hsmarks into the array marks and then sorts it
    nummarks = len(hsmarks)
    marks = list(hsmarks)
    marks.sort()

    # generates a three dimensional array with each chromosome, the lengths of each chr divided by binsize, and the number of marks (always 1 in complete)
    # grid is a regular np array with Dim 2 as the length of chromosome 1/binsize (longest) 
    grid = np.zeros((nchroms, lengths[0] // nbinsize, nummarks), dtype=int)

    # set variables for control data
    # gridcontrol - generates a three dimensional array with each chromosome, the lengths of the largest chr divided by binsize, and the number of marks (always 1 in complete)
    if options.control:
        gridcontrol = np.ones((nchroms, lengths[0] // nbinsize, nummarks), dtype=int)
        sumgridcontrol = np.zeros((nchroms, lengths[0] // nbinsize, nummarks), dtype=int)
        bpresentcontrol = [False] * nchroms
        hmfilescontrol = {f"{szcell}\t{szmark}": ["subsampled.ctrl.bed"]}

    # ----------------------------------------

    for szcell in hscells:
        options.info(f"{szcell}: start grid")
        bpresent = [False] * nchroms

        # loading data for the cell type
        grid, bpresent = cload_grid(
            options.info,
            grid,
            bpresent,
            marks,
            options.nshift,
            nbinsize,
            options.noffsetleft,
            options.noffsetright,
            hmfiles,
            szcell,
            options.szmarkdir,
            chmchrom,
            0,
            options.control,
            lengths,
        )

        # ensure bpresent is list of bools
        bpresent = list(map(bool, bpresent))

        if options.control:
            # we have control data loading cell type data for that
            gridcontrol, bpresentcontrol = cload_grid(
                options.info,
                gridcontrol,
                bpresentcontrol,
                marks,
                options.nshift,
                nbinsize,
                options.noffsetleft,
                options.noffsetright,
                hmfilescontrol,
                szcell,
                options.szcontroldir,
                chmchrom,
                options.npseudocountcontrol,
                options.control,
                lengths,
            )


    nummarks_m1 = nummarks - 1

    count, total = 0, 0  # number of marks & total bins in each chr
    if options.control:
        # binarization will be based on control data
        sumgridcontrol = window_sum_grid(gridcontrol, 
                                        sumgridcontrol, 
                                        options.nflankwidthcontrol, 
                                        lengths, 
                                        nbinsize)

        # determiming thresholds for each mark and background depth
        thresholds = determine_mark_thresholds_from_binned_data_array_against_control(
            grid,
            sumgridcontrol,
            bpresent,
            bpresentcontrol,
            options.dpoissonthresh,
            options.dfoldthresh,
            options.bcontainsthresh,
            options.dcountthresh,
            lengths,
            nbinsize,
        )

        # extract number of bins in grid that are greater than threshold
        # sumgridcontrol is used as an index of threshold
        # create bool array for chr present in mark and control
        bpresentboth = [a==1 and b==1 for a, b in zip(bpresent, bpresentcontrol)]

        # reduce grid and sumgridcontrol to only those chrm present for both mark and control
        grid = grid[bpresentboth]
        sumgridcontrol = sumgridcontrol[bpresentboth]
        filt_lengths = list(compress(lengths, bpresentboth))

        # create a grid where each value is the threshold for that position (from sumgridcontrol indices)
        flat_indices = sumgridcontrol.flatten()
        flat_values = np.asarray(thresholds[0])[flat_indices]
        threshold_grid = flat_values.reshape(sumgridcontrol.shape)

        # true for all values where grid is higher than threshold: 
        mask = threshold_grid <= grid

        # write in false to all values beyond chromosome lengths
        for nchrom in range(sum(bpresentboth)):
           mask[nchrom][filt_lengths[nchrom] // nbinsize : ] = False
           
        # count all and all true
        count = np.sum(mask)
        total = sum([x // nbinsize for x in filt_lengths])

    else:  ## if no control file
        total, thresholds = determine_mark_thresholds_from_binned_data_array(
            grid,
            bpresent,
            options.dpoissonthresh,
            options.dfoldthresh,
            options.bcontainsthresh,
            options.dcountthresh,
            lengths, 
            nbinsize,
        )
        # extract number of bins in grid which are greater than threshold
        ## this replaces all of the old functionality writing to a file
        ## total now comes from determine_mark_thres
        count = np.sum(thresholds[nummarks_m1] <= grid)

    return count, total


# --------------------------------------


def cload_grid(
    info,
    grid,
    bpresent,
    marks,
    nshift,
    nbinsize,
    noffsetleft,
    noffsetright,
    hmfiles,
    szcell,
    szmarkdir,
    hmchrom,
    ninitval,
    bcontrol,
    lengths,
):
    """Converts read level data for a cell into integer count information.

    This function populates the provided grid with data obtained from mark files. It iterates through mark files for each
    cell type and mark, calculates the bin index for each data point, and updates the corresponding grid cell.

    :param grid: The array to load the read counts into.
    :type grid: numpy.ndarray
    :param bpresent: Indicates if there is a read for a chromosome with an index in hmchrom.
    :type bpresent: list
    :param marks: Contains the names of the header marks.
    :type marks: list
    :param nshift: The number of bases a read should be shifted in the 5' to 3' direction of a read.
    :type nshift: int
    :param nbinsize: The number of base pairs in a bin.
    :type nbinsize: int
    :param noffsetleft: The amount that should be subtracted from the left coordinate so it is 0-based inclusive.
    :type noffsetleft: int
    :param noffsetright: The amount that should be subtracted from the right coordinate so it is 0-based inclusive.
    :type noffsetright: int
    :param hmfiles: Maps cell and mark to an actual read file.
    :type hmfiles: dict
    :param szcell: The cell we are interested in.
    :type szcell: str
    :param szmarkdir: The directory with bedfiles to read.
    :type szmarkdir: str
    :param hmchrom: Maps chromosome names to an index.
    :type hmchrom: dict
    :param ninitval: The value that the data should be initialized to.
    :type ninitval: int
    :param bcontrol: if True, data is control data
    :type bcontrol: bool
    :param lengths: lengths of all chroms
    :type lengths: list

    :return: grid (numpy.ndarray): The updated grid with the read counts.

    Calculation Details:
        - Initializes necessary variables and data structures.
        - Iterates over each mark file for each cell type.
        - Checks for the presence of control data or missing data and handles accordingly.
        - Retrieves mark file data, extracts relevant information, and calculates the bin index.
        - Updates the corresponding grid cell with the calculated index and mark information.
        - Sets the presence flags for chromosomes and marks based on the loaded data.

    Raises:
        ValueError: If the column index exceeds the maximum index for the mark file data.
        ValueError: If an invalid strand is encountered in the mark file data.

    """
    nummarks = len(grid[0][0])

    # columns
    nchromcol, nbegincol, nendcol, nstrandcol, nmaxindex = 0, 1, 2, -1, -1

    # going through all the mark files in each cell type
    for nmark in range(nummarks):
        alfiles = hmfiles.get(f"{szcell}\t{marks[nmark]}")
        info(f"alfiles is {alfiles}")
        if alfiles is None:
            if bcontrol:
                print(
                    f"Warning did not find control data for {szcell} {marks[nmark]} treating as missing"
                )
            else:
                print(
                    f"Warning did not find data for {szcell} {marks[nmark]} treating as missing"
                )

            if not bcontrol:
                # slight efficiency improvement here in v1.04
                for nchrom in range(len(grid)):
                    for nbin in range(lengths[nchrom]//nbinsize):
                        grid[nchrom, nbin, nmark] = -1
        else:
            for szfile in alfiles:
                grid, bpresent = chromTools.c_io.read_to_grid(
                    os.path.join(szmarkdir, szfile),
                    hmchrom,
                    nchromcol,
                    nstrandcol,
                    nbegincol,
                    nendcol,
                    noffsetleft,
                    noffsetright,
                    nshift,
                    nmark,
                    nbinsize,
                    grid,
                    bpresent,
                )
    return grid, bpresent


def determine_mark_thresholds_from_binned_data_array(
    grid, bpresent, dpoissonthresh, dfoldthresh, bcontainsthresh, dcountthresh, lengths, nbinsize
):
    """Determines the Poisson cutoffs based on the provided data.

    This function calculates the thresholds for each mark based on the provided binned data array. It uses the Poisson
    distribution and other parameters to determine the cutoffs.

    :param grid: The integer data values from which to determine the Poisson cutoffs.
    :type grid: numpy.ndarray
    :param bpresent: A vector indicating which indices of 'grid' to include in the analysis.
    :type bpresent: numpy.ndarray
    :param dpoissonthresh: The tail probability threshold on the Poisson distribution.
    :type dpoissonthresh: float
    :param dfoldthresh: The fold threshold required for a present call.
    :type dfoldthresh: float
    :param bcontainsthresh: If True, the Poisson cutoff should be the highest value that still contains the 'dpoissonthresh' probability.
                            If False, it requires strictly greater.
    :type bcontainsthresh: boolean
    :param dcountthresh: The absolute signal threshold for a present call.
    :type dcountthresh: float

    :return:  thresholds (list): Each element in this list represents the Poisson cutoff for a specific mark.

    Calculation Details:
        - Initializes variables and data structures.
        - Computes the sum of data values for each mark across all relevant bins and chromosomes.
        - Calculates the total number of bins considered for threshold determination.
        - Computes the threshold for each mark based on the Poisson distribution.
        - Adjusts the threshold based on the bcontainsthresh parameter.
        - Sets the final threshold for each mark as the maximum value among various calculations.

    Note:
        - The thresholds are computed based on the Poisson distribution and various parameters provided.

    """
    dcumthreshold = 1 - dpoissonthresh
    nummarks = len(grid[0][0])
    ntotallocs = 0
    sumtags = np.array([0] * nummarks)
    thresholds = [0] * nummarks

    for nchrom in range(len(grid)):
        if bpresent[nchrom]:
            sumtags[0] += np.sum(
                grid[nchrom]
            )  ## would need to be adapted if nummarks > 1
            ntotallocs += lengths[nchrom] // nbinsize

    for nj in range(len(sumtags)):
        dlambda = sumtags[nj] / ntotallocs
        dcum, nthresh, dlogfactorial = 0, 0, 0

        while dcum <= dcumthreshold:
            dprob = math.exp(math.log(dlambda) * nthresh - dlambda - dlogfactorial)
            dcum += dprob
            nthresh += 1
            
            dlogfactorial += math.log(nthresh)

        if bcontainsthresh:
            nthresh -= 1

        thresholds[nj] = max(int(dfoldthresh * dlambda), nthresh, int(dcountthresh))
    return ntotallocs, thresholds


def determine_mark_thresholds_from_binned_data_array_against_control(
    grid,
    sumgridcontrol,
    bpresent,
    bpresentcontrol,
    dpoissonthresh,
    dfoldthresh,
    bcontainsthresh,
    dcountthresh,
    lengths,
    nbinsize
):
    """Determines the Poisson cutoffs based on the provided data.

    This function calculates the thresholds for each mark by comparing the binned data array against the control data.
    It uses the Poisson distribution to determine the thresholds based on various parameters.

    :param grid: The integer data values from which to determine the Poisson cutoffs.
    :type grid: numpy.ndarray
    :param sumgridcontrol: The control data to which the thresholds will be relative.
    :type sumgridcontrol: numpy.ndarray
    :param bpresent: A vector indicating which indices of 'grid' to include in the analysis.
    :type bpresent: numpy.ndarray
    :param bpresentcontrol: A vector indicating which indices of 'gridcontrol' to include in the analysis.
    :type bpresentcontrol: numpy.ndarray
    :param dpoissonthresh: The tail probability threshold on the Poisson distribution.
    :type dpoissonthresh: float
    :param dfoldthresh: The fold threshold required for a present call.
    :type dfoldthresh: float
    :param bcontainsthresh: If True, the Poisson cutoff should be the highest value that still contains the 'dpoissonthresh' probability.
                            If False, it requires strictly greater.
    :type bcontainsthresh: boolean
    :param dcountthresh: The absolute signal threshold for a present call.
    :type dcountthresh: float

    :return:  thresholds (list): Each element in this list represents the Poisson cutoff for a specific mark.

    Calculation Details:
        - Initializes variables and data structures.
        - Computes the total number of reads for each mark and its matched control.
        - Determines the maximum control value found for each mark.
        - Determines the background control values encountered for each mark.
        - Calculates thresholds for each mark based on the observed data and control values.

    Note:
        - The thresholds are computed based on the Poisson distribution and various parameters provided.

    """
    dcumthreshold = 1 - dpoissonthresh

    nummarks = len(grid[0][0])
    sumtags = np.array([0] * nummarks)
    sumtagscontrol = np.array([0] * nummarks)
    thresholds = [0] * nummarks
    maxcontrol = np.array([0] * nummarks)

    hscontrol_np = np.empty(0, dtype=int)
    for nchrom in range(len(grid)):
        if bpresent[nchrom] and bpresentcontrol[nchrom]:
            sumtags[0] += np.sum(grid[nchrom])
            sumtagscontrol[0] += np.sum(sumgridcontrol[nchrom][0:lengths[nchrom] // nbinsize])
            hscontrol_np = np.concatenate(
                (hscontrol_np, sumgridcontrol[nchrom][0:lengths[nchrom] // nbinsize]), axis=None
            )
            hscontrol_np = np.unique(hscontrol_np)
            maxcontrol[0] = max(hscontrol_np)
            # --------------
    hscontrol = []
    hscontrol.append(set(hscontrol_np))
    print(f"sumtgas is {sumtags}")
    print(f"sumtgascontrol is {sumtagscontrol}")
    print(f"maxcontrol is: {maxcontrol}")


    for nmark in range(len(sumtags)):
        thresholds[nmark] = [0] * (maxcontrol[nmark] + 1)
        davgratio = sumtags[nmark] / sumtagscontrol[nmark]
        thresholds[nmark][0] = max(int(dcountthresh), 1)
        for nbackground in range(1, maxcontrol[nmark] + 1):
            # bug fixed in 1.14 that changes less than to less than equal
            if nbackground in hscontrol[nmark]:
                dlambda = davgratio * nbackground
                dcum, nthresh, dlogfactorial = 0, 0, 0
                while dcum <= dcumthreshold:
                    dprob = math.exp(
                        math.log(dlambda) * nthresh - dlambda - dlogfactorial
                    )
                    dcum += dprob
                    nthresh += 1
                    dlogfactorial += math.log(nthresh)

                if bcontainsthresh:
                    # decreasing to include the dpoissonthreshold probability
                    nthresh -= 1

                thresholds[nmark][nbackground] = max(
                    int(dfoldthresh * dlambda), nthresh, int(dcountthresh)
                )

    return thresholds


@nb.njit
def window_sum_grid(gridcontrol, sumgridcontrol, nflankwidthcontrol, lengths, nbinsize):
        """Calculates the windowed sum of values in the control grid.

    This function iterates over chromosomes, bin positions, and marks in the control grid. For each bin position,
    it sums the values within a window defined by the flank width. The resulting sum is stored in the corresponding
    position of the sum grid.

    :param gridcontrol: The control grid containing integer data values.
    :type gridcontrol: numpy.ndarray
    :param sumgridcontrol: The sum grid for storing the calculated sums.
    :type sumgridcontrol: numpy.ndarray
    :param nflankwidthcontrol: The flank width for the window (inclusive).
    :type nflankwidthcontrol: int

    :return: A tuple containing the updated control grid and sum grid.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]

    Iteration Details:
        - Iterates over each chromosome in the control grid.
        - For each chromosome, iterates over bin positions and marks.
        - Calculates the windowed sum for each bin position and mark.
        - Updates the sum grid with the calculated sums.

    Note:
        - The flank width determines the number of neighboring positions on each side of a bin to include in the sum.
          For example, if nflankwidthcontrol is set to 5, it includes 5 positions to the left and 5 positions to the
          right of the current bin (including the current bin itself).
    """
        for nchrom in range(len(gridcontrol)):
            bin_len = lengths[nchrom] // nbinsize
            for nbin in range(0, bin_len):
                nstart = max(0, nbin - nflankwidthcontrol)
                nend = min(nbin + nflankwidthcontrol, bin_len - 1)
                for nmark in range(len(sumgridcontrol[nchrom][nbin])):
                    nsum = 0
                    for nrow in range(nstart, nend + 1):
                        nval = gridcontrol[nchrom][nrow][nmark]
                        if nval > 0:
                            nsum += nval
                    sumgridcontrol[nchrom][nbin] = nsum
        return sumgridcontrol
