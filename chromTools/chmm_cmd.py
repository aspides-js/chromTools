#!/usr/bin/python


# ------------------------------------
# modules
# ------------------------------------

import math
import os
import pathlib
import time

import numba as nb
import numpy as np

# ------------------------------------
# Main function
# ------------------------------------


def make_control_grid(options):
    """
    Generate control grid and sum grid for binning control data.

    This function initializes starting options, reads the chromosome length information file, loads the marks,
    generates a control grid and sum grid, and loads control data into the grids.

    :param options: The options for grid generation.
    :type options: Namespace
    :return: The control grid, sum grid, and presence flags for control marks.
    :rtype: tuple[np.ndarray, np.ndarray, list[bool]]
    :raises ValueError: Invalid line found in the chromosome length file if fewer than two columns are found.
    """

    ## set starting options
    nbinsize = 200
    szcell = "control"
    szmark = "mark"
    hscells = set([szcell])
    hsmarks = set([szmark])
    hscellcontrol = hscells
    hscellnocontrol = set()
    hmfilescontrol = {f"{szcell}\t{szmark}": ["downsampled.ctrl.bed"]}

    ## reads in the chromosome length information file
    # the first column of this file is the chromosome and the second is the chromsome length
    with open(options.szchromlengthfile) as brchrom:
        allines = brchrom.readlines()
        chroms = []  # stores the chromosome name
        lengths = []  # stores the chromosome length
        hmchrom = {}  # stores a mapping of chromosome to chromosome index
        for szLine in allines:
            st = szLine.strip().split()
            if len(st) < 2:
                raise ValueError(f"Invalid line found in {options.szchromlengthfile}")
            chrom = st[0]
            length = int(st[1])
            chroms.append(chrom)
            hmchrom[chrom] = len(chroms) - 1
            lengths.append(length)

    # loads all the marks in hsmarks into the array marks and then sorts it
    nummarks = len(hsmarks)
    numcontrolmarks = -1
    marks = list(hsmarks)
    marks.sort()

    # generates a three dimensional array with each chromosome, the lengths of each chr divided by binsize, and the number of marks (always 1 in complete)
    gridcontrol = np.empty((len(chroms),), dtype=np.ndarray)
    sumgridcontrol = np.empty((len(chroms),), dtype=np.ndarray)
    bpresentcontrol = [False] * len(chroms)
    bpresentmarkscontrol = [False] * nummarks

    # ----------------------------------------

    for szcell in hscells:
        bmissing = szcell in hscellnocontrol
        if hscellcontrol is None:
            # no control data for this cell type
            bcontrolfile = False
        else:
            bcontrolfile = True
            if len(hscellcontrol) == 1 and not bmissing:
                # we have one control for all marks
                numcontrolmarks = 1
            else:
                # will allocate the full memory for all marks
                numcontrolmarks = nummarks

        if bcontrolfile:
            if gridcontrol[0] is None or len(gridcontrol[0][0]) != numcontrolmarks:
                # reallocate if changing array size
                # allowed to go between single and matched
                for ni in range(len(chroms)):
                    gridcontrol[ni] = np.ones(
                        (lengths[ni] // nbinsize, numcontrolmarks), dtype=int
                    )
                    sumgridcontrol[ni] = np.zeros(
                        (lengths[ni] // nbinsize, numcontrolmarks), dtype=int
                    )

            # we have control data loading cell type data for that
            print("LOAD GRID CONTROL")
            gridcontrol, bpresentcontrol, bpresentmarkscontrol = load_grid(
                gridcontrol,
                bpresentcontrol,
                bpresentmarkscontrol,
                marks,
                options.nshift,
                nbinsize,
                options.noffsetleft,
                options.noffsetright,
                hmfilescontrol,
                szcell,
                options.szcontroldir,
                hmchrom,
                options.npseudocountcontrol,
                True,
            )
    return gridcontrol, sumgridcontrol, bpresentcontrol


def make_binary_data_from_bed(
    n, gridcontrol, sumgridcontrol, bpresentcontrol, options
):
    """Binarize BED data, both directly and with control.

    :param n: File number.
    :type n: int
    :param options: Namespace object.
    :type options: object
    :param control: If not False (default), control contains a list of control files.
    :type control: bool or list
    :param szchromfile: The name of the file containing chromosome information. It is a two-column file with the first column representing the chromosome and the second column representing the chromosome length.
    :type szchromfile: str
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
    :param bcenterinterval: If True, the center of the read is used instead of shifting if the read has already been extended.
    :type bcenterinterval: bool
    :param noffsetleft: The amount that should be subtracted from the left coordinate to make it 0-based inclusive.
    :type noffsetleft: int
    :param noffsetright: The amount that should be subtracted from the right coordinate to make it 0-based inclusive.
    :type noffsetright: int
    :param szoutputsignaldir: If not None, the intermediate signal data will be printed to this directory.
    :type szoutputsignaldir: str
    :param szoutputbinarydir: The directory where the binarized data will be printed.
    :type szoutputbinarydir: str
    :param szoutputcontroldir: If not None, the intermediate control signal data will be printed to this directory.
    :type szoutputcontroldir: str
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
    :param szcolfields: A comma-delimited string indicating the 0-based columns of the chromosome, start, end,
                        and optionally strand position. If set to None, the default values are used (0, 1, 2) for
                        chromosome, start, and end, with the strand as the sixth column or the last column if fewer
                        columns are present.
    :type szcolfields: str
    :param dcountthresh: The absolute signal threshold for a present call.
    :type dcountthresh: float
    :param bbinarizebam: If True, reads files as BAM files; otherwise, reads files as BED files.
    :type bbinarizebam: bool

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
    hmfiles = {f"{szcell}\t{szmark}": [f"downsampled.{n}.bed"]}

    if not options.control:
        hscellnocontrol = hscells
        hscellcontrol = None
        bcontrol = False
    else:
        hscellcontrol = hscells
        hscellnocontrol = set()
        hmfilescontrol = {f"{szcell}\t{szmark}": ["downsampled.ctrl.bed"]}
        bcontrol = True

    ## reads in the chromosome length information file
    # the first column of this file is the chromosome and the second is the chromsome length
    with open(options.szchromlengthfile) as brchrom:
        allines = brchrom.readlines()
        chroms = []  # stores the chromosome name
        lengths = []  # stores the chromosome length
        hmchrom = {}  # stores a mapping of chromosome to chromosome index
        for szLine in allines:
            st = szLine.strip().split()
            if len(st) < 2:
                raise ValueError(f"Invalid line found in {options.szchromlengthfile}")
            chrom = st[0]
            length = int(st[1])
            chroms.append(chrom)
            hmchrom[chrom] = len(chroms) - 1
            lengths.append(length)

    # loads all the marks in hsmarks into the array marks and then sorts it
    nummarks = len(hsmarks)
    numcontrolmarks = -1
    marks = list(hsmarks)
    marks.sort()

    # generates a three dimensional array with each chromosome, the lengths of each chr divided by binsize, and the number of marks (always 1 in complete)
    grid = np.empty((len(chroms),), dtype=np.ndarray)
    for ni in range(len(chroms)):
        grid[ni] = np.zeros((lengths[ni] // nbinsize, nummarks), dtype=int)

    bpresentmarks = [False] * nummarks

    # ----------------------------------------

    for szcell in hscells:
        bpresent = [False] * len(chroms)
        # going through each declared cell type
        # hscellcontrol = hmfilescellcontrol.get(szcell)
        bmissing = szcell in hscellnocontrol
        if hscellcontrol is None:
            # no control data for this cell type
            bcontrolfile = False
        else:
            bcontrolfile = True
            if len(hscellcontrol) == 1 and not bmissing:
                # we have one control for all marks
                numcontrolmarks = 1
            else:
                # will allocate the full memory for all marks
                numcontrolmarks = nummarks

        # loading data for the cell type
        print(time.time() - start_time)
        print("LOAD GRID")
        grid, bpresent, bpresentmarks = load_grid(
            grid,
            bpresent,
            bpresentmarks,
            marks,
            options.nshift,
            nbinsize,
            options.noffsetleft,
            options.noffsetright,
            hmfiles,
            szcell,
            options.szmarkdir,
            hmchrom,
            0,
            bcontrol,
        )

        print(time.time() - start_time)
    nummarks_m1 = nummarks - 1

    count, total = 0, 0  # number of marks & total bins in each chr
    if bcontrolfile:
        # binarization will be based on control data
        print(time.time() - start_time)
        print("windowsumgrid_numba")
        for nchrom in range(len(gridcontrol)):
            sumgridcontrol_nchrom = sumgridcontrol[nchrom]
            for nbin in range(len(sumgridcontrol_nchrom)):
                (
                    gridcontrol[nchrom],
                    sumgridcontrol[nchrom][nbin],
                ) = window_sum_grid_numba(
                    gridcontrol[nchrom],
                    sumgridcontrol_nchrom[nbin],
                    nbin,
                    options.nflankwidthcontrol,
                )
        # determiming thresholds for each mark and background depth
        options.info("DETERMINE FROM BINNED CTL_numba")
        thresholds = (
            determine_mark_thresholds_from_binned_data_array_against_control(
                grid,
                sumgridcontrol,
                bpresent,
                bpresentcontrol,
                options.dpoissonthresh,
                options.dfoldthresh,
                options.bcontainsthresh,
                options.dcountthresh,
            )
        )
        print(time.time() - start_time)
        for nchrom in range(len(chroms)):
            if bpresent[nchrom] and bpresentcontrol[nchrom]:
                # we have both primary and control data for the mark
                szfile = pathlib.Path(
                    options.szoutputbinarydir,
                    f"{szcell}_{chroms[nchrom]}_binary.txt",
                )
                # print("Writing to file", szfile)
                with open(szfile, "w") as pw:
                    pw.write(f"{szcell}\t{chroms[nchrom]}\n")
                    for nmark in range(nummarks_m1):
                        pw.write(f"{marks[nmark]}\t")
                    pw.write(f"{marks[nummarks_m1]}\n")

                    for nbin in range(len(grid[nchrom])):
                        for nmark in range(nummarks_m1):
                            if numcontrolmarks == 1:
                                ncontrolval = sumgridcontrol[nchrom][nbin][0]
                            else:
                                ncontrolval = sumgridcontrol[nchrom][nbin][nmark]

                            # printing one if count exceeds background threshold
                            if not bpresentmarks[nmark]:
                                pw.write("2\t")
                            elif (
                                thresholds[nmark][ncontrolval]
                                <= grid[nchrom][nbin][nmark]
                            ):
                                pw.write("1\t")
                                count += 1
                                total += 1
                            else:
                                pw.write("0\t")
                                total += 1

                        if numcontrolmarks == 1:
                            ncontrolval = sumgridcontrol[nchrom][nbin][0]
                        else:
                            ncontrolval = sumgridcontrol[nchrom][nbin][nummarks_m1]

                        if not bpresentmarks[nummarks_m1]:
                            pw.write("2\n")
                        elif (
                            thresholds[nummarks_m1][ncontrolval]
                            <= grid[nchrom][nbin][nummarks_m1]
                        ):
                            pw.write("1\n")
                            count += 1
                            total += 1
                        else:
                            pw.write("0\n")
                            total += 1
        print(time.time() - start_time)
    else:  ## if no control file
        print("Treating as no control")
        options.info("DETERMINE FROM BINNED_numba")
        thresholds = determine_mark_thresholds_from_binned_data_array(
            grid,
            bpresent,
            options.dpoissonthresh,
            options.dfoldthresh,
            options.bcontainsthresh,
            options.dcountthresh,
        )
        print(time.time() - start_time)
        for nchrom in range(len(chroms)):
            if bpresent[nchrom]:
                szfile = os.path.join(
                    options.szoutputbinarydir,
                    szcell + "_" + chroms[nchrom] + "_binary.txt",
                )
                # print("Writing to file " + szfile)
                with open(szfile, "w") as pw:
                    pw.write(f"{szcell}\t{chroms[nchrom]}\n")
                    for nmark in range(len(marks) - 1):
                        pw.write(f"{marks[nmark]}\t")
                    pw.write(f"{marks[nummarks_m1]}\n")
                    for nbin in range(len(grid[nchrom])):
                        for nmark in range(nummarks_m1):
                            if not bpresentmarks[nmark]:
                                pw.write("2\t")
                            elif thresholds[nmark] <= grid[nchrom][nbin][nmark]:
                                pw.write("1\t")
                            else:
                                pw.write("0\t")
                        if not bpresentmarks[nummarks_m1]:
                            pw.write("2\n")
                        elif thresholds[nummarks_m1] <= grid[nchrom][nbin][nummarks_m1]:
                            pw.write("1\n")
                            count += 1
                            total += 1
                        else:
                            pw.write("0\n")
                            total += 1
        print(time.time() - start_time)
    return count, total


# --------------------------------------
def load_grid(
    grid,
    bpresent,
    bpresentmarks,
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
):
    """Converts read level data for a cell into integer count information.

    This function populates the provided grid with data obtained from mark files. It iterates through mark files for each
    cell type and mark, calculates the bin index for each data point, and updates the corresponding grid cell.

    :param grid: The array to load the read counts into.
    :type grid: numpy.ndarray
    :param bpresent: Indicates if there is a read for a chromosome with an index in hmchrom.
    :type bpresent: list
    :param bpresentmark: Indicates if the mark is present in the cell type.
    :type bpresentmark: list
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
    :param bcontrol: if True, data is control dataz
    :type bcontrol: bool

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

        if alfiles is None:
            if bcontrol:
                print(
                    f"Warning did not find control data for {szcell} {marks[nmark]} treating as missing"
                )
            else:
                print(
                    f"Warning did not find data for {szcell} {marks[nmark]} treating as missing"
                )

            bpresentmarks[nmark] = False
            if not bcontrol:
                # slight efficiency improvement here in v1.04
                for nchrom in range(len(grid)):
                    numbins = len(grid[nchrom])
                    for nbin in range(numbins):
                        grid[nchrom][nbin][nmark] = -1
        else:
            bpresentmarks[nmark] = True
            for nfile in range(len(alfiles)):
                szfile = alfiles[nfile]
                with open(os.path.join(szmarkdir, szfile), "r") as brbed:
                    for szLine in brbed:
                        szLineA = szLine.split()
                        szchrom = szLineA[nchromcol]
                        objInt = hmchrom.get(szchrom)

                        if nmaxindex >= len(szLineA):
                            raise ValueError(
                                f"Column index {nmaxindex} exceeds maximum index {len(szLineA) - 1} indices are 0-based"
                            )

                        if objInt is not None:
                            nchrom = objInt
                            szstrand = szLineA[nstrandcol]
                            if szstrand == "+":
                                nbin = (
                                    int(szLineA[nbegincol]) - noffsetleft + nshift
                                ) // nbinsize
                            elif szstrand == "-":
                                nbin = (
                                    int(szLineA[nendcol]) - noffsetright - nshift
                                ) // nbinsize
                            else:
                                raise ValueError(f"{szstrand} is an invalid strand!")
                            if nbin >= 0 and nbin < len(grid[nchrom]):
                                grid[nchrom][nbin][nmark] += 1
                                bpresent[nchrom] = True
    return grid, bpresent, bpresentmarks


def determine_mark_thresholds_from_binned_data_array(
    grid, bpresent, dpoissonthresh, dfoldthresh, bcontainsthresh, dcountthresh
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
            grid_nchrom = grid[nchrom]
            ntotallocs, sumtags = determine_sumtags_numba(
                grid_nchrom,
                nummarks,
                sumtags,
                ntotallocs,
            )
    print(f"sumtags is: {sumtags}")
    print(f"hscontrol is: {ntotallocs}")
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

    return thresholds


def determine_mark_thresholds_from_binned_data_array_against_control(
    grid,
    gridcontrol,
    bpresent,
    bpresentcontrol,
    dpoissonthresh,
    dfoldthresh,
    bcontainsthresh,
    dcountthresh,
):
    """Determines the Poisson cutoffs based on the provided data.

    This function calculates the thresholds for each mark by comparing the binned data array against the control data.
    It uses the Poisson distribution to determine the thresholds based on various parameters.

    :param grid: The integer data values from which to determine the Poisson cutoffs.
    :type grid: numpy.ndarray
    :param gridcontrol: The control data to which the thresholds will be relative.
    :type gridcontrol: numpy.ndarray
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
    numcontrolmarks = len(gridcontrol[0][0])
    sumtags = np.array([0] * nummarks)
    sumtagscontrol = np.array([0] * nummarks)
    thresholds = [0] * nummarks
    maxcontrol = np.array([0] * nummarks)
    x = 0
    for nchrom in range(len(grid)):
        x += len(grid[nchrom])

    hscontrol_np = np.empty([len(range(nummarks)), x], dtype=int)

    x = -1
    for nchrom in range(len(grid)):
        if bpresent[nchrom] and bpresentcontrol[nchrom]:
            grid_nchrom = grid[nchrom]
            gridcontrol_nchrom = gridcontrol[nchrom]
            (
                hscontrol_np,
                sumtags,
                sumtagscontrol,
                maxcontrol,
                x,
            ) = determine_sumtags_ctrl_numba(
                grid_nchrom,
                nummarks,
                sumtags,
                gridcontrol_nchrom,
                numcontrolmarks,
                maxcontrol,
                hscontrol_np,
                sumtagscontrol,
                x,
            )
    hscontrol = []
    for i in hscontrol_np:
        hscontrol.append(set(i))

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
def determine_sumtags_numba(
    grid_nchrom,
    nummarks,
    sumtags,
    ntotallocs,
):
    """
    Calculate the sum of tags for each mark in a specific chromosome.

    This function calculates the sum of tags for each mark in a specific chromosome by iterating over the bins and
    marks in the grid.

    :param grid_nchrom: The grid for a specific chromosome.
    :type grid_nchrom: np.ndarray
    :param nummarks: The number of marks.
    :type nummarks: int
    :param sumtags: The array to store the sum of tags for each mark.
    :type sumtags: np.ndarray
    :param ntotallocs: The total number of locations.
    :type ntotallocs: int
    :return: The updated total number of locations and sum of tags for each mark.
    :rtype: tuple[int, np.ndarray]
    """
    for nbin in range(len(grid_nchrom)):
        for nmark in range(nummarks):
            sumtags[nmark] += grid_nchrom[nbin][nmark]
    ntotallocs += len(grid_nchrom)
    return ntotallocs, sumtags


@nb.njit
def determine_sumtags_ctrl_numba(
    grid_nchrom,
    nummarks,
    sumtags,
    gridcontrol_nchrom,
    numcontrolmarks,
    maxcontrol,
    hscontrol,
    sumtagscontrol,
    nchrom_nbin,
):
    """
    Calculate the sum of tags for marks and control marks in a specific chromosome.

    This function calculates the sum of tags for marks and control marks in a specific chromosome by iterating over the
    bins and marks in the grids. It updates the maximum control value, the control set, and the sum of tags for each
    mark and control mark.

    :param grid_nchrom: The grid for a specific chromosome.
    :type grid_nchrom: np.ndarray
    :param nummarks: The number of marks.
    :type nummarks: int
    :param sumtags: The array to store the sum of tags for each mark.
    :type sumtags: np.ndarray
    :param gridcontrol_nchrom: The control grid for a specific chromosome.
    :type gridcontrol_nchrom: np.ndarray
    :param numcontrolmarks: The number of control marks.
    :type numcontrolmarks: int
    :param maxcontrol: The array to store the maximum control value for each mark.
    :type maxcontrol: np.ndarray
    :param hscontrol: The set to store the control values.
    :type hscontrol: list[set]
    :param sumtagscontrol: The array to store the sum of control tags for each mark.
    :type sumtagscontrol: np.ndarray
    :param nchrom_nbin: The count of the bins in the chromosome.
    :type nchrom_nbin: int
    :return: The updated control set, sum of tags for marks, sum of control tags, maximum control values, and count of bins.
    :rtype: tuple[list[set], np.ndarray, np.ndarray, np.ndarray, int]
    """    
    for nbin in range(len(grid_nchrom)):
        nchrom_nbin += 1
        for nmark in range(nummarks):
            # int nval = grid_nchrom_nbin[nmark];
            if numcontrolmarks == 1:
                ncontrolval = gridcontrol_nchrom[nbin][0]
            else:
                ncontrolval = gridcontrol_nchrom[nbin][nmark]
            if ncontrolval > maxcontrol[nmark]:
                maxcontrol[nmark] = ncontrolval
            hscontrol[nmark][nchrom_nbin] = ncontrolval
            sumtags[nmark] += grid_nchrom[nbin][nmark]
            sumtagscontrol[nmark] += ncontrolval
    return hscontrol, sumtags, sumtagscontrol, maxcontrol, nchrom_nbin


@nb.njit
def window_sum_grid_numba(
    gridcontrol_nchrom, sumgridcontrol_chrom_nbin, nbin, nflankwidthcontrol
):
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
    nstart = max(0, nbin - nflankwidthcontrol)
    nend = min(nbin + nflankwidthcontrol, len(gridcontrol_nchrom) - 1)
    for nmark in range(len(sumgridcontrol_chrom_nbin)):
        nsum = 0
        for nrow in range(nstart, nend + 1):
            nval = gridcontrol_nchrom[nrow][nmark]
            if nval > 0:
                nsum += nval
        sumgridcontrol_chrom_nbin[nmark] = nsum
    return gridcontrol_nchrom, sumgridcontrol_chrom_nbin
