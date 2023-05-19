#!/usr/bin/python


# ------------------------------------
# python modules
# ------------------------------------

import os
import numpy as np
import math

# ------------------------------------
# Main function
# ------------------------------------

def make_binary_data_from_bed( n, options ):
    """_summary_

    Args:
        options (_type_): _description_

    Raises:
        ValueError: _description_
        ValueError: _description_
        ValueError: _description_

    Returns:
        _type_: _description_
    """ 
    ## check all options back to options
    szchromlengthfile = options.szchromlengthfile
    szoutputbinarydir = options.szoutputbinarydir
    szmarkdir = options.szmarkdir
    szcontroldir = options.szcontroldir

    nshift = options.nshift
    bcenterinterval = options.bcenterinterval
    noffsetleft = options.noffsetleft
    noffsetright = options.noffsetright
    npseudocountcontrol = options.npseudocountcontrol

    szcolfields = options.szcolfields

    dpoissonthresh = options.dpoissonthresh
    dfoldthresh = options.dfoldthresh
    bcontainsthresh = options.bcontainsthresh
    dcountthresh = options.dcountthresh
    nflankwidthcontrol = options.nflankwidthcontrol

    nbinsize = 200
    bcontrol = False
    szcontrolfile = ''

    control = False
    szcell = f'{n}'
    szmark = 'mark'
    hscells = set()
    hsmarks = set()
    hscells.add(szcell)
    hsmarks.add(szmark)
    hmfiles = {}    
    hmfilescontrol = {}
    hmfiles[f'{szcell}\t{szmark}'] = [f'downsampled.{n}.bed']
    if control:
        hscellcontrol = hscells
        hmfilescontrol[f'{szcell}\t{szmark}'] = ['downsampled.ctrl.bed']
        bcontrol = True
    else:
        hscellnocontrol = hscells
        hmfilescellcontrol = {}

    ## reads in the chromosome length information file
    # the first column of this file is the chromosome and the second is the chromsome length
    with open(szchromlengthfile) as brchrom:
        allines = brchrom.readlines()

    chroms = []  # stores the chromosome name
    lengths = []  # stores the chromosome length
    hmchrom = {}  # stores a mapping of chromosome to chromosome index
    for szLine in allines:
        st = szLine.strip().split()
        if len(st) < 2:
            raise ValueError("Invalid line found in " + szchromlengthfile)
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
        grid[ni] = np.zeros((lengths[ni]//nbinsize, nummarks), dtype=int)

    gridcontrol = None
    sumgridcontrol = None
    bpresentcontrol = None
    bpresentmarkscontrol = None

    if bcontrol:
        gridcontrol = np.empty((len(chroms),), dtype=np.ndarray)
        sumgridcontrol = np.empty((len(chroms),), dtype=np.ndarray)
        bpresentcontrol = [False for _ in range(len(chroms))]
        bpresentmarkscontrol = [False for _ in range(nummarks)]

    bpresentmarks = [False for _ in range(nummarks)]
    #----------------------------------------


    for szcell in hscells:
        # added in v1.18 to consistently reset chroms considered
        bpresent = [False for _ in range(len(chroms))]
        # going through each declared cell type
        hscellcontrol = hmfilescellcontrol.get(szcell)
        bmissing = szcell in hscellnocontrol
        if hscellcontrol is None:
            # no control data for this cell type
            bcontrolfile = False
        else:
            bcontrolfile = True
            if len(hscellcontrol) == 1 and not bmissing:
                # update in v1.14
                # we have one control for all marks
                numcontrolmarks = 1
            else:
                # will allocate the full memory for all marks
                numcontrolmarks = nummarks

        # loading data for the cell type
        print("loadgrid")
        grid = load_grid(grid, bpresent, bpresentmarks, marks, nshift, nbinsize, bcenterinterval, noffsetleft,
                noffsetright, hmfiles, szcell, szmarkdir, hmchrom, 0, szcolfields, False, False)

        if bcontrolfile:
            if gridcontrol[0] is None or len(gridcontrol[0][0]) != numcontrolmarks:
                # reallocate if changing array size
                # allowed to go between single and matched
                for ni in range(len(chroms)):
                    gridcontrol[ni] = np.ones((lengths[ni]//nbinsize, numcontrolmarks), dtype=int)
                    sumgridcontrol[ni] = np.zeros((lengths[ni]//nbinsize, numcontrolmarks), dtype=int)

            # we have control data loading cell type data for that
            gridcontrol = load_grid(gridcontrol, bpresentcontrol, bpresentmarkscontrol, marks, nshift, nbinsize, bcenterinterval, noffsetleft, noffsetright,
                    hmfilescontrol, szcell, szcontroldir, hmchrom, npseudocountcontrol, szcolfields, False, True)



    nummarks_m1 = nummarks - 1


    if bcontrolfile:
        # binarization will be based on control data

        # smoothing control data
        window_sum_grid(gridcontrol, sumgridcontrol, nflankwidthcontrol)

        # determiming thresholds for each mark and background depth
        thresholds = determine_mark_thresholds_from_binned_data_array_against_control(
                grid, sumgridcontrol, bpresent, bpresentcontrol, dpoissonthresh, dfoldthresh, bcontainsthresh, dcountthresh
            )
        for nchrom in range(len(chroms)):
            if bpresent[nchrom] and bpresentcontrol[nchrom]:
                szfile = os.path.join(szoutputbinarydir, szcell + "_" + chroms[nchrom] + "_binary.txt")
                print("Writing to file", szfile)
                with open(szfile, 'w') as pw:
                    # we have both primary and control data for the mark
                    pw.write(szcell + "\t" + chroms[nchrom] + "\n")
                    for nmark in range(nummarks_m1):
                        pw.write(marks[nmark] + "\t")
                    pw.write(marks[nummarks_m1] + "\n")

                    grid_nchrom = grid[nchrom]
                    sumgridcontrol_nchrom = sumgridcontrol[nchrom]
                    for nbin in range(len(grid_nchrom)):
                        grid_nchrom_nbin = grid_nchrom[nbin]
                        sumgrid_nchrom_nbin = sumgridcontrol_nchrom[nbin]

                        for nmark in range(nummarks_m1):
                            if numcontrolmarks == 1:
                                ncontrolval = sumgrid_nchrom_nbin[0]
                            else:
                                ncontrolval = sumgrid_nchrom_nbin[nmark]

                            # printing one if count exceeds background threshold
                            if not bpresentmarks[nmark]:
                                pw.write("2\t")
                            elif (thresholds[nmark][ncontrolval] <= grid_nchrom_nbin[nmark]):
                                pw.write("1\t")
                            else:
                                pw.write("0\t")

                        if numcontrolmarks == 1:
                            ncontrolval = sumgrid_nchrom_nbin[0]
                        else:
                            ncontrolval = sumgrid_nchrom_nbin[nummarks_m1]

                        if not bpresentmarks[nummarks_m1]:
                            pw.write("2\n")
                        elif (thresholds[nummarks_m1][ncontrolval] <= grid_nchrom_nbin[nummarks_m1]):
                            pw.write("1\n")
                        else:
                            pw.write("0\n")

    else: ## if no control file
        thresholds = determine_mark_thresholds_from_binned_data_array(grid, bpresent, dpoissonthresh, dfoldthresh, bcontainsthresh, dcountthresh)
        count=0
        total=0
        for nchrom in range(len(chroms)):
            if bpresent[nchrom]:
                szfile = os.path.join(szoutputbinarydir, szcell + "_" + chroms[nchrom] + "_binary.txt")
                print("Writing to file " + szfile)
                with open(szfile, 'w') as pw:
                    pw.write(szcell + "\t" + chroms[nchrom] + "\n")
                    for nmark in range(len(marks)-1):
                        pw.write(str(marks[nmark]) + "\t")
                    pw.write(str(marks[nummarks_m1]) + "\n")
                    grid_nchrom = grid[nchrom]
                    for nbin in range(len(grid_nchrom)):
                        grid_nchrom_nbin = grid_nchrom[nbin]
                        for nmark in range(nummarks_m1):
                            if not bpresentmarks[nmark]:
                                pw.write("2\t")
                            elif (thresholds[nmark] <= grid_nchrom_nbin[nmark]):
                                pw.write("1\t")
                            else:
                                pw.write("0\t")
                        if not bpresentmarks[nummarks_m1]:
                            pw.write("2\n")
                        elif (thresholds[nummarks_m1] <= grid_nchrom_nbin[nummarks_m1]):
                            pw.write("1\n")
                            count+=1
                            total+=1
                        else:
                            pw.write("0\n")
                            total+=1
    return count, total

#--------------------------------------

def load_grid(grid, bpresent, bpresentmarks, marks, nshift, nbinsize, bcenterinterval, noffsetleft, noffsetright, hmfiles, szcell, szmarkdir, hmchrom, ninitval, szcolfields, bpeaks, bcontrol):
    """     * Converts read level data for a cell into integer count information
     * grid is the array for which the read counts will be loaded into"""
    nummarks = len(grid[0][0])
    # columns
    nchromcol = 0
    nbegincol = 1
    nendcol = nbegincol+1
    nstrandcol = -1
    nmaxindex = -1
    # going through all the mark files in each cell type
    for nmark in range(nummarks):
        alfiles = hmfiles.get(f"{szcell}\t{marks[nmark]}")

        if (nummarks == 1) and bcontrol:
            # this code was added in version 1.04 to handle the situation in which there is a missing mark
            # in the first position, but control data can be listed elsewhere
            # trying to find a listing where control data is available.
            nfilemark = 1
            while (alfiles is None) and (nfilemark < len(marks)):
                # only one control for the cell looking for a valid one
                alfiles = hmfiles.get(f"{szcell}\t{marks[nfilemark]}")
                nfilemark += 1

        if alfiles is None:
            if bcontrol:
                print(f"Warning did not find control data for {szcell} {marks[nmark]} treating as missing")
            else:
                print(f"Warning did not find data for {szcell} {marks[nmark]} treating as missing")

            bpresentmarks[nmark] = False
            if not bcontrol:
                # slight efficiency improvement here in v1.04
                for nchrom in range(len(grid)):
                    grid_nchrom = grid[nchrom]
                    numbins = len(grid_nchrom)
                    for nbin in range(numbins):
                        grid_nchrom[nbin][nmark] = -1
        else:
            bpresentmarks[nmark] = True
            print(len(alfiles))
            for nfile in range(len(alfiles)):
                szfile = alfiles[nfile]
                with open(os.path.join(szmarkdir, szfile), 'r') as brbed:
                    for szLine in brbed:
                        szLineA = szLine.split()
                        szchrom = szLineA[nchromcol]
                        objInt = hmchrom.get(szchrom)

                        if nmaxindex >= len(szLineA):
                            raise ValueError("Column index "+nmaxindex+" exceeds maximum index "+str(len(szLineA)-1)+" indices are 0-based")

                        if objInt is not None:
                            nchrom = objInt
                            szstrand = szLineA[nstrandcol]
                            if szstrand == "+":
                                nbin = (int(szLineA[nbegincol])-noffsetleft+nshift)//nbinsize
                            elif szstrand == "-":
                                nbin = (int(szLineA[nendcol])-noffsetright-nshift)//nbinsize
                            else:
                                raise ValueError(szstrand+" is an invalid strand!")
                            if nbin>=0 and nbin<len(grid[nchrom]):
                                grid[nchrom][nbin][nmark] += 1
                                bpresent[nchrom] = True
    return grid


def determine_mark_thresholds_from_binned_data_array(grid, bpresent, dpoissonthresh, dfoldthresh, bcontainsthresh, dcountthresh):
    dcumthreshold = 1 - dpoissonthresh
    nummarks = len(grid[0][0])
    ntotallocs = 0
    sumtags = [0] * nummarks
    thresholds = [0] * nummarks

    for nchrom in range(len(grid)):
        if bpresent[nchrom]:
            grid_nchrom = grid[nchrom]
            for nbin in range(len(grid_nchrom)):
                grid_nchrom_nbin = grid[nchrom][nbin]
                for nmark in range(nummarks):
                    sumtags[nmark] += grid_nchrom_nbin[nmark]
            ntotallocs += len(grid_nchrom)

    for nj in range(len(sumtags)):
        dlambda = sumtags[nj] / ntotallocs
        dcum = 0
        nthresh = 0
        dlogfactorial = 0
        
        while dcum <= dcumthreshold:
            dprob = math.exp(math.log(dlambda)*nthresh - dlambda - dlogfactorial)
            dcum += dprob
            nthresh += 1
            dlogfactorial += math.log(nthresh)

        if bcontainsthresh:
            nthresh -= 1
        
        thresholds[nj] = max(int(math.ceil(dfoldthresh*dlambda)), nthresh, 1, int(math.ceil(dcountthresh)))
    print(thresholds)    
    return thresholds


def determine_mark_thresholds_from_binned_data_array_against_control(grid, gridcontrol, bpresent, bpresentcontrol, dpoissonthresh, dfoldthresh, bcontainsthresh, dcountthresh):
    dcumthreshold = 1 - dpoissonthresh

    nummarks = len(grid[0][0])
    numcontrolmarks = len(gridcontrol[0][0])

    # stores the the total number of reads for each mark and its matched control
    sumtags = [0] * nummarks
    sumtagscontrol = [0] * nummarks

    # stores the thresholds for each mark and background value
    thresholds = [[] for _ in range(nummarks)]

    # stores the maximum control value found for each mark
    maxcontrol = [0] * nummarks

    # stores which background control values have been found for each mark
    hscontrol = [set() for _ in range(nummarks)]

    for nchrom in range(len(grid)):
        if bpresent[nchrom] and bpresentcontrol[nchrom]:
            grid_nchrom = grid[nchrom]
            gridcontrol_nchrom = gridcontrol[nchrom]

            for nbin in range(len(grid_nchrom)):
                grid_nchrom_nbin = grid_nchrom[nbin]
                gridcontrol_nchrom_nbin = gridcontrol_nchrom[nbin]
                for nmark in range(nummarks):
                    # int nval = grid_nchrom_nbin[nmark];
                    if numcontrolmarks == 1:
                        ncontrolval = gridcontrol_nchrom_nbin[0]
                    else:
                        ncontrolval = gridcontrol_nchrom_nbin[nmark]

                    if ncontrolval > maxcontrol[nmark]:
                        maxcontrol[nmark] = ncontrolval
                    hscontrol[nmark].add(ncontrolval)
                    sumtags[nmark] += grid_nchrom_nbin[nmark]
                    sumtagscontrol[nmark] += ncontrolval

    for nmark in range(len(sumtags)):
        # computing threshold for each mark
        thresholds[nmark] = [0] * (maxcontrol[nmark] + 1)
        thresholds_nmark = thresholds[nmark]

        # determine the relative enrichment for real reads versus the local expected
        davgratio = sumtags[nmark] / sumtagscontrol[nmark]

        # sets a background of 0 threshold to 1
        thresholds_nmark[0] = max(int(dcountthresh), 1)

        # going through each background value
        for nbackground in range(1, maxcontrol[nmark] + 1):
            # bug fixed in 1.14 that changes less than to less than equal
            if nbackground in hscontrol[nmark]:
                # only compute the background threshold for values we observed

                # expected number of reads is the local background number of reads times the global
                # read depth enrichment for sumtags
                dlambda = davgratio * nbackground

                dcum = 0
                nthresh = 0
                dlogfactorial = 0

                while dcum <= dcumthreshold:
                    dprob = math.exp(math.log(dlambda) * nthresh - dlambda - dlogfactorial)
                    dcum += dprob
                    nthresh += 1
                    dlogfactorial += math.log(nthresh)

                if bcontainsthresh:
                    # decreasing to include the dpoissonthreshold probability
                    nthresh -= 1

                thresholds_nmark[nbackground] = max(int(dfoldthresh * dlambda), nthresh, int(dcountthresh))
    return thresholds




def window_sum_grid(gridcontrol, sumgridcontrol, nflankwidthcontrol):
    # iterates over chromosome, bin position, and mark
    for nchrom in range(len(gridcontrol)):
        gridcontrol_nchrom = gridcontrol[nchrom]
        sumgridcontrol_nchrom = sumgridcontrol[nchrom]

        for nbin in range(len(sumgridcontrol_nchrom)):
            sumgridcontrol_nchrom_nbin = sumgridcontrol_nchrom[nbin]
            nstart = max(0, nbin - nflankwidthcontrol)
            nend = min(nbin + nflankwidthcontrol, len(gridcontrol_nchrom) - 1)

            for nmark in range(len(sumgridcontrol_nchrom_nbin)):
                nsum = 0
                for nrow in range(nstart, nend+1):
                    nval = gridcontrol_nchrom[nrow][nmark]
                    if nval > 0:
                        nsum += nval

                sumgridcontrol_nchrom[nbin][nmark] = nsum