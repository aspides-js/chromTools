#!/usr/bin/env python

"""
Description: chromTools main executable.

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).

"""

# ------------------------------------
# Modules
# ------------------------------------

import argparse as ap
import sys

from chromTools import cmd
from chromTools.constants import COMPLETE_VERSION
from chromTools.validate import args_validator


# ------------------------------------
# Main function
# ------------------------------------
def main():
    """
    The Main function/pipeline for chromTools.

    """

    description = (
        "chromTools -- By read count analysis of dataset completeness"  # %(prog)s
    )

    # Parse arguments:
    argparser = ap.ArgumentParser(description=description)

    argparser.add_argument(
        "--version", action="version", version="%(prog)s " + COMPLETE_VERSION
    )

    subparsers = argparser.add_subparsers(dest="subcommand", required=True)

    # command for 'complete'
    add_complete_parser(subparsers)

    args = argparser.parse_args()

    ## Validate arguments
    options = args_validator(args)

    ## Run
    if args.subcommand == "complete":
        cmd.run(options)
        options.info("Complete")


def add_complete_parser(subparsers):
    complete_parser = subparsers.add_parser(
        "complete", help="By read count analysis of dataset completeness"
    )
    complete_parser.add_argument(
        "-f",
        "--files",
        type=str,
        nargs="+",
        required=True,
        help="Input BED files. Should include full path to file",
    )
    complete_parser.add_argument(
        "-c",
        "--control",
        type=str,
        nargs="+",
        default=False,
        dest="control",
        help="Control BED files. Should include full path to file",
    )
    complete_parser.add_argument(
        "-i",
        "--increment",
        type=int,
        default=50_000_000,
        help="Read number to increase increment by. Default is 50_000_000",
    )
    complete_parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        type=str,
        default="",
        help="Directory where output files are written to. If unspecified files will be written to the current working directory",
    )
    complete_parser.add_argument(
        "-g",
        "--genome",
        type=str,
        default="hg38",
        help="Path to chromosome and chromosome length file. Genome assemblies hg18, hg19, hg38, mm9, mm10, rn5, rn6, danRer7, danRer10, dm3, dm6, ce6, and ce10 can be accessed with their genome assembly name",
    )
    complete_parser.add_argument(
        "--gsize",
        default=False,
        help="Required integer if specifying own genome chromosome length file. Default is FALSE",
    )
    complete_parser.add_argument(
        "-s",
        "--seed",
        type=int,
        help="User can set random seen. If unspecified this will be generated randomly.",
    )
    complete_parser.add_argument(
        "--paired",
        dest="paired",
        default=False,
        action="store_true",
        help="Specify whether the input data is paired or unpaired. Default is FALSE",
    )
    complete_parser.add_argument(
        "--force-overwrite",
        dest="force",
        default=False,
        action="store_true",
        help="If TRUE, files and directories in OUTDIR will be overwritten. Default is FALSE",
    )
    return


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted\n")
