# this script merges the results from featureCounts in to a single gene X cell counts matrix

"""
Merge count tables for genes for multiple samples into a gene X single cell/sample counts matrix

Input files as a list of comma-separated file paths

Default output is stdout, redirect using the --output flag
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os
import re

# set up the command line parser
parser = argparse.ArgumentParser(description='Merge gene count matrices together')
parser.add_argument("--output", dest="output", type="string",
                    default=sys.stdout,
                    help="Output file path")

args = parser.parse_args()
input_files = sys.argv[-1].split(",")

# this WILL take a long time for many thousands of files
# sequentially open and merge files based on the following columns:
# Geneid, Chr, Start, End, Strand, Length
# parse the file name input filename

n_files = 0
for cfile in input_files:
    fname = cfile.split("/")[-1]
    samp_name = fname.split(".")[:-2]
    if n_files == 0:
        count_df = pd.read_table(cfile, sep="\t", header=0, index_col=None, 
                                 comment='#')
        columns = list(count_df.columns[:-1])
        columns.append(samp_name[0])
        count_df.columns = columns
        n_files += 1
    else:
        _df = pd.read_table(cfile, sep="\t", header=0, index_col=None, 
                                 comment='#')
        columns = list(_df.columns[:-1])
        columns.append(samp_name[0])
        _df.columns = columns
        
        # merge
        count_df = pd.merge(count_df, _df,
                            left_on=["Geneid", "Chr", "Start", "End", "Strand", "Length"],
                            right_on=["Geneid", "Chr", "Start", "End", "Strand", "Length"])

count_df.to_csv(args.output, sep="\t",
                index=None)
