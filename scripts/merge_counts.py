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
import logging

# setup the logger
logging.basicConfig(level=logging.INFO,
                        filename="/nfs/research1/marioni/mdmorgan/Thymus/logs/Merge_quant.log")

# set up the command line parser
parser = argparse.ArgumentParser(description='Merge gene count matrices together')
parser.add_argument("--output", dest="output", type=str,
                    default=sys.stdout,
                    help="Output file path")

parser.add_argument("--input-directory", dest="input_dir", type=str,
                    help="The directory containing the input files")

parser.add_argument("--file-regex", dest="file_regex", type=str,
                    help="Regular expression to match to input files. Accepts "
                    "normal regular expression syntax")

args = parser.parse_args()

# can't take thousands of input files because posix has a limiti to the
# character limit of an argument to a script
# need to pass an input directory and a glob/regex instead
# could be a problem for snakemake to handle?
reg_compile = re.compile(args.file_regex)
found_files = [ft for ft in os.listdir(args.input_dir) if re.search(reg_compile, ft)]
input_files = [os.path.join(args.input_dir, fx) for fx in found_files]
#input_files = sys.argv[-1].split(",")

logging.info("Found {} count files to merge".format(len(input_files)))

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
        n_files += 1
        # merge
        count_df = pd.merge(count_df, _df,
                            left_on=["Geneid", "Chr", "Start", "End", "Strand", "Length"],
                            right_on=["Geneid", "Chr", "Start", "End", "Strand", "Length"])
        if not n_files % 100:
            logging.info("Processed {} count files".format(n_files))

logging.info("Writing merged gene counts to file: {}".format(args.output))
count_df.to_csv(args.output, sep="\t",
                index=None)
