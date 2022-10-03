#!/usr/bin/env python

"""
Calculate the fraction of reads given two cutoffs.

"""

import sys
import os

sample_name = sys.argv[1]
csv = sys.argv[2]
outdir = sys.argv[3]
cutoff = sys.argv[4]

count_below, count_above, count_all = 0, 0, 0
cutoff_below, cutoff_above = cutoff.strip().split()
cutoff_below, cutoff_above = int(cutoff_below), int(cutoff_above)

with open(csv, "r") as f:
    for line in f:
        length, c = line.strip().split(",")
        if not length in ["plus", "problem"]:
            length, c = int(length), int(c)
            count_all = count_all + c
            if length < cutoff_below:
                count_below += c
            if length > cutoff_above:
                count_above += c
        elif length == "plus":
            count_all += int(c)
            count_above += int(c)

outfile = os.path.join(outdir, sample_name + "_frac_" + str(cutoff_below) + "_" + str(cutoff_above) + ".csv")
with open(outfile, "w") as f:
    f.write(sample_name + "," + str(count_below/count_all) + "," + str(count_above/count_all) + "\n")
