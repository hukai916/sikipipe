#!/usr/bin/env python

"""
To visualize the read length distribution by UMI group.

"""

import sys
import os
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import pandas as pd

csv = sys.argv[1]
sample_name = sys.argv[2]
cutoff      = sys.argv[3]
outdir_stat = sys.argv[4]
output_plot = sys.argv[5]
group_num   = sys.argv[6]
bin_number  = sys.argv[7]

# read csv into dict
d_umi = {}
for line in open(csv):
    umi, length = line.strip().split(",")
    umi = umi.split("_")[-1]
    if not umi in d_umi:
        d_umi[umi] = [length]
    else:
        d_umi[umi].append(length)

d_cutoff = {}
for x in d_umi:
    if len(d_umi[x]) == int(cutoff):
        if not len(d_cutoff) > int(group_num):
            d_cutoff[x] = d_umi[x]
        else:
            break

# save to stat
output_stat = os.path.join(outdir_stat, sample_name + ".csv")
os.makedirs(os.path.dirname(output_stat), exist_ok=True)
with open(output_stat, "w") as f:
    for k in d_cutoff:
        for c in d_cutoff[k]:
            f.write(k + "," + str(c) + "\n")

# save to plot
if os.path.getsize(output_stat) == 0: # create empty png for empty csv file
    _output_plot = os.path.join(output_plot, sample_name + "_empty.png")
    os.makedirs(os.path.dirname(_output_plot), exist_ok=True)
    with open(_output_plot, "w") as f:
        pass
else:
    df = pd.read_csv(output_stat, sep = ",", header = None)
    for group_id, group_df in df.groupby(df.iloc[:,0]):
        span_n = group_df.iloc[:, 1][group_df.iloc[:,1].apply(lambda x: str(x).isdigit())].tolist()
        span_n = [int(x) for x in span_n]
        span_s = group_df.iloc[:, 1][group_df.iloc[:,1].apply(lambda x: not str(x).isdigit())].tolist()

        if not len(span_n) == 0:
            span = max(span_n) - min(span_n)
        else:
            span = 0

        if not len(span_n) == 0:
            bins = max(1 + len(span_s), max(span_n) + len(span_s))
        else:
            bins = len(span_s)

        _tem = group_df.iloc[:,1].tolist()
        _tem = [int(x) for x in _tem if str(x).isdigit()] # exclude "plus" and "problem"

        _output_plot_tem = os.path.join(output_plot, sample_name + "_" + group_id + ".txt")
        with open(_output_plot_tem, "w") as f:
            f.write(",".join([str(x) for x in sorted(_tem)]))
            f.write("\n")
            f.write(str(bins))

        if bin_number == "auto":
            bin_n = bins
        else:
            bin_n = int(bin_number)
        # plt.hist(sorted(_tem), bins = bin_n, range = (0, bin_n)) # range controls the width of hist columns
        plt.hist(sorted(_tem))
        # plt.hist(group_df.iloc[:,1], bins = bins)

        # sparse the x-axis ticks:
        # N = 100  # 1 tick every 5
        # myticks = [i for i in range(bin_n + N  + 1) if not i%N]
        # newlabels = [i for i in range(bin_n + N + 1) if not i%N]
        #
        # plt.xticks(myticks, newlabels, rotation = 85)
        # plt.tick_params(labelsize=8)

        plt.xlabel("read length")
        plt.ylabel("count")
        plt.title("Read Length Distribution within UMI Group: " + group_df.iloc[0,0])

        _output_plot = os.path.join(output_plot, sample_name + "_" + group_id + ".png")
        os.makedirs(os.path.dirname(_output_plot), exist_ok=True)
        plt.savefig(_output_plot, dpi = 300)

        plt.close()
