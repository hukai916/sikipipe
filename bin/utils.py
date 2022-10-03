import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import os
import pandas as pd

oligo = {
    "univ_3p":    "TGGATCACTTGTGCAAGCATCACATCGTAG",
    "univ_3p_rc": "CTACGATGTGATGCTTGCACAAGTGATCCA",
    "univ_5p" :   "GCAGTCGAACATGTAGCTGACTCAGGTCAC",
    "univ_5p_rc": "GTGACCTGAGTCAGCTACATGTTCGACTGC"
}

def hamming_distance(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def plot_repeat_dist(csv, output_file, sample_name, N, bin_number = 250):
    """
    csv:
        col1 col2
        repeat_length   count
    N = 5  # 1 tick every 5
    """
    # output plot:
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # read in stat csv: ensure 'plus' and 'problem' are there
    df = pd.read_csv(csv, sep = ",", header = None)
    x, y = [], []
    _dict = {"plus": 0, "problem": 0} # ensure that 'plus' and 'problem' are there
    for i in range(len(df.iloc[:,0])):
        _dict[str(df.iloc[i, 0])] = df.iloc[i, 1]

    # determine max count:
    max_count = max([int(i) for i in df.iloc[:,0] if str(i).isdigit()])
    # plus, problem = 0, 0
    # for i, v in enumerate(df.iloc[:, 0]):
    #     if v == "plus":
    #         plus = df.iloc[i, 1]
    #     elif v == "problem":
    #         problem = df.iloc[i, 1]

    # determine bin_n considering 'auto' option:
    if bin_number == 'auto': # if 'auto', use max_count as bin_number
        bin_n = max_count # include "plus" and "problem"
    else:
        bin_n = int(bin_number)
        assert max_count <= bin_n, "X-axis scale too small!"

    # determine x and weight (y), add N number of gap between count and "plus" and "problem" for better visualization
    x, y = [], []
    for i, v in enumerate(list(range(bin_n + N)) + ["plus"] + [0] * N +  ["problem"]):
        x = x + [str(i)]
        if not v in ["plus", "problem"]:
            if not str(i) in _dict:
                y = y + [0]
            else:
                y = y + [_dict[str(i)]]
        elif v == "plus":
            y = y + [_dict["plus"]]
        elif v == "problem":
            y = y + [_dict["problem"]]

    # sparse the ticks:
    myticks = [i for i in range(bin_n) if not i%N]
    mylabels = [i for i in range(bin_n) if not i%N]
    if not "plus" in mylabels:
        mylabels = mylabels + ["plus"]
        myticks = myticks + [bin_n + N]
    if not "problem" in mylabels:
        mylabels = mylabels + ["problem"]
        myticks = myticks + [bin_n + 2* N]

    # plot:
    n, bins, patches = plt.hist(x, weights = y, bins = bin_n + 2 * N, range = (0, bin_n + 2 * N))
    plt.xticks(myticks, mylabels, rotation = 85)
    plt.tick_params(labelsize = 8)

    plt.savefig(output_file, dpi = 300, bbox_inches = "tight")
    plt.close()

def print_test():
    return("inside here")
