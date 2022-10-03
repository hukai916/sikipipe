#!/usr/bin/env python

"""
Plot UMI pattern: "reads per UMI" vs "log2(count)"
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

csv = sys.argv[1]
sample_name = sys.argv[2]
outfile = sys.argv[3]
bin_number = sys.argv[4]

data=pd.read_csv(csv, sep=',')
x = data.iloc[:,0].to_frame()
weight = np.log2(data[["count"]])

# n, bins, patches = axs[_i][_j].hist(x, weights = weight, bins=len(x), edgecolor='black')
if bin_number == "auto":
    bin_n = int(x.iloc[-1])
else:
    bin_n = int(bin_number)
n, bins, patches = plt.hist(x, weights = weight, bins = bin_n, range = (0, bin_n))

plt.title(sample_name)
plt.xlabel("reads per UMI")
plt.ylabel("log2(count)")
plt.savefig(outfile, dpi = 300)
plt.close()
