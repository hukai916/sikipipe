#!/usr/bin/env python

"""
To count the repeat length by measuring distance.

Dev notes:
1. Nextflow won't resume if utils changes.

"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os
from collections import Counter
import regex
import gzip
from mimetypes import guess_type
from functools import partial
import numpy as np
import pandas as pd
from utils import plot_repeat_dist

r1 = sys.argv[1]
r2 = sys.argv[2]
sample_name = sys.argv[3]
output_dir  = sys.argv[4]
r1_flanking = sys.argv[5]
r2_flanking = sys.argv[6]
mismatch = sys.argv[7]
bin_number = sys.argv[8]

# r1_flanking = "TCGAGTCCCTCAAGTCCTTC"
# r2_flanking = "CCGCCACCGCCGCCGCCGCC"
# mismatch = 2
encoding = guess_type(r1)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

r1_dna = Seq(r1_flanking)
r2_dna = Seq(r2_flanking)
r1_flanking_rc = str(r1_dna.reverse_complement())
r2_flanking_rc = str(r2_dna.reverse_complement())

dict_repeat = {"problem": 0, "plus": 0}
dict_count  = {}

with _open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        r1_search = regex.search("(" + r1_flanking + ")" + "{s<=" + str(mismatch) + "}", str(record.seq))
        r2_search = regex.search("(" + r2_flanking_rc + ")" + "{s<=" + str(mismatch) + "}", str(record.seq))
        if not r1_search:
            dict_repeat["problem"] += 1 # should not happen
            dict_count[record.name] = "problem"
        elif not r2_search:
            dict_repeat["plus"] += 1
            dict_count[record.name] = "plus"
        else:
            r1_match_start = r1_search.start()
            r1_match_length = len(r1_search.captures()[0])
            r2_match_end = r2_search.end()
            r2_match_length = len(r2_search.captures()[0])
            repeat_length = r2_match_end - r1_match_start - r1_match_length - r2_match_length
            dict_count[record.name] = repeat_length

            if not repeat_length in dict_repeat:
                dict_repeat[repeat_length] = 1
            else:
                dict_repeat[repeat_length] += 1
keys = [x for x in dict_repeat.keys() if not x == "problem" and not x == "plus"]

# output stat:
output_file = os.path.join(output_dir, "stat_r1", sample_name + ".csv")
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as f:
    for i in sorted(keys):
        f.write(str(i) + "," + str(dict_repeat[i]) + "\n")
    f.write("plus" + "," + str(dict_repeat["plus"]) + "\n")
    f.write("problem" + "," + str(dict_repeat["problem"]) + "\n")

# output plot:
N = 5

output_plot = os.path.join(output_dir, "plot_r1", sample_name + ".png")
os.makedirs(os.path.dirname(output_plot), exist_ok=True)
plot_repeat_dist(output_file, output_plot, sample_name, N, bin_number)

# ouput raw count:
output_count = os.path.join(output_dir, "count_r1", sample_name + ".csv")
os.makedirs(os.path.dirname(output_count), exist_ok=True)
with open(output_count, "w") as f:
    for x in dict_count:
        f.write(x + "," + str(dict_count[x]) + "\n")

## Now use R2:
dict_repeat = {"problem": 0, "plus": 0}
dict_count  = {}

with _open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        r2_search = regex.search("(" + r2_flanking + ")" + "{s<=" + str(mismatch) + "}", str(record.seq))
        r1_search = regex.search("(" + r1_flanking_rc + ")" + "{s<=" + str(mismatch) + "}", str(record.seq))
        if not r2_search:
            dict_repeat["problem"] += 1 # should not happen
            dict_count[record.name] = "problem"
        elif not r1_search:
            dict_repeat["plus"] += 1
            dict_count[record.name] = "plus"
        else:
            r2_match_start = r2_search.start()
            r2_match_length = len(r2_search.captures()[0])
            r1_match_end = r1_search.end()
            r1_match_length = len(r1_search.captures()[0])
            repeat_length = r1_match_end - r2_match_start - r2_match_length - r1_match_length
            dict_count[record.name] = repeat_length

            if not repeat_length in dict_repeat:
                dict_repeat[repeat_length] = 1
            else:
                dict_repeat[repeat_length] += 1
keys = [x for x in dict_repeat.keys() if not x == "problem" and not x == "plus"]

# output stat:
output_file = os.path.join(output_dir, "stat_r2", sample_name + ".csv")
os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as f:
    for i in sorted(keys):
        f.write(str(i) + "," + str(dict_repeat[i]) + "\n")
    f.write("plus" + "," + str(dict_repeat["plus"]) + "\n")
    f.write("problem" + "," + str(dict_repeat["problem"]) + "\n")

# output plot:
output_plot = os.path.join(output_dir, "plot_r2", sample_name + ".png")
os.makedirs(os.path.dirname(output_plot), exist_ok=True)
plot_repeat_dist(output_file, output_plot, sample_name, N, bin_number)

# ouput raw count:
output_count = os.path.join(output_dir, "count_r2", sample_name + ".csv")
os.makedirs(os.path.dirname(output_count), exist_ok=True)
with open(output_count, "w") as f:
    for x in dict_count:
        f.write(x + "," + str(dict_count[x]) + "\n")
