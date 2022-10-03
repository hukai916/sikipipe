#!/usr/bin/env python

"""
UMI collapse by using mode.

"""

import sys
import os
from collections import Counter
from statistics import mean
import gzip
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO
from collections import Counter

csv         = sys.argv[1]
sample_name = sys.argv[2]
read        = sys.argv[3]
outfile     = sys.argv[4]
cutoff      = sys.argv[5]

cutoff      = int(cutoff)
read_length_mode = {}

encoding = guess_type(read)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

# read in umi and its mode
with open(csv, "r") as f:
    for line in f:
        umi, count, mean, mode = line.split(",")
        count, mean, mode = int(float(count)), int(float(mean)), int(mode)
        if count >= cutoff:
            if not umi in read_length_mode:
                read_length_mode[umi] = mode
# read in fastq:
dict_umi = {} # record that with count >= cutoff
with _open(read) as f:
    for record in SeqIO.parse(f, 'fastq'):
        read_length = len(record.seq)
        umi = record.name.split("_")[-1]
        if umi in read_length_mode:
            if not umi in dict_umi:
                dict_umi[umi] = [record]
            else:
                dict_umi[umi].append(record)

# for each umi (with coutn >= cutoff), find the most common read seq and output it:
records = []
for umi in dict_umi:
    most_common_read = Counter([str(record.seq) for record in dict_umi[umi]]).most_common(1)[0][0]
    for record in dict_umi[umi]:
        if str(record.seq) == most_common_read:
            records.append(record)
            break
SeqIO.write(records, outfile, "fastq")
