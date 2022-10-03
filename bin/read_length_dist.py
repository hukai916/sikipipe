#!/usr/bin/env python

"""
To count the read length.

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

r = sys.argv[1]
sample_name = sys.argv[2]
output_dir  = sys.argv[3]

encoding = guess_type(r)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

dict_count  = {}

with _open(r) as f:
    for record in SeqIO.parse(f, 'fastq'):
        read_length = len(record.seq)
        dict_count[record.name] = read_length

# ouput raw count:
output_count = os.path.join(output_dir, "length", sample_name + ".csv")
os.makedirs(os.path.dirname(output_count), exist_ok=True)
with open(output_count, "w") as f:
    for x in dict_count:
        f.write(x + "," + str(dict_count[x]) + "\n")
