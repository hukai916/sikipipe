#!/usr/bin/env python

"""
Preprocessing data:
1. if first 30 nt > 10 distance away from expected universal nt: send to problematic set.
2. if first 30 nt closer to 3p-univ: TGGATCACTTGTGCAAGCATCACATCGTAG, reverse complement it.
3. if read length outside of 1000 - 5000: abnormal


Usage:
python prep_data.py reads.fastq(.gz) outfile_normal.fastq.gz outfile_abnormal.fastq.gz outfile_stat.csv

References:
1. parse fastq + gz: https://stackoverflow.com/questions/42757283/seqio-parse-on-a-fasta-gz
"""

import sys
import gzip
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO
from Bio.Seq import Seq
from utils import hamming_distance
import os
import numpy as np

filename        = sys.argv[1]
sample_name     = sys.argv[2]
outfile_normal  = sys.argv[3]
outfile_abnormal = sys.argv[4]
outfile_stat     = sys.argv[5]
univ_5p          = sys.argv[6]
univ_3p          = sys.argv[7]

oligo = {}
oligo["univ_3p"] = univ_3p
oligo["univ_5p"] = univ_5p
oligo["univ_3p_rc"] = str(Seq(univ_3p).reverse_complement())
oligo["univ_5p_rc"] = str(Seq(univ_5p).reverse_complement())

encoding = guess_type(filename)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open # must be rt for SeqIO

_path, _basename = os.path.dirname(filename), os.path.basename(filename)
_basename = _basename + ".gz" if not _basename.endswith(".gz") else _basename
# _basename = _basename.replace(".fastq", ".fasta")
f_out = gzip.open(outfile_normal, 'wb')
f_out_abnormal = gzip.open(outfile_abnormal, 'wb')

count_normal, count_abnormal = 0, 0
with _open(filename) as f:
    for record in SeqIO.parse(f, 'fastq'):
        # print(record.name)
        # exit()
        direction = "forward" # if read direction same as primer, do nothing
        dist = [hamming_distance(oligo["univ_5p"], record.seq[:30]), hamming_distance(oligo["univ_3p"], record.seq[:30])]
        dist_rc = [hamming_distance(oligo["univ_5p_rc"], record.seq[-30:]), hamming_distance(oligo["univ_3p_rc"], record.seq[-30:])]
        direction = "forward" if np.argmin(dist) == 0 else "reverse"
        if min(dist) > 10 or min(dist_rc) > 10 or len(record.seq) < 1000 or len(record.seq) > 5000:
            """
            If either universal sequence is missing, treat the read as abnormal.
            """
            count_abnormal += 1
            f_out_abnormal.write(bytes('@' + record.id + '\n', 'ascii'))
            f_out_abnormal.write(bytes(str(record.seq) + '\n', 'ascii')) # gzip require byte format
            qual = ''.join([chr(x + 33) for x in record.letter_annotations['phred_quality']])
            f_out_abnormal.write(bytes('+\n' + qual + '\n', 'ascii'))
        else:
            if direction == "forward":
                count_normal += 1
                f_out.write(bytes('@' + record.id + '\n', 'ascii'))
                f_out.write(bytes(str(record.seq) + '\n', 'ascii'))
                qual = ''.join([chr(x + 33) for x in record.letter_annotations['phred_quality']])
                f_out.write(bytes('+\n' + qual + '\n', 'ascii'))
            elif direction == "reverse":
                count_normal += 1
                f_out.write(bytes('@' + record.id + '\n', 'ascii'))
                f_out.write(bytes(str(record.reverse_complement().seq) + '\n', 'ascii'))
                qual = ''.join([chr(x + 33) for x in record.letter_annotations['phred_quality'][::-1]])
                f_out.write(bytes('+\n' + qual + '\n', 'ascii'))

f_out.close()
f_out_abnormal.close()

with open(outfile_stat, "w") as f:
    f.write(sample_name + "," + str(count_normal) + "," + str(count_abnormal) + "\n")
