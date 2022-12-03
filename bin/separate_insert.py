#!/usr/bin/env python

"""
To separate reads into insert-containing and non-insert-contaning reads.

Usage:
    separate_insert.py xxx.bam insert_fasta out_insert.bam out_non_insert.bam 15
"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import pysam
import regex

bam             = sys.argv[1]
insert_fasta    = sys.argv[2]
out_bam_insert  = sys.argv[3]
out_bam_non_insert = sys.argv[4]
insert_frac_size   = int(sys.argv[5])
mismatch = int(int(insert_frac_size) * 0.1)

# break insert into tiles
insert_list = []
for record in SeqIO.parse(insert_fasta, "fasta"):
    ref = str(record.seq)
    insert_list = [ref[i:i+insert_frac_size] for i in range(len(ref) - insert_frac_size)]
    break

samfile = pysam.AlignmentFile(bam, "rb")
outfile_insert = pysam.AlignmentFile(out_bam_insert, "wb", template = samfile)
outfile_non_insert = pysam.AlignmentFile(out_bam_non_insert, "wb", template = samfile)

for alignment in samfile:
    query = alignment.query_sequence
    contain_insert = 0
    for frac in insert_list:
        if regex.search("(?i)(" + frac + ")" + "{e<=" + str(mismatch) + "}", query):
            contain_insert = 1
            break
    if contain_insert:
        outfile_insert.write(alignment)
    else:
        outfile_non_insert.write(alignment)
