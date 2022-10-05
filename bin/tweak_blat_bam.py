#!/usr/bin/env python

"""
Tweak BLAT bam by adding seq and qual:

Version1: awk solution, can't reverse complement seq and qual
# add seq and qual:
zcat < $reads | awk -v ORS="" '{if (NR%4 == 0) {print \$0"\\n"} else {print \$0"\t"} }' > ${outdir}/blat_bam_tuned/${prefix}.fastq.tsv

samtools view -h ${outdir}/blat_bam_tuned/${prefix}.soft.bam | awk 'BEGIN{FS=OFS="\t"} NR==FNR { seq_dict[substr(\$1, 2)] = \$2; qual_dict[substr(\$1, 2)] = \$4; next } { if (FNR < 4) {print \$0} else { if (\$1 in seq_dict) { \$10 = seq_dict[\$1]; \$11 = qual_dict[\$1]; print \$0 } } } ' ${outdir}/blat_bam_tuned/${prefix}.fastq.tsv - | samtools view -o ${outdir}/blat_bam_tuned/${prefix}.tuned.bam
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
import pysam

bam     = sys.argv[1]
fastq   = sys.argv[2]
outfile = sys.argv[3]

encoding = guess_type(fastq)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

samfile = pysam.AlignmentFile(bam, "rb")
outfile = pysam.AlignmentFile(outfile, "w", template = samfile)

encoding = guess_type(fastq)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

# read fastq into dictionary:
dict_fastq = {} # keep track of read name, and read seq in rc, and qual in reverse order
with _open(fastq) as f:
    for record in SeqIO.parse(f, 'fastq'):
        tem = record.format("fastq").split("\n")
        id, qual = str(tem[0])[1:], tem[3][::-1]
        seq_rc = str(record.reverse_complement().seq)
        dict_fastq[id] = [seq_rc, qual]

for read in samfile:
    if (read.is_reverse):
        read.query_sequence = dict_fastq[read.qname][0]
        read.qual = dict_fastq[read.qname][1]
        outfile.write(read)
    else:
        outfile.write(read)

#
# # if seq and qual already in bam:
# for read in samfile:
#     if (read.is_reverse):
#         outfile1.write(read)
#         query_sequence = read.query_sequence
#         qual = read.qual
#         read.query_sequence = read.get_forward_sequence()
#         read.qual = qual[::-1] # also reverse qualities.
#         outfile.write(read)
#     else:
#         outfile.write(read)
