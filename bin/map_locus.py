#!/usr/bin/env python

"""
To annotate reads into htt, misprimed, and problem reads:
    if both ends match: htt_locus;
    if no end matches: problem_reads
Usage:
    python map_locus.py r1.fastq r2.fastq dir_to_store_htt dir_to_store_misprimed dir_to_store_problem output_stat_file r1_flanking r2_flanking
"""

from Bio import SeqIO
import sys
import os
from collections import Counter
import gzip
from mimetypes import guess_type
from functools import partial
import regex

r1 = sys.argv[1]
r2 = sys.argv[2]
htt_locus_dir = sys.argv[3]
misprimed_locus_dir = sys.argv[4]
problem_reads_dir = sys.argv[5]
output_stat_file  = sys.argv[6]
r1_flanking = sys.argv[7]
r2_flanking = sys.argv[8]
error       = sys.argv[9]

r1_match = {}
r2_match = {}
r_match  = {}

encoding = guess_type(r1)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

with _open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if regex.search("(" + r1_flanking + ")" + "{e<=" + str(error) + "}", str(record.seq[:len(r1_flanking)])):
            r1_match[record.name] = 1
        else:
            r1_match[record.name] = 0
with _open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if regex.search("(" + r2_flanking + ")" + "{e<=" + str(error) + "}", str(record.seq[:len(r2_flanking)])):
            r2_match[record.name] = 1
        else:
            r2_match[record.name] = 0

for k in r1_match:
    assert k in r2_match, "not equal records!"
    r_match[k] = r1_match[k] + r2_match[k]

# output:
out_htt_r1 = os.path.join(htt_locus_dir, os.path.basename(r1).split(".gz")[0])
out_htt_r2 = os.path.join(htt_locus_dir, os.path.basename(r2).split(".gz")[0])
out_mis_r1 = os.path.join(misprimed_locus_dir, os.path.basename(r1).split(".gz")[0])
out_mis_r2 = os.path.join(misprimed_locus_dir, os.path.basename(r2).split(".gz")[0])
out_problem_r1 = os.path.join(problem_reads_dir, os.path.basename(r1).split(".gz")[0])
out_problem_r2 = os.path.join(problem_reads_dir, os.path.basename(r2).split(".gz")[0])

os.makedirs(os.path.dirname(out_htt_r1), exist_ok=True)
os.makedirs(os.path.dirname(out_htt_r2), exist_ok=True)
os.makedirs(os.path.dirname(out_mis_r1), exist_ok=True)
os.makedirs(os.path.dirname(out_mis_r2), exist_ok=True)
os.makedirs(os.path.dirname(out_problem_r1), exist_ok=True)
os.makedirs(os.path.dirname(out_problem_r2), exist_ok=True)

r1_htt = []
r1_mis = []
r1_problem = []
with _open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if r_match[record.name] == 2:
            r1_htt.append(record)
        elif r_match[record.name] == 1:
            r1_mis.append(record)
        elif r_match[record.name] == 0:
            r1_problem.append(record)

SeqIO.write(r1_htt, out_htt_r1, "fastq")
SeqIO.write(r1_mis, out_mis_r1, "fastq")
SeqIO.write(r1_problem, out_problem_r1, "fastq")

r2_htt = []
r2_mis = []
r2_problem = []
with _open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if r_match[record.name] == 2:
            r2_htt.append(record)
        elif r_match[record.name] == 1:
            r2_mis.append(record)
        elif r_match[record.name] == 0:
            r2_problem.append(record)
SeqIO.write(r2_htt, out_htt_r2, "fastq")
SeqIO.write(r2_mis, out_mis_r2, "fastq")
SeqIO.write(r2_problem, out_problem_r2, "fastq")

# print some stats:
res = Counter(r_match.values())
with open(output_stat_file, "w") as f:
    name = os.path.basename(r1).split("_1.fastq.gz")[0]
    p2 = str(res[2]/(sum([res[2], res[1], res[0]])))
    p1 = str(res[1]/(sum([res[2], res[1], res[0]])))
    p0 = str(res[0]/(sum([res[2], res[1], res[0]])))

    f.write(name + "\t" + str(res[2]) + "\t" + p2 + "\t" + str(res[1]) + "\t" + p1 + "\t" + str(res[0]) + "\t" + p0 + "\n")
