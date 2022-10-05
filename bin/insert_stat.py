#!/usr/bin/env python

"""
Parse out large insertion information from BAM file

Dev notes:
        consume_query consume_reference
    M 0 yes yes
    I 1 yes no
    D 2  no yes
    N 3  no yes
    S 4 yes no
    H 5  no no
    P 6  no no
    = 7  yes yes
    X 8  yes yes
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

bam = sys.argv[1]
outfile = sys.argv[2]

samfile = pysam.AlignmentFile(bam, "r")

def consume_seq(query_s, ref_s, cigar_tuple):
    # update query_s:
    if (cigar_tuple[0]) in [0, 1, 4, 7, 8]:
        query_s = query_s + cigar_tuple[1]
    # update ref_s:
    if (cigar_tuple[0]) in [0, 2, 3, 7, 8]:
        ref_s = ref_s + cigar_tuple[1]
    return([query_s, ref_s])

# Note, supplementary reads don't have query (None)
insert = [] #: [query_name, CIGAR, insert_size, ref_s, insert_seq, read_seq]
with open(outfile, "w") as f:
    f.write("read_id,pos_insert,size_insert,cigar,seq_insert\n")
    for read in samfile:
        query_n = read.qname
        cigar   = read.cigarstring
        query_s = 0
        query_seq = read.query_sequence
        ref_s   = read.reference_start
        for m in read.cigar:
            if m[0] == 1: # if insert
                insert_size = m[1]
                if query_seq == None:
                    insert_seq = "None"
                    f.write(",".join([query_n, str(ref_s), str(insert_size), cigar, str(insert_seq)]) + "\n")
                else:
                    insert_seq = query_seq[query_s:query_s + insert_size]
                    f.write(",".join([query_n, str(ref_s), str(insert_size), cigar, insert_seq]) + "\n")

            query_s, ref_s = consume_seq(query_s, ref_s, m)
