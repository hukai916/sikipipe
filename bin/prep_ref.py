#!/usr/bin/env python

"""
Simply rc ref
"""

import sys
import gzip
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO
from Bio.Seq import Seq


filename   = sys.argv[1]
outfile    = sys.argv[2]

encoding = guess_type(filename)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open # must be rt for SeqIO

records = []
with _open(filename) as f:
    for record in SeqIO.parse(f, 'fasta'):
        tem = record
        tem.seq = record.reverse_complement().seq
        records.append(tem)
SeqIO.write(records, outfile, "fasta")
