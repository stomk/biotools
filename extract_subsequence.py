#!/usr/bin/env python

# --.py fasta begin end

import sys
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()

parser.add_argument('seq_file')
parser.add_argument('-s', '--seq', default = None)
parser.add_argument('-b', '--bgn', type=int)
parser.add_argument('-e', '--end', type=int)
parser.add_argument('-p', '--prefix')

args = parser.parse_args()
seq_file = args.seq_file
seq_name = args.seq
bgn = args.bgn
end = args.end
prefix = args.prefix

records = SeqIO.parse(seq_file, 'fasta')

if seq_name:
    for record in records:
        if record.id == seq_name:
            target = record
            break
    if not target:
        print "Sequence %s was not found." % seq_name
        sys.exit()
else:
    target = next(records)

prefix = prefix if prefix else target.id
subseq = SeqRecord(target.seq[bgn:end], id=(prefix + '_' + str(bgn) + '_' + str(end)), description='')

SeqIO.write(subseq, subseq.id + '.fa', 'fasta')


