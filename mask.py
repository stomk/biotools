#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument('fasta_file')
parser.add_argument('interval_file')
args = parser.parse_args()

def mask(seq, interval):
    bgn, end = map(lambda x:x-1, interval)  # TRF output is 1-origin, while seq object is 0-origin
    return seq[:bgn] + 'N'*(end - bgn + 1) + seq[end+1:]  
    

interval_table = {}
with open(args.interval_file, 'r') as interval_f:
    for line in interval_f:
        line = line.rstrip('\n').split('\t')
        name = line[0]
        bgn, end = map(int, line[1:3])
        interval_table[name] = interval_table.setdefault(name, []) + [(bgn, end)]
   

out_records = []
for seq_record in SeqIO.parse(args.fasta_file, 'fasta'):
    name = seq_record.id
    if not name in interval_table:
        continue
    seq = seq_record.seq
    intervals = interval_table[name]
    for interval in intervals:
        seq = mask(seq, interval)
    out_records.append(SeqRecord(seq, id=name, description=''))

out_file = args.fasta_file + '.masked'
SeqIO.write(out_records, out_file, 'fasta')
