#!/usr/bin/env python

import sys
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()

parser.add_argument('seq_file')
parser.add_argument('-l', '--trim_len', type=int)
parser.add_argument('-p', '--prefix')

args = parser.parse_args()
seq_file = args.seq_file
trim_len = args.trim_len
prefix = args.prefix

def write_subseqence_fasta(seq_record, bgn, end, prefix=''):
    prefix = prefix + '_' if prefix else ''
    subseq = SeqRecord(seq_record.seq[bgn : end], id=(prefix + seq_record.id + '_' + str(bgn) + '_' + str(end)), description='')
    SeqIO.write(subseq, subseq.id + '.fa', 'fasta')

for seq_record in SeqIO.parse(seq_file, 'fasta'):
    seq_len = len(seq_record.seq)
    write_subseqence_fasta(seq_record, 0, trim_len, prefix)
    write_subseqence_fasta(seq_record, (seq_len - trim_len), seq_len, prefix)



