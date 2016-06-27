#!/usr/bin/env python

import sys
import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq

from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Generate random sequence of specified length in FASTA format.')
parser.add_argument('-l', '--length', type=int, default=100, help='length of the sequence')
parser.add_argument('-o', '--out', dest='out_file', default='random.fasta', help='output file name')
parser.add_argument('-n', '--name', dest='header_name', default='random seq', help='header name')
parser.add_argument('--with_n', action='store_true', help='include N')
parser.add_argument('--lower', dest='lower_option', action='store_true', help='Output sequence in lowercase')

args = parser.parse_args()

chars = 'ACGT'
if args.with_n:
    chars += 'N'
if args.lower_option:
    chars = chars.lower()

seq = ""
for i in range(args.length):
    seq += random.choice(chars)

SeqIO.write(SeqRecord(Seq(seq), id=args.header_name, description=''), args.out_file, 'fasta')
