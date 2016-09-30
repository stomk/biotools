#!/usr/bin/env python

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Split multi-fasta into individual fasta file')
parser.add_argument('fasta_file', help='')
parser.add_argument('-p', '--prefix', default='', help='Prefix for file name and fasta header')
parser.add_argument('-d', '--outdir', help='Output directory')

args = parser.parse_args()

if args.outdir:
    outdir = args.outdir
    if not os.path.exists(outdir):
        os.mkdir(outdir)
else:
    outdir = ''

for r in SeqIO.parse(args.fasta_file, 'fasta'):
    name = args.prefix + '_' + r.id if args.prefix else r.id
    r.id = name
    r.description = ''
    SeqIO.write(r, os.path.join(outdir, name + '.fasta'), 'fasta')
