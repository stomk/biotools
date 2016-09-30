#!/usr/bin/env python

import sys
import os.path
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Produce reverse complement of a given sequence in FASTA format')
parser.add_argument('fasta_file', help='')

args = parser.parse_args()

record = SeqIO.read(args.fasta_file, 'fasta')

rev_seq = record.seq.reverse_complement()
rev_rec = SeqRecord(rev_seq, id=record.id+'_rev', description='')

out_file = os.path.splitext(args.fasta_file)[0] + '.reverse.fasta'
SeqIO.write(rev_rec, out_file, 'fasta')
