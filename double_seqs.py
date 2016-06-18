#!/usr/bin/env python

import sys
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument('fasta_file')

args = parser.parse_args()

doubled_seqs = []

for seq_record in SeqIO.parse(args.fasta_file, 'fasta'):
    seq = Seq(str(seq_record.seq) * 2)
    
    id = seq_record.id + '_double'
    doubled_seqs.append(SeqRecord(seq, id, description=''))

out_file = args.fasta_file + '.doubled'
SeqIO.write(doubled_seqs, out_file, 'fasta')
                        
