#!/usr/bin/env python

import sys
from Bio import SeqIO

seq_file = sys.argv[1]
out_file = open(seq_file + '.lcr', 'w')

reads = SeqIO.parse(seq_file, 'fasta')
for read in reads:
	seq = read.seq
	if len(seq) == 0:
		correct_rate = 0.0
	else:
		num_corrected_letters = 0
		num_corrected_letters += seq.count('A') + seq.count('C') + seq.count('G') + seq.count('T')
		correct_rate = float(num_corrected_letters) / len(seq)
	out_file.write('%s\t%d\t%.3f\n' % (read.id, len(seq), correct_rate))


