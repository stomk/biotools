#!/usr/bin/env python

#[0]query_name  [1]target_name  [2]%identity  [3]alignment_length  [4]#mismatch  [5]#gap_open  [6]query_start  [7]query_end  [8]target_start  [9]target_end  [10]e-value  [11]bit_score
# vertebrate_telomere	24	100.00	30	0	0	1	30	19064738	19064767	5e-08	56.5

import sys
import argparse

def calc_average_identity(blast_file):
    with open(blast_file, 'r') as f:
        hit_table = {}

        for line in f:
            line = line.strip().split('\t')
            target, identity, hit_len = (line[1], float(line[2]), int(line[3]))
            if hit_len > 140:
                continue
            hit_table[target] = hit_table.setdefault(target, []) + [identity]

        for target, identities in hit_table.items():
            print '%s\t%.2f' % (target, sum(identities)/len(identities))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('blast_file')
    args = parser.parse_args()
    calc_average_identity(args.blast_file)
