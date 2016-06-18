#!/usr/bin/env python

#[0]query_name  [1]target_name  [2]%identity  [3]alignment_length  [4]#mismatch  [5]#gap_open  [6]query_start  [7]query_end  [8]target_start  [9]target_end  [10]e-value  [11]bit_score
# vertebrate_telomere	24	100.00	30	0	0	1	30	19064738	19064767	5e-08	56.5

import sys
import argparse

def calc_total_match_len_per_hit(blast_file, min_total_match_len):
    def print_output(target, match_len):
        if match_len > min_total_match_len:
            print '%s\t%d' % (target, match_len)

    with open(blast_file, 'r') as f:
        prev_target = ""
        total_match_len = 0
        for line in f:
            line = line.strip().split('\t')
            target, aln_len = (line[1], int(line[3]))
            if target == prev_target:
                total_match_len += aln_len
            else:
                print_output(prev_target, total_match_len)
                total_match_len = aln_len
                prev_target = target
        print_output(prev_target, total_match_len)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('blast_file')
    parser.add_argument('-m', '--min_total_match_len', type=int, default=0)
    args = parser.parse_args()
    calc_total_match_len_per_hit(args.blast_file, args.min_total_match_len)
