#!/usr/bin/env python

#[0]query_name  [1]target_name  [2]%identity  [3]alignment_length  [4]#mismatch  [5]#gap_open  [6]query_start  [7]query_end  [8]target_start  [9]target_end  [10]e-value  [11]bit_score
# vertebrate_telomere	24	100.00	30	0	0	1	30	19064738	19064767	5e-08	56.5

import sys
import argparse

def estimate_span(blast_file, max_gap_len, min_span_len, min_occupancy, delimiter):
    def should_print(bgn, end, occupied_bp):
        span_len = end - bgn + 1
        occupancy = 1.0 * occupied_bp / span_len
        return (span_len > min_span_len and occupancy > min_occupancy)

    def print_output(target, bgn, end, occupied_bp, dist):
        span_len = end - bgn + 1
        occupancy = 1.0 * occupied_bp / span_len
        print delimiter.join(map(str, [target, bgn, end, end-bgn+1, '%.3f' % occupancy, dist]))

    with open(blast_file, 'r') as f:
        hit_table = {}

        for line in f:
            line = line.strip().split('\t')
            target, t_bgn, t_end = (line[1], int(line[8]), int(line[9]))
            if t_bgn > t_end:
                tmp = t_bgn
                t_bgn = t_end
                t_end = tmp
            if target in hit_table:
                hit_table[target] = hit_table[target] + [(t_bgn, t_end)]
            else:
                hit_table[target] = [(t_bgn, t_end)]

        for target, intervals in hit_table.items():
            intervals.sort(key=lambda x: x[0])
            span_bgn, prev_end = intervals[0]
            occupied_bp = prev_end - span_bgn + 1
            prev_printed_end = 0
            for bgn, end in intervals[1:]:
                if bgn > prev_end + max_gap_len:
                    if should_print(span_bgn, prev_end, occupied_bp):
                        dist_from_prev = span_bgn - prev_printed_end if prev_printed_end != 0 else '-'
                        print_output(target, span_bgn, prev_end, occupied_bp, dist_from_prev)
                        prev_printed_end = prev_end
                    span_bgn = bgn
                    occupied_bp = 0
                prev_end = end
                occupied_bp += end - bgn + 1
            if should_print(span_bgn, prev_end, occupied_bp):
                dist_from_prev = span_bgn - prev_printed_end if prev_printed_end != 0 else '-'
                print_output(target, span_bgn, prev_end, occupied_bp, dist_from_prev)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Parse BLAST output file (outfmt=6) and estimate continual intervals of hits')
    parser.add_argument('blast_file')
    parser.add_argument('-g', '--max_gap_len',   type=int, default=0)
    parser.add_argument('-s', '--min_span_len',  type=int, default=0)
    parser.add_argument('-p', '--min_occupancy', type=float, default=0.0)
    parser.add_argument('-d', '--delimiter', default='\t')
    args = parser.parse_args()
    estimate_span(args.blast_file, args.max_gap_len, args.min_span_len, args.min_occupancy, args.delimiter)
