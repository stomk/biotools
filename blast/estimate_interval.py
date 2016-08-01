#!/usr/bin/env python

#[0]query_name  [1]target_name  [2]%identity  [3]alignment_length  [4]#mismatch  [5]#gap_open  [6]query_start  [7]query_end  [8]target_start  [9]target_end  [10]e-value  [11]bit_score
# vertebrate_telomere	24	100.00	30	0	0	1	30	19064738	19064767	5e-08	56.5

import sys
import argparse

def estimate_span(blast_file, contig_file, max_gap_len, min_span_len, min_occupancy, delimiter):
    def should_print(bgn, end, occupied_bp):
        span_len = end - bgn + 1
        occupancy = 1.0 * occupied_bp / span_len
        return (span_len > min_span_len and occupancy > min_occupancy)

    def print_output(target, bgn, end, occupied_bp, dist):
        span_len = end - bgn + 1
        occupancy = 1.0 * occupied_bp / span_len
        print delimiter.join(map(str, [target, bgn, end, end-bgn+1, '%.3f' % occupancy, dist]))

    def within_same_contig(contig_table, chrom, pos_a, pos_b):
        if not contig_table:
            return True
        for pos in contig_table[chrom]:
            if pos_a > pos and pos_b > pos:
                continue
            elif pos_a <= pos and pos_b <= pos:
                return True
            else:
                return False
        return True  # Both pos exist in the last contig


    if contig_file:
        with open(contig_file) as f:
            contig_table = {}
            for line in f:
                fields = line.strip().split()
                chrom = fields[0]
                bgn = int(fields[1])
                contig_table[chrom] = contig_table.setdefault(chrom, []) + [bgn]
    else:
        contig_table = None

    with open(blast_file, 'r') as f:
        hit_table = {}

        for line in f:
            line = line.strip().split('\t')
            target, t_bgn, t_end = (line[1], int(line[8]), int(line[9]))
            if t_bgn > t_end:
                t_bgn, t_end = t_end, t_bgn
            hit_table[target] = hit_table.setdefault(target, []) + [(t_bgn, t_end)]

        for target, intervals in hit_table.items():
            intervals.sort(key=lambda x: x[0])
            span_bgn, prev_end = intervals[0]
            occupied_bp = prev_end - span_bgn + 1
            prev_printed_end = 0
            for bgn, end in intervals[1:]:
                if bgn > prev_end + max_gap_len or not within_same_contig(contig_table, target, prev_end, bgn):
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
    parser.add_argument('-f', '--contig_file', default=None, help='File of contig intervals')
    parser.add_argument('-g', '--max_gap_len',   type=int, default=0, help='Allow gap at maximum of this value')
    parser.add_argument('-s', '--min_span_len',  type=int, default=0, help='Output interval with span length being at least this value')
    parser.add_argument('-p', '--min_occupancy', type=float, default=0.0, help='Output interval with occupancy being at least this value')
    parser.add_argument('-d', '--delimiter', default='\t', help='Specify delimiter of output')
    args = parser.parse_args()
    estimate_span(args.blast_file, args.contig_file, args.max_gap_len, args.min_span_len, args.min_occupancy, args.delimiter)
