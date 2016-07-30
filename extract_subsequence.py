#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='')

parser.add_argument('seq_file', help='Target FASTA file')
parser.add_argument('-f', '--file', dest='interval_file', help='File containing list of seq, bgn, end ')
parser.add_argument('-o', '--out', dest='out_file', default='extracted.fasta', help='Output file')
parser.add_argument('-s', '--seq', dest='seq_name', default=None, help='Specify single sequence name to extract')
parser.add_argument('-b', '--bgn', type=int, help='Specify beginning position')
parser.add_argument('-e', '--end', type=int, help='Specify end position')

args = parser.parse_args()


def gen_subseq(records, name, bgn, end):
    if not name in records:
        print >> sys.stderr, "Sequence '%s' does not exist in target file" % name
        return
    target_record = records[name]
    subseq = target_record.seq[bgn-1:end]
    record = SeqRecord(subseq, id='_'.join(map(str, [name, bgn, end])), description='')
    return record


if args.interval_file:
    with open(args.interval_file) as f:
        records = SeqIO.index(args.seq_file, 'fasta')
        extracted_records = []
        for line in f:
            name, bgn, end = line.strip().split()[:3]
            bgn, end = int(bgn), int(end)
            record = gen_subseq(records, name, bgn, end)
            if record:
                extracted_records.append(record)

    SeqIO.write(extracted_records, args.out_file, 'fasta')

elif args.seq_name:
    records = SeqIO.index(args.seq_file, 'fasta')
    record = gen_subseq(records, args.seq_name, args.bgn, args.end)
    SeqIO.write(record, args.out_file, 'fasta')

else:
    print >> sys.stderr, "[Error] Specify target sequence name by using either '-f' or '-s' option"


