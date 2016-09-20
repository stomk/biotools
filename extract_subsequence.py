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
parser.add_argument('-n', '--seq', dest='seq_name', default=None, help='Specify single sequence name to extract')
parser.add_argument('-s', '--single', action='store_true', help='Process FASTA file with a single sequence in it')
parser.add_argument('-b', '--bgn', type=int, help='Specify beginning position')
parser.add_argument('-e', '--end', type=int, help='Specify end position')
parser.add_argument('-r', '--orientation', action='store_true')

args = parser.parse_args()


def gen_subseq_from_record(record, bgn, end, reverse=False):
    subseq = record.seq[bgn-1:end]
    if reverse:
        subseq = subseq.reverse_complement()
    return SeqRecord(subseq, id='_'.join(map(str, [record.id, bgn, end])), description='')

def gen_subseq_from_records(records, name, bgn, end, reverse=False):
    if not name in records:
        print >> sys.stderr, "Sequence '%s' does not exist in target file" % name
        return
    target_record = records[name]
    return gen_subseq_from_record(target_record, bgn, end, reverse)


if args.interval_file:
    with open(args.interval_file) as f:
        extracted_records = []
        if args.single:
            record = SeqIO.read(args.seq_file, 'fasta')
            for line in f:
                fields = line.strip().split()
                bgn, end = map(int, fields[0:2])
                reverse = True if args.orientation and fields[4] == '-' else False
                s_record = gen_subseq_from_record(record, bgn, end, reverse)
                if s_record:
                    extracted_records.append(s_record)
        else:
            records = SeqIO.index(args.seq_file, 'fasta')
            for line in f:
                fields = line.strip().split()
                name = fields[0]
                bgn, end = map(int, fields[1:3])
                if args.orientation:
                    reverse = True if fields[3] == '-' else False
                else:
                    reverse = False
                if args.single:
                    record = gen_subseq_from_record(record, bgn, end, reverse)
                else:
                    record = gen_subseq_from_records(records, name, bgn, end, reverse)
                if record:
                    extracted_records.append(record)
    SeqIO.write(extracted_records, args.out_file, 'fasta')

elif args.seq_name:
    records = SeqIO.index(args.seq_file, 'fasta')
    record = gen_subseq_from_records(records, args.seq_name, args.bgn, args.end)
    SeqIO.write(record, args.out_file, 'fasta')

else:
    print >> sys.stderr, "Sequence name was not specified. Assuming the fasta file has only one sequence"
    record = SeqIO.read(args.seq_file, 'fasta')
    extracted_record = gen_subseq_from_record(record, args.bgn, args.end)
    SeqIO.write(extracted_record, args.out_file, 'fasta')

    

