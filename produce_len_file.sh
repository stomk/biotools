#!/bin/bash

FASTA=$1

grep ">" $FASTA | sed 's/>//' > $FASTA.temp1
fatt len $FASTA > $FASTA.temp2
paste $FASTA.temp1 $FASTA.temp2 > $FASTA.len
rm $FASTA.temp1 $FASTA.temp2

