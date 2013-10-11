#!/bin/bash

for file in ./*.fasta
do
  ~/apps/bioawk/bioawk -c fastx '{SUM+=length($seq)} END {printf("Showing statistics for sequence file: "FILENAME"\nNumber of sequences: %d\nTotal length: %d\nMean length: %d\n",NR,SUM,SUM/NR)}' $file &
done
