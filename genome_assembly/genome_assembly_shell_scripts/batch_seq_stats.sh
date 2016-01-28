#!/bin/bash

## The output is indentical to genometools seq (i.e., gt seq -stat file.fastq)
## However, genometools does not interpret the filename correctly when reading
## from a pipe, and this is no good for logging purposes.

set -e
set -u
set -o pipefail

for file in ./*.gz
do
    bioawk -c fastx '{SUM+=length($seq)} END {printf("Showing statistics for sequence file: "FILENAME"\nNumber of sequences: %d\nTotal length: %d\nMean length: %d\n",NR,SUM,SUM/NR)}' $file
done
