#!/bin/bash

set -euo pipefile

## Three methods for converting Fastq to Fasta. It appears that seqtk is the easist to 
## use and the fastest method.

# emboss
#seqret -sequence $1 -sformat fastq -osformat fasta -stdout -auto

# awk
#awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' input.fastq > output.fasta

# seqtk
for file in ./*.fastq
do
    f=$(echo ${file%.*})
    fa=${f}.fasta
    seqtk seq -A $file > $fa &
done
