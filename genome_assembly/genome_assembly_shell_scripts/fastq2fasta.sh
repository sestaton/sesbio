#!/bin/bash

seqret -sequence $1 -sformat fastq -osformat fasta -stdout -auto

#awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' input.fastq > output.fasta