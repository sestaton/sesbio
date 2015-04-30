#!/bin/bash

cd `pwd`

db=/home/statonse/db/Ha412v1r1_genome_no_cp-mt-rd.fasta
dbbase=$(echo ${db%.*})
gff=${dbbase}_tirvish.gff3
gff_sort=${dbbase}_tirvish_sort.gff3
index=${db}.index
hmm=/home/statonse/db/Pfam-A.hmm
gt=/home/statonse/apps/genometools-unstable/bin/gt

## index genome
## NB: the '-mirrored' option is required for tirvish(to search both strands)
##     Also, the database is incompatible with ltrharvest, so see the ltrdigest.sh 
##     script for details on the usage of that program
time $gt suffixerator -db $db -indexname $index -tis -suf -lcp -des -ssp -dna -mirrored -v

## run tirvish
time $gt tirvish -index $index -hmms $hmm > $gff
time $gt gff3 -sort $gff > $gff_sort