#!/bin/bash

## NB: This script is for testing fgclust, which is now RepeatExplorer.
##     It is advised to use Transposome for analysis, and this script for any benchmarking of RepeatExplorer.

cd `pwd`

hitsort2clusters=/home/jmblab/statonse/apps/seqgrapher/fgclust/hitsort2clusters.R
mgblast=/home/jmblab/statonse/apps/gicl/mgblast 

$mgblast \
-i ha412ho-short_IDs.fna.clean_nochloro_unique.fasta.1 \
-d ha412ho_1_clean \
-F "m D" \
-D 4 \
-p 85 \
-W18 \
-UT \
-X40 \
-KT \
-JF \
-v90000000 \
-b90000000 \
-C80 \
-H 320 \
-a 4 \
-o ha412ho_1_self_mgblast_out.txt

perl ~/github/transposon_annotation/parse_mgblast.pl -i ha412ho_1_self_mgblast_out.txt -o ha412ho_1_self_mgblast_out_parsed.txt -id 90.00 -cov 0.15

$hitsort2clusters -i ha412ho_1_self_mgblast_out_parsed.txt