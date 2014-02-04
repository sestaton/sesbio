#!/bin/bash

## This code runs the clustering method of Novak et al. 2009, and was used for the analysis of 
## the sunflower genome in Staton et al. 2012, TPJ.
##
## NB: Don't use this method. It's prohibitively slow  and does not make use of paired-end reads.
##     Use Transposome (https://github.com/sestaton/Transposome) instead.

cd /scratch/statonse/fgclust_ha412ho

#/home/jmblab/statonse/apps/gicl/mgblast \
#-i ha412ho-short_IDs.fna.clean_nochloro_unique.fasta.1.1 \
#-d ha412ho_1_1_unique \
#-F "m D" \
#-D 4 \
#-p 85 \
#-W18 \
#-UT \
#-X40 \
#-KT \
#-JF \
#-v90000000 \
#-b90000000 \
#-C80 \
#-H 320 \
#-a 4 \
#-o ha412ho_1_1_self_mgblast_out.txt

#perl ~/ePerl/parse_mgblast.pl -i ha412ho_1_1_self_mgblast_out.txt -o ha412ho_1_1_self_mgblast_out_parsed.txt -id 90.00 -cov 0.15

/home/jmblab/statonse/apps/seqgrapher/fgclust/clusters2graphs.R -H ha412ho_1_self_mgblast_out_parsed.txt -c ha412ho_1_self_mgblast_out_parsed.txt.cls -p 4 -s ha412ho_1_clusters