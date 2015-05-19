#!/bin/bash

gt=/home/statonse/apps/genometools-unstable/bin/gt
db=/home/statonse/db/Ha412v1r1_genome_no_cp-mt-rd.fasta
dbbase=$(echo ${db%.*})
gff=${dbbase}_ltrharvest99.gff3
gff_sort=${dbbase}_ltrharvest99_sort.gff3
gff_h=${dbbase}_ltrdigest99.gff3
index=${db}.index
hmm=/home/statonse/db/Pfam-A.hmm
trnas=/home/statonse/db/plant_tRNAs.fasta
###################################################
#LTRharvest
###################################################
# create suffix array for ltrharvest
#time $gt suffixerator -db $db -indexname $index -tis -suf -lcp -ssp -sds -des -dna -v

time $gt ltrharvest \
-longoutput \
-mintsd 4 \
-maxtsd 6 \
-minlenltr 100 \
-maxlenltr 6000 \
-mindistltr 1500 \
-maxdistltr 25000 \
-motif tgca \
-similar 99 \
-vic 10 \
-index $index \
-out pred-all_Ha412v1.1_99 \
-outinner pred-inner_Ha412v1.1_99 \
-gff3 $gff

###################################################
#LTRdigest
###################################################
# sort the gff3 file(s) prior to running ltrdigest
$gt gff3 -sort $gff > $gff_sort

time $gt ltrdigest \
-trnas $trnas \
-hmms $hmm \
-aliout yes \
-aaout yes \
-seqfile $db \
-matchdescstart \
-seqnamelen 50 \
-outfileprefix ltrdigest_Ha412v1.1_99_red $gff_sort $index \
-o $gff_h
