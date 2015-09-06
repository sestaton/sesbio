#!/bin/bash

cd genometools-unstable
gt=`pwd`/bin/gt
cd ..

#gt=/usr/local/bioinfo/genometools/genometools-unstable/bin/gt
db=Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta
dbbase=$(echo ${db%.*})
gff=${dbbase}_trim_ltrharvest85.gff3
gff_sort=${dbbase}_trim_ltrharvest85_sort.gff3
gff_h=${dbbase}_trim_ltrdigest85.gff3
index=${db}.index
#hmm=/db/Pfam-A.hmm
hmm=/db/transposable+element_hmms/transposable+element-4.hmm
trnas=/db/plant_tRNAs.fasta
###################################################
#LTRharvest
###################################################
# create suffix array for ltrharvest
#time $gt suffixerator -db $db -indexname $index -tis -suf -lcp -ssp -sds -des -dna -v
# reuse db from LTR analysis

time $gt ltrharvest \
-longoutput yes \
-seqids yes \
-mintsd 4 \
-maxtsd 6 \
-minlenltr 70 \
-maxlenltr 500 \
-mindistltr 280 \
-maxdistltr 1500 \
-similar 85 \
-vic 10 \
-index $index \
-out pred-all_Ha412v1.1_trim85 \
-outinner pred-inner_Ha412v1.1_trim85 \
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
-matchdescstart yes \
-seqnamelen 50 \
-o $gff_h \
-outfileprefix ltrdigest_Ha412v1.1_trim85 $gff_sort

