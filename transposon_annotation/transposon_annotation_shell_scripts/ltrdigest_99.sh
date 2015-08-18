#!/bin/bash

## get genometools
#wget -O- http://genometools.org/pub/genometools-unstable.tar.gz | tar xzf -
cd genometools-unstable
#&& make -j 12 64bit=yes
gt=`pwd`/bin/gt
cd ..

#gt=/home/statonse/apps/genometools-unstable/bin/gt
db=Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta
dbbase=$(echo ${db%.*})
gff=${dbbase}_ltrharvest99.gff3
gff_sort=${dbbase}_ltrharvest99_sort.gff3
gff_h=${dbbase}_ltrdigest99.gff3
index=${db}.index
hmm=/db/transposable+element_hmms/transposable+element.hmm
trnas=/db/eukaryotic-tRNAs.fas
###################################################
#LTRharvest
###################################################
# create suffix array for ltrharvest
#time $gt suffixerator -db $db -indexname $index -tis -suf -lcp -ssp -sds -des -dna

# ltrharvest
time $gt ltrharvest \
-longoutput yes \
-seqids yes \
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
-out Ha412v1r1_ltrharvest99_pred-all \
-outinner Ha412v1r1_ltrharvest99_pred-inner \
-gff3 $gff

###################################################
#LTRdigest
###################################################
# sort the gff3 file(s) prior to running ltrdigest
time $gt gff3 -sort $gff > $gff_sort

# ltrdigest
time $gt ltrdigest \
-trnas $trnas \
-hmms $hmm \
-aliout yes \
-aaout yes \
-seqfile $db \
-matchdescstart yes \
-seqnamelen 50 \
-o $gff_h \
-outfileprefix Ha412v1r1_ltrdigest99 $gff_sort
