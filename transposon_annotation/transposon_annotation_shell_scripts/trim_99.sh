#!/bin/bash

#wget -O- http://genometools.org/pub/genometools-unstable.tar.gz | tar xzf -
cd genometools-unstable
#&& make -j 12 64bit=yes
gt=`pwd`/bin/gt
cd ..

#gt=/home/statonse/apps/genometools-unstable/bin/gt
db=Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta
dbbase=$(echo ${db%.*})
gff=${dbbase}_trim_ltrharvest99.gff3
gff_sort=${dbbase}_trim_ltrharvest99_sort.gff3
gff_h=${dbbase}_trim_ltrdigest99.gff3
index=${db}.index
hmm=/db/transposable+element_hmms/transposable+element-3.hmm
trnas=/db/eukaryotic-tRNAs.fas

###################################################
#LTRharvest
###################################################
# create suffix array for ltrharvest
#time $gt suffixerator -db $db -indexname $index -tis -suf -lcp -ssp -sds -des -dna -v

time $gt ltrharvest \
-seqids yes \
-longoutput yes \
-mintsd 4 \
-maxtsd 6 \
-minlenltr 70 \
-maxlenltr 500 \
-mindistltr 280 \
-maxdistltr 1500 \
-motif tgca \
-similar 99 \
-vic 10 \
-index $index \
-out pred-all_Ha412v1.1_trim99 \
-outinner pred-inner_Ha412v1.1_trim99 \
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
-seqnamelen 50 \
-matchdescstart yes \
-o $gff_h \
-outfileprefix ltrdigest_Ha412v1.1_trim99 $gff_sort

