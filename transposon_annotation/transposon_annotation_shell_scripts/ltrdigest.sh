#!/bin/bash

##NB: 1) Having the database file in the working directory is perferable
##    so that all the results files are in that directory.
##    2) The suffixerator command uses a lot of memory, so be aware (~20 Gb RAM
##    for a plant genome).
##    3) The index commands for ltrharvest are different than those for tirvish.

gt=/home/statonse/apps/genometools-unstable/bin/gt
db=Ha412v1r1_genome_no_cp-mt-rd.fasta
dbbase=$(echo ${db%.*})
gff=${dbbase}_ltrharvest.gff3
gff_sort=${dbbase}_ltrharvest_sort.gff3
gff_h=${dbbase}_ltrdigest.gff3
index=${db}.index
hmm=/home/statonse/db/Pfam-A.hmm
trnas=/home/statonse/db/plant_tRNAs.fasta
###################################################
#LTRharvest
###################################################
# create suffix array for ltrharvest
time $gt suffixerator -db $db -indexname $index -tis -suf -lcp -ssp -sds -des -dna -v

# run ltrharvest
time $gt ltrharvest \
-seqids \
-longoutput \
-mintsd 3 \
-minlenltr 100 \
-maxlenltr 6000 \
-maxdistltr 25000 \
-motif tgca \
-similar 99 \
-vic 10 \
-index $index \
-out pred-all_Ha412v1.1 \
-outinner pred-inner_Ha412v1.1 \
-gff3 $gff

###################################################
#LTRdigest
###################################################
# sort the gff3 file(s) prior to running ltrdigest
$gt gff3 -sort $gff > $gff_sort

# run ltrdigest
time $gt ltrdigest \
-trnas $trnas \
-hmms $hmm \
-aliout yes \
-aaout yes \
-seqfile $db \
-matchdescstart \
-seqnamelen 50 \
-o $gff_h \
-outfileprefix ltrdigest_Ha412v1.1 $gff_sort $index
