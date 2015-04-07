#!/bin/bash

fsamp=$1
rsamp=$2
fbase=$(echo ${fsamp%.*})
rbase=$(echo ${rsamp%.*})
fsamp_scr=${fbase}_scr.fastq
rsamp_scr=${rbase}_scr.fastq
db=/home/jmblab/statonse/db/ScreenDB_Hannus-cpDNA_UniVec-5-2_Plant-mtDNA_LSU-SSU-DB.fasta
out=out_%.fq

time ~/apps/bbmap/bbsplit.sh -Xmx8g in1=$fsamp in2=$rsamp ref=$db basename=$out outu1=$fsamp_scr outu2=$rsamp_scr
