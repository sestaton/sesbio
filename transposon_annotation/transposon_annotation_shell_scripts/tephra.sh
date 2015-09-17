#!/bin/bash

cd `pwd`

genome=Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta
trnadb=plant_tRNAs.fasta
hmmdb=transposable+element.hmm

time tephra findltrs -g $genome -t $trnadb -d $hmmdb -c
