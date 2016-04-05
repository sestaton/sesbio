#!/bin/bash

cd `pwd`

export OMP_NUM_THREADS=3
export OMP_THREAD_LIMIT=4

velvetopt=/home/jmblab/statonse/apps/VelvetOptimiser-2.2.0/VelvetOptimiser.pl 

for file in ./*cpDNA_seqs.fasta
do
  species=$(echo $file | perl -e 'my $line = <>; $line =~ s/\_.*//; $line =~ s/\.\///; print $line;')
  qryFile=$(echo ${file%.*})
  fsFile=${qryFile}_1_s.fasta
  rsFile=${qryFile}_2_s.fasta
  iFile=${qryFile}_interl.fasta

  perl $velvetopt -s 79 -e 99 -f "-fasta -shortPaired $iFile -fasta -short $fsFile -fasta -short $rsFile" \
-t 8 --optFuncKmer 'n50' -p ${species}_VelvetOpt_k79-99 -d ${species}_VelvetOpt_k79-99
done