#!/bin/bash

cd `pwd`

for file in ./*cpDNA_seqs.fasta
do
  qryFile=$(echo ${file%.*})
  fFile=${qryFile}_1.fasta
  rFile=${qryFile}_2.fasta
  fpFile=${qryFile}_1_p.fasta
  rpFile=${qryFile}_2_p.fasta
  fsFile=${qryFile}_1_s.fasta
  rsFile=${qryFile}_2_s.fasta
  iFile=${qryFile}_interl.fasta

  # split the read pairs into separate files
  perl ~/ePerl/interleaved2pairs.pl -i $file -f $fFile -r $rFile

  # pair the reads and put singletons into separate files
  perl ~/ePerl/pairfq.pl -f $fFile -r $rFile -fp $fpFile -rp $rpFile -fs $fsFile -rs $rsFile -im

  # interleave the read pairs for velvet
  ~/apps/bin/shuffleSeqs-fast $fpFile $rpFile > $iFile
done