#!/bin/bash

cd `pwd`

hmmscan=/usr/local/hmmer/latest/bin/hmmscan

for file in ./*faa
do
  filebase=$(echo ${file%.*})
  fileext=$(echo ${file##*.})
  outfile=$filebase"_hmmscan-pfamA.out"
  domtblout=$filebase"_hmmscan-pfamA.domtblout"
  tblout=$filebase"_hmmscan-pfamA.tblout"

  $hmmscan -o $outfile --tblout $tblout --domtblout $domtblout --acc --noali --cpu 4 /db/pfam/latest/Pfam-A.hmm $file
 done