#!/bin/bash

set -e
set -u
set -o pipefail

##TODO: take db and cpus as options

cd `pwd`

for prog in hmmscan
do
  hash $prog &> /dev/null
  if [ $? -eq 1 ]; then
      echo >&2 "$prog is required but it's not installed. Exiting."
      exit 1
  fi
done

cpu=4

for file in ./*.faa
do
  filebase=$(echo ${file%.*})
  fileext=$(echo ${file##*.})
  outfile=$filebase"_hmmscan-pfamA.out"
  domtblout=$filebase"_hmmscan-pfamA.domtblout"
  tblout=$filebase"_hmmscan-pfamA.tblout"

  hmmscan -o $outfile --tblout $tblout --domtblout $domtblout \
      --acc --noali --cpu $cpu /db/pfam/latest/Pfam-A.hmm $file
done
