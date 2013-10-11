#!/bin/bash

FILES=/scratch/statonse/for_AAARF_split_sm/all_builds_blasthits/blast_all_4_repeat_dbs/AAARF_split10/*fasta

for f in $FILES
do
  cat $f
done