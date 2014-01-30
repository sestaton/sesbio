#!/bin/bash

OF=ha412ho_500k_split1_clusteredlist.txt

for file in ./*fasta
do
  grep ">" $file | sed 's/>//' >> $OF
done