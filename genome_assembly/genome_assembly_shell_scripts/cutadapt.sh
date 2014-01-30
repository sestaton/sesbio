#!/bin/bash

for file in ./*trimmed.fasta
do
    f=$(echo ${file%.*})
    fo=${f}_noadapter.fasta
    fl=${f}_cutadapt.log
    python /usr/local/bioinfo/cutadapt/cutadapt-1.1/bin/cutadapt --times 3 -f fasta -g CAAGTCGT -g TCACCTAG -g GACACAGT -m 50 $file -o $fo 2>&1 > $fl 
done
