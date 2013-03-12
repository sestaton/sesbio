#/bin/bash

for file in Phoeb_ATCACG_L005_R1_001_trimmed.fastq Phoeb_ATCACG_L005_R2_001_trimmed.fastq
do
    f=$(echo ${file%.*})
    fa=${f}.fasta
    seqtk seq -A $file > $fa &
done
