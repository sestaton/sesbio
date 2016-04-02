#!/bin/bash

## NB: Remove the "-phred64" option for recently generated data.
##     This was needed to processes some Illumina data generated in Nov. 2010.

cd `pwd`

fastqc=$HOME/apps/FastQC/fastqc

files=($(ls ./*[12].raw.fastq.gz))
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    echo "Working on ${files[i]}" "${files[i+1]}"
    fq_base=$(echo ${files[i]} | sed 's/\.f.*$//' -)
    fq_base2=$(echo ${files[i+1]} | sed 's/\.f.*$//' -)
    out_base=$(echo ${fq_base} | sed 's/\.1\.raw//' -)
    paired=${fq_base}_trimmed_p_Bclip.fq.gz
    paired2=${fq_base2}_trimmed_p_Bclip.fq.gz
    unpaired=${fq_base}_trimmed_s_Bclip.fq.gz
    unpaired2=${fq_base2}_trimmed_s_Bclip.fq.gz
    trimmed=${out_base}_all_trimmed.fq.gz
    fqc_out=${out_base}_fastqc
    mkdir $fqc_out

    cat $paired $paired2 $unpaired $unpaired2 > $trimmed

    $fastqc -t 4 -o $fqc_out $trimmed
    rm $trimmed
done
