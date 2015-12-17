#!/bin/bash

## NB: Remove the "-phred64" option for recently generated data.
##     This was needed to processes some Illumina data generated in Nov. 2010.

cd `pwd`

trimfile=$HOME/apps/Trimmomatic-0.35/adapters/TruSeq2-PE.fa
trimmomatic=$HOME/apps/Trimmomatic-0.35/trimmomatic-0.35.jar
Btrim=$HOME/github/sesbio/genome_assembly/remove_trailing_Bs.pl 

files=($(ls ./*[12].raw.fastq.gz))
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    echo "Working on ${files[i]}" "${files[i+1]}"
    fq_base=$(echo ${files[i]} | sed 's/\.f.*$//' -)
    fq_base2=$(echo ${files[i+1]} | sed 's/\.f.*$//' -)
    log=${fq_base}_trimmomatic.log
    paired=${fq_base}_trimmed_p.fq.gz
    paired2=${fq_base2}_trimmed_p.fq.gz
    unpaired=${fq_base}_trimmed_s.fq.gz
    unpaired2=${fq_base2}_trimmed_s.fq.gz

    goodnoBp=${fq_base}_trimmed_p_Bclip.fq.gz
    goodnoBp2=${fq_base2}_trimmed_p_Bclip.fq.gz
    goodnoBs=${fq_base}_trimmed_s_Bclip.fq.gz
    goodnoBs2=${fq_base2}_trimmed_s_Bclip.fq.gz

    ## quality trim
    java -jar $trimmomatic PE -phred64 \
	-threads 12 -trimlog $log \
	${files[i]} ${files[i+1]} \
	$paired $unpaired \
	$paired2 $unpaired2 \
	ILLUMINACLIP:$trimfile:2:30:10 LEADING:3 TRAILING:3 \
	SLIDINGWINDOW:4:15 MINLEN:40

    ## B-trim
    perl $Btrim -i $paired -l 40 -o STDOUT | gzip > $goodnoBp
    perl $Btrim -i $paired2 -l 40 -o STDOUT | gzip > $goodnoBp2
    perl $Btrim -i $unpaired -l 40 -o STDOUT | gzip > $goodnoBs
    perl $Btrim -i $unpaired2 -l 40 -o STDOUT | gzip > $goodnoBs2

    rm $paired $paired2 $unpaired $unpaired2
done
