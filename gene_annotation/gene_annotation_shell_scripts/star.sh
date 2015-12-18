#!/bin/bash

set -euo pipefail

script=$(basename $0)

function usage() {
cat <<EOF

USAGE: $script <subj_seq_file>

subj_seq_file   :   A (nucleotide) genome FASTA file to index.
gtf_file        :   The GTF file for the genome.

EOF
}

function print_error() {
cat <<ERR

ERROR: Command line not parsed correctly. Check input.

ERR
}


cd `pwd`

star=$HOME/github/STAR/bin/Linux_x86_64/STAR
samtools=$HOME/github/samtools/samtools

if [ $# -lt 2 ]; then
    print_error
    usage
    exit 1
fi

subjSeq=$1
gtfFile=$2
genome=$(basename $subjSeq)
subjseqFile=$(echo ${genome%.*})
base=$(echo ${subjseqFile%.*})

## create the index                                                                                          
$star --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir /home/statonse/db \
    --genomeFastaFiles $subjSeq \
    --sjdbGTFfile $gtfFile

files=($(ls ./*[12].raw.fastq.gz))
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    #echo "Working on ${files[i]}" "${files[i+1]}"
    fq_base=$(echo ${files[i]} | sed 's/\.f.*$//' -)
    fq_base2=$(echo ${files[i+1]} | sed 's/\.f.*$//' -)
    out_base=$(echo ${fq_base} | sed 's/\.1\.raw//' -)
    paired=${fq_base}_trimmed_p_Bclip.fq.gz
    #paired2=${fq_base2}_trimmed_p_Bclip.fq.gz
    paired2=${fq_base}_trimmed_s_Bclip.fq.gz

    bamsort=${out_base}_${base}_sort.bam
    echo "Working on $paired $paired2"
    echo "bamsort: $bamsort"

    ## run star
    $star --runMode alignReads \
	--runThreadN 4 \
	--genomeDir /home/statonse/db \
	--readFilesIn $paired $paired2 \
	--outFileNamePrefix ${out_base}_aligned \
	--outSAMtype BAM SortedByCoordinate \
	--readFilesCommand zcat \
	--outStd BAM_SortedByCoordinate > $bamsort

    echo -e "star done...\n"
done
