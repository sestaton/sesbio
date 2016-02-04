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

##TODO: make this an option as it takes a long time
## create the index
#$star --runThreadN 4 \
#    --runMode genomeGenerate \
#    --genomeDir /home/statonse/db \
#    --genomeFastaFiles $subjSeq \
#    --sjdbGTFfile $gtfFile

files=($(ls ./*[12].raw_trimmed_p.fq.gz))
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    out_base=$(echo ${files[i]} | sed 's/\.1\.raw.*//' -)
    paired=${files[i]}
    paired2=${files[i+1]}

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

    ## NB: may have to manually increase buffer size in __init__.py of HTSeq code
    ## or it will die
    out=${out_base}_counts-exon.txt
    samtools view $bamsort | htseq-count -f sam --order pos --minaqual 1 -s no -t exon -i gene_id - $gtfFile > $out
    echo -e "htseq-count done...\n"
done
