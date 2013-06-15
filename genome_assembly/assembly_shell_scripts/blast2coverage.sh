#!/bin/bash

function usage() {
cat <<EOF

USAGE: $0 <subj_seq_file> <qry_seq_file> 

subj_seq_file   :   A (nucleotide) Fasta file to index.
qry_seq_file    :   A (nucleotide) Fasta file to query the index.

EOF
}

function print_error() {
cat <<ERR

ERROR: Command line not parsed correctly. Check input.

ERR
}

if [ $# -lt 2 ]; then
    print_error
    usage
    exit 1
fi

for prog in formatdb blastall samtools 
do
  hash $prog &> /dev/null
  if [ $? -eq 1 ]; then
      echo >&2 "$prog is required but it's not installed. Exiting."
      exit 1
  fi
done

subjSeq=$1
qrySeq=$2
subjseqFile=$(echo ${subjSeq%.*})
qryseqFile=$(echo ${qrySeq%.*})
db=${subjSeq}_blastdb
blastout=${qryseqFile}_${subjseqFile}.blastn
sam=${qryseqFile}_${subjseqFile}.sam
bam=${qryseqFile}_${subjseqFile}.bam
bamsort=${qryseqFile}_${subjseqFile}_sort

## create the database
formatdb -p F -i $subjSeq -n $db

## run blast
blastall -p blastn -e 1e-10 -a 8 -i $qrySeq -d $db -o $blastout 2> /dev/null

echo -e "blast done...\n"
## convert blast to sam
perl blast2sam.pl $blastout > $sam

echo -e "blast2sam conversion done...\n"
## convert sam to bam, using reference to generate sam headers
samtools view -T $subjSeq -bS $sam > $bam

echo -e "sam to bam conversion done...\n"
## sort the bam
samtools sort $bam $bamsort

echo -e "sam sort done...\n"
## calculate the depth of coverage
samtools depth ${bamsort}.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'

echo -e "samtools depth done...\n"
## clean up
rm ${db}.*