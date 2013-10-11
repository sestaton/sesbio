#!/bin/bash

function usage() {
cat <<EOF

USAGE: $0 <seq_file> <kmer_len> 

seq_file   :   A (nucleotide) Fasta file to analyze.
kmer_len   :   An integer value to use in the analysis (e.g., 20).

N.B. The number of threads to use and threshold for k-mers to write to a sequence file are 
hard coded (lines 43 and 44). Set as appropriate.

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

for prog in meryl
do
  hash $prog &> /dev/null
  if [ $? -eq 1 ]; then
      echo >&2 "$prog is required but it's not installed. Exiting."
      exit 1
  fi
done

#
# Change these values to what is appropriate for your:
#
threads=2    # computational resources
thresh=100   # biological question

seq=$1
merLen=$2
seqFile=$(echo ${seq%.*})
db=$seqFile"_"$merLen"mer_meryldb"
merCt=$seqFile"_"$merLen"mer_counts.txt"
merHist=$seqFile"_"$merLen"mer_hist.txt"
merSeq=$seqFile"_"$merLen"mer_seqs.fasta"

# create the database
meryl -threads $threads -v -B -m $merLen -s $seq -o $db

# dump mer counts
meryl -threads $threads -s $db -Dc > $merCt

# dump a histogram of counts
meryl -threads $threads -s $db -Dh > $merHist

# dump all $merLen mers above a threshhold 
meryl -threads $threads -s $db -Dt -n $thresh > $merSeq