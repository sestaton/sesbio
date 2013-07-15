#!/bin/bash

# transall.sh 
# 11/30/11 SES
# Updated 12/5/11

usage() {
cat << EOF

USAGE: $0 <fasta_file> <split_size>

fasta_file   : file of reads to translate
split_size   : integer representing the number of splits 
               to create for faster processing. Every ORF
               will be rejoined in a single file upon 
               completion.
 
EOF
}

print_error() {
cat << ERR

ERROR: Command line not parsed correctly. Check input.

ERR
}

if [ $# -lt 1 ]; then
    print_error
    usage
    exit 1
fi

# gt is required for this script to work
# the below assumes you are running under the bash shell
hash gt &> /dev/null
if [ $? -eq 1 ]; then
    echo >&2 "genometools is required but it's not installed. Exiting."
    exit 1
fi

#
#
fasta=$(echo ${1%.*})
trans=$fasta"_ORFs.faa"

gt splitfasta -numfiles $2 $1

for i in $(seq $2)    # for older bash
#for i in {1..33}       # for bash 3.0+
do
  mkdir $i $i"_tr"
  cd $i
  perl ~/ePerl/external/bp_seqretsplit.pl ../$1"."$i
  cd ../
  perl batch_sixpack3.pl -i $i -o $i"_tr" --clean
  rm -rf $i
  rm $1"."$i
  # Move the cat $file part here, then remove the ORF dir
  # keeping everything tidy
  for file in ./$i"_tr"/*
  do
    cat $file >> $trans
  done
  rm -rf $i"_tr"
done

d=$(date)
echo "$0 finished translating $1 at $d"
