#!/bin/bash

set -euo pipefail

function usage() {
cat <<EOF

USAGE: $0 <seq_file> <kmer_len> 

seq_file   :   A (nucleotide) Fasta file to analyze.
kmer_len   :   An integer value to use in the analysis (e.g., 20).

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

for prog in gt
do
  hash $prog &> /dev/null
  if [ $? -eq 1 ]; then
      echo >&2 "genometools is required but it's not installed. Exiting."
      exit 1
  fi
done

seq=$1
merLen=$2
seqFile=$(echo ${seq%.*})
db=$seqFile"_"$merLen"mer_tallymer_suffixarr"
merCt=$seqFile"_"$merLen"mer_tallymer_counts.txt"
merOcc=$seqFile"_"$merLen"mer_tallymer_occratio.txt"
#merSeq=$seqFile"_"$merLen"mer_seqs.fasta"

#-------------------------
# build the suffix array
#-------------------------
gt suffixerator -dna -pl -tis -suf -lcp -v -parts 4 -db $seq -indexname $db

#------------------------------------
# call occratio; get range of k-mers
#------------------------------------
#gt tallymer occratio -output unique nonunique -minmersize 10 -maxmersize 20 -esa reads

#------------------------------------------------
# call occratio; get ratio relative to the total
#------------------------------------------------
gt tallymer occratio -output unique relative -minmersize $merLen -maxmersize 180 -esa $db > $merOcc

#---------------------------------------------
# -mersizes option restricts the calculations
#---------------------------------------------
#gt tallymer occratio -output unique nonunique -mersizes 10 13 17 -esa reads

#------------------------------------------------------------------------------
# while occratio runs for a range of mer sizes, mkindex runs on a fixed length
#------------------------------------------------------------------------------
gt tallymer mkindex -mersize $merLen -minocc 40 -esa $db > $merCt

#-------------------------------------------------------------------------------------------------------------
# The program search now uses the index reads and matches all 19-mers of the input sequence U89959 against it
#-------------------------------------------------------------------------------------------------------------
#gt tallymer search -output qseqnum qpos counts sequence -tyr tyr-reads -q U89959.fna
