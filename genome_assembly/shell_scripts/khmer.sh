#!/bin/bash

script=$(basename $0)

function usage() {
cat <<EOF
USAGE: $script <seq> <mer>

seq   : fasta/q file to index and summarize/normalize for k-mer counts
mer   : k-mer length to use in querying <seq> data

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

## input/output
qrySeq=$1
merLen=$2
qryFile=$(echo ${qrySeq%.*})
basename=${qryFile}_${merLen}mer
khmerHash=${basename}_khmer-hash
khmerHist=${basename}_abund.hist
khmerLog=${basename}_khmer.log
khmerHashNorm=${basename}_khmer-hash_norm
khmerSeqNorm=${qrySeq}.keep
khmerHistNorm=${basename}_abund.hist.norm

#debug
#echo $qrySeq $merLen $qryFile $basename $khmerHash $khmerHist $khmerLog $khmerHashNorm $khmerSeqNorm $khmerHistNorm
#exit 0

## khmer PATH
python=/usr/local/python/2.7.2/bin/python
khmer_scripts=/home/jmblab/statonse/apps/khmer/scripts

## load-into-counting.py
$python $khmer_scripts/load-into-counting.py \
-N 4 \
-x 2e9 \
-k $merLen \
$khmerHash \
$qrySeq

## abundance-dist.py
$python $khmer_scripts/abundance-dist.py \
-z $khmerHash \
$qrySeq \
$khmerHist

## filter-abund.py
#$python $khmer_scripts/filter-abund.py \                     # filter k-mers below a threshold
#-C 1 \                                                       # theshold set to 1 in this case
#$khmerHash \                                                 # hash from which to read k-mer counts
#$qrySeq                                                      # sequence used to construct the hash

## normalize-by-median.py
$python $khmer_scripts/normalize-by-median.py \
-k $merLen \
-N 4 \
-x 2e9 \
-C 20 \
-l $khmerHash \
-s $khmerHashNorm \
-R $khmerLog \
$qrySeq

## filter-abund.py
$python $khmer_scripts/filter-abund.py \
-C 1 \
$khmerHashNorm \
$khmerSeqNorm

## abundance-dist.py
$python $khmer_scripts/abundance-dist.py \
$khmerHashNorm \
$khmerSeqNorm \
$khmerHistNorm