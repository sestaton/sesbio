#!/bin/bash

script=$(basename $0)

function usage() {
cat <<EOF
USAGE: $script <query> <target> <output>

query   : A file of ESTs to query against a set of genomic sequences.
target  : A set of genomic DNA sequences to be used as the target.
output  : A file to place the alignment results. The format is (tab-delimited):

query_id target_id perc_id query_len target_len num_mismatches query_aln_beg query_aln_end target_aln_beg target_aln_end aln_score 

EOF
}

function print_error() {
cat <<ERR

ERROR: Command line not parsed correctly. Check input.

ERR
}

if [ $# -lt 3 ]; then
    print_error
    usage
    exit 1
fi

query=$1
target=$2
output=$3

exonerate=/usr/local/exonerate/latest/bin/exonerate
$exonerate --ryo "%qi\t%ti\t%pi\t%qal\t%ql\t%tl\t%em\t%qab\t%qae\t%tab\t%tae\t%s\n" \
--verbose 0 \
--showvulgar 0 \
--showalignment 0 \
--model est2genome \
$query $target > $output 
