#!/bin/bash

script=$(basename $0)

function usage() {
cat <<EOF

USAGE: $script <query_file> 

query_file    :   A FASTA/Q file to use for calculating basic statistics.

EOF
}

function print_error() {
cat <<ERR

ERROR: Command line not parsed correctly. Check input.

ERR
}

if [ ! $# == 1 ]; then
    print_error
    usage
    exit 1
fi

bioawk -c fastx '{SUM+=length($seq)} END {printf("Showing statistics for sequence file: "FILENAME"\nNumber of sequences: %d\nTotal length: %d\nMean length: %d\n",NR,SUM,SUM/NR)}' $1


