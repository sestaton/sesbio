#!/bin/bash

script=$(basename $0)

function usage() {
cat <<EOF 

USAGE: $script forward reverse screendb

forward      :     The file of forward reads to screen
reverse      :     The file of reverse reads to screen
screendb     :     The FASTA database of sequences to filter out

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

fsamp=$1
rsamp=$2
db=$3
fbase=$(echo ${fsamp%.*})
rbase=$(echo ${rsamp%.*})
fsamp_scr=${fbase}_scr.fastq
rsamp_scr=${rbase}_scr.fastq
out=out_%.fq

time ~/apps/bbmap/bbsplit.sh -Xmx8g \
threads=2 \
in1=$fsamp \
in2=$rsamp \
ref=$db \
basename=$out \
outu1=$fsamp_scr \
outu2=$rsamp_scr
