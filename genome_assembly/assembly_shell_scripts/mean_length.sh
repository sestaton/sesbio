#!/bin/bash

function usage() {
cat <<EOF

USAGE: $0 <query_file> 

query_file    :   A Fasta/q file to use for calculating read lengths.

EOF
}

function print_error() {
cat <<ERR

ERROR: Command line not parsed correctly. Check input.

ERR
}

if [ $# -lt 1 ]; then
    print_error
    usage
    exit 1
fi

function timer() {
    if [[ $# -eq 0 ]]; then
        echo $(date '+%s')
    else
        local  stime=$1
        etime=$(date '+%s')

        if [[ -z "$stime" ]]; then stime=$etime; fi

        dt=$((etime - stime))
        ds=$((dt % 60))
        dm=$(((dt / 60) % 60))
        dh=$((dt / 3600))
        printf '%d:%02d:%02d' $ddh $dm $ds
    fi
}

~/apps/bioawk/bioawk -c fastx '{print length($seq)}' $1 | awk '{TOTAL+=$1} END{printf("COUNT:%d, TOTAL:%d, MEAN:%d\n",NR,TOTAL,TOTAL/NR)}'