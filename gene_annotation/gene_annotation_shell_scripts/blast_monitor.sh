#!/bin/bash

function usage() {
cat <<EOF

USAGE: $0 <blast> <seq_file> 

blast      :   A tab-delimited BLAST report (only tested with blastall style).
seq_file   :   A (nucleotide) Fasta file to analyze.

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

blast=$1
query=$2

if [ ! -e $blast ]; then
    echo "ERROR: Could not find file: $blast"
    usage
    exit 1
fi

if [ ! -e $query ]; then
    echo "ERROR: Could not find file: $query"
    usage
    exit 1
fi

echo "The BLAST report is: $blast"
echo "The Fasta query is : $query"
echo ""
curquery=$(tail -2 $blast | head -1 | cut -f 1)
curline=$(grep ">" $query | grep -n $curquery | cut -f 1 -d ':')
nblines=$(grep -c ">" $query)
percent=$(echo "($curline/$nblines) *100" | bc -l | cut -c 1-4)
echo "The BLAST job is $percent% done..."
