#!/bin/bash

set -e
set -u
set -o pipefail

script=$(basename $0)

function usage() {
cat <<EOF

USAGE: $script <blast> <seq_file> 

blast      :   A tab-delimited BLAST report (legacy "-m 8" or blast+ "-outfmt 6").
seq_file   :   The FASTA query file.

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
echo "The FASTA query is : $query"
echo ""
curquery=$(tail -2 $blast | head -1 | cut -f 1)
curline=$(grep ">" $query | grep -Fwn $curquery | cut -f 1 -d ':')
nblines=$(grep -c ">" $query)
percent=$(echo "($curline/$nblines) *100" | bc -l | cut -c 1-4)
echo "The BLAST job is $percent% done..."
