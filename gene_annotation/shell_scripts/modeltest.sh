#!/bin/bash

set -e
set -u
set -o pipefail

script=$(basename $0)

function usage() {
cat <<EOF

USAGE: $script <seq_file>  

seq_file   :   A (nucleotide) Fasta alignment file to test.

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

seq=$1

# "modeltest" is an alias on my system for the variable below; used for simplicity
modeltest="java -jar /usr/local/bioinfo/jmodeltest/jmodeltest-2.1.1/jModelTest.jar"
$modeltest -d $seq -g 4 -i -f -AIC -a 2> /dev/null | tail -1 | cut -f2
