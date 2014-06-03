#!/bin/bash

function usage() {
cat <<EOF
USAGE: $0 <fastq_file>

fastq_file   : file or reads to filter
         
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

java -Xmx1024m -classpath /home/jmblab/statonse/apps/FastQC uk.ac.bbsrc.babraham.FastQC.FastQCApplication $1