#!/bin/bash

set -euo pipefail

script=$(basename $0)

function usage() {
cat <<EOF

USAGE: $script <gff>
gff   :  A GFF file of gene predictions for training snap


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

maker2zff=$HOME/apps/maker/bin/maker2zff
snapbin=$HOME/apps/maker/exe/snap

gff=$1
$maker2zff $gff
$snapbin/fathom -categorize 1000 genome.ann genome.dna
$snapbin/fathom -export 1000 -plus uni.ann uni.dna
$snapbin/forge export.ann export.dna
$snapbin/hmm-assembler.pl ha412 . > ha412.hmm
