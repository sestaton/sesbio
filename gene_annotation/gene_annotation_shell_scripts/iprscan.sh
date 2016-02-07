#!/bin/bash

set -euo pipefail

if [ $# -lt 1 ]; then
    echo -e "No input found. Exiting.\n"
    echo -e "USAGE: $0 seqs_prot.faa\n"
    exit 1
fi

iprscan=$HOME/apps/interproscan-5.16-55.0/interproscan.sh

time $iprscan -i $1 \
    --goterms \
    --iprlookup \
    -b canon_cds_pep_iprscan \
    --pathways UniPathway,KEGG,MetaCyc,Reactome \
    -f tsv,xml,gff3
