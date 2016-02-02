#!/bin/bash

set -euo pipefail

cd `pwd`

iprscan=$HOME/apps/interproscan-5.16-55.0/interproscan.sh

time $iprscan -i $1 \
    --goterms \
    --iprlookup \
    -b canon_cds_pep_iprscan \
    --pathways UniPathway,KEGG,MetaCyc,Reactome \
    -f tsv,xml,gff3
