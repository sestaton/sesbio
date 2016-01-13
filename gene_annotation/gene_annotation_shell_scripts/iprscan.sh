#!/bin/bash

cd `pwd`

iprscan=/home/statonse/apps/iprscan/interproscan-5.14-53.0/interproscan.sh

$iprscan -i seqs.fas \
	 --goterms \
	 --iprlookup \
	 -b seqs_iprscan \
	 --pathways UniPathway,KEGG,MetaCyc,Reactome \
	 -f tsv,xml,gff3 > ipr.out

ipr_update_gff genome.all.gff ipr.out > genome.all.ipr.gff
