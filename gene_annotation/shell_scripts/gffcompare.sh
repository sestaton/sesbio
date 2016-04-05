#!/bin/bash

set -euo pipefail

cd `pwd`

gffcompare=$HOME/github/gffcompare/gffcompare

$gffcompare -r Ha412v1r1_genes.gtf \
-s $HOME/db/Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta \
-i stringtie_refgtf_list.txt

