#!/bin/bash

set -euo pipefail

cd `pwd`

stringtie=$HOME/github/stringtie/stringtie
gtf=$1

for file in ./*bam
do
    base=$(echo ${file%.*})
    stringassem=${base}_stringtie_assembly.gtf
    stringabund=${base}_stringtie_gene_abund.txt
    ballgowndir=${base}_balldown_tables

    $stringtie $file -G $gtf -p 12 -o $stringassem -A $stringabund -b $ballgowndir -v
done
