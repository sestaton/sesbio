#!/bin/bash

set -euo pipefail

gt=$HOME/.tephra/gt/bin/gt
gff=$1
base=$(echo ${gff%.*})

gff_sort=${base}_sort.gff3
glen_stat=${base}_sort_genelengths.txt
elen_stat=${base}_sort_exonlengths.txt
enum_stat=${base}_sort_exon-nums.txt
ilen_stat=${base}_sort_intronlengths.txt
clen_stat=${base}_sort_cdslengths.txt

# some identifying prefix for the output
prefix=bronze-maker-round3

time $gt gff3 -sort -tidy $gff > $gff_sort

time $gt stat -genelengthdistri -addintrons -o ${prefix}-genelengthdist.txt $gff_sort > $glen_stat
time $gt stat -exonlengthdistri -o ${prefix}-exonlengthdist.txt $gff_sort > $elen_stat
time $gt stat -exonnumberdistri -o ${prefix}-exonnumberdist.txt $gff_sort > $enum_stat
time $gt stat -intronlengthdistri -addintrons -o ${prefix}-intronlengthdist.txt $gff_sort > $ilen_stat
time $gt stat -cdslengthdistri -o ${prefix}-cdslengthdist.txt $gff_sort > $clen_stat
