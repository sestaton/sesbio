#!/bin/bash

set -euo pipefail

cd `pwd`

samtools=$HOME/github/samtools/samtools
bed=eps.bed
dir=$1
out=$2

mkdir $out

for bam in ./$dir/*bam
do
    file=$(basename $bam)
    bamfile=$(echo ${file%.*})
    covfile=${bamfile}_cov.tsv
    $samtools depth -b $bed $bam > $out/$covfile
done

perl sum_eps_cov.pl $out $out/all_cov_summaries.tsv

echo "coverage est. done";
