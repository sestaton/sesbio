#!/bin/bash

set -euo pipefail

## NB: rmblastn is set to use 4 CPUs by default, so be aware if running
##     multiple threads
cd `pwd`

rm=$HOME/apps/maker/exe/RepeatMasker/RepeatMasker
lib=Ha412v1r1_complete_transposons_8-16.fas
genome=Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta

$rm -excln -pa 24 -no_is -gff -lib $lib $genome
