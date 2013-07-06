#!/bin/bash

cd `pwd`

abacus=/usr/local/bioinfo/pagit/PAGIT/ABACUS/abacus.pl
image=/usr/local/bioinfo/pagit/PAGIT/IMAGE/image.pl 
image_sum=/usr/local/bioinfo/pagit/PAGIT/IMAGE/image_run_summary.pl

## run abacus.pl
perl $abacus -r ref.seq \
-q assembly.fas \
-p nucmer \
-m -b -o species_target_abacus

## run image.pl
perl $image -scaffolds Calyc_chloroplast_genome_contigs_t.fasta \
-prefix Calyc_reads \
-iteration 1 \
-all_iteration 9 \
-dir_prefix ite \
-kmer 55 > image.out

## run image_run_summary.pl
perl $image_sum ite >> image.out
