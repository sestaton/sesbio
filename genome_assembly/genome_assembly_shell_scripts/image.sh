#!/bin/bash

## NB: The abacus command below is run prior to running image,
##     though the command you see is for a different run (which is
##     why filenames don't match). The command is there for reference.

cd `pwd`

abacus=/usr/local/bioinfo/pagit/PAGIT/ABACUS/abacus.pl
image=/usr/local/bioinfo/pagit/PAGIT/IMAGE/image.pl 
image_sum=/usr/local/bioinfo/pagit/PAGIT/IMAGE/image_run_summary.pl
contigs2scaffs=/usr/local/bioinfo/pagit/PAGIT/IMAGE/contigs2scaffolds.pl

## run abacus.pl
#perl $abacus -r ref.seq \
#-q assembly.fas \
#-p nucmer \
#-m -b -o species_target_abacus

## run image.pl
perl $image -scaffolds Ager_HA383cp_abacus.MULTIFASTA.fa \
-prefix Ager_reads \
-iteration 1 \
-all_iteration 9 \
-dir_prefix ite \
-kmer 55 

## run image_run_summary.pl
perl $image_sum ite 

## convert contigs from final iteration to scaffolds
cd ite9
perl $contigs2scaffs new.fa new.read.placed 300 0 scaffolds
