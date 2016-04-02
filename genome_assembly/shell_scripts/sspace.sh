#!/bin/bash

cd `pwd`

perl ~/apps/SSPACE-1.2_linux-x86_64/SSPACE_v1-2.pl \
-l libraries.txt \
-s Ray_31mer_full_RHA801.Contigs_over200_nocp.fasta \
-k 5 \
-a 0.7 \
-x 1 \
-m 30 \
-o 20 \
-b Rayk31_scaffolds_extension