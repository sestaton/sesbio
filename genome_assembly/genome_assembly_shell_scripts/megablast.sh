#!/bin/bash

## NB: This file was for testing fgclust against Transposome.

cd /scratch/statonse/fgclust_ha412ho

megablast -i ha412ho-short_IDs.fna.clean_nochloro_unique.fasta \
-d ha412ho_unique_clean \
-F "m D" \
-D 2 \
-p 85 \
-W18 \
-UT \
-X40 \
-JF \
-v 200000000 \
-b 0 \
-a 4 \
-o ha412_self_mgblast_out.txt
