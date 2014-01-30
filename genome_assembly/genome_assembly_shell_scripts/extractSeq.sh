#!/bin/bash
cd /scratch/statonse/tmp_fgclust
perl ~/ePerl/extractSeq.pl -i ha412ho-short_IDs.fna.clean_nochloro_unique.fasta.1.3 -o ha412_1_3_CL3.fasta -b CL3_idlist.txt -f fasta --extracthits 