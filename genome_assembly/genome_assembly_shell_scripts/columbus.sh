#!/bin/bash

cd `pwd`

/usr/local/velvet/latest/velveth tks_columbus_75 75 -reference -fasta ~/db/HA383_chloroplast_seq.fasta -shortPaired -sam TKS_HA383_bowtie2_sort.sam

/usr/local/velvet/latest/velvetg tks_columbus_75 -exp_cov auto -ins_length 400 -min_contig_lgth 200 -cov_cutoff 5