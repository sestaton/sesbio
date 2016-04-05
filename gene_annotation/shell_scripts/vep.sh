#!/bin/bash

set -euo pipefail

VEP=/home/statonse/apps/ensembl-tools-release-84/scripts/variant_effect_predictor 
vep_script=${VEP}/variant_effect_predictor.pl

## build custom database
perl -I$VEP $vep_script \
-i /moonriseNFS/HA412/annotation/ubc_annotation/genes/Ha412v1r1_genes.gff3.gz \
-f ~/db/Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta \
-d 84 \
-s sunflower 2>&1 > gff2vep.out 

## run vep
perl -I$VEP $vep_script -offline \
-t Ha412v1r1_genome_no_cp-mt-rd_chr-q_gatk_rnaseq_vars/Ha412v1r1_genome_no_cp-mt-rd_chr-q_gatk_rnaseq_vars_merged_vars.vcf \
--species sunflower \
-o Ha412v1r1_genome_no_cp-mt-rd_chr-q_gatk_rnaseq_vars/Ha412v1r1_genome_no_cp-mt-rd_chr-q_gatk_rnaseq_vars_merged_vars_vep.txt 2>&1 > vep_gbautes_rnaseq.out
