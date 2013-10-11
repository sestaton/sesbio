#!/bin/bash

cd `pwd`

ulimit -v 500000000

## Senecio
perl repeat_explorer_v0.13.pl \
-i Senecio_TTAGGC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast.bln \
-o Senecio_repeat_explorer_results_500k \
-id 90.00 \
-cov 0.55 \
-f Senecio_TTAGGC_prinseq_trimmed_clean_desc_paired_scr_500k_interl.fasta \
-r Senecio_TTAGGC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1.txt \
-m 500 \
-d ~/db/RepBase1801_Hannuus_LTR-RT_families \
-j ~/db/repbase1801_full.json \
--in_memory

## Ageratina
perl repeat_explorer_v0.13.pl \
-i Ageratina_CAGATC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast.bln \
-o Ageratina_repeat_explorer_results_500k \
-id 90.00 \
-cov 0.55 \
-f Ageratina_CAGATC_prinseq_trimmed_clean_desc_paired_scr_500k_interl.fasta \
-r Ageratina_CAGATC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1.txt \
-m 500 \
-d ~/db/RepBase1801_Hannuus_LTR-RT_families \
-j ~/db/repbase1801_full.json \
--in_memory

## Calyc
perl repeat_explorer_v0.13.pl \
-i Calyc_GATCAG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast.bln \
-o Calyc_repeat_explorer_results_500k \
-id 90.00 \
-cov 0.55 \
-f Calyc_GATCAG_prinseq_trimmed_clean_desc_paired_scr_500k_interl.fasta \
-r Calyc_GATCAG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1.txt \
-m 500 \
-d ~/db/RepBase1801_Hannuus_LTR-RT_families \
-j ~/db/repbase1801_full.json \
--in_memory

## CP
perl repeat_explorer_v0.13.pl \
-i CP_AGTCAA_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast.bln \
-o CP_repeat_explorer_results_500k \
-id 90.00 \
-cov 0.55 \
-f CP_AGTCAA_prinseq_trimmed_clean_desc_paired_scr_500k_interl.fasta \
-r CP_AGTCAA_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1.txt \
-m 500 \
-d ~/db/RepBase1801_Hannuus_LTR-RT_families \
-j ~/db/repbase1801_full.json \
--in_memory

## Dasyphyllum
perl repeat_explorer_v0.13.pl \
-i Dasyphllum_ATCACG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast.bln \
-o Dasyphyllum_repeat_explorer_results_500k \
-id 90.00 \
-cov 0.55 \
-f Dasyphllum_ATCACG_prinseq_trimmed_clean_desc_paired_scr_500k_interl.fasta \
-r Dasyphllum_ATCACG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1.txt \
-m 500 \
-d ~/db/RepBase1801_Hannuus_LTR-RT_families \
-j ~/db/repbase1801_full.json \
--in_memory

## Gerb
perl repeat_explorer_v0.13.pl \
-i Gerb_GGCTAC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast.bln \
-o Gerb_repeat_explorer_results_500k \
-id 90.00 \
-cov 0.55 \
-f Gerb_GGCTAC_prinseq_trimmed_clean_desc_paired_scr_500k_interl.fasta \
-r Gerb_GGCTAC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1.txt \
-m 500 \
-d ~/db/RepBase1801_Hannuus_LTR-RT_families \
-j ~/db/repbase1801_full.json \
--in_memory

## Gnaph
perl repeat_explorer_v0.13.pl \
-i Gnaph_ACAGTG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast.bln \
-o Gnaph_repeat_explorer_results_500k \
-id 90.00 \
-cov 0.55 \
-f Gnaph_ACAGTG_prinseq_trimmed_clean_desc_paired_scr_500k_interl.fasta \ 
-r Gnaph_ACAGTG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1.txt \
-m 500 \
-d ~/db/RepBase1801_Hannuus_LTR-RT_families \
-j ~/db/repbase1801_full.json \
--in_memory

## Saff
perl repeat_explorer_v0.13.pl \
-i Saff_GATCAG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast.bln \
-o Saff_repeat_explorer_results_500k \
-id 90.00 \
-cov 0.55 \
-f Saff_GATCAG_prinseq_trimmed_clean_desc_paired_scr_500k_interl.fasta \
-r Saff_GATCAG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1.txt \
-m 500 \
-d ~/db/RepBase1801_Hannuus_LTR-RT_families \
-j ~/db/repbase1801_full.json \
--in_memory

## TKS
perl repeat_explorer_v0.13.pl \
-i TKS_CAGATC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast.bln \
-o TKS_repeat_explorer_results_500k \
-id 90.00 \
-cov 0.55 \
-f TKS_CAGATC_prinseq_trimmed_clean_desc_paired_scr_500k_interl.fasta \
-r TKS_CAGATC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1.txt \
-m 500 \
-d ~/db/RepBase1801_Hannuus_LTR-RT_families \
-j ~/db/repbase1801_full.json \
--in_memory


