#!/bin/bash

## NB: This script deprecated; it preceeded the existence of Transposome (https://github.com/sestaton/Transposome),
##     which will generate genome coverage summaries from clustering reports automatically. Use Transposome,
##     it's far more convenient and accurate. 10/19/13 SES


## Ageratina
perl annotation_summary_to_gcov.pl -i \
Ageratina_CAGATC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1_annotations_summary.tsv \
-sfo Ageratina_500k_superfamily_summary.tsv \
-fo Ageratina_500k_family_summary.tsv

## Ann1238
perl annotation_summary_to_gcov.pl -i \
Ann1238_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-3_annotations_summary.tsv \
-sfo Ann1238_500k_superfamily_summary.tsv \
-fo Ann1238_500k_family_summary.tsv

## Calyc
perl annotation_summary_to_gcov.pl -i \
Calyc_GATCAG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1_annotations_summary.tsv \
-sfo Calyc_500k_superfamily_summary.tsv \
-fo Calyc_500k_family_summary.tsv

## CP
perl annotation_summary_to_gcov.pl -i \
CP_AGTCAA_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1_annotations_summary.tsv \
-sfo CP_500k_superfamily_summary.tsv \
-fo CP_500k_family_summary.tsv

## Dasy
perl annotation_summary_to_gcov.pl -i \
Dasyphllum_ATCACG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1_annotations_summary.tsv \
-sfo Dasyphyllum_500k_superfamily_summary.tsv \
-fo Dasyphyllum_500k_family_summary.tsv

## Gerb
perl annotation_summary_to_gcov.pl -i \
Gerb_GGCTAC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1_annotations_summary.tsv \
-sfo Gerb_500k_superfamily_summary.tsv \
-fo Gerb_500k_family_summary.tsv

## Gnaph
perl annotation_summary_to_gcov.pl -i \
Gnaph_ACAGTG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-3_annotations_summary.tsv \
-sfo Gnaph_500k_superfamily_summary.tsv \
-fo Gnaph_500k_family_summary.tsv

## Harg
perl annotation_summary_to_gcov.pl -i \
Harg_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-3_annotations_summary.tsv \
-sfo Harg_500k_superfamily_summary.tsv \
-fo Harg_500k_family_summary.tsv

## Hport
perl annotation_summary_to_gcov.pl -i \
Hport_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-3_annotations_summary.tsv \
-sfo Hport_500k_superfamily_summary.tsv \
-fo Hport_500k_family_summary.tsv 

## Hteph
perl annotation_summary_to_gcov.pl -i \
Hteph_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-3_annotations_summary.tsv \
-sfo Hteph_500k_superfamily_summary.tsv \
-fo Hteph_500k_family_summary.tsv

## Hvert
perl annotation_summary_to_gcov.pl -i \
Hvert_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-3_annotations_summary.tsv \
-sfo Hvert_500k_superfamily_summary.tsv \
-fo Hvert_500k_family_summary.tsv

## Phoeb
perl annotation_summary_to_gcov.pl -i \
Phoeb_330k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-3_annotations_summary.tsv \
-sfo Phoeb_330k_superfamily_summary.tsv \
-fo Phoeb_330k_family_summary.tsv

## Saff
perl annotation_summary_to_gcov.pl -i \
Saff_GATCAG_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1_annotations_summary.tsv \
-sfo Saff_500k_superfamily_summary.tsv \
-fo Saff_500k_family_summary.tsv

## Senecio
perl annotation_summary_to_gcov.pl -i \
Senecio_TTAGGC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1_annotations_summary.tsv \
-sfo Senecio_500k_superfamily_summary.tsv \
-fo Senecio_500k_family_summary.tsv

## TKS
perl annotation_summary_to_gcov.pl -i \
TKS_CAGATC_prinseq_trimmed_clean_desc_paired_scr_500k_interl_mgblast_90PID_55PCOV_merged_cluster_report_4-1_annotations_summary.tsv \
-sfo TKS_500k_superfamily_summary.tsv \
-fo TKS_500k_family_summary.tsv

## mk outdirs
mkdir all_family_summaries all_superfamily_summaries
cp *_family_summary.tsv all_family_summaries
cp *_superfamily_summary.tsv all_superfamily_summaries
