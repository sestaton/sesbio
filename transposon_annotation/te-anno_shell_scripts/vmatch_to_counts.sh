#!/bin/bash

perl vmatch_to_counts_v0.01.pl repbase_idlist.txt Ageratina_CAGATC_prinseq_trimmed_clean_desc_vmatch_IDs_ct.txt | sort -nrk 2 > Ageratina_superfamily_ct.tsv &
perl vmatch_to_counts_v0.01.pl repbase_idlist.txt Calyc_GATCAG_prinseq_trimmed_clean_desc_vmatch_IDs_ct.txt | sort -nrk 2 > Calyc_superfamily_ct.tsv &
perl vmatch_to_counts_v0.01.pl repbase_idlist.txt CP_AGTCAA_prinseq_trimmed_clean_desc_vmatch_IDs_ct.txt | sort -nrk 2 > CP_superfamily_ct.tsv &
perl vmatch_to_counts_v0.01.pl repbase_idlist.txt Dasyphllum_ATCACG_prinseq_trimmed_clean_desc_vmatch_IDs_ct.txt | sort -nrk 2 > Dasyphyllum_superfamily_ct.tsv &
perl vmatch_to_counts_v0.01.pl repbase_idlist.txt Gerb_GGCTAC_prinseq_trimmed_clean_desc_vmatch_IDs_ct.txt | sort -nrk 2 > Gerb_superfamily_ct.tsv &
perl vmatch_to_counts_v0.01.pl repbase_idlist.txt Gnaph_ACAGTG_prinseq_trimmed_clean_desc_vmatch_IDs_ct.txt | sort -nrk 2 > Gnaph_superfamily_ct.tsv &
perl vmatch_to_counts_v0.01.pl repbase_idlist.txt PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_vmatch_IDs_ct.txt | sort -nrk 2 > PI603989_superfamily_ct.tsv &
perl vmatch_to_counts_v0.01.pl repbase_idlist.txt Saff_GATCAG_prinseq_trimmed_clean_desc_vmatch_IDs_ct.txt | sort -nrk 2 > Saff_superfamily_ct.tsv &
perl vmatch_to_counts_v0.01.pl repbase_idlist.txt Senecio_TTAGGC_prinseq_trimmed_clean_desc_vmatch_IDs_ct.txt | sort -nrk 2 > Senecio_superfamily_ct.tsv &
perl vmatch_to_counts_v0.01.pl repbase_idlist.txt TKS_CAGATC_prinseq_trimmed_clean_desc_vmatch_IDs_ct.txt | sort -nrk 2 > TKS_superfamily_ct.tsv &
