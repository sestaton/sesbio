#!/bin/bash

set -euo pipefail

genome=Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta
base=$(echo ${genome%.*})
hmms=$HOME/db/transposable+element_hmms/transposable+element.hmm
trnas=$HOME/db/plant_tRNAs.fasta
repdb=$HOME/db/RepBase1801.fasta

## LTRs
time tephra findltrs \
-g $genome \
-d $hmms \
-t $trnas \
-c tephra_ltr_config.yml \
-o ${base}_tephra_ltrs.gff3

time tephra classifyltrs \
-g $genome \
-d $repdb \
-t 12 \
-f ${base}_tephra_ltrs.gff3 \
-o ${base}_classified_ltrs

time tephra maskref \
-g $genome \
-d ${base}_classified_ltrs/${base}_combined_LTR_families.fasta \
-o ${base}_masked.fas

## solo-LTRs
time tephra sololtr -i ${base}_classified_ltrs/${base}_tephra_ltrs_copia \
-g ${base}_masked.fas \
-o ${base}_masked_copia_sololtrs.gff3 \
-r ${base}_masked_copia_sololtr_rep.tsv \
-s ${base}_masked_copia_sololtr_seqs.fas

time tephra sololtr -i ${base}_classified_ltrs/${base}_tephra_ltrs_gypsy \
-g ${base}_masked.fas \
-o ${base}_masked_copia_sololtrs.gff3 \
-r ${base}_masked_gypsy_sololtr_rep.tsv \
-s ${base}_masked_gypsy_sololtr_seqs.fas 

## ltrage
time tephra ltrage \
-g ${base}_masked.fas \
-t 12 \
-o ${base}_gypsy_ltrages \
-f ${base}_tephra_ltrs_gypsy.gff3 \
--clean

time tephra ltrage \
-g ${base}_masked.fas \
-t 12 \
-o ${base}_copia_ltrages \
-f ${base}_tephra_ltrs_gypsy.gff3 \ \
--clean

## illrecomb
time tephra illrecomb -i ${base}_classified_ltrs/${base}_tephra_ltrs_copia \
-o ${base}_masked_copia_illrecomb.fas \
-r ${base}_masked_copia_illrecomb_rep.tsv \
-s ${base}_masked_copia_illrecomb_stats.tsv \
-t 12

time tephra illrecomb -i ${base}_classified_ltrs/${base}_tephra_ltrs_gypsy \
-o ${base}_masked_gypsy_illrecomb.fas \
-r ${base}_masked_gypsy_illrecomb_rep.tsv \
-s ${base}_masked_gypsy_illrecomb_stats.tsv \
-t 12

## TRIMs
time tephra findtrims \
-g ${base}_masked.fas \
-d $hmms \
-t $trnas

time tephra maskref \
-g ${base}_masked.fas \
-d ${base}_masked_trim_ltrdigest85_combined_filtered.fasta \
-o ${base}_masked2.fas

## Helitrons
time tephra findhelitrons \
-g ${base}_masked2.fas \
-o ${base}_masked2_helitrons.gff3

time tephra maskref \
-g ${base}_masked2.fas \
-d ${base}_masked2_tephra_hscan_helitrons.hel.fa \
-o ${base}_masked3.fas

## TIR elements
time tephra findtirs \
-g ${base}_masked3.fas \
-d $hmms \
-o ${base}_masked3_tirs.gff3
 
time tephra classifytirs \
-g ${base}_masked3.fas \
-f ${base}_masked3_tirs_filtered.gff3

time tephra maskref \
-d ${base}_masked3_tirs.fasta \
-g ${base}_masked3.fas \
-o ${base}_masked4.fas

## non-LTRs
time tephra findnonltrs \
-g ${base}_masked4.fas

time tephra maskref \
-g ${base}_masked4.fas \
-d nonLTRs_out/nonLTRs_out_tephra_nonltr.fasta \
-o ${base}_masked5.fas

