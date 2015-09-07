#!/bin/bash

##NB: this is for human data with prebuilt database.

cd `pwd`

vep_dir=$HOME/apps/ensembl/ensembl-tools-release-75/scripts/variant_effect_predictor

perl $vep_dir/variant_effect_predictor.pl \
--offline \
--no_stats \
--everything \
--xref_refseq \
--check_existing \
--total_length \
--allele_number \
--no_escape \
--fork 2 \
--fasta \
~/.vep \
--input_file example.vcf \
--output_file example.vep.txt
