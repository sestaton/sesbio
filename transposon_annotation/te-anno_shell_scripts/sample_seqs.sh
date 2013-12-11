#!/bin/bash

cd `pwd`

time perl -MTransposome::SeqUtil -e 'Transposome::SeqUtil->new( file => shift, sample_size => 500_000, no_store => 1)->sample_seq' \
TKS_CAGATC_prinseq_trimmed_clean_desc_paired_scr_interl.fasta \
> zsample_500k.fasta