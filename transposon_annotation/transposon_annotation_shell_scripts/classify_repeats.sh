#!/bin/bash

##NB: The idea here is to compare many databases along with
##    domain-finding methods for repeat consensi

#cd /scratch/statonse/for_RepeatScout/RM_20404.TueNov300822312010

# ncbi-latest
#blastall -p blastn \
#-d /db/ncbiblast-latest/nt \
#-i MID21_contig00001.fasta \
#-o ~/assemb_contigs/MID_21/MID21_BLASThits/MID21_contig00001.bln 

# Repbase
#blastall -p blastn \
#-e 1e-5 \
#-i consensi.fa.classified \
#-d ~/db/Repbase15.06 \
#-a 2 \
#-o CGP-repeats_RMod_Rebase1506.bln \
#-m 8

# classify with RepeatMasker
#RepeatMasker -lib consensi.fa.classified \
#-e wublast \
#CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.fasta

# translate for HMMER
#perl ~/ePerl/external/bp_translate_seq.pl \
#consensi.fa.classified \
#> consensi.fa.classified.faa

# HMMER3
#/usr/local/hmmer-latest/bin/hmmscan -o consensi.fa.classified_hmmer3.out \
#--cpu 4 \
#--tblout=CGP_BACs_RMout-hittable_pfam-a_out.txt \
#--domtblout=CGP_BACs_RMout-domain_table_pfam-a_out.txt \
#/db/pfam/latest/Pfam-A.hmm \
#consensi.fa.classified.faa

# parse HMMER3n
#perl ~/ePerl/parse_hmmer3.pl -i consensi.fa.classified_hmmer3.out \
#-o consensi.fa.classified_hmmer3.out-parsed.tab
