#!/bin/bash

# create index
~/apps/bwa/bwa index HA383_chloroplast_seq.fasta

# aln reads to the index
~/apps/bwa/bwa mem -M -t 8 -v 1 -a -p HA383_chloroplast_seq.fasta \
Dasyphyllum_cpDNA_seqs_interl.fasta | sed 's/^\[.*//g;/^ *$/d' \
> Dasyphyllum_cpDNA_seqs_interl_HA383cp_bwa-mem_aln.sam

# add header to sam
samtools view -bT HA383_chloroplast_seq.fasta \
Dasyphyllum_cpDNA_seqs_interl_HA383cp_bwa-mem_aln.sam \
> Dasyphyllum_cpDNA_seqs_interl_HA383cp_bwa-mem_aln.bam

# run calmd to add MD tag to sam
samtools calmd -bS Dasyphyllum_cpDNA_seqs_interl_HA383cp_bwa-mem_aln.sam \
HA383_chloroplast_seq.fasta \
> Dasyphyllum_cpDNA_seqs_interl_HA383cp_bwa-mem_aln_md.bam

# run samstats to generate alignment report
~/apps/samstat/src/samstat Dasyphyllum_cpDNA_seqs_interl_HA383cp_bwa-mem_aln_md.bam