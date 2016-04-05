#!/bin/bash

cd `pwd`

gt=/home/jmblab/statonse/apps/readjoiner/readjoiner-1.1/bin/gt

$gt readjoiner prefilter -readset safflower \
-db s_1_1_prinseq_trimmed_filtered_orphaned_reads.fasta \
s_1_2_prinseq_trimmed_filtered_orphaned_reads.fasta \
s_1_interleaved.fasta \
s_2_1_prinseq_trimmed_filtered_orphaned_reads.fasta \
s_2_2_prinseq_trimmed_filtered_orphaned_reads.fasta \
s_2_interleaved.fasta \
s_3_1_prinseq_trimmed_filtered_orphaned_reads.fasta \
s_3_2_prinseq_trimmed_filtered_orphaned_reads.fasta \
s_3_interleaved.fasta \
s_4_1_prinseq_trimmed_filtered_orphaned_reads.fasta \
s_4_2_prinseq_trimmed_filtered_orphaned_reads.fasta \
s_4_interleaved.fasta \
s_5_1_prinseq_trimmed_filtered_orphaned_reads.fasta \
s_5_2_prinseq_trimmed_filtered_orphaned_reads.fasta \
s_5_interleaved.fasta

$gt readjoiner overlap -readset safflower -l 35

$gt readjoiner assembly -readset safflower -l 35