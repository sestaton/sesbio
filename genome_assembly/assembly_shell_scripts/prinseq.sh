nn#!/bin/bash

#cd /iob_scratch/statonse/RHA801_assemb/prinseq_trimmed

perl ~/apps/prinseq-lite-0.17.1/prinseq-lite.pl -fastq SRR350226_1.fastq -out_format 1 -out_good SRR350226_1_prinseq_trimmed -log SRR350226_1_prinseq.log -min_len 40 -noniupac -min_qual_mean 15 -lc_method entropy -lc_threshold 60 -trim_ns_right 10 -ns_max_p 20

#perl ~/apps/prinseq-lite-0.15/prinseq-lite.pl -fastq s_7_2_sequence.fastq -stats_info -stats_len
#perl ~/apps/prinseq-lite-0.15/prinseq-lite.pl -fastq s_7_2_sequence.fastq -stats_info -stats_len