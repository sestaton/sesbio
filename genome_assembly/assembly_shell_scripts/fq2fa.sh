#!/bin/bash

perl ~/ePerl/fq2fas.pl -i HA0001_61U0AAAXX_1_1_concat_qseq.fastq -o HA0001_61U0AAAXX_1_1_concat_qseq.fasta 
perl ~/ePerl/fq2fas.pl -i HA0001_61U0AAAXX_1_2_concat_qseq.fastq -o HA0001_61U0AAAXX_1_2_concat_qseq.fasta
perl ~/ePerl/fq2fas.pl -i HA0001_61U0AAAXX_2_1_concat_qseq.fastq -o HA0001_61U0AAAXX_2_1_concat_qseq.fasta
perl ~/ePerl/fq2fas.pl -i HA0001_61U0AAAXX_2_2_concat_qseq.fastq -o HA0001_61U0AAAXX_2_2_concat_qseq.fasta

