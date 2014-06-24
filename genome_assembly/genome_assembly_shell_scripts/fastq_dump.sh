#!/bin/bash

## Download the SRA file for a Zea mays Illumina paired-end run,
## and split the archive into separte fastq files.

cd `pwd`

curl -O ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX142/SRX142106/SRR486236/SRR486236.sra

/usr/local/sra/latest/bin/fastq-dump -A SRR486236 --split-3 SRR486236.sra