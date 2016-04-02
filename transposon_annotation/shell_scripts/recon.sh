#!/bin/bash

set -euo pipefail

#NB: RECON with generate a segmenation fault if there are periods (".") in
#    the sequence names.
#    RECON will also generate a segmentation fault with more than ~10k sequences,
#    and it is very slow. Usage not recommended, but below is a guide on how to make it faster
#    if the usage is necessary.
#
cd `pwd`

# link to script: https://github.com/sestaton/sesbio/blob/master/gene_annotation/parallel_blast.pl
perl ~/github/sesbio/gene_annotation/parallel_blast.pl \
-i SRR486236_trimmed_filteredall_10k_interl.fasta \
-p blastn \
-d maize_10k_blastdb \
-n 1000 \
-t 10 \
-o SRR486236_trimmed_filteredall_10k_interl_allvall.bln

# link to script: https://github.com/sestaton/sesbio/blob/master/transposon_annotation/blast2MSP.pl
perl ~/github/sesbio/transposon_annotation/blast2MSP.pl \
SRR486236_trimmed_filteredall_10k_interl_allvall.bln \
> 10k_MSP.txt

script=RECON-1.08/scripts/recon.pl

time perl $script \
10k_seqnames.txt \
10k_MSP.txt \
1