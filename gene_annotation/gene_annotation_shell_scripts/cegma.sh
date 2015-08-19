#!/bin/bash

export PATH=/usr/local/bioinfo/cegma/CEGMA_v2/wise2.2.3-rc7/src/bin:$PATH
export CEGMA=/usr/local/bioinfo/cegma/CEGMA_v2
export PATH=$PATH:$CEGMA/bin
export PATH=$PATH:$CEGMA/geneid/bin
export CEGMATMP=/tmp
export PERL5LIB="$PERL5LIB:$CEGMA/lib"

cd `pwd`

genome=Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta
base=$(echo ${genome%.*})
cegma_out=${base}_cegma

cegma -v --ext -o $cegma_out --genome $genome -threads 12
