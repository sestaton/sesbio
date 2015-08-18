#!/bin/bash

export CEGMA=/usr/local/bioinfo/cegma/CEGMA_v2
export PATH=$PATH:$CEGMA/bin
export PATH=$PATH:$CEGMA/geneid/bin
export CEGMATMP=/tmp
export PERL5LIB="$PERL5LIB:$CEGMA/lib"

cd `pwd`

genome=Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta

cegma --genome $genome -threads 12
