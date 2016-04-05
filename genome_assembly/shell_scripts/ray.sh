#!/bin/bash

cd `pwd`

mpiexec=/usr/local/openmpi/1.4.4/gcc412/bin/mpiexec
export LD_LIBRARY_PATH=/usr/local/openmpi/1.4.4/gcc412/lib:${LD_LIBRARY_PATH}

$mpiexec -n 24 ~/apps/ray/Ray -k17 \
-p PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr_1_p.fasta \
PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr_2_p.fasta \
-o k17-RayOutput