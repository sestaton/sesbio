#!/bin/bash

## This script is for testing my own version VelvetOptimiser
## https://github.com/sestaton/VelvetOptimiser.git

export OMP_THREAD_LIMIT=10
export OMP_NUM_THREADS=2

cd `pwd`

perl ~/github/VelvetOptimiser/VelvetOptimiser.pl -s 77 \
-e 89 \
-f '-fasta -shortPaired TKS_500k_interl.fasta' \
-t 8 \
--optFuncKmer 'n50' \
-p VelvetOpt_k77-89_all \
-d VelvetOpt_k77-89_all