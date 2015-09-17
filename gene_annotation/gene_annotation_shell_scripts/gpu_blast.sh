#!/bin/bash

##NB: This is mostly for understanding how to get gpu_blast working. It does
##    not actually provide any real advantage in speed up (if you have multiple CPUs),
##    but this may help for testing.

## path to cuda libs
export LD_LIBRARY_PATH=/usr/local/cuda/4.0.17/cuda/lib64:${LD_LIBRARY_PATH}

cd `pwd`

## use the gpu_blast executables
## NB: 'latest' is a symlink: /usr/local/gpu-blast/1.1/cuda-4.0.17/bin/blastp
makeblastdb=/usr/local/gpu-blast/latest/bin/makeblastdb
blastp=/usr/local/gpu-blast/latest/bin/blastp

$makeblastdb \
-in goodProteins.fasta \
-out goodProteins_gpu \
-sort_volumes \
-max_file_sz 500MB \
-dbtype prot

#$blastp \
#-db goodProteins_gpu \
#-query goodProteins.fasta \
#-gpu t \
#-method 2 \
#-gpu_blocks 256 \
#-gpu_threads 32

$blastp \
-db goodProteins_gpu \
-query goodProteins.fasta \
-gpu t \
-num_threads 12 \
-outfmt 6 \
-num_descriptions 100000 \
-num_alignments 100000 \
-seg yes \
-out goodProteins_allvall_gpu.bln


