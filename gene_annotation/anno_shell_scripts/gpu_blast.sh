#!/bin/bash

export LD_LIBRARY_PATH=/usr/local/cuda/4.0.17/cuda/lib64:${LD_LIBRARY_PATH}

cd /escratch/statonse_Feb_07

#/usr/local/gpu-blast/latest/bin/makeblastdb -in goodProteins.fasta -out goodProteins_gpu -sort_volumes -max_file_sz 500MB -dbtype prot

#/usr/local/gpu-blast/1.1/cuda-4.0.17/bin/blastp -db goodProteins_gpu -query goodProteins.fasta -gpu t -method 2 -gpu_blocks 256 -gpu_threads 32

/usr/local/gpu-blast/1.1/cuda-4.0.17/bin/blastp -db goodProteins_gpu -query goodProteins.fasta -gpu t -num_threads 12 -outfmt 6 -num_descriptions 100000 -num_alignments 100000 -seg yes -out goodProteins_allvall_gpu.bln


