#!/bin/bash

cd screened_454_reads

ssaha2 -output ssaha2 -kmer 8 -skip 1 -seeds 1 -score 12 -cmatch 9 -ckmer 6 -454 1 pIndigo_BAC536.fasta MID17_in.fasta > MID17_ssaha2_vectorscreen_in.txt