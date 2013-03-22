#!/bin/bash

khmer_scripts=/home/statonse/khmer_test/khmer/scripts

$khmer_scripts/load-into-counting.py \                                           # create the hash of k-mer counts
-N 4 \                                                                           # create this many hashes
-x 2e9 \                                                                         # create hash at least this size (set larger if you have >16 GB RAM)
-k 20 \                                                                          # k-mer size to hash
PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr_20mer_khmer-hash \    # hash name 
PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr.fasta                 # sequence name

$khmer_scripts/abundance-dist.py \                                               # create a histogram of k-mer coverage
-z PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr_20mer_khmer-hash \ # hash to read k-mer counts
PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr.fasta \               # sequence used to construct hash
PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr_abund.hist            # name of the histogram of counts

$khmer_scripts/filter-abund.py \                                                 # filter k-mers below a threshold
-C 1 \                                                                           # theshold set to 1 in this case
PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr_20mer_khmer-hash \    # hash from which to read k-mer counts
PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr.fasta                 # sequence used to construct the hash

$khmer_scripts/normalize-by-median.py \                                               # remove redundant k-mers from the graph
-k 17 \                                                                               # use this k-mer size
-N 4 \                                                                                # construct this many hashes
-x 2e9 \                                                                              # create hash of at least this size
-C 5 \                                                                                # coverage cutoff
-l PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr_20mer_khmer-hash \      # hash   
-s PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr_20mer_khmer-hash_norm \ # normalized hash
-R PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr_diginorm.txt \          # log of results
PI603989C03YUABXX_s1_prinseq_trimmed_clean_desc_paired_scr.fasta                      # input sequence
