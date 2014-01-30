#!/bin/bash
cd /scratch/statonse/for_cleaned
perl ~/ePerl/filterSeqbyHit3.pl -i ha412ho-short_IDs.fna.clean_nochloro_unique.fasta -o ha412ho-short_IDs.fna.clean_nochloro_nomt_unique.fasta -b ha412_mt.bln -f fasta --removehits