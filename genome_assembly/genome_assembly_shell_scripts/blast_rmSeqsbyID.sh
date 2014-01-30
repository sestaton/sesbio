#!/bin/bash
cd /scratch/statonse/for_cleaning
#blastall -p blastn -d /db/ncbiblast-latest/nt -i MID21_contig00001.fasta -o ~/assemb_contigs/MID_21/MID21_BLASThits/MID21_contig00001.bln 
#blastall -p blastn -e 1e-10 -i ha412ho-short_IDs.fna.clean_nochloro_nomt_unique.fasta -d ~/db/LSU_SSU_rRNA_Ath -a 2 -o ha412_rRNAs.bln -m 8

perl ~/ePerl/filterSeqbyHit3.pl -i ha412ho-short_IDs.fna.clean_nochloro_nomt_unique.fasta -o ha412ho-short_IDs.fna.clean_nochloro_nomt_norRNA_unique.fasta -b ha412_rRNAs.bln -f fasta --removehits

