#!/bin/bash
cd /scratch/statonse/for_cleaning
#blastall -p blastn -d /db/ncbiblast-latest/nt -i MID21_contig00001.fasta -o ~/assemb_contigs/MID_21/MID21_BLASThits/MID21_contig00001.bln 
blastall -p blastn -e 1e-10 -i ha412ho-short_IDs.fna.clean_nochloro_unique.fasta -d ~/db/At_mtgenome -a 2 -o ha412_mt.bln -m 8
