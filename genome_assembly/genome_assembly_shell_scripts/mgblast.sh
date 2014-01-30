#!/bin/bash
cd /scratch/statonse/fgclust_ha412ho
/home/jmblab/statonse/apps/gicl/mgblast -i ha412ho-short_IDs.fna.clean_nochloro_nomt_norRNA_unique.fasta -d ha412ho_clean -F "m D" -D 4 -p 85 -W18 -UT -X40 -KT -JF -v90000000 -b90000000 -C80 -H 320 -a 8 -o ha412ho_self_mgblast_out.txt

perl ~/ePerl/parse_mgblast.pl -i ha412ho_self_mgblast_out.txt -o ha412ho_self_mgblast_out_parsed.txt -id 90.00 -cov 0.15