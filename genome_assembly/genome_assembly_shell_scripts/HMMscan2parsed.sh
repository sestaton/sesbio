#!/bin/bash
#--------------------------------
# take RS out and run HMMscan
#--------------------------------
#
#
#
#

cd /scratch/statonse/for_RepeatScout
perl ~/ePerl/external/bp_translate_seq.pl CGP-repeats-filtered.2 > CGP_JMBLAB_MCU_sunflowerBACs_11-09_RS-out.faa

/usr/local/hmmer-latest/bin/hmmscan -o CGP_BACs_RS-out_hmmscan_pfam-a_out.txt --cpu 4 --tblout=CGP_BACs_RS-out-hittable_pfam-a_out.txt --domtblout=CGP_BACs_RS-out-domain_table_pfam-a_out.txt /db/pfam/latest/Pfam-A.hmm CGP_JMBLAB_MCU_sunflowerBACs_11-09_RS-out.faa

perl ~/ePerl/parse_hmmer3.pl -i CGP_BACs_RS-out_hmmscan_pfam-a_out.txt -o CGP_BACs_RS-out_hmmscan_pfam-a_out-parsed.tab