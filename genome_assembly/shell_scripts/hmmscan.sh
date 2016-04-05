#!/bin/bash

#/usr/local/hmmer-latest/bin/hmmscan -o CGP_BACs_hmmscan_pfam-b_out.txt \
#--cpu 4 \
#--tblout=hittable_pfam-b_out.txt \
#--domtblout=domain_table_pfam-b_out.txt \
#/db/pfam/latest/Pfam-B.hmm \
#CGP_JMBLAB_MCU_sunflowerBACs_11-09.faa

usr/local/hmmer-latest/bin/hmmscan -o test_json_out --cpu 4 --json /db/pfam/latest/Pfam-A.hmm CGP_JMBLAB_MCU_sunflowerBACs_11-09.faa