#!/bin/bash

###################################################
#LTRharvest
###################################################
# create suffix array for ltrharvest

gt suffixerator -db MCU_07A15.fasta -indexname MCU_BAC -tis -suf -lcp -ssp -sds -des -dna -v

gt ltrharvest -longoutput -mintsd 3 -maxlenltr 4000 -index MCU_BAC -out pred-all_MCU_BAC -outinner pred-inner_MCU_BAC -gff3 MCU_BAC.gff3

###################################################
#LTRdigest
###################################################
# sort the gff3 file(s) prior to running ltrdigest
gt gff3 -sort MCU_BAC.gff3 > MCU_BAC_sorted.gff3

gt ltrdigest -trnas /iob_home/jmblab/statonse/db/tRNAdb/plant_tRNAs.fasta -hmms /iob_home/jmblab/statonse/db/HMMs/*.hmm -aliout yes -aaout yes -outfileprefix ltrdigest_MCU_BAC MCU_BAC_sorted.gff3 MCU_BAC > MCU_BAC_ltrdigest.gff3
