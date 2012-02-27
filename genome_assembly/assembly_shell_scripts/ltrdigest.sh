#!/bin/bash

export PATH=$PATH:/home/jmblab/statonse/apps/genometools-1.3.4/bin

###################################################
#LTRharvest
###################################################
# create suffix array for ltrharvest

cd /home/jmblab/statonse/draft_jmblab_BACs/for_ltrdigest/MCU_BAC
#do
gt suffixerator -db MCU_07A15.fasta -indexname MCU_BAC -tis -suf -lcp -ssp -sds -des -dna -v
#done
#do
gt ltrharvest -longoutput -maxlenltr 4000 -index MCU_BAC -out pred-all_MCU_BAC -outinner pred-inner_MCU_BAC -gff3 MCU_BAC.gff3
#done
###################################################
#LTRdigest
###################################################
# create suffix array for ltrdigest
#gt suffixerator -tis -des -dna -ssp -db all_CGP_Burke_sunflower_BACs.fasta -indexname all_Sunflower_BACs_ltrdigest

# sort the gff3 file(s) prior to running ltrdigest
gt gff3 -sort MCU_BAC.gff3 > MCU_BAC_sorted.gff3

gt ltrdigest -trnas /home/jmblab/statonse/db/tRNAdb/plant_tRNAs.fasta -hmms /home/jmblab/statonse/db/HMMs/*.hmm -aliout yes -aaout yes -outfileprefix ltrdigest_MCU_BAC MCU_BAC_sorted.gff3 all_Sunflower_BACs > MCU_BAC_ltrdigest.gff3
