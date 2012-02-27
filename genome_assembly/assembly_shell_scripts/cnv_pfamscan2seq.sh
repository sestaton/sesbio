#!/bin/bash

# create 3 files: a tab-delimited report of hits, a file of the sequences, a file of the sequence ID and hmm name
perl ~/ePerl/cnv_pfamscan2tab.pl -i CGP_JMBLAB_MCU_sunflowerBACs_11-09_pfamscan-out.txt -o CGP_JMBLAB_MCU_sunflowerBACs_11-09_pfamscan-out.tab --hittype --seq

# convert the pfamscan report to a 3 column tab report
#
# seqname hmm_name sequence
paste CGP_JMBLAB_MCU_sunflowerBACs_11-09_pfamscan-out.tab_HMM-name CGP_JMBLAB_MCU_sunflowerBACs_11-09_pfamscan-out.tab_SEQ > CGP_JMBLAB_MCU_sunflowerBACs_11-09_HMM-name_SEQ.txt

# get the reverse transcriptase sequences
grep "RVT_1" CGP_JMBLAB_MCU_sunflowerBACs_11-09_HMM-name_SEQ.txt > CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_1.txt
grep "RVT_2" CGP_JMBLAB_MCU_sunflowerBACs_11-09_HMM-name_SEQ.txt > CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_2.txt

# create a tab-delimited file of the sequence ID and sequence
cut -f 1,3 CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_1.txt > CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_1.seq
cut -f 1,3 CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_2.txt > CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_2.seq

# convert the tab-delimited seq files to fasta
perl ~/apps/repminer/scripts/cnv_tab2fasta.pl -i CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_1.seq -o CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_1.fasta
perl ~/apps/repminer/scripts/cnv_tab2fasta.pl -i CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_2.seq -o CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_2.fasta

# align the sequences
muscle -in CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_1.fasta -out CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_1_muscle.aln
muscle -in CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_2.fasta -out CGP_JMBLAB_MCU_sunflowerBACs_11-09_RVT_2_muscle.aln

# clean up
rm *seq *HMM* *tab *1.txt *2.txt *SEQ