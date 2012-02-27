#!/bin/bash
# 
# TODO: compare the converions with blast->btab>gff3 to bp_seach2gff3.pl (blast->gff3)

cd /scratch/statonse/latest_sunflowerBACS-11-13/TransposonPSI
#perl ~/apps/TransposonPSI_08222010/transposonPSIcreate/BPbtab < CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX_TEpep.blx > CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX_TEpep.btab

perl ~/apps/TransposonPSI_08222010/scripts/TPSI_btab_to_gff3.pl CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX_TEpep.btab > CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX_TEpep_btab.gff3
