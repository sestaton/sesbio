#!/bin/bash
cd /scratch/statonse/latest_sunflowerBACS-11-13/TransposonPSI

# NCBI-BLAST
#formatdb -p T -i ~/apps/TransposonPSI_08222010/transposon_ORF_lib/transposon_db.pep -n transposon_db_pep
#tblastn -i CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.faa -d transposon_db_pep -o CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX_TEpep.tbln

# WUBLAST
#xdformat -p -o transposon_db_pep ~/apps/TransposonPSI_08222010/transposon_ORF_lib/transposon_db.pep
#blastx transposon_db_pep CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.fasta > CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX_TEpep.blx
blastp transposon_db_pep CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.faa > CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX_TEpep.blp