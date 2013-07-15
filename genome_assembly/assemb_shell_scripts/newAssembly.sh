#!/bin/bash
cd /home/jmblab/statonse/
/usr/local/454/bin/newAssembly dir BAC_17_fasta

cd /home/jmblab/statonse/
/usr/local/454/bin/addRun dir BAC_17_fasta sfffile /home/jmblab/statonse/454/sff/BAC_Region_1_MID.MID17.sff 

cd /home/jmblab/statonse/
/usr/local/454/bin/addRun dir BAC_17_fasta readfastafile BAC_17_LrgContigs.fasta