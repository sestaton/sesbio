#!/bin/bash

cd `pwd`

perl /usr/local/pfam-scan/latest/pfam_scan.pl \
-fasta CGP_JMBLAB_MCU_sunflowerBACs_11-09.faa \
-align \
-cpu 4 \
-dir /db/pfam/latest/ \
-outfile CGP_JMBLAB_MCU_sunflowerBACs_11-09_pfamscan-out.txt
