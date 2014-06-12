#!/bin/bash

cd `pwd`

##NB: vt is the vector file

runMapping -o mir_ref -vt ~/db/pIndigo_BAC536.fasta -ref BAC17_mira_out_over500.fasta_default.unpadded.fasta -read sub_1-BAC_Region_1_MID.MID17.sff
