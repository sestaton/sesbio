#!/bin/bash

cd /scratch/statonse/ngs_HA412HO/reads/cleaned

perl ~/ePerl/bp_removeSeqbyID2.pl -i ha412ho.fna.clean -o ha412ho.fna.clean_nochloro -b ha412_cp.bls -f fasta