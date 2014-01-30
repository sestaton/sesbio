#!/bin/bash
cd /scratch/statonse/ngs_HA412HO/reads/cleaned
perl ~/ePerl/filterSeqbyHit.pl -i ha412ho.qual.clean -o ha412ho.qual.clean_nochloro -b ha412_cp.bls -f qual --removehits