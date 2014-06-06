#!/bin/bash

## NB: This 'seqio.py' script is part of the "backbone" pipeline. Not sure
##     if it is maintained anymore.

cd `pwd`

python2.6 /usr/local/Python-2.6.2/bin/seqio.py \
-s lb_ha412.pl_454.sm_412ho.fasta \
-q lb_ha412.pl_454.sm_412ho.qual \
-f fasta \
-o lb_ha412.pl_454.sm_412ho.sfastq \
-t sfastq