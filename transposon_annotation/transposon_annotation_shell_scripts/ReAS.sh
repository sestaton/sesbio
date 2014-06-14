#!/bin/bash

## NB: This program has bugs and never finshes, don't use it.

cd `pwd`

export PATH=$PATH:/home/jmblab/statonse/apps/ReAS_2.02/bin

reas=/home/jmblab/statonse/apps/ReAS_2.02/bin/reas_all.pl 

perl $reas -read ha412ho.fna.11 -k 20 -pa 4 -output consensus_test.fa -log ReAS_log