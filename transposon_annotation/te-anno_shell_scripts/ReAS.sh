#!/bin/bash

export PATH=$PATH:/home/jmblab/statonse/apps/ReAS_2.02/bin

cd /scratch/statonse/for_ReAS

perl /home/jmblab/statonse/apps/ReAS_2.02/bin/reas_all.pl -read ha412ho.fna.11 -k 20 -pa 4 -output consensus_test.fa -log ReAS_log