#!/bin/bash

cd `pwd`

#/usr/local/RepeatModeler-latest/BuildDatabase -name CGP-repeats-BACs CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.fasta
xdformat -n -I -o CGP-repeats-BACs CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.fasta

time /usr/local/RepeatModeler-latest/RepeatModeler -database CGP-repeats-BACs >& run.out &