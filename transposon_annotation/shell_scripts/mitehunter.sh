#!/bin/bash

cd `pwd`

db=/home/statonse/db/Ha412v1r1_genome_no_cp-mt-rd.fasta

perl ~/apps/MITE_Hunter/MITE_Hunter_manager.pl -i $db -g HAmites -n 8 -S 12345678