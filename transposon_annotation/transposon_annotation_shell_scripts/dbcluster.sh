#!/bin/bash

## NB: self-clustering using vmatch

mkvtree -db ha412ho.fna.1 -dna -indexname ha412_part1 -tis -ois -suf -bwt -bck -v -pl

vmatch -dbcluster 80 7 cluster -pred-all -p -d -seedlength 50 -l 1101 -exdrop 9 ha412_part1