#!/bin/bash
# ///////////////////////////////////////////////////////////////////////////////////////
# be sure to create a symlink to the input file or copy it to the analysis directory
# otherwise, the mummer output files will be written to the directory of the input file 
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
cd /scratch/statonse/for_ReRep/run_ReRep

time perl run_ReRep.pl -i ha412ho.fna.1 -l 400 -r 3 --prefix part1-5 --clean