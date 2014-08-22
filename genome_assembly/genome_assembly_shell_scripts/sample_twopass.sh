#!/bin/bash

##New option with 'seqtk sample' see: https://www.biostars.org/p/110107/#110248

#Downsample a fixed number (2-pass mode):

seqtk sample -2 read1.fa.gz 20000 > sub1.fa
seqtk sample -2 read2.fa.gz 20000 > sub2.fa

##It reads the input twice (so twice as slow). 
##In the first pass, it finds the sampled read indices. 
##In the second pass, it outputs reads at the stored indices. 
##The peak RAM is about the number of sampled reads multiplied by 24, again independent of the input. 
##You need the latest seqtk for this 2-pass mode.