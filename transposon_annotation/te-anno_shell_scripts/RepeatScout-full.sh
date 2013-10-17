#!/bin/bash
cd /scratch/statonse/for_RepeatScout/

#----------------------------------------------------------
# RepeatScout-full.sh - simple shell script to run RS
#----------------------------------------------------------
#
# Author:   S. Evan Staton
# Date:     11/29/10
# Contact:  statonse at uga dot edu
# 
# NOTES:
# 
# ENV variables for nseg and trf have to be set in shell
#
# export PATH=$PATH:$HOME/apps/bin
# alias trf='trf400.linuxAMD64.exe'                               # RepeatScout filter-stage-2.prl uses alias 'trf'
# export TRF_COMMAND=$HOME/apps/bin/trf400.linuxAMD64.exe         # RepeatScout filter-stage-2.prl uses ENV var TRF_COMMAND
# 
# About setting the size of l from the manual:
#
# The default value of l, which is the "length of l-mer to consider", is set
# to be ceil(log_4(L)+1) where:
#   ceil(x) = smallest integer greater than x
#   log_4(x) = log base 4 of x
#   L is the length of the input sequence
# This value can be adjusted by giving the "-l" parameter, but it is essential
# that the same value of -l be given to both build_lmer_table and RepeatScout.
# It is not clear that values of l other than the default are sensible, but
# the options are there if you need them.
#-----------------------------------------------------------
# TODO: could run chained jobs to control stopping but the
# errors are diagnosable in verbose+ mode
#
# classify results at runtime (easy) - this script will be just for RepeatScout

#------------------------------
# the actual business is below 
#------------------------------      
# build the lmer table with the default as recommended
# IMPORTANT: v1.0.5 will read a multifasta but X's are not allowed in the sequence
time build_lmer_table -sequence CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.fasta -freq CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.freq

# run RepeatScout in verbose mode for debug, -stopafter sets the number of columns to search after no good alignment is found. 
# Default is 100 to decrease running time, but 500 was used in the paper and may give a better result
time RepeatScout -sequence CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.fasta -output CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.rs -freq CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.freq -stopafter 500 -vvvv

# filter-stage-1 removes tandem and low-complexity repeats that are unlikely to be mobile elements
time perl /usr/local/RepeatScout-latest/filter-stage-1.prl CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.rs > CGP-repeats-filtered.1

# create the masked file and lib from filtered repeats 
time RepeatMasker CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.fasta -e wublast -lib CGP-repeats-filtered.1

# select only the repeats above a certain threshold, small number of BACs so --thresh=2 (default --thresh=10 for genome)
time cat CGP-repeats-filtered.1 | /usr/local/RepeatScout-latest/filter-stage-2.prl --cat=CGP_JMBLAB_MCU_sunflowerBACs_11-09_noX.fasta.out --thresh=2 > CGP-repeats-filtered.2

