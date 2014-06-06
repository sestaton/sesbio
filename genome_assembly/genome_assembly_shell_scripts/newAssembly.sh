#!/bin/bash

## Description - Create a new assembly project for Newbler from various sources.

newAssembly=/usr/local/454/bin/newAssembly 
addRun=/usr/local/454/bin/addRun

$newAssembly dir BAC_17_fasta

$addRun dir BAC_17_fasta sfffile /home/jmblab/statonse/454/sff/BAC_Region_1_MID.MID17.sff 

$addRun dir BAC_17_fasta readfastafile BAC_17_LrgContigs.fasta