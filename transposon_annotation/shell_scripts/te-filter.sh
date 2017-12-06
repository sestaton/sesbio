#!/bin/bash

## Filter TE-related protein/domain matches from a repeat set to identify gene fragments
## USAGE: bash te-filter.sh some_annotation_table.tsv > filtered.tsv

# rules for filtering matches
grep -v -i "transposon" $1 |
grep -v -i "helicase" |  
grep -v -i "hypothetical protein" | 
grep -v -i "ribonuclease H-like" | 
grep -v -i "transposase" | 
grep -v -i "chromo" | 
grep -v -i "zinc finger" | 
grep -v -i "tyrosine-protein kinase" | 
grep -v -i "dna binding" | 
grep -v -i "reverse transcriptase" | 
grep -v -i "aspartic peptidase" | 
grep -v -i "gag-pre" | 
grep -v -i "integrase" | 
grep -v -i "dna-binding" | 
grep -v -i "gag domain" | 
grep -v -i "gag-pol" | 
grep -v -i "gag protein" |
#grep -v -i "arabidopsis retrotransposon" | 
grep -v -i "copia protein" |
grep -v -i "zinc knuckle" |
grep -v -i "polymerase"


