#!/bin/bash
# take the raw output of mgblast which is 12 fields
# and get the query, subject, %identity, bitscore
cut -f 1,5,9,10 mgblast_out.txt > cut_fields.txt

# print each line if the %identity is above 90.00
awk '$3>90.00' cut_fields.txt > parsed_cut.txt

# now prepare for fgclust by only taking the query, subject, and score
cut -f 1,2,4 parsed_cut.txt > fgclust_in.txt

## try: cut -f 1,5,9,10 | awk | cut 