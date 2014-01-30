#!/bin/bash

#remove duplicate entries : http://www.unix.com/shell-programming-scripting/127159-remove-duplicate-line-detail-based-column-one-data.html
#awk 'n!=$1{print;n=$1}' infile

perl ../ePerl/top_blast_hits3.pl -i bld14_bt1an2kb_repbase1506.blast -f blast -o out | awk 'n!=$1{print;n=$1}' out > parsed
