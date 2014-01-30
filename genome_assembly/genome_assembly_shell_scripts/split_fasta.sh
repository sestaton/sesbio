#!/bin/bash

N=$1
IN=$2
OUT=$3

awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%$N==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $2 > $3