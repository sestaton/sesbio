#!/bin/bash

awk 'NR>1 {if(/>/){print c"\t"l;c=$1}else {l=$1}}' 454AlignmentInfo.tsv