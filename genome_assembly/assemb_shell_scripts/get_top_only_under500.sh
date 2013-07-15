#!/bin/bash
cd /scratch/statonse/for_AAARF_split_sm/build_30_cleaned/bld30_blasthits

#remove duplicate entries : http://www.unix.com/shell-programming-scripting/127159-remove-duplicate-line-detail-based-column-one-data.html
#awk 'n!=$1{print;n=$1}' infile > outfile_top

blastall -p blastn -e 1e-5 -i build_30_under500bp.fasta -d ~/db/Repbase15.06 -o bld30_under500bp_repbase1506.bl7 -m7

perl ~/ePerl/parse_blast.pl -i bld30_under500bp_repbase1506.bl7 -f blastxml -o bld30_under500bp_repbase1506_parsed.txt 

./get_top_only_bw500and1.sh
./get_top_only_bw1and2.sh
./get_top_only_bw2and4.sh
./get_top_only_over4.sh