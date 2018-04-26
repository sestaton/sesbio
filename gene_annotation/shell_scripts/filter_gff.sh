#!/bin/bash

set -euo pipefail

list=nonTE_gene_list_80cov_90pid.txt

function split_chrs {
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr01 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr1.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr02 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr2.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr03 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr3.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr04 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr4.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr05 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr5.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr06 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr6.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr07 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr7.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr08 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr8.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr09 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr9.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr10 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr10.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr11 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr11.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr12 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr12.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr13 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr13.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr14 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr14.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr15 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr15.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr16 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr16.gff3
    cat <(head -1 HanXRQr1.0-20151230.fa.all_sort.gff3) <(grep Chr17 HanXRQr1.0-20151230.fa.all_sort.gff3) > chr17.gff3
    egrep -v "Chr01|Chr02|Chr03|Chr04|Chr05|Chr06|Chr07|Chr08|Chr09|Chr10|Chr11|Chr12|Chr13|Chr14|Chr15|Chr16|Chr17" \
	HanXRQr1.0-20151230.fa.all_sort.gff3 > chrX.gff3
}

#split_chrs()

for file in ./chr{1..17}.gff3 ./chrX.gff3
do
    echo "=====> working on $file"
    base=$(echo ${file%.*})
    before=$(grep -v "^#" $file | awk '$3 == "gene"' | wc -l)
    perl gff_filter_gene_list.pl -g $file -l $list > ${base}_filtered.gff3
    after=$(grep -v "^#" ${base}_filtered.gff3 | awk '$3 == "gene"' | wc -l)
    echo "===> Before filtering: $before"
    echo "===> After filtering : $after"
done

echo "=====> sorting final gff3"
gt gff3 -sort *filtered.gff3 > all_chrs_sorted.gff3
listct=$(wc -l $list)
finalct=$(grep -v "^#" all_chrs_sorted.gff3 | awk '$3 == "gene"' | wc -l)
rm *filtered.gff3
echo "===> Gene count in input list: $listct"
echo "===> Gene count in final GFF3: $finalct"
echo "=====> Done."
