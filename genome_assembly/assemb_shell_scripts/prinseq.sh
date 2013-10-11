#!/bin/bash

## NB: Remove the "-phred64" option for recently generated data.
##     This was needed to processes some Illumina data generated in Nov. 2010.

cd `pwd`

prinseq_lite=/usr/local/bioinfo/prinseq/latest/prinseq-lite.pl
prinseq_graphs=/usr/local/bioinfo/prinseq/latest/prinseq-graphs.pl

for file in ./*valid.fastq
do
    fq_base=$(echo $file | sed 's/\.f.*$//' -)
    good=${fq_base}_trimmed
    psl_log=${fq_base}_prinseq.log
    psg_log=${fq_base}_prinseq-graphs.log
    graph_data=${fq_base}_graph_data.txt
    out=${fq_base}_prinseq.out
    graphs_out=${fq_base}_prinseq_graphs

    perl $prinseq_lite -phred64 -fastq $file -out_format 3 -out_good $good -log $psl_log \
	-min_len 40 -noniupac -min_qual_mean 15 -lc_method entropy -no_qual_header \
	-lc_threshold 60 -trim_ns_right 10 -ns_max_p 20 -graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn \
	-graph_data $graph_data 2>&1 /dev/null

    perl $prinseq_graphs -i $graph_data -o $graphs_out -html_all -log $psg_log
done
