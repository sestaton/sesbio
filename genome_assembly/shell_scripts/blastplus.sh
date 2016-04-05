#!/bin/bash

blastdir=/usr/local/ncbi-blast+-latest/bin

#$blastdir/makeblastdb -in -dbtype nucl -out -parse_seqids

$blastdir/blastn -query contigs.fa -db /usr/local/ncbiblast-latest/nt -outfmt 6 -out contigs_nt.bln -num_threads 8