#!/bin/bash

export PATH=/home/statonse/apps/blat:$PATH
scipio=/home/statonse/apps/scipio/scipio-1.4/scipio.1.4.1.pl
yaml2gff=/home/statonse/apps/scipio/scipio-1.4/yaml2gff.1.4.pl
scipiogff2gff=/home/statonse/apps/maker/exe/augustus/scripts/scipiogff2gff.pl
yaml2log=/home/statonse/apps/maker/scipio/scipio-1.4/yaml2log.1.4.pl
gff2gbSmallDNA=/home/statonse/apps/maker/exe/augustus/scripts/gff2gbSmallDNA.pl
etraining=/home/statonse/apps/maker/exe/augustus/bin/etraining
filterGenes=/home/statonse/apps/maker/exe/augustus/scripts/filterGenes.pl

$scipio --blat_output=prot.vs.genome.psl genome.fa proteins.aa > scipio.yaml
cat scipio.yaml | $yaml2gff > scipio.scipiogff
$scipiogff2gff --in=scipio.scipiogff --out=scipio.gff
cat scipio.yaml | $yaml2log > scipio.log
$gff2gbSmallDNA scipio.gff genome.fa 1000 genes.raw.gb
$etraining --species=generic --stopCodonExcludedFromCDS=true genes.raw.gb 2> train.err

cat train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > badgenes.lst
$filterGenes badgenes.lst genes.raw.gb > genes.gb
