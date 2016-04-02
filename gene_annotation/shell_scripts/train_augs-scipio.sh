#!/bin/bash

set -euo pipefail

AUGUSTUS=$HOME/apps/augustus-3.2.1
export PATH=$PATH:$AUGUSTUS/bin
export PATH=$HOME/apps/blat:$PATH
export AUGUSTUS_CONFIG_PATH=$AUGUSTUS/config

scipio=$HOME/apps/scipio/scipio-1.4/scipio.1.4.1.pl
yaml2gff=$HOME/apps/scipio/scipio-1.4/yaml2gff.1.4.pl
scipiogff2gff=$AUGUSTUS/scripts/scipiogff2gff.pl
yaml2log=$HOME/apps/scipio/scipio-1.4/yaml2log.1.4.pl
gff2gbSmallDNA=$AUGUSTUS/scripts/gff2gbSmallDNA.pl
etraining=$AUGUSTUS/bin/etraining
augustus=$AUGUSTUS/bin/augustus
filterGenes=$AUGUSTUS/scripts/filterGenes.pl
randomSplit=$AUGUSTUS/scripts/randomSplit.pl
newSpecies=$AUGUSTUS/scripts/new_species.pl
optimizeAugustus=$AUGUSTUS/scripts/optimize_augustus.pl

genome=$HOME/db/Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta.masked

## NB: This first step is the longest as it involves alignment with blat
echo "running scipio"
$scipio --blat_output=prot.vs.genome.psl $genome proteins.aa > scipio.yaml
echo "Done running scipio."
cat scipio.yaml | $yaml2gff > scipio.scipiogff
echo "Converting to gff now..."
$scipiogff2gff --in=scipio.scipiogff --out=scipio.gff
echo "to gff"
cat scipio.yaml | $yaml2log > scipio.log
echo "writing log.."
$gff2gbSmallDNA scipio.gff genome.fa 1000 genes.raw.gb
echo "converting scipio gff to smalldna"
$etraining --species=generic --stopCodonExcludedFromCDS=true genes.raw.gb 2> train.err
echo "running etraining.."

cat train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > badgenes.lst
$filterGenes badgenes.lst genes.raw.gb > genes.gb
echo "filtering genes.."

## Now create a new species and make the training set
echo "making random split.."
$randomSplit genes.gb 100
echo "making new species..."
$newSpecies --species=sunflower
echo "running etraining on training set..."
$etraining --species=sunflower genes.gb.train
echo "running augustus on test set..."
$augustus --species=sunflower genes.gb.test | tee firsttest.out
grep -A 22 Evaluation firsttest.out
echo "running optimize procedure..."
$optimizeAugustus --species=sunflower genes.gb.train
echo "running training round 2..."
$etraining --species=sunflower genes.gb.train
echo "accuracy prediction by running augustus round 2..."
$augustus --species=sunflower genes.test.gb
