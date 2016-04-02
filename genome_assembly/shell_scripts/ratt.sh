#!/bin/bash

## NB: This was for aligning chloroplast scaffolds relative to a referece chloroplast genome.
##     In this case, I was aligning Senecio to Helianthus.

export PAGIT_HOME=/usr/local/bioinfo/pagit/PAGIT
export RATT_HOME=$PAGIT_HOME/RATT

$RATT_HOME/start.ratt.sh ../Senecio Senecio_scaffolds.fasta Senecio_annotations Species Helianthus.fasta 
