#!/bin/bash

export PAGIT_HOME=/usr/local/bioinfo/pagit/PAGIT
export RATT_HOME=$PAGIT_HOME/RATT

$RATT_HOME/start.ratt.sh ../Senecio Senecio_scaffolds.fasta Senecio_annotations Species Helianthus.fasta 
