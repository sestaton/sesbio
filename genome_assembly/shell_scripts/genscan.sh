#!/bin/bash

export DP_GENSCAN_BIN='/home/jmblab/statonse/apps/bin/genscan'
export DP_GENSCAN_LIB='/home/jmblab/statonse/apps/genscan/Arabidopsis.smat'

perl ~/apps/dawgpaws/scripts/batch_genscan.pl -i all_HA_BACs/ -o all_HA_BACs_genscan --gff-ver 3 
