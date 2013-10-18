#!/bin/bash

#export PATH=$PATH:/home/jmblab/statonse/apps/dawgpaws/scripts
export TRNA_DB='/home/jmblab/statonse/apps/ltr_finder/tRNAdb/Os-tRNAs.fa'
export PROSITE_DIR='/home/jmblab/statonse/apps/ltr_finder/ps_scan'
export LTR_FINDER='/home/jmblab/statonse/apps/ltr_finder/ltr_finder'

cd `pwd`

perl ~/apps/dawgpaws/scripts/batch_ltrfinder.pl -i all_CGP_jmblab_BACs \
-o all_CGP_jmblab_BACs_bltrf_6-14 \
-c /home/jmblab/statonse/apps/dawgpaws/scripts/config/batch_ltrfinder.jcfg \
-f --verbose