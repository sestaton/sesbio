#!/bin/bash

# 1/27/11 - DAWGPAWS batch_blast.pl does not work for me so, I wrote my own batch_blast.pl
# -- SES

cd /scratch/statonse/for_fgclust_pipe/ha412ho_500k_split1/ha412ho-500k-split1_self_mgblast_out_parsed_cls-cluster_fasta_files

#------------------------------------------------
#export DP_BLAST_DIR='/home/jmblab/statonse/db/'

#perl ~/apps/dawgpaws/scripts/batch_blast.pl -i cluster_assemblies-Wed_Jan_26_22-07-07_EST_2011 -o /scratch/statonse/for_fgclust_pipe/ha412ho_500k_split1/ha412ho-500k-split1_self_mgblast_out_parsed_cls-cluster_fasta_files/cluster_assemblies-Wed_Jan_26_22-07-07_EST_2011 -d ~/db/ -c ~/apps/dawgpaws/scripts/config/batch_blast_tab.ecfg
#------------------------------------------------

perl ~/ePerl/batch_blast.pl -i cluster_assemblies-Wed_Jan_26_22-07-07_EST_2011/renamed_assembly_contigs