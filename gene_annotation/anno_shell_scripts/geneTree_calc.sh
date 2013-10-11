#!/bin/bash

cd `pwd`

perl orthoMCLgroups2GeneFamilies_trees_per_group_v0.10.pl -g groupcopy -gb goodProteins_group_representatives_Vvinifera_peptide_top_reciprocal_hits_edit.txt -onf Vvinifera_145_cds_edit.fa -opf Vvinifera_145_peptide_edit.fa -nf goodGenes_nt.fasta -pf goodProteins.fasta -os v10stats -ts ortho_group_tree_stats.tsv