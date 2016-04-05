#!/bin/bash

set -euo pipefail

if [ $# -lt 1 ]; then
    echo -e "No input found. Exiting.\n"
    echo -e "USAGE: $0 seqs_prot.faa\n"
    exit 1
fi

iprscan=$HOME/apps/interproscan-5.16-55.0/interproscan.sh

time $iprscan -i $1 \
    -appl Phobius,SUPERFAMILY,Gene3D,Coils,ProSiteProfiles,TIGRFAM,PRINTS,SMART,PIRSF,ProSitePatterns,Pfam,SignalP_EUK,ProDom \
    --goterms \
    --iprlookup \
    -b maker3_cds_pep_iprscan \
    --pathways UniPathway,KEGG,MetaCyc,Reactome \
    -f tsv,xml,gff3

#Phobius (1.01) : A combined transmembrane topology and signal peptide predictor
#        SignalP_GRAM_NEGATIVE (4.1) : SignalP (organism type gram-negative prokaryotes) predicts the presence and location of signal peptide cleavage sites in amino acid sequences for gram-negative prokaryotes.
#                  SUPERFAMILY (1.75) : SUPERFAMILY is a database of structural and functional annotation for all proteins and genomes.
#                       Gene3D (3.5.0) : Structural assignment for whole genes and genomes using the CATH domain structure database
#                        Hamap (201511.02) : High-quality Automated and Manual Annotation of Microbial Proteomes
#                        Coils (2.2.1) : Prediction of Coiled Coil Regions in Proteins
#              ProSiteProfiles (20.119) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them
#                      TIGRFAM (15.0) : TIGRFAMs are protein families based on Hidden Markov Models or HMMs
#                       PRINTS (42.0) : A fingerprint is a group of conserved motifs used to characterise a protein family
#                        SMART (6.2) : SMART allows the identification and analysis of domain architectures based on Hidden Markov Models or HMMs
#                        PIRSF (3.01) : The PIRSF concept is being used as a guiding principle to provide comprehensive and non-overlapping clustering of UniProtKB sequences into a hierarchical order to reflect their evolutionary relationships.
#              ProSitePatterns (20.119) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them
#                         Pfam (28.0) : A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs)
#                  SignalP_EUK (4.1) : SignalP (organism type eukaryotes) predicts the presence and location of signal peptide cleavage sites in amino acid sequences for eukaryotes.
#                       ProDom (2006.1) : ProDom is a comprehensive set of protein domain families automatically generated from the UniProt Knowledge Database.
#        SignalP_GRAM_POSITIVE (4.1) : SignalP (organism type gram-positive prokaryotes) predicts the presence and location of signal peptide cleavage sites in amino acid sequences for gram-positive prokaryotes.

#Deactivated analyses:
#                      PANTHER (10.0) : Analysis Panther is deactivated, because the resources expected at the following paths do not exist: data/panther/10.0/model
#                        TMHMM (2.0c) : Analysis TMHMM is deactivated, because the resources expected at the following paths do not exist: bin/tmhmm/2.0c/decodeanhmm, data/tmhmm/2.0c/TMHMM2.0c.model
