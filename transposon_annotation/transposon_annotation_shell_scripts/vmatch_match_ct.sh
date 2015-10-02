#!/bin/bash

script=$(basename $0)

function usage() {
cat <<EOF

USAGE: $script <subj_seq_file> <qry_seq_file> <kmer_len> 

subj_seq_file   :   A (nucleotide) Fasta file to index.
qry_seq_file    :   Fasta file to search the index.
kmer_len        :   An integer value to use in the analysis (e.g., 20).

EOF
}

function print_error() {
cat <<ERR

ERROR: Command line not parsed correctly. Check input.

ERR
}

if [ $# -lt 3 ]; then
    print_error
    usage
    exit 1
fi

for prog in mkvtree vmatch vseqselect 
do
  hash $prog &> /dev/null
  if [ $? -eq 1 ]; then
      echo >&2 "$prog is required but it's not installed. Exiting."
      exit 1
  fi
done

#
# We want the results in the current working directory, not in /home
# because the indices and map files are enormous.
#
cwd=$(pwd)

subjSeq=$1
qrySeq=$2
merLen=$3
seqFile=$(echo ${subjSeq%.*})
qryFile=$(echo ${qrySeq%.*})
basename=${cwd}/${seqFile}_${merLen}mer
db=${basename}_mkvtreedb
vmerSearchFull=${basename}_${qryFile}_vmatch_full.out
vmerSearchIDsCt=${basename}_${qryFile}_vmatch_IDs_ct.txt
vmerSearchSupfamMatchCt=${basename}_${qryFile}_superfamily_vmatch_ct.tsv
#vmatch_to_counts=$HOME/github/sesbio/transposon_annotation/vmatch_to_counts.pl

# construct the persistent index
#
# For option "-pl n," n is the prefix length.
# Option "-allout" is major convenience for constructing all the index tables automatically
# instead of specifying each one individually on the command line (which is impossible to remember).
mkvtree -db $subjSeq -dna -indexname $db -allout -v -pl 

# run Vmatch for some query
#
# Option "-showdesc 0" maps the result ID back to the original IDs.
# Option "-l n" sets the match length where n is the length.
vmatch -v -showdesc 0 -q $qrySeq -l $merLen -identity 100 $db |\
 grep -v "^#" |\
 sed -e 's/^[ \t]*//g;/^ *$/d' |\
 egrep -v "Sbjct|Query" |\
 perl -lane 'print join("\t",@F)' |\
 cut -f2 |\
 sort |\
 uniq -c |\
 sort -bnr > $vmerSearchIDsCt

#perl $vmatch_to_counts all_vmatch_20mer_repbase_mapped_files/repbase_idlist.txt $vmerSearchIDsCt |\
# sort -nrk 2 > $vmerSearchSupfamMatchCt

# Clean up?
rm ${db}.*
