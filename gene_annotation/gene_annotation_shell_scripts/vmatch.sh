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

subjSeq=$1
qrySeq=$2
merLen=$3
seqFile=$(echo ${subjSeq%.*})
qryFile=$(echo ${qrySeq%.*})
db=$seqFile"_"$merLen"mer_mkvtreedb"
vmerSearch=$seqFile"_"$merLen"mer_"$qryFile"_vmatch.out"
vmerSearchFull=$seqFile"_"$merLen"mer_"$qryFile"_vmatch_full.out"
vmerSearchIDs=$seqFile"_"$merLen"mer_"$qryFile"_vmatch_IDs.txt"
vmerSearchSeqnum=$seqFile"_"$merLen"mer_"$qryFile"_vmatch_seqnum.txt"
vmerSeq=$seqFile"_"$merLen"mer_"$qryFile"_vmatch_seqs.fasta"
vmerSubSeq=$seqFile"_"$merLen"mer_"$qryFile"_vmatch_subject_seqs.fasta"
vmerSubSeqCt=$seqFile"_"$merLen"mer_"$qryFile"_vmatch_subject_seqs_ct.txt"
#vmerSubSeqCtsort=$seqFile"_"$merLen"mer_"$qryFile"_vmatch_subject_seqs_ct_sort.txt"

# construct the persistent index
#
# For option "-pl n," n is the prefix length.
# Option "-allout" is major convenience for constructing all the index tables automatically
# instead of specifying each one individually on the command line (which is impossible to remember).
mkvtree -db $subjSeq -dna -indexname $db -allout -v -pl $merLen

# run Vmatch for some query
#
# Option "-showdesc 0" maps the result ID back to the original IDs.
# Option "-l n" sets the match length where n is the length.
vmatch -s -v -showdesc 0 -q $qrySeq -l $merLen $db > $vmerSearchFull
vmatch -q $qrySeq -l $merLen $db > $vmerSearch

# create an ID list of query sequences matching the index
grep -v "^#" $vmerSearchFull | sed -e 's/^[ \t]*//g;/^ *$/d' | egrep -v "Sbjct|Query" | perl -lane 'print join("\t",@F)' | cut -f2 > $vmerSearchIDs
grep -v "^#" $vmerSearch | sed -e 's/^[ \t]*//g' | perl -lane 'print $F[1]' > $vmerSearchSeqnum

# select the matching subject sequences from the index
vseqselect -seqnum $vmerSearchSeqnum $db > $vmerSeq # this requires the sequence numbers not the IDs

# grab the matching subject strings and summarize unique string counts
grep "Sbjct" $vmerSearchFull | perl -lane 'print $F[1]' | perl ~/ePerl/get_unique_word_counts.pl | sort -nrk 2 > $vmerSubSeqCt

# summarize unique string counts
#perl ~/ePerl/get_unique_word_counts.pl -i $vmerSubSeq -o $vmerSubSeqCt 
#sort -nrk 2 $vmerSubSeqCt > $vmerSubSeqCtsort

# clean up
#mv $vmerSubSeqCtsort $vmerSubSeqCt 
#rm $vmerSubSeqCtsort
rm ${db}.*
