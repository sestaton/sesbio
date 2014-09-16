#!/bin/bash


# -- screening searches: large contaminant hit expected; that is, 
#    at least 60% of the query should get a hit of at least 96% identity
#
# megablast -d ~/db/screendb -i SRR349656_prinseq_trimmed.fasta \
# -p 96 -W18 -JF -F "m D" -X30 -D3 -a 8 | grep -v "^#" | sort -k1,1 -k12,12 -nr > $tabout


usage() {
cat <<EOF

USAGE: $0 <fasta_file> <split_size>

fasta_file   : file of reads to translate
split_size   : integer representing the number of splits 
               to create for faster processing. Every ORF
               will be rejoined in a single file upon 
               completion.
 
EOF
}

print_error() {
cat <<ERR

ERROR: Command line not parsed correctly. Check input.

ERR
}

if [ $# -lt 1 ]; then
    print_error
    usage
    exit 1
fi

# gt is required for this script to work
# the below assumes you are running under the bash shell
hash gt &> /dev/null
if [ $? -eq 1 ]; then
    echo >&2 "genometools is required but it's not installed. Exiting."
    exit 1
fi

fastact=$(grep -c ">" $1)
gt splitfasta -numfiles $2 $1

fasta=$(echo ${1%.*})
fastaext=$(echo ${1##*.})

bl=$fasta"_bl_screendb_tmp.mbln"
filtered=$fasta"_filteredall.fas"
fh=$HOME/github/sesbio/genome_assembly/genome_assembly_shell_scripts/filter_hits.sh

#for i in {1..$2}
for i in $(seq $2)
do
  mkdir $i"_bl" && cd $i"_bl"
  mv ../$1"."$i .
  tmpbl=$bl"."$i
  t=$(date)
  echo "Running megablast on $1"."$i at $t"
  megablast -d ~/db/screendb -i $1"."$i \
      -p 96 -W18 -JF -F "m D" -X30 -D3 -a 8 | grep -v "^#" | sort -k1,1 -k12,12 -nr > $tmpbl
  bash $fh $1"."$i $tmpbl
  cd ../
  cat $i"_bl"/*filtered.$i >> $filtered
  rm -rf $i"_bl"
done

filteredct=$(grep -c ">" $filtered)

d=$(date)
echo "$0 finished blasting $1 at $d"
