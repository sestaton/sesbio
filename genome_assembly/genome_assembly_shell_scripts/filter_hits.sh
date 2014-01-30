#!/bin/bash
####################################################
# filter_hits.sh
####################################################
# Description: match sequences to an reference
#              and remove them 
# 
# input is a fasta file of reads to filter and a 
# reference genome or sequences
#
#
# Author: S. Evan Staton
# Date: 
# Contact: statonse at gmail dot com
#
# Notes: Some of the functions implemented in this 
# script were taken from the public domain. 
# See end of script for details.
#
# 
####################################################
# TODO: 

function usage() {
cat << EOF
USAGE: $0 <fasta_file> <blast_file>

fasta_file   : file or reads to filter
blast_file   : tab-delimited blast output 
         
EOF
}

function print_error() {
cat << ERR

ERROR: Command line not parsed correctly. Check input.

ERR
}

if [ $# -lt 2 ]; then
    print_error
    usage
    exit 1
fi

# cdbfasta is required for this script to work
# the below assumes you are running under the bash shell
hash cdbfasta &> /dev/null
if [ $? -eq 1 ]; then
    echo >&2 "cdbfasta is required but it's not installed. Exiting."
    exit 1
fi

hash cdbyank &> /dev/null
if [ $? -eq 1 ]; then
    echo >&2 "cdbyank is required but it's not installed. Exiting."
    exit 1
fi

function timer() {
    if [[ $# -eq 0 ]]; then
        echo $(date '+%s')
    else
        local  stime=$1
        etime=$(date '+%s')

        if [[ -z "$stime" ]]; then stime=$etime; fi

        dt=$((etime - stime))
        ds=$((dt % 60))
        dm=$(((dt / 60) % 60))
        dh=$((dt / 3600))
        printf '%d:%02d:%02d' $ddh $dm $ds
    fi
}

function initlog () {
    LOGFILE=$1
    echo '' > ${LOGFILE}
}

function log () {
    echo $* >> ${LOGFILE}
}

date=$(date | sed 's/\ //g;s/\:/\_/g')
initlog "filter_hits_"$date".log"

log =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
log $0 executed on `date`
log Sequence file: $1
log BLAST file   : $2
log =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
log `printf "\n"`

# set the timer
tmr=$(timer)

fasta=$(echo ${1%.*})
blast=$(echo ${2%.*})

fastaext=$(echo ${1##*.})
blastext=$(echo ${2##*.})

# pull each read id and sort 
sorted_fasta_ids=$fasta"_sorted_ids_"$date
sorted_blast_ids=$blast"_sorted_ids_"$date
    
echo "Capturing IDs for $1 ..."
grep ">" $1 | sed 's/>//;s/\s.*//g' | sort > $sorted_fasta_ids  # modified 12/8/11 to shorten IDs SES
  
echo "Capturing IDs for $2 ..."
cut -f 1 $2 | sed 's/\s.*//g' | sort -u > $sorted_blast_ids     # modified 12/8/11 SES
  
# capture counts
fastact=$(awk 'END{print NR}' $sorted_fasta_ids)
blastct=$(awk 'END{print NR}' $sorted_blast_ids)
    
# get the complement of the two sorted sets
comp_id_list="all_"$fasta"_"$blast"_seq_ids_"$date

echo "Filtering ID lists for $1 ..."

comm -1 -3 $sorted_blast_ids $sorted_fasta_ids > $comp_id_list

filteredct=$(awk 'END{print NR}' $comp_id_list)
sumct=$(($fastact+$blastct))
rm $sorted_fasta_ids $sorted_blast_ids
    
# create a list of ids to pull from the index
#sorted_idlist_to_yank=$fasta"_sorted_idlist_"$date
    
# create each index
echo "Creating index for $1 ..."
cdbfasta $1 -o $1.index &> /dev/null

# do we want to perform this action on a fastq?
#cdbfasta $1 -Q -o $1.index &> /dev/null
    
# pull the mated seqs for each index
if [ -e "$1.index" ]; then
    echo "Pulling reads for $1 ..."
    if [ -e "$fasta_filtered.$fastaext" ]; then
	fasta_filtered=$fasta"_filtered_"$date"."$fastaext
    else
	fasta_filtered=$fasta"_filtered."$fastaext
    fi
    cat $comp_id_list | cdbyank $1.index > $fasta_filtered
else
    echo >&2 "Something went wrong with $1.index index creation. Exiting."
    exit 1
fi

rm $comp_id_list $1.index

echo "Done."
echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

log Total num reads: $sumct
log Total discarded reads: $blastct
log Total filtered reads: $filteredct

echo "Total num reads: $sumct"
echo "Total discarded reads: $blastct"
echo "Total filtered reads: $filteredct"

log `printf 'Total elapsed time: %s\n' $(timer $tmr)`
printf 'Total elapsed time: %s\n' $(timer $tmr)
log =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

############################################################################################
# Acknowledgements:
#
# Hashing method was recommeded on stack overflow:
# http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
#
# Timer function taken from Linux Journal:
# http://www.linuxjournal.com/content/use-date-command-measure-elapsed-time 
# 
# Logging method taken from:
# http://www.cv.nrao.edu/~jmalone/talks/bash.pdf 