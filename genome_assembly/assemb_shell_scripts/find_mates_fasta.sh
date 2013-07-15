#!/bin/bash
####################################################
# find_mates.sh
####################################################
# Description: mate paired reads in two fastq files
# and discared orphaned reads
# 
# input is two fastq files that have been trimmed 
# output is the paired sequences
#
#
# Author: S. Evan Staton
# Date: 6/23/11
# Contact: statonse at gmail dot com
#
# Notes: Some of the functions implemented in this 
# script were taken from the public domain. 
# See end of script for details.
#
# To my knowledge the fastq format for the Illumina 1.8
# pipeline has changed and that means this script won't work
# for these data.
####################################################
# TODO: 

function usage() {
cat << EOF
USAGE: $0 <left_pair> <right_pair> <rm-desc>

left_pair   : left mate in fastq format with header ending with /1
right_pair  : right mate in fastq format with header ending with /2
rm-desc     : remove the comments in the fastq file. Type "rm-desc"
              without quotes after <right_pair> name to minimize the 
              comment line to "+" thereby reducing the file size.

NB: This script works for Illumina 1.3+ but won't work for 1.8+ as
the file format has changed. This will be incorporated in a later 
version of this script.           
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
initlog "find_mates_"$date".log"

#wd=`pwd`

log =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
log $0 executed on `date`
log Left pair: $1
log Right pair: $2
log =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
log `printf "\n"`

# set the timer
tmr=$(timer)

leftmate=$(echo ${1%.*})
rightmate=$(echo ${2%.*})

leftext=$(echo ${1##*.})
rightext=$(echo ${2##*.})

# get the instrument name for safe matching ids in the index
# this assumes each pair came from the same instrument (of course, this should be true)
read -r fq_header1 < $1
instr1_name=$(echo $fq_header1 | grep -oE '>\w+.?\w+') # could check the source lane '^\@\w+\:\w'

read -r fq_header2 < $2
instr2_name=$(echo $fq_header2 | grep -oE '>\w+.?\w+')

# test to make sure they are from the same instrument
if [ $instr1_name != $instr2_name ]; then
    echo >&2 "$1 and $2 do not appear to be mates. They are from different instruments. Exiting"
    exit 1
fi

# pull each read id and sort 
sorted_left_ids=$leftmate"_sorted_ids_"$date
sorted_right_ids=$rightmate"_sorted_ids_"$date
    
echo "Capturing IDs for $1 ..."
#echo `grep "$instr_name" $1 | sed 's/\/1$//;s/^\@//' | sort -T $wd > $sorted_left_ids`
grep "$instr1_name" $1 | sed 's/\/1$//;s/>//' | sort > $sorted_left_ids

echo "Capturing IDs for $2 ..."
#echo `grep "$instr_name" $2 | sed 's/\/2$//;s/^\@//' | sort -T $wd > $sorted_right_ids`
grep "$instr1_name" $2 | sed 's/\/2$//;s/>//' | sort > $sorted_right_ids
    
# capture counts
leftct=$(awk 'END{print NR}' $sorted_left_ids)
rightct=$(awk 'END{print NR}' $sorted_right_ids)
    
# get the intersection of the two sorted sets
paired_id_list="all_"$leftmate"_"$rightmate"_paired_seq_ids_"$date
echo "Filtering ID lists for mates ..."
comm -1 -2 $sorted_left_ids $sorted_right_ids > $paired_id_list
paired=$(awk 'END{print NR}' $paired_id_list)
pairedct=$(($paired*2))
unpairedct=$(($leftct+$rightct))
orphanedct=$(($unpairedct-$pairedct))
rm $sorted_left_ids $sorted_right_ids
    
# create a list of ids to pull from the index
sorted_left_idlist_to_yank=$leftmate"_sorted_idlist_"$date
sorted_right_idlist_to_yank=$rightmate"_sorted_idlist_"$date
sed 's/$/\/1/' $paired_id_list > $sorted_left_idlist_to_yank
sed 's/$/\/2/' $paired_id_list > $sorted_right_idlist_to_yank
rm $paired_id_list
    
# create each index
echo "Creating indices for $1 and $2 ..."
cdbfasta $1 -o $1.index &> /dev/null
cdbfasta $2 -o $2.index &> /dev/null
    
# pull the mated seqs for each index
if [ -e "$1.index" ]; then
    echo "Pulling mated reads for $1 ..."
    if [ -e "$leftmate_mated.$leftext" ]; then
	first_mates=$leftmate"_mated_"$date"."$leftext
    else
	first_mates=$leftmate"_mated."$leftext
    fi
    cat $sorted_left_idlist_to_yank | cdbyank $1.index > $first_mates
else
    echo >&2 "Something went wrong with $1.index index creation. Exiting."
    exit 1
fi

rm $sorted_left_idlist_to_yank $1.index

if [ -e "$2.index" ]; then
    echo "Pulling mated reads for $2 ..."
    if [ -e "$rightmate_mated.$rightext" ]; then
	second_mates=$rightmate"_mated_"$date"."$rightext
    else
	second_mates=$rightmate"_mated."$rightext
    fi
    cat $sorted_right_idlist_to_yank | cdbyank $2.index > $second_mates
else
    echo >&2 "Something went wrong with $2.index index creation. Exiting."
    exit 1
fi

rm $sorted_right_idlist_to_yank $2.index

# remove the description
if [ "$3" == "rm-desc" ]; then
    echo "Removing comments from paired reads ..."
    first_mates_trim=$leftmate"_mated_"$date"_nodesc."$leftext
    second_mates_trim=$rightmate"_mated_"$date"_nodesc."$rightext
    sed 's/^\+.*/+/' $first_mates > $first_mates_trim
    sed 's/^\+.*/+/' $second_mates > $second_mates_trim

    rm $first_mates $second_mates
fi

echo "Done."
echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

log Total left paired reads: $leftct
log Total right paired reads: $rightct
log Total unpaired reads: $unpairedct
log Total paired reads: $pairedct
log Total orphaned reads: $orphanedct

echo "Total left paired reads: $leftct"
echo "Total right paired reads: $rightct"
echo "Total unpaired reads: $unpairedct"
echo "Total paired reads: $pairedct"
echo "Total orphaned reads: $orphanedct"

log `printf 'Total elapsed time: %s\n' $(timer $tmr)`
printf 'Total elapsed time: %s\n' $(timer $tmr)

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