#!/bin/bash

#=============================================================
# Name: create_GO-mapping_files.sh
#
# Purpose: This is a pipeline script to take HMMER2GO
#          output and generate a GO term association file
#          for each gene. Also, a population and study set
#          of gene IDs is created for use with Ontologizer.
#          It is possible to also run Ontologizer from the
#          command line, but the results are different
#          from the GUI version.
#
# Author: S. Evan Staton
# Date: 4/2/12
# Updated: 8/13/14
#
# NB: This has been updated to use HMMER2GO (https://github.com/sestaton/HMMER2GO)
#     but it is untested. Should only take minor changes if there are
#     bugs...be warned though.
#==============================================================
# TODO: 

function usage() {
cat <<EOF 

USAGE: $0 input_study_directory_name input_popn_dir_name output_directory_name species_name

input_study_directory_name  : Name of directory holding ONLY the study HMMscan output files (output of option --tblout).
input_popn_dir_name         : Name of directory holding population HMMscan output file.
output_directory_name       : Name of directory to place the GO mapping files.
species_name                : A species name must be given for use in the association file.

EOF
}

function print_error() {
cat <<ERR

ERROR: Command line not parsed correctly. Check input.

ERR
}

if [ $# -lt 4 ]; then
    print_error
    usage
    exit 1
fi

studydir=$1
popndir=$2
outdir=$3
species=$4

if [ -d "$outdir" ]; then
    if [ ! -L "$outdir" ]; then
    echo -e "$outdir exists already. Exiting.\n"
    usage
    exit 1
    else
	echo -d "$outdir appears to be a symlink. Exiting.\n";
	usage 
	exit 1
    fi
fi

if [ ! -a "$outdir" ]; then    # make sure $dir is not a file
    mkdir $outdir              # add some granularity (bitmask)
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
initlog "create_GO-mapping_files_"$date".log"

log =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
log $0 executed at `date`
log `printf "\n"`

# set the timer
tmr=$(timer)

#
# Get Pfam -> GO mapping file
#
hmmer2go fetchmap -o pfam2go

# 
# Generate custom gene -> GO term mapping file from HMMscan report
#
cd $studydir
for file in * 
do
  filebase=$(echo ${file%.*})
  outfile=$filebase"_parsed.txt"
  studyIDs=$filebase"_studyIDs.txt"

  #perl ~/ePerl/hmm2go.pl -i $file -p ~/db/pfam2go -o $outfile
  hmmer2go mapterms -i $file -p pfam2go -o $outfile

  cut -f1 $outfile > $studyIDs
 
  #
  # Send the results to the output dir
  #
  mv $outfile ../$outdir
  mv $studyIDs ../$outdir
done
cd ../$popndir

# This is the full genome, which is needed to create the 'population' set of IDs. 
# The other files will be used as the 'study' sets. 
#
for file in *
do
  filebase=$(echo ${file%.*})
  outfile=$filebase"_parsed.txt"
  gomapfile=$filebase"_parsed_GOterm_mapping.txt"
  gaffile=$filebase"_parsed_GOterm_mapping.gaf"
  popnIDs=$filebase"_parsed_populationIDs.txt"

  #perl ~/ePerl/hmm2go.pl -i $file -p ~/db/pfam2go -o $outfile --map
  hmmer2go mapterms -i $file -p pfam2go -o $outfile --map

  #
  # Generate gaf file to use with Ontologizer
  #
  #perl ~/ePerl/generate_GO_association.pl -i $gomapfile -o $gaffile -s $species -g ~/db/GO.terms_alt_ids
  hmmer2go map2gaf -i $gomapfile -o $gafffile -s $species

  #
  # Generate population gene list
  #
  cut -f1 $outfile > $popnIDs

  #
  # Send the results to the output dir
  #
  mv $gomapfile ../$outdir
  mv $outfile ../$outdir
  mv $popnIDs ../$outdir
  mv $gaffile ../$outdir
done
cd ..

#
# Print time to completion
#
log `printf 'Total elapsed time in minutes: %s\n' $(timer $tmr)`
log =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
