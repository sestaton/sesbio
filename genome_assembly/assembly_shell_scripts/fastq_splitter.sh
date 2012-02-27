#!/bin/bash
# 
# splits a FASTQ file into pieces, with each read being added to the files in a
#  round-robin fashion
# don't forget to change the file name being catted
#
#
# changed from http://chadburrus.com/tag/file-splitter/

# to verify this works right, diff the original against the output of the
#  following, possibly ignoring newline differences
#
#    ls | grep "\-split.fa" | sort -n | xargs -i{} cat {} > merged

function usage() {
cat << EOF
USAGE: $0 <fastq> <target_num>

fastq       : fastq file to split
target_num  : number of files to split the fastq into

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

fastq=$(echo ${1%.*})
fastqext=$(echo ${1##*.})

awk '
BEGIN {
  # Change these to the number of files you want and the number of lines to
  #  print in each file before moving on to the next file
  NUM_FILES = $2;
  NUM_LINES = 10000000;
}
{
  lines += 1;
  if (lines == 1)
  {
    files += 1
  };
  print $1 > files "-split.fq";
  if (lines == NUM_LINES)
  {
    lines = 0;
    if (files == NUM_FILES)
    {
      files = 0
    }
  }
}' $1