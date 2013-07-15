#!/bin/bash
# Usage: sh remove_duplicates.sh < infile > outfile
#http://www.unix.com/shell-programming-scripting/97827-how-can-i-remove-those-duplicate-sequence-unix-what-command-line-i-should-type.html

awk '/^>/{if(hdr=="") hdr=$0}/^[^>]/{x[$0]++;if(x[$0]==1) {print hdr;print} hdr=""}' 