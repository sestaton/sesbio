#!/bin/bash

cd /scratch/statone/for_AAARF_split_sm/all_builds_blasthits/parsed_hits
#foreach file *.txt
sed '/^#/d' all_AAARF_builds_repbase1506_parsed.txt > all_AAARF_builds_repbase1506_parsed_nocomments.txt
#end