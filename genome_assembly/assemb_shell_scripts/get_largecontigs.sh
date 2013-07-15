#!/bin/bash

#run the perl script to get large contigs and a list of contig names
perl ~/Desktop/ePerl/get_largecontigs.pl -i 454AllContigs.fna -o test_out.txt

#strip the > from every line
sed 's/^>//g' test_out.txt > test2.txt