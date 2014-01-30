#!/bin/bash

perl ~/ePerl/parse_blast_more.pl -i bld30_under500bp_repbase1506.bl7 -f blastxml -o bld30_under500bp_repbase1506_parsed_more.txt
perl ~/ePerl/parse_blast_more.pl -i bld30_between500bpand1kb_repbase1506.bl7 -f blastxml -o bld30_between500bpand1kb_repbase1506_parsed_more.txt
perl ~/ePerl/parse_blast_more.pl -i bld30_between1kband2kb_repbase1506.bl7 -f blastxml -o bld30_between1kband2kb_repbase1506_parsed_more.txt
perl ~/ePerl/parse_blast_more.pl -i bld30_between2kband4kb_repbase1506.bl7 -f blastxml -o bld30_between2kband4kb_repbase1506_parsed_more.txt
perl ~/ePerl/parse_blast_more.pl -i bld30_over4kb_repbase1506.bl7 -f blastxml -o bld30_over4kb_repbase1506_parsed_more.txt