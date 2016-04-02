#!/bin/bash

set -euo pipefail

dir=`pwd`

db=ref.fas
db_base=$(echo ${db%.*})
g_headlcvs=${db_base}_hscan_head.lcvs
g_taillcvs=${db_base}_hscan_tail.lcvs
g_paired=${db_base}_hscan_paired.txt

jar=/home/statonse/apps/HelitronScanner/HelitronScanner.jar 
lcvs=/home/statonse/apps/HelitronScanner/TrainingSet/head.lcvs
rcvs=/home/statonse/apps/HelitronScanner/TrainingSet/tail.lcvs

java -jar $jar scanHead -g $db -lf $lcvs -o $g_headlcvs --rc -tl 10 -buffer_size 1000000
java -jar $jar scanTail -g $db -lf $rcvs -o $g_taillcvs --rc -tl 10 -buffer_size 1000000
java -jar $jar pairends -hs $g_headlcvs -ts $g_taillcvs --rc -o $g_paired
java -jar $jar draw -p $g_paired -g $db -o bronze_helitrons_hscan --pure --flanking --ext -ext5 100 -ext3 100
