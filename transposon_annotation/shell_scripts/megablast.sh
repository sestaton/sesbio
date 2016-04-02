#!/bin/bash

query=test_data/t_reads.fas
db=t_readsdb

megablast -i $query \
-d $db \
-F T \
-D 4 \
-m 0 \
-V T \
-p 0 \
-W 18 \
-U T \
-X 40 \
-J F \
-a 2
