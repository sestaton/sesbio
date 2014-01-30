#!/bin/bash

cat HA0002_620J0AAXX_2_3_concat_qseq.txt | awk '$11 == 1' > HA0002_620J0AAXX_2_3_concat_qseq_filtered.txt
