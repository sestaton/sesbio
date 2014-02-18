#!/bin/bash

## The read pairs have to be in the same order,
## otherwise overlaps will be found between the wrong pairs

## Specifying the mode is important, without it, the program fails silently
## -m 0 == simple overlap mode

 ~/apps/cope-v1-1-3/src/cope \
-a PAM_S2_L001_R1_001_trimmed_p.fasta \
-b PAM_S2_L001_R2_001_trimmed_p.fasta \
-o PAM_connected_paired.fasta \
-2 pam_1_noncon_paired.fasta \
-3 pam_2_noncon_paired.fasta \
-m 0 \
> pam_cope_paired.log 2> pam_cope_paired.err