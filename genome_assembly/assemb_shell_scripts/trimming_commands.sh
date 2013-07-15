#!/bin/bash

#TODO: Add quality trimming as Step 2 and shift everything else down the list.

# 1. Run FastQC to analyze statistics.

# 2. Trim the first 15 bases (NB: the actual need for this step and lengths will be data-dependent).
awk 'NR % 2 == 0 { print substr($1, 16, 150 - 15) } NR % 2 == 1' PF_S3_L001_R2_001_trimmed.fastq > PF_S3_L001_R2_001_trimmed_cutfirst15bp.fastq

# 3. Take a peek at the reads for the presence of adapters.
#grep -v "@" PAM_S2_L001_R1_001_trimmed.fastq | grep -v "+" | grep -v "[^ATCG]" | grep AGGTCCTGTCTCTTATACACATCTCCGAGCCCACGAGACCTAGGTGA
ack AGGTCCTGTCTCTTATACACATCTCCGAGCCCACGAGACCTAGGTGA PAM_S2_L001_R1_001_trimmed.fastq

# 4. Trim the adapter.
cutadapt -a AGGTCCTGTCTCTTATACACATCTCCGAGCCCACGAGACCTAGGTGA -q 20 -m 50 -e 0.3 PF_S3_L001_R2_001_trimmed_cutfirst15bp.fastq > PF_S3_L001_R2_001_trimmed_cutfirst15bp_cutadapt.fastq

# 5. Run FastQC again and compare to the report generated in Step 1.
