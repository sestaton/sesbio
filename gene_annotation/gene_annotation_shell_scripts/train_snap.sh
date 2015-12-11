#!/bin/bash

set -euo pipefail

maker2zff=$HOME/apps/maker/bin/maker2zff
snapbin=$HOME/apps/maker/exe/snap

$maker2zff $1
$snapbin/fathom -categorize 1000 genome.ann genome.dna
$snapbin/fathom -export 1000 -plus uni.ann uni.dna
$snapbin/forge export.ann export.dna
$snapbin/hmm-assembler.pl ha412 . > ha412.hmm
