#!/bin/bash

set -euo pipefail

cd `pwd`

runAssembly=/usr/local/454/bin/runAssembly 
out=newbler_out
fas=ha412.fna
qual=ha412.qual

time $runAssembly -o $out -consed $fas $qual
