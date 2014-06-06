#!/bin/bash

cd `pwd`

runAssembly=/usr/local/454/bin/runAssembly 

time $runAssembly -o HA412_newbler -consed ha412ho.fna ha412ho.qual
