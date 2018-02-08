#!/bin/bash

set -euo pipefail

cd `pwd`

##TODO: check if we want to run with mpi and how many cpus to use
export LD_PRELOAD=/usr/lib64/openmpi/lib/libmpi.so

maker=$HOME/apps/maker
mpiexec=/usr/lib64/openmpi/bin/mpiexec
export PERL5LIB=$maker/perl/lib

time $mpiexec -n 24 $maker/bin/maker -RM_off maker_opts.ctl maker_bopts.ctl maker_exe.ctl
#$maker --debug -RM_off maker_opts.ctl maker_bopts.ctl maker_exe.ctl
