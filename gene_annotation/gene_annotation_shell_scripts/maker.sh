#!/bin/bash

cd `pwd`

##TODO: check if we want to run with mpi and how many cpus to use
export LD_PRELOAD=/usr/lib64/openmpi/lib/libmpi.so

maker=/home/statonse/apps/maker/bin/maker
mpiexec=/usr/lib64/openmpi/bin/mpiexec

$mpiexec -n 24 $maker maker_opts.ctl maker_bopts.ctl maker_exe.ctl
