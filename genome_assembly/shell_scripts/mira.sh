#!/bin/bash

cd `pwd`

time mira -project=RL11_sub1_30x -fasta -job=denovo,genome,draft,454 -highlyrepetitive 

## MPI code below
#echo $LSB_HOSTS
#cat /dev/null > mlist.$$
#for variable in $LSB_HOSTS; do
#echo $variable >> mlist.$$
#done
#mpirun -np 4 -machinefile mlist.$$ time mira -project=HA412_draft -fasta -job=denovo,genome,draft,454 -highlyrepetitive -noclipping=yes 
#rm -f mlist.$$
