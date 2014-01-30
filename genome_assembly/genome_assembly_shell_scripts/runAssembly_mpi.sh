#!/bin/bash
cd /iob_home/jlblab/statonse/454/
echo $LSB_HOSTS
cat /dev/null > mlist.$$
for variable in $LSB_HOSTS; do
   echo $variable >> mlist.$$
done
mpirun -np 4 -machinefile mlist.$$ /usr/local/454/bin/runAssembly -consed BAC_Region_1_MID.MID17.sff.fna
rm -f mlist.$$