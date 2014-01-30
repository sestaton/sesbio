#!/bin/bash
cd /home/jmblab/statonse/HA412HO_UBC/for_mira_assemb
echo $LSB_HOSTS
cat /dev/null > mlist.$$
for variable in $LSB_HOSTS; do
echo $variable >> mlist.$$
done
mpirun -np 4 -machinefile mlist.$$ time mira -project=HA412_draft -fasta -job=denovo,genome,draft,454 -highlyrepetitive -noclipping=yes 
rm -f mlist.$$