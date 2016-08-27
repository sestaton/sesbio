#!/bin/bash

set -euo pipefail

orthoMCLdir=$HOME/apps/orthomcl/orthomclSoftware-v2.0.9/bin
mcl=$HOME/apps/mcl-14-137/bin/mcl
out=orthomcl 

# log start
start=$(date)
echo "=====> Starting OrthoMCL at $start."

## orthomclInstallSchema
#mysql> CREATE DATABASE orthomcl;
#mysql> GRANT SELECT,INSERT,UPDATE,DELETE,CREATE VIEW,CREATE, INDEX, DROP on orthomcl.* TO evan@localhost;
#echo "Starting orthomclInstallSchema..."
#time perl $orthoMCLdir/orthomclInstallSchema \
#    orthomcl.config $out/install_schema.log

#echo "Done with orthomclInstallSchema."
## orthomclAdjustFasta
## orthomclFilterFasta
## --- BLAST ---

## orthomclBlastParser
echo "Starting orthomclBlastParser..."
time perl $orthoMCLdir/orthomclBlastParser \
    goodProteins_allvall.bln \
    $out > similarSequences.txt 2> blastparsed.out

echo "Done with orthomclBlastParser."
## orthomclLoadBlast
echo "Starting orthomclLoadBlast..."
time perl $orthoMCLdir/orthomclLoadBlast \
    orthomcl.config \
    similarSequences.txt 2>&1 > loadblast.out

echo "Done with orthomclLoadBlast."
## orthomclPairs
echo "Starting orthomclPairs..."
time perl $orthoMCLdir/orthomclPairs \
    orthomcl.config \
    $out/pairs.log cleanup=no 2>&1 > pairs.out

echo "Done with orthomclPairs."
## orthomclDumpPairsFiles
echo "Starting orthomclDumpPairsFiles..."
time perl $orthoMCLdir/orthomclDumpPairsFiles \
    orthomcl.config 2>&1 > dumppairs.out 

echo "Done with orthomclDumpPairsFiles."
## --- mcl ---
echo "Starting mcl..."
time $mcl mclInput --abc -I 1.5 -o $out/mclOutput

echo "Done with mcl."
## orthomclMclToGroups
echo "Starting orthomclMclToGroups..."
time perl $orthoMCLdir/orthomclMclToGroups \
    orthoGroup \
    1000 < $out/mclOutput > $out/groups.txt

echo "Done with orthomclMclToGroups."

# log end
end=$(date)
echo "=====> Done with OrthoMCL at $end."
