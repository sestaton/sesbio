#!/bin/bash

set -e
set -u
set -o pipefail

export GTHDATADIR=/home/statonse/apps/genomethreader/gth-1.6.5-Linux_x86_64-64bit/bin/gthdata
export BSSMDIR=/home/statonse/apps/genomethreader/gth-1.6.5-Linux_x86_64-64bit/bin/bssm

cd `pwd`

gthbin=/home/statonse/apps/genomethreader/gth-1.6.5-Linux_x86_64-64bit/bin
gth=$gthbin/gth
gthbssmtrain=$gthbin/gthbssmtrain
gthbssmbuild=$gthbin/gthbssmbuild

for file in ./Ha{1..17}.fa
do
    filebase=$(echo ${file%.*})
    seqfile=$file
    cdna=HA412_trinity.fa.gz
    protein=plant_refseq_non-hypothetical.faa
    gff=${filebase}_trinity_gth.gff3
    log=${filebase}_gth_phase1.out
    #gff_bssm=${filebase}_trinity_gth_bssm.gff3
    #datadir=ha412_bssm
    
    # gth
    time $gth -v \
	 -genomic $seqfile \
	 -cdna $cdna \
	 -maskpolyatails \
	 -gff3out \
	 -skipalignmentout \
	 -md5ids \
	 -o $gff \
	 -introndistri \
	 -exondistri \
	 -refseqcovdistri > $log
    
    # gthbssmtrain
    #time $gthbssmtrain -seqfile $seqfile \
    #	      -outdir $datadir $gff
    
    # gthbssmbuild
    #time $gthbssmbuild -gtdonor \
    #	      -gcdonor \
    #	      -agacceptor \
    #	      -datapath $datadir \
    #	      -bssmfile ha412.bssm

    # gth
    #time $gth -bssm ha412 \
    #     -genomic $seqfile \
    #     -cdna $cdna \
    #     -protein $protein \
    #     -o $gff_bssm \
    #     -paralogs \
    #     -introndistri \
    #     -exondistri \
    #     -refseqcovdistri
done
