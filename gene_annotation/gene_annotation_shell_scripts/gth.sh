#!/bin/bash

set -euo pipefail

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
done

comb_gff=HA412_trinity_gth.gff3
comb_gff_sort=HA412_trinity_gth_sort.gff3
gff_bssm=HA412_trinity_gth_sort_bssm.gff3
datadir=ha412_bssm

## merge GFF, then train and build models
# script link: https://github.com/sestaton/sesbio/blob/master/gene_annotation/merge_gth_gffs.pl
perl ~/github/sesbio/gene_annotation/merge_gth_gffs.pl > $comb_gff

## sort features
gt -gff3 -sort -addintrons $comb_gff > $comb_gff_sort

# gthbssmtrain
time $gthbssmtrain -v -usedesc -seqfile $seqfile -outdir $datadir $comb_gff_sort

# gthbssmbuild
time $gthbssmbuild \
     -gtdonor \
     -gcdonor \
     -agacceptor \
     -datapath $datadir \
     -bssmfile ha412.bssm

# gth
time $gth -v \
     -bssm ha412 \
     -genomic $seqfile \
     -cdna $cdna \
     -protein $protein \
     -o $gff_bssm \
     -paralogs \
     -introndistri \
     -exondistri \
     -refseqcovdistri
