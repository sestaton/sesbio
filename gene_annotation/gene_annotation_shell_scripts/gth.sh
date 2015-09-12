#!/bin/bash

export GTHDATADIR=/home/statonse/apps/genomethreader/gth-1.6.5-Linux_x86_64-64bit/bin/gthdata
export BSSMDIR=/home/statonse/apps/genomethreader/gth-1.6.5-Linux_x86_64-64bit/bin/bssm

cd `pwd`

gthbin=/home/statonse/apps/genomethreader/gth-1.6.5-Linux_x86_64-64bit/bin
gth=$gthbin/gth
gthbssmtrain=$gthbin/gthbssmtrain
gthbssmbuild=$gthbin/gthbssmbuild

seqfile=Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta.masked
cdna=HA412_trinity.fa.gz
protein=plant_refseq_all.faa.gz
gff=H412v1r1_trinity_gth.gff3
gff_bssm=H412v1r1_trinity_gth_bssm.gff3
datadir=ha412_bssm

# gth
time $gth -genomic $seqfile \
     -cdna $cdna \
     -maskpolyatails \
     -gff3out \
     -skipalignmentout \
     -md5ids \
     -o $gff \
    -introndistri \
    -exondistri \
    -refseqcovdistri

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
