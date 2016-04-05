#!/bin/bash

cd TPSI_gff3_files*

for $filename in ./*.gff3
do
  sed -n '1,2 p' $filename > top2

  sed '/^#/d' $filename > nocomments

  cat top2 nocomments > $filename.valid

  rm top2 nocomments
done

cd ..