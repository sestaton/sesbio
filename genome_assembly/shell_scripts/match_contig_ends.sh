#!/bin/csh
# match_contig_ends
##### Match ends of contigs for gap closure by performing blast searches
##### on contig ends to look for matching split genes.  Also matches
##### contig ends against entire contigs to look for false joins.
##### Can optionally screen ends of contigs for repeats before matching.
#
# Written by: James D. White, University of Oklahoma, Advanced Center for
#   Genome Technology
#
# 20070822 JDW - Add blast_pgm set variable and document additional
#                required programs.
#
# 20070618 JDW - Addded -n and -p options to use custom blastn and blastp
#                databases.
# 20061205 JDW - Changed to use blastall.  Removed all2many calls used to
#                split out contigs.  Now blasting all contigs together.
# 20020516 JDW - Removed pipe to merge_bin_contigs for merging duplicate
#                contigs.
#Date Written: Mar 25, 1999

# Other Programs Required:
#
# From ftp://ftp.ncbi.nih.gov/blast/
#   blastall (comes with the blast package)
#   formatdb (comes with the blast package)
#
# From http://repeatmasker.genome.washington.edu
#   repeatmasker (optional)
#
# From http://www.genome.ou.edu/informatics.html
#   get_contig_ends
#   match_contigs
#   Blast2table (Requires BioPerl for parsing blast output)
#
# From http://www.bioperl.org
#   BioPerl


set Date_Last_Modified="August 22, 2007"

set blast_pgm="blastall"

set blastn_db="/hd1/src/trna_work/rna.fa_new"
set blastx_db="nr"
# get values of -n and -x options in either order
if ($#argv > 1 && "$1" == "-n") then
  shift
  set blastn_db="$1"
  echo "blastn_db=$blastn_db"
  shift
  if ($#argv > 1 && "$1" == "-x") then
    shift
    set blastx_db="$1"
    echo "blastx_db=$blastx_db"
    shift
  endif
else
  if ($#argv > 1 && "$1" == "-x") then
    shift
    set blastx_db="$1"
    echo "blastx_db=$blastx_db"
    shift
    if ($#argv > 1 && "$1" == "-n") then
      shift
      set blastn_db="$1"
      echo "blastn_db=$blastn_db"
      shift
    endif
  endif
endif

if ($#argv < 2) then
  cat <<ENDHELP

Match ends of contigs for gap closure by performing blast searches
on contig ends to look for matching split genes.  Also matches
contig ends against entire contigs to look for false joins.  Creates
output files 'fasta_file'.endmatches and 'fasta_file'.contigmatches.
Can optionally screen ends of contigs for repeats before matching.

Usage: $0 [-n blastn_db] [-x blastx_db] end_len \\
           fasta_file [goto_label/repeatmasker [repeatmasker_options]]

       blastn_db  is the name of the blastn database to use.  The
       		  default is "/hd1/src/trna_work/rna.fa_new".

       blastx_db  is the name of the blastx database to use.  The
       		  default is "nr".

       end_len    is the length of the ends of the contigs to use
                  for blast searches

       fasta_file is the name of the fasta file to search

       goto_label is used to restart the program at a checkpoint
                  if it was interrupted before completion.  goto_label
                  may be omitted, or must be one of:
                  blastn_start, blastn_continue, blastn_finish,
                  blastx_start, blastx_continue, blastx_finish,
                  append_data, output_data,
                  blastc_start, blastc_continue, blastc_finish,
                  or output_data2.  For blastn_continue, blastx_continue,
                  or blastc_continue, you may edit the file 'filesn',
                  'filesx', or 'filesc', respectively, to delete the
                  names of files for which the corresponding blast search
                  has completed.

       repeatmasker is the word 'repeatmasker', which may be followed
                  by options to be passed to the repeatmasker program.
                  If repeatmasker is specified, then the program
                  repeatmasker is run to mask out common repeat sequences
                  from the ends of the contigs before looking for matches.


Date_Last_Modified: $Date_Last_Modified

ENDHELP

  exit
endif

  set endlen="$1"
  shift
  set fafile="$1"
  shift

##### Use third argument to go to label for partial continued run
  if ($#argv == 1 && "$1" != "repeatmasker") then
    set gotolabel="$1"
    shift
    goto $gotolabel
  endif

##### extract the ends of the contigs
get_ends:
  echo "extract ends of contigs"
  if (! -d blast_end_dir) mkdir blast_end_dir
  \rm blast_end_dir/* >& /dev/null
  get_contig_ends -x -s -e $endlen $fafile > blast_end_dir/${fafile}.ends
#  get_contig_ends -x -s -e $endlen $fafile | merge_bin_contigs -i -j : -  blast_end_dir/${fafile}.ends


##### extract the ends of the contigs
repeatmasker:
  if ($#argv >= 1 && "$1" == "repeatmasker") then
    echo "mask repeats"
    shift
    cd blast_end_dir
    cp ${fafile}.ends ${fafile}.ends.masking
    repeatmasker $* ${fafile}.ends.masking
    if ($status) then
      echo "repeatmasker error, exiting"
      exit $status
    endif
    if (! -e ${fafile}.ends.masking.masked) then
      echo "error: cannot find repeatmasker .masked file, exiting"
      exit 1
    endif
    cp ${fafile}.ends.masking.masked ${fafile}.ends.masked
    cd ..
    echo "end of mask repeats"
  else
    cp blast_end_dir/${fafile}.ends blast_end_dir/${fafile}.ends.masked
  endif
  

##### run blastn searches for tRNAs
blastn_start:
  echo "Prepare for blastn searches"
  cd blast_end_dir
  tr "xX" "NN" < ${fafile}.ends.masked > ${fafile}.ends.N
  cd ..

blastn_continue:
  echo "Begin blastn searches"
  cd blast_end_dir
  $blast_pgm -p blastn -i ${fafile}.ends.N -o ${fafile}.blastn -d $blastn_db -e 0.5 >& blastn.errors
  if ($status) then
    echo "Error: file=${fafile}.ends.N, status=$status, blastn failed" | tee -a ${fafile}.error 
  endif
  echo "Finished blastn search"
  cd ..

blastn_finish:
  echo "Converting blastn output files to table format"
  cd blast_end_dir
  Blast2table -format 6 ${fafile}.blastn > ${fafile}.tablen
  echo "End of Blastn searches"
  cd ..

##### run blastx searches
blastx_start:

blastx_continue:
  echo "Begin blastx searches"
  cd blast_end_dir
  $blast_pgm -p blastx -i ${fafile}.ends.N -o ${fafile}.blastx -d $blastx_db >& blastx.errors
  if ($status) then
    echo "Error: file=${fafile}.ends.N, status=$status, blastx failed" | tee -a ${fafile}.error 
  endif
  echo "Finished blastx search"
  cd ..

blastx_finish:
  echo "Converting blastx output files to table format"
  cd blast_end_dir
  Blast2table -format 6 ${fafile}.blastx  > ${fafile}.tablex
  echo "End of Blastx searches"
  cd ..

##### append blast search files
append_data:
  echo "Append blastn and blastx search files"
  cd blast_end_dir
  cat ${fafile}.tablex ${fafile}.tablen > ${fafile}.table
  cd ..

##### match blast entries and output data
output_data:
  echo "${fafile}.endmatches" > ${fafile}.endmatches
  date | echo >> ${fafile}.endmatches
  match_contigs -m 5 blast_end_dir/${fafile}.table >> ${fafile}.endmatches

##### blast contig ends against full contigs file
blastc_start:
  echo "Prepare to blast ends against full contigs"
  perl -pe "tr/Xx/NN/;s/^>.*\.Contig/>Contig/" < $fafile > blast_end_dir/CONTIGS

  cd blast_end_dir
  echo "Making Temporary Blast DB"
  formatdb -p F -i CONTIGS
  cd ..

blastc_continue:
  echo "Begin blast ends against full contigs"
  cd blast_end_dir
  $blast_pgm -p blastn -i ${fafile}.ends.N -o ${fafile}.blastc -d CONTIGS -e 0.5 >& blastc.errors
  if ($status) then
    echo "Error: file=${fafile}.ends.N, status=$status, blastc failed" | tee -a ${fafile}.error
  endif
  echo "Finished blast ends against full contigs"
  cd ..

blastc_finish:
  echo "Converting blast against full contigs output files to table format"
  cd blast_end_dir
  Blast2table -format 6 ${fafile}.blastc > ${fafile}.tablec
  echo "End of blast ends against full contigs"
  cd ..

##### match blast entries and output data
output_data2:
  echo "${fafile}.contigmatches" > ${fafile}.contigmatches
  date | echo >> ${fafile}.contigmatches
  match_contigs -m 5 blast_end_dir/${fafile}.tablec >> ${fafile}.contigmatches

##### we are done
  echo "End of blast_contig_ends"
