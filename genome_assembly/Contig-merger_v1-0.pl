#!/usr/bin/perl

######################################################################
#Marten Boetzer 24-08-2010                                           #
#Contig-merger_v1-0.pl                                                   #
#This script;                                                        #
#  -multiple contig sets with -c option                              #
#  -minimal overlap between the seed and the contigs with -o option  #
#  -minimal required contigs having overlap -o with -m option        #
######################################################################


#USAGE:
#perl Contig-merger_v1-0.pl -c <contigfile> -c <contigfile> -m <k-mer> -o <min_cov> -b <output-basename>

# -c = Set of contig sequence in .fasta format. Multiple -c parameters can be used.
# -o = Contigs are cut into fragments of this size. (default = 50)
# -m = Minimum number of contigs having overlap -o required to extend an assembly (default = number of input contig sets)
# -b = base name for the output files
#-------------------------------------------------LOAD LIBRARIES
  use strict;
  use Storable;
  #require "getopts.pl";
use Getopt::Std;
  use File::Path;
  use File::Basename;
#-------------------------------------------------MAKE PARAMETERS
  use vars qw($opt_o $opt_s $opt_b $opt_c $opt_m $opt_l);
  getopts('o:s:b:c:m:l:');

  my ($min_overlap,$outfile_base)=(50,"standard_output");
  my $seplines = ("-" x 80)."\n";

  my $verbose = 0;

#-------------------------------------------------GET AND STORE INPUT PARAMETERS
  my @confiles;
  my $contigfiles = $opt_c if($opt_c);
  @confiles =split(" ", $contigfiles);
  if(!$opt_c || $#confiles < 1){
    die "At least two sets of contigs should be inserted with the -c option\n";
  }
  my $min_find = ($#confiles+1);

  $min_overlap = $opt_o if($opt_o);
  my $min_ctg_len = ($min_overlap+1);
  $min_ctg_len = $opt_l if($opt_l);
  $outfile_base = $opt_b if($opt_b);
  $min_find = $opt_m if($opt_m);
  die "Minimal coverage should be <= ".($#confiles+1)." -- fatal\n" if($min_find >($#confiles+1));

  my @confiles =split(" ", $contigfiles);

  print "Your inputs are;\n\t-o = $min_overlap\n";
  print "\t-c = $contigfiles\n";
  print "\t-m = $min_find\n";
  print "\t-l = $min_ctg_len\n";
  print "\t-b = $outfile_base\n\n";

#-------------------------------------------------PARAMETER CHECKS
  die "ERROR: Minimal overlap option (-o) should be higher than 1. Your insert is $min_overlap. Exiting...\n" if(!($min_overlap>=1) || !($min_overlap * 1 eq $min_overlap));


  my $bin;
  my $num_input = 0;
  print "shredding contig files into $min_overlap"."bp fragments\n";
  foreach my $file (@confiles){
    $num_input++;
    die "ERROR: Invalid contig file $file ...Exiting.\n" if(! -e $file);
    print "\t$file...";
    breakContigs($file);
    print "done!\n";
  }
  my $outfile = "$outfile_base.mergedcontig.fa";
  open (TIG, ">$outfile") || die "Can't write to $outfile -- fatal\n";
  print "\nusing ".keys ( %$bin)." sequences to merge contigs\n\n";

  print "assembling...";
  my $ctg = 0;
  foreach my $seq (sort {$bin->{$b}<=>$bin->{$a}} keys %$bin){#cycle through the input reads
    if(defined $bin->{$seq}){
      last if($bin->{$seq} < $min_find);
      deleteData($seq);

      my $orig_mer = length($seq);
      my $start_sequence = uc($seq);
      ($seq) = doExtension("3", $orig_mer, $seq, $min_overlap);

      my $seqrc = reverseComplement($seq);
      ($seqrc) = doExtension("5", $orig_mer, $seqrc, $min_overlap);

      if(length($seqrc) > $min_ctg_len){
        $ctg++;
        print TIG ">contig$ctg|".length($seqrc)."\n$seqrc\n";
      }
    }
  }
  close TIG;
  print "done!\n\n";

  my $time = (time - $^T);
  my $minutes = int ($time / 60);
  $time = $time % 60;
  my $summary = writesummaryfiles($outfile);
  print "\nProcess run succesfully in $minutes minutes and $time seconds\n\n";



### FUNCTION TO READ THE ADDITIONAL CONTIGS TO A HASH
sub breakContigs{
  my ($file) = @_;
  my ($counter, $totalNotFound, $seq) = (0,0,'');
  open(IN,$file);
  while(<IN>){
    s/\r\n/\n/;
    chomp;
    my $line = $_;
    $seq.= uc($line) if(eof(IN));
    if (/\>(\S+)/ || eof(IN)){
      if ($seq ne '' && length($seq) >= $min_overlap){
        my $subct = 0;
        while($subct < length($seq)-$min_overlap){
          my $subnor = substr($seq, $subct, $min_overlap+1);
          if(index($subnor, "N") < 0){
            if(defined $bin->{$subnor}){
              $bin->{$subnor}++;
            }else{
              my $subrv =reverseComplement($subnor);
              $bin->{$subrv}++;
            }
          }
          $subct++;
        }
      }
      $seq='';
    }else{
      $seq.=uc($line);
    }
  }
  close IN;
}

###Do a de novo assembly
sub doExtension{
  my ($direction, $orig_mer, $seq, $min_overlap) = @_;
  my ($previous,$extended) = ($seq,1);
  CONSENSUS:
  while($extended){
    my $subseq = substr($seq, -$min_overlap);
    my $revseq = reverseComplement($subseq);
    my $overhang;
    $overhang->{'A'} = $bin->{$subseq."A"}+$bin->{"T$revseq"};
    $overhang->{'C'} = $bin->{$subseq."C"}+$bin->{"G$revseq"};
    $overhang->{'G'} = $bin->{$subseq."G"}+$bin->{"C$revseq"};
    $overhang->{'T'} = $bin->{$subseq."T"}+$bin->{"A$revseq"};
  
    my $coverage = $overhang->{'A'}+$overhang->{'C'}+$overhang->{'G'}+$overhang->{'T'};
    if ($coverage < $min_find){
      last CONSENSUS;
    }
    my ($ct_dna, $previous_bz) = (0, "");
    my $extend_nuc = "";
    BASE:
    foreach my $bz (sort {$overhang->{$b}<=>$overhang->{$a}} keys %$overhang){
      if($ct_dna == 1){## the two most abundant bases at that position
        if($previous_bz ne "" && ($overhang->{$previous_bz} >= $min_find) && $overhang->{$previous_bz} > $overhang->{$bz}){### a simple consensus btw top 2

          $extend_nuc = "$previous_bz";
          last BASE;
        }else{
          last CONSENSUS;
        }
      }
      $previous_bz = $bz;
      $ct_dna++;
    }
    my $checkseq = $seq . $extend_nuc;#lc($consensus);
    deleteData($subseq."$extend_nuc");
    $seq = $checkseq;
    $extended = 1;
  }###while get the OK for extension
  return $seq;
}

sub writesummaryfiles{
  my ($input_file) = @_;

  open (INFILE, $input_file) || die "Can't open input file $input_file.\n";
  my $counter=0;
  my @line;
  my $seq = "";
  my $name = "";
  my @lengths;
  my $sum = 0;
  my $totalNcount = 0;
  my $string;
  while (<INFILE>) { 
    chomp;
    $seq.=$_ if(eof(INFILE));
    if ($_ =~ /^[>]/ || eof(INFILE)) { 
      if($seq ne ""){
        $counter++;
         push(@lengths, length($seq));
         if(length($seq) < 30){
           print "counter = $seq\n";
         }
  
         $sum+= length($seq);
  
         my $Ncount = () = $seq =~ /[Nn]/g;
             # my $Ncount = $seq =~ tr/[N]//; # count the number of commas. this does not modify the string in any way
         $totalNcount += $Ncount;
         $name = "";
         $seq = "";
      }
  
      $name = $_;
    }
    else {
       $seq .= $_;
    }               
  }
  
  my $half_length = $sum/2;
  my $N25 = $half_length/2;
  my $N75 = $half_length/2+$half_length;
  
  my @lengths2 = reverse sort { $a <=> $b } @lengths;
  
  my $sumN50 = 0;
  my ($foundN50, $foundN25, $foundN75) = (0,0,0);
  for(my $i = 0; $i <= $#lengths; $i++)
  {
    $sumN50 += @lengths2[$i];
    if($sumN50 >= $half_length && $foundN50 == 0){
      $foundN50 = @lengths2[$i];
    }
    if($sumN50 >= $N25 && $foundN25 == 0){
      $foundN25 = @lengths2[$i];
    }
    if($sumN50 >= $N75 && $foundN75 == 0){
      $foundN75 = @lengths2[$i];
    }
  }
  print "Results summary: \n";
  print "Total number of contigs =", $counter, "\n";
  print "Sum (bp) = ", $sum, "\n";
  print "Max contig size = ", @lengths2[0],"\n";
  
  print "Min contig size = ". @lengths2[$#lengths]."\n";
  print "Average contig size = ".int($sum/$counter)."\n";
  print "N25 = ", $foundN25, "\n";
  print "N50 = ", $foundN50, "\n";
  print "N75 = ", $foundN75, "\n";
  
  close (INFILE);
}


###DELETE READ DATA IF IT HAS BEEN USED FOR EXTENDING A CONTIG
sub deleteData {
   my ($sequence) = @_;
   my $comp_seq = reverseComplement($sequence);
   delete $bin->{$comp_seq};
   delete $bin->{$sequence};
}
###PRINTS A COUNTER ON THE SCREEN AND OVERWRITES PREVIOUS LIN
sub CounterPrint{
  my $countingMessager = shift;
  print "\r$countingMessager";
  $|++;
}

###FUNCTION TO REVERSE COMPLEMENT A SEQUENCE
sub reverseComplement{
   $_ = shift;
   tr/ATGCatgc/TACGtacg/;
   return (reverse());
}
