#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;
#use Getopt::Long;  # get options such as: --align, --combine, --report
                    # default would be: -l -r or something similar


my $usage = "\nUSAGE: LTRseqs2align2.pl 5primeseqs.fasta 3primeseqs.fasta

\tThe files must be formatted with cng_header2BACname2.pl at this time to
\twork with clustalw2.\n";
my $left = $ARGV[0] || die "\n$usage\n";
my $right = $ARGV[1] || die "\n$usage\n";


my %seq_in = ( 
               '5prime' => Bio::SeqIO->new('-file'   => "<$left",
                                           '-format' => 'fasta'),
               '3prime' => Bio::SeqIO->new('-file'   => "<$right",
					   '-format' => 'fasta'),
             );

#my $count;
#my @fiveIDarr;

while ( my $fiveltr = $seq_in{'5prime'}->next_seq() ) {

     my $fiveid = $fiveltr->id;
     #if ($fiveid =~ m/^5prime_/) {         # my format created by cng_header2BACname2.pl
	 $fiveid =~ s/5prime_//;

	 #push(@fiveIDarr,$fiveid);
     #} 
     #elsif ($fiveid =~ m/^LTR5\'_/) {      # Gydb format
	 #print "\n$fiveid\n";
	 #$fiveid =~ s/LTR5\'_//;
	
	 #push(@fiveIDarr,$fiveid);
     #}
     
     #foreach my $id (@fiveIDarr) {
       
	 #print $id,"\n";              # WORKS AS INTENDED HERE

	 my $LTRseqs_name = $fiveid;
         $LTRseqs_name .= "_LTRseqs.fasta";
	 my $seq_out = Bio::SeqIO->new(-format=>'fasta', -file=>">$LTRseqs_name");     # WORKS AS INTENDED HERE; but writing only first set of LTRs

	 while ( my $threeltr = $seq_in{'3prime'}->next_seq() ) {
    
	     my $threeid = $threeltr->id;
	     #if ($threeid =~ m/^3prime_/) {
		 $threeid =~ s/3prime_//;
	     #}
	     #elsif ($threeid =~ m/^LTR3\'_/) {
		 #$threeid =~ s/LTR3\'_//;
		 #print "\n$threeid\n";                       WORKS HERE 2/2/11
	     #}

	     if ($fiveid =~ /$threeid/) {		 
		 
		 #print "\n$id\n";
		 $seq_out->write_seq($fiveltr,$threeltr);
		 last;

	     }
	 
	 }
     
     #}

 }
exit;
