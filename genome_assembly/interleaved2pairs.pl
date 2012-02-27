#!/usr/bin/perl -w


# TODO: Split fastq or fasta into separate files

use strict;

my $in = $ARGV[0];

if (!$in) {
    die "\n$0 infasta\n\n";
}

my $pair1 = $in;
$pair1 =~ s/\.[^\.]*$//;    # http://www.perlmonks.org/?node_id=729477
$pair1 .= "_pair1.fasta";
my $pair2 = $in;
$pair2 =~ s/\.[^\.]*$//;    
$pair2 .= "_pair2.fasta";

open(my $fas, '<', $in) or die "\nERROR: Could not open file: $in\n";
open(my $faspair1, '>', $pair1) or die "\nERROR: Could not open file: $pair1\n";
open(my $faspair2, '>', $pair2) or die "\nERROR: Could not open file: $pair2\n";

while(my $header = <$fas>) {
    chomp $header;
    my $seq = <$fas>;
    chomp $seq;
    print $faspair1 $header."\n".$seq."\n" if $header =~ m/1$/;
    print $faspair2 $header."\n".$seq."\n" if $header =~ m/2$/;
}

close($fas);
close($faspair1);
close($faspair2);
