#!/usr/bin/env perl

# This script formats a repeat database for input to 
# Transposome: https://github.com/sestaton/Transposome

use 5.012;
use strict;
use warnings;
use autodie qw(open);
use Transposome::SeqIO;
use Getopt::Long;

my $usage = "\n$0 -i infile -o outfile\n";
my $infile;
my $outfile;

GetOptions( 'i|infile=s' => \$infile, 'o|outfile=s' => \$oufile );

say $usage and exit(1) if !$infile or !$outfile;
open my $out, '>', $outfile;

my $seqio = Transposome::SeqIO->new( file => $infile );

while (my $seq = $seqio->next_seq) {
    my $nt = $seq->get_seq;
    $nt =~ s/.{60}\K/\n/g;
    if ($seq->get_id =~ /^RLG/) {
        say $out join "\n", ">".$seq->get_id."\t"."Gypsy"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^RLC/) {
        say $out join "\n", ">".$seq->get_id."\t"."Copia"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^DHH/) {
	say $out join "\n", ">".$seq->get_id."\t"."Helitron"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^DTA/) {
	say $out join "\n", ">".$seq->get_id."\t"."hAT"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^DTC/) {
	say $out join "\n", ">".$seq->get_id."\t"."CACTA"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^DTH/) {
	say $out join "\n", ">".$seq->get_id."\t"."PIF/Harbinger"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^DTM/) {
	say $out join "\n", ">".$seq->get_id."\t"."Mutator"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^DTT/) {
	say $out join "\n", ">".$seq->get_id."\t"."Tc1/Mariner"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^PPP/) {
	say $out join "\n", ">".$seq->get_id."\t"."Penelope"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^RIL/) {
	say $out join "\n", ">".$seq->get_id."\t"."L1"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^RIT/) {
	say $out join "\n", ">".$seq->get_id."\t"."RTE"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^RLX/) {
	say $out join "\n", ">".$seq->get_id."\t"."Unknown_LTR"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^RST/) {
	say $out join "\n", ">".$seq->get_id."\t"."tRNA"."\t"."Zea mays", $nt;
    }
    else {
        # should never get here, but the data may be malformed
        say STDERR "\n[ERROR]: ",$seq->get_id," does not seem to match known TE superfamilies";
    }
}
close $out;
