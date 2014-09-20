#!/usr/bin/env perl

# This script formats a repeat database for input to 
# Transposome: https://github.com/sestaton/Transposome

use 5.012;
use strict;
use warnings;
use Transposome::SeqIO;

my $usage = "\n$0 infile > outfile\n";
my $infile = shift or die $usage;

my $seqio = Transposome::SeqIO->new( file => $infile );

while (my $seq = $seqio->next_seq) {
    my $nt = $seq->get_seq;
    $nt =~ s/.{60}\K/\n/g;
    if ($seq->get_id =~ /^RLG/) {
        say join "\n", ">".$seq->get_id."\t"."Gypsy"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^RLC/) {
        say join "\n", ">".$seq->get_id."\t"."Copia"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^DHH/) {
	say join "\n", ">".$seq->get_id."\t"."Helitron"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^DTA/) {
	say join "\n", ">".$seq->get_id."\t"."hAT"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^DTC/) {
	say join "\n", ">".$seq->get_id."\t"."CACTA"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^DTH/) {
	say join "\n", ">".$seq->get_id."\t"."PIF/Harbinger"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^DTM/) {
	say join "\n", ">".$seq->get_id."\t"."Mutator"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^DTT/) {
	say join "\n", ">".$seq->get_id."\t"."Tc1/Mariner"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^PPP/) {
	say join "\n", ">".$seq->get_id."\t"."Penelope"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^RIL/) {
	say join "\n", ">".$seq->get_id."\t"."L1"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^RIT/) {
	say join "\n", ">".$seq->get_id."\t"."RTE"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^RLX/) {
	say join "\n", ">".$seq->get_id."\t"."Unknown_LTR"."\t"."Zea mays", $nt;
    }
    elsif ($seq->get_id =~ /^RST/) {
	say join "\n", ">".$seq->get_id."\t"."tRNA"."\t"."Zea mays", $nt;
    }
    else {
        # should never get here, but the data may be malformed
        say STDERR "\n[ERROR]: ",$seq->get_id," does not seem to match known TE superfamilies";
    }
}
