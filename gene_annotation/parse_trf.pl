#!/usr/bin/env perl

##TODO: join steps: 1) mask simple repeats, 
##                  2) hard mask, 
##                  3) filter N percent, 
##                  4) filter by length (optionally),
##                  5) run TRF, 
##                  6) parse (here)
##NB: running trf without the html output results in one output file

use 5.010;
use strict;
use warnings;
use Bio::Tools::TandemRepeatsFinder;
use Data::Dump;

my $usage  = "USAGE: perl $0 <trf.out>\n";
my $infile = shift or die $usage;
my $monomer_len;
$monomer_len //= 50;

my $trf = Bio::Tools::TandemRepeatsFinder->new( -file => $infile );

my (@copy_num, @period_size, %monomers);

while ( my $feat = $trf->next_result ) {

    ## Doesn't work, perhaps it did with a previous version of TRF?
    #my @tags = $feat->get_tag_values();
    #say join " ", @tags;

    my ($id) = $feat->seq_id();
    my ($period_size) = $feat->get_tag_values('period_size');
    my ($copy_number) = $feat->get_tag_values('copy_number');
    my ($consensus_size) = $feat->get_tag_values('consensus_size');
    my ($percent_matches) = $feat->get_tag_values('percent_matches');
    my ($percent_indels) = $feat->get_tag_values('percent_indels');
    my ($percent_a) = $feat->get_tag_values('percent_a');
    my ($percent_c) = $feat->get_tag_values('percent_c');
    my ($percent_g) = $feat->get_tag_values('percent_g');
    my ($percent_t) = $feat->get_tag_values('percent_t');
    my ($entropy) = $feat->get_tag_values('entropy');
    my ($consensus_sequence) = $feat->get_tag_values('consensus_sequence');
    my ($repeat_sequence) = $feat->get_tag_values('repeat_sequence');

    # these values generate no results with the latest TRF
    #my ($run_parameters) = $feat->get_tag_values('run_parameters');
    #my ($sequence_description) = $feat->get_tag_values('sequence_description');

    if ($period_size >= $monomer_len) {
	#say "\n# ID:                  $id";
	#say "#  period_size:          $period_size"; 
	#say "#  copy_number:          $copy_number";
	#say "#  consensus_size:       $consensus_size";
	#say "#  percent_matches:      $percent_matches";
	#say "#  percent_indels:       $percent_indels";
	#say "#  percent_a:            $percent_a";
	#say "#  percent_c:            $percent_c";
	#say "#  percent_g:            $percent_g";
	#say "#  percent_t:            $percent_t";
	#say "#  entropy:              $entropy";
	#say "#  consensus_sequence:   $consensus_sequence";
	#say "#  repeat_sequence:      $repeat_sequence\n";
	
	#push @copy_num, $copy_number;
	#push @period_size, $period_size;
	#$monomers{$consensus_sequence}++;
	say join "\t", $period_size, $copy_number;
    }
}

#dd \%monomers;
#say join "\n", @period_size;
#say join "\n", @copy_num;
