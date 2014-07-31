#!/usr/bin/env perl

##TODO: join steps: 1) mask simple repeats, 2) hard mask, 3) filter N percent, 4) filter by length (optionally),
##                  5), run TRF, 6) parse (here)
use 5.010;
use strict;
use warnings;
use Bio::Tools::TandemRepeatsFinder;

my $usage  = "USAGE: perl $0 <trf.out>\n";
my $infile = shift or die $usage;

my $trf = Bio::Tools::TandemRepeatsFinder->new( -file => $infile );

while ( my $feat = $trf->next_result ) {

    ## available tags:

    #  period_size
    #  copy_number
    #  consensus_size
    #  percent_matches
    #  percent_indels
    #  percent_a
    #  percent_c
    #  percent_g
    #  percent_t
    #  entropy
    #  consensus_sequence
    #  repeat_sequence
    #  run_parameters
    #  sequence_description

    my @tags = $feat->get_tag_values();
    say join " ", @tags;

    # print the source sequence id, start, end, percent matches, and the consensus sequence
    #my ($percent_matches)    = $feat->get_tag_values('percent_matches');
    #my ($consensus_sequence) = $feat->get_tag_values('consensus_sequence');
    #print $feat->seq_id()."\t".$feat->start()."\t".$feat->end()."\t$percent_matches\t$consensus_sequence\n";
}
