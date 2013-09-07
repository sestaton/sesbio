#!/usr/bin/env perl

use strict;
use warnings;
use Bio::AlignIO;

my $usage = "USAGE: cnv_2phylip.pl aln_in aln_out\n";

my $aln_file = shift or die $usage;
my $noninter = shift or die $usage;

my $aln_in = Bio::AlignIO->new(-file   => $aln_file,
			       -format => 'clustalw');

my $aln_out = Bio::AlignIO->new(-file   => ">$noninter",
				-format => 'phylip',
			        -interleave => 1);
    
while (my $aln = $aln_in->next_aln()) {
    $aln_out->write_aln($aln);
}
