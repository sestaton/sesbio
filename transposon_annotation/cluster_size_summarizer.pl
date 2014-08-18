#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);

my $usage = "$0 infile\n";
my $infile = shift or die $usage;

summarize_clusters($infile);

exit;


sub summarize_clusters {
    my $cls_file = shift;

    my $seqtot;

    {
        local $/ = '>';
        
        open my $in, '<', $cls_file;   
        while (my $line = <$in>) {
            $line =~ s/>//g;
            next if !length($line);
            my ($clsid, $seqids) = split /\n/, $line;
            my ($id, $seqnum) = split /\s+/, $clsid;
	    $seqtot += $seqnum;
	    say join "\t", $id, $seqnum;
        }
        close $in;
	say join "\t", "Total clustered reads in $cls_file: ", $seqtot;
    }
}
