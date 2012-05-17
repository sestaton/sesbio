#!/usr/bin/perl -w

use strict;

my $usage = "$0 infile > out\n";
my $infile = shift or die $usage;

open( my $in, '<', $infile) or die "\nERROR: Could not open file: $infile\n";

while(my $line = <$in>) {
    chomp $line;
    my ($header, $seq) = split(/\t/,$line);
    $seq =~ s/(.{60})/$1\n/gs;  
    print join("\n",">".$header,$seq), "\n";
}

close($in);

exit;
