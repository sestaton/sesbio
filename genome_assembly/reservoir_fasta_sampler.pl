#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use autodie qw(open);

my $usage = "$0 -i infile -o outfile -n num\n";
my $infile;
my $outfile;
my $k;

GetOptions(
	   'i|infile=s'   => \$infile,
	   'o|outfile=s'  => \$outfile,
	   'n|numseqs=i'  => \$k,
	   );

die $usage if !$infile or !$outfile or !$k;

open(my $in, '<', $infile);
open(my $out, '>', $outfile);

# Use reservoir sampling to select $k random sequences:
my @samples;
my $n = 0;  # total number of sequences read
my $i;      # index of current sequence
while (<$in>) {
    if (/^\>/) {
        # Select a random sequence from 0 to $n-1 to replace:
        $i = int rand ++$n;
        # Save all samples until we've accumulated $k of them:
        $samples[$n-1] = $samples[$i] if $n <= $k;
        # Only actually store the new sequence if it's one of the $k first ones:
        $samples[$i] = [] if $i < $k;
    }
    push @{ $samples[$i] }, $_ if $i < $k;
}
close($in);

warn "Only read $n < $k sequences, selected all.\n" if $n < $k;
warn "Selected $k out of $n sequences (", 100 * $k / $n, "%).\n" if $n >= $k;

# Print sampled sequences:
print $out @$_ for @samples;
close($out);
