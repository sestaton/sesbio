#!/usr/bin/env perl

##NB: The assumption is that you have one large scaffold,
##    or pseudomolecule in FASTA format that you want to break into pieces
##    of a defined size.
##    The code is modifed from: http://stackoverflow.com/a/16660010/1543853

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;

my $usage = "$0 -i seqfile -l length -c count";
my $infile;
my $outfile;
my $size;
my $count;
my $mark  = 'X';
my %substrings;

GetOptions(
           'i|infile=s'  => \$infile,
           'o|outfile=s' => \$outfile,
           'l|length=i'  => \$size,
           'c|count=i'   => \$count,
    );

die $usage if !$infile or !$outfile;

$size //= 20;
$count //= 20;
my ($seqid, $seq, @seqparts);

open my $in, '<', $infile;
open my $out, '>', $outfile;

{
    local $/ = '>';

    while (my $line = <$in>) {
	chomp $line;
	($seqid, @seqparts) = split /\n/, $line;
	$seq = join '', @seqparts;
	next unless defined $seqid && defined $seq;
	$seqid =~ s/,//g;
	$seqid =~ s/\s+/_/g;
    }
}
close $in;

if (2*$size*$count-$size-$count >= length($seq)) {
    die "selection may not complete; choose a shorter length or fewer substrings, ".
	"or provide a longer input string\n";
}

my @substrings;

while (@substrings < $count) {
    my $pos = int rand(length($seq)-$size+1);
    my $end = $pos+$size-1;
    my $id = $seqid."_".$pos."_".$end;
    push @substrings, { $id => substr($seq, $pos, $size, $mark x $size) }
	if substr($seq, $pos, $size) !~ /\Q$mark/;
}

for my $seqh (@substrings) {
    for my $seq (keys %$seqh) {
	say $out join "\n", ">".$seq, $seqh->{$seq};
    }
}
close $out;


