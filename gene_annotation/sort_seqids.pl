#!/usr/bin/env perl

##TODO: don't store sequences

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

my $infile;
my $outfile;
my $usage = "$0 -i in -o sorted\n";

GetOptions(
	   'i|infile=s'   => \$infile,
	   'o|outfile=s'  => \$outfile,
	   );

if (!$infile or !$outfile) {
    say $usage,"\nERROR: No input was given. Exiting.";
    exit(1);
}

open my $in, '<', $infile;
open my $out, '>', $outfile;

my @aux = undef;
my ($name, $comm, $seq, $qual);
my %seqhash;
my $seqct = 0;
my $t0 = gettimeofday();

while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
    $seqhash{$name} = $seq;
    $seqct++;
}

for my $k (sort keys %seqhash) {
    say $out join "\n", ">".$k,$seqhash{$k};
}

close $in;
close $out;

my $t1 = gettimeofday();
my $elapsed = $t1 - $t0;
my $time = sprintf("%.2f",$elapsed);

say "========== Sorted $seqct sequences in $time seconds.";

sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!@$aux);
    return if ($aux->[1]);
    if (!defined($aux->[0])) {
        while (<$fh>) {
            chomp;
            if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                $aux->[0] = $_;
                last;
            }
        }
        if (!defined($aux->[0])) {
            $aux->[1] = 1;
            return;
        }
    }
    my ($name, $comm);
    defined $_ && do {
        ($name, $comm) = /^.(\S+)(?:\s+)(\S+)/ ? ($1, $2) : 
	                 /^.(\S+)/ ? ($1, '') : ('', '');
    };
    my $seq = '';
    my $c;
    $aux->[0] = undef;
    while (<$fh>) {
        chomp;
        $c = substr($_, 0, 1);
        last if ($c eq '>' || $c eq '@' || $c eq '+');
        $seq .= $_;
    }
    $aux->[0] = $_;
    $aux->[1] = 1 if (!defined($aux->[0]));
    return ($name, $comm, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $comm, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
}
