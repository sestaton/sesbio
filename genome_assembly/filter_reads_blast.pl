#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;
use File::Basename;

# lexical vars
my $infile;
my $blast;
my $outfile;
my $matchlen;
my $identity;
my $help;

GetOptions(
	   'i|infile=s'             => \$infile,
	   'b|blastfile=s'          => \$blast,
	   'o|outfile=s'            => \$outfile,
	   'l|match_len=i'          => \$matchlen,
	   'pid|percent_identity=i' => \$identity,
	   'h|help'                 => \$help,
	   );

# help?
usage() and exit(0) if $help;

# check @ARGV
if (!$infile || !$blast || !$outfile) {
    usage();
    exit(1);
}

# set defaults
my $totalseqs = 0;
my $filteredseqs = 0;
$matchlen //= 50;
$identity //= 50;
my @aux = undef;
my ($name, $comm, $seq, $qual);

open my $fh, '<', $infile;
open my $out, '>', $outfile;
my $matches = filter_hits($blast, $matchlen, $identity);

while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
    $totalseqs++;
    unless (exists $matches->{$name}) {
	$filteredseqs++;
	say $out join "\n", ">".$name, $seq;
    }
}
close $fh;
close $out;

my $perc_filtered = sprintf("%0.2f", $totalseqs/$filteredseqs);

say "======== ${perc_filtered}% percent ($totalseqs".
    "/"."$filteredseqs) of the reads in $infile were filtered from $blast.";

#
# methods
# 
sub filter_hits {
    my ($blast, $matchlen, $identity) = @_;

    open my $in, '<', $blast;

    my %matches;

    while (<$in>) {
	chomp;
	my @fields = split;
	if ($fields[2] >= $identity && $fields[3] >= $matchlen) {
	    $matches{$fields[0]}++;
	}
    }
    close $in;

    return \%matches;
}
	
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

sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script -i infile -o outfile -b blasttable -l matchlen [-h] [-pid]

Required:
 -i|infile             :       A multifasta to screen for contamination.
 -b|blastfile          :       The blast table file of the input to a database of contaminants.
 -o|outfile            :       A file to put the screened sequences.

Options:
 -l|matchlen           :       Minimum length of the query to keep (Default: 50 bp).
 -pid|percent_identity :       Set the minimum percent identity (integer) threshold (Default: 50%).
 -h|help               :       Print a usage statement.

END
}
