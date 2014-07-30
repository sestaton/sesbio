#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $infile;
my $outfile;
my $help;
my $man;

GetOptions(
	   'i|infile=s'  => \$infile,
	   'o|outfile=s' => \$outfile,
	   'h|help'      => \$help,
	   'm|man'       => \$man,
	   );

usage() and exit(0) if $help;

if (!$infile || !$outfile) {
    say "\nERROR: Command line not parsed correctly. Check input.\n";
    usage();
    exit(1);
}

open my $in, '<', $infile or die "\nERROR: Could not open file: $infile\n";
open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

my @aux = undef;
my ($name, $comm, $seq, $qual);

while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
    while ($seq =~ /([atcg])/) {
	$seq =~ s/$1/N/g;
    }
    $seq =~ s/.{60}\K/\n/g;
    say $out join "\n", ">".$name, $seq;
}

#
# methods
#
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
    print STDERR<<EOF
USAGE: $script [-h] [-m] [--version]

Required:
    -i|infile         :      The softmasked fasta file.
    -o|outfile        :      A file to place the hardmasked fasta sequences.

Options:
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation (NOT IMPLEMENTED).

EOF
}
