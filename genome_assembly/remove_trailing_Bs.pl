#!/usr/bin/env perl

##TODO: Add POD, help/man options.

use 5.010;
use strict;
use warnings;
use autodie;
use File::Basename;
use Getopt::Long;
use Bio::Kseq;

#
# lexical vars
#
my $infile;
my $outfile;
my $length;
my $help;
my $man;

GetOptions(
           'i|infile=s'    => \$infile,
           'o|outfile=s'   => \$outfile,
           'l|length=s'    => \$length,
           'h|help=s'      => \$help,
           'm|man=s'       => \$man,
          );

#
# check input
#
usage() and exit(0) if $help;
usage() and exit(1) if !$infile or !$outfile;

$length //= '50';

# create Bio::Kseq objects
my $knseq = Bio::Kseq->new($infile);
my $nt_it = $knseq->iterator;

open my $out, '>', $outfile;

while (my $nseq = $nt_it->next_seq) {
    my $seq = $nseq->{seq};
    my $qual = $nseq->{qual};
    my $qual_len = length($qual);
    if ($qual =~ /(B+)$/) {
	my $b_len = length($1);
	my $good_len = $qual_len - $b_len;
	my $no_b_qual = substr($qual,0,$good_len);
	my $no_b_seq  = substr($seq,0,$good_len);
	my $no_b_seq_len = length($no_b_seq);
	if ($no_b_seq_len >= $length) {
	    say $out "@".$nseq->{name};
	    say $out $no_b_seq;
	    say $out "+";
	    say $out $no_b_qual;
	}
    }
    else {
	say $out "@".$nseq->{name};
	say $out $seq;
	say $out "+";
	say $out $qual;
    }
}
close $out;

#
# subs
#
sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-i] [-o] [-l] [-h] [-m]

Required:
     -i|infile         :      A fastq file.
     -o|outfile        :      A fastq file with trailing Bs as quality encodings removed.

Options:
    -l|length         :       Minimum length threshold (Default: 50bp).
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation (NB: not yet implemented).

EOF
}
