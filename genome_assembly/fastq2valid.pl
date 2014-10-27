#!/usr/bin/env perl

##NB: This is to fix sequences for prinseq. This issue is that prinseq won't 
##    accept sequences that start with a dot character.

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use File::Basename;
use Getopt::Long;
use Bio::Kseq;

#
# lexical vars
#
my $infile;
my $outfile;

GetOptions(
           'i|infile=s'    => \$infile,
           'o|outfile=s'   => \$outfile,
          );

#
# check input
#
if (!$infile || !$outfile) {
    usage();
    exit(1);
}

my $knseq = Bio::Kseq->new($infile);
my $nt_it = $knseq->iterator;

open my $out, '>', $outfile;

while (my $nseq = $nt_it->next_seq) {
    my $seq = $nseq->{seq};
    my $qual = $nseq->{qual};
    if ($seq =~ /^\./) {
	$seq =~ s/^.//;
	$qual =~ s/^.//;
    }
    say $out "@".$nseq->{name};
    say $out $seq;
    say $out "+";
    say $out $qual;
}
close $out;

sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-i] [-o] 

Required:
     -i|infile         :      A fastq file.
     -o|outfile        :      A prinseq-compliant fastq file.

Options:
    -h|help           :       Print a usage statement (Not implemented).
    -m|man            :       Print the full documentation (Not implemented).

EOF
}
