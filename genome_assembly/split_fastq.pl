#!/usr/bin/env perl

##TODO: add POD
##NB: For 4 line records, head -N is preferred for simplicity

use 5.014;
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
my $numreads;

GetOptions(
           'i|infile=s'    => \$infile,
           'o|outfile=s'   => \$outfile,
           'n|numreads=s'  => \$numreads,
          );

#
# check input
#
if (!$infile || !$outfile) {
    usage();
    exit(1);
}

$numreads //= '20000000';
my $num = 0;

my $knseq = Bio::Kseq->new($infile);
my $nt_it = $knseq->iterator;

open my $out, '>', $outfile;

while (my $nseq = $nt_it->next_seq) {
    $num++;
    if ($num <= $numreads) {
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
}
close $out;


#
#
#
sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-i] [-o] 

Required:
     -i|infile         :      A fastq file.
     -o|outfile        :      A prinseq-compliant fastq file.

Options:
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}
