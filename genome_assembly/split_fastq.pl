#!/usr/bin/env perl

##TODO: add POD

use v5.14;
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

my $knseq = Bio::Kseq->new($infile);
my $nt_it = $knseq->iterator;

open(my $out, '>', $outfile);

while (my $nseq = $nt_it->next_seq) {
    my $seq = $nseq->{seq};
    my $qual = $nseq->{qual};
    if ($seq =~ /^\./) {
	$seq =~ s/^.//;
	$qual =~ s/^.//;
    }
    print $out "@".$nseq->{name},"\n";
    print $out $seq,"\n";
    print $out "+\n";
    print $out $qual,"\n";
}
close($out);


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
