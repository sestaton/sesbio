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

open(my $out, '>', $outfile);

while (my $nseq = $nt_it->next_seq) {
    my $seq = $nseq->{seq};
    my $qual = $nseq->{qual};
    my $qual_len = length($qual);
    if ($qual =~ /(B*)$/) {
	my $b_len = length($1);
	my $good_len = $qual_len - $b_len;
	my $no_b_qual = substr($qual,0,$good_len);
	my $no_b_seq = substr($seq,0,$good_len);
	print $out "@".$nseq->{name},"\n";
	print $out $no_b_seq,"\n";
	print $out "+\n";
	print $out $no_b_qual,"\n";
    }
    else {
	print $out "@".$nseq->{name},"\n";
	print $out $seq,"\n";
	print $out "+\n";
	print $out $qual,"\n";
    }
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
     -o|outfile        :      A fastq file with trailing Bs as quality encodings removed.

Options:
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}
