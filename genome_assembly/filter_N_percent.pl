#!/usr/bin/env perl

##TODO: Add POD, help/man options.

use 5.010;
use strict;
use warnings;
use autodie;
use File::Basename;
use Getopt::Long;
use Bio::DB::HTS::Kseq;

#
# lexical vars
#
my $infile;
my $outfile;
my $percent;
my $help;
my $man;

GetOptions(
           'i|infile=s'    => \$infile,
           'o|outfile=s'   => \$outfile,
           'p|percent=i'   => \$percent,
           'h|help=s'      => \$help,
           'm|man=s'       => \$man,
          );

#
# check input
#
usage() and exit(0) if $help;
usage() and exit(1) if !$infile or !$outfile;

$percent //= '50';

my $knseq = Bio::DB::HTS::Kseq->new($infile);
my $nt_it = $knseq->iterator;

my $out = get_outfh($outfile);

while (my $nseq = $nt_it->next_seq) {
    my $seq = $nseq->{seq};
    my $seq_len = length($seq);
    my $n_count = ($seq =~ tr/Nn//);
    my $n_perc  = sprintf("%.2f",($n_count/$seq_len)*100);
    if ($n_perc <= $percent) {
	say $out join "\n", "@".$nseq->{name}, $seq, "+", $nseq->{qual};
    }
}
close $out;

# methods
sub get_outfh {
    my ($file) = @_;

    my $fh;
    if ($file =~ /^-$|STDOUT/i) {
	open $fh, '>&', \*STDOUT or die "\nERROR: Could not open STDOUT\n";
    }
    else {
	open $fh, '>', $file or die "\nERROR: Could not open file: $file\n";
    }

    return $fh;
}

sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-i] [-o] [-l] [-h] [-m]

Required:
     -i|infile         :      A fastq file.
     -o|outfile        :      A fastq file with trailing Bs as quality encodings removed.

Options:
    -p|percent        :       The threshold of N% above which reads will be discarded (Default: 50).
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation (NB: not yet implemented).

EOF
}
