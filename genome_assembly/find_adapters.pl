
#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
#use Data::Dump;

#
# lexical vars
#
my $infile;
my $five_prime_outfile;
my $three_prime_outfile;
my $threshold;
my $merlength;
my $help;

GetOptions(
           'i|infile=s'              => \$infile,
           'a|five_prime_outfile=s'  => \$five_prime_outfile,
	   'b|three_prime_outfile=s' => \$three_prime_outfile,
           't|threshold=i'           => \$threshold,
           'l|merlength=i'           => \$merlength,
           'h|help'                  => \$help,
           );

#
# Check @ARGV
#pod2usage( -verbose => 2 ) if $man;
usage() and exit(0) if $help;

if (!$infile  || !$five_prime_outfile || !$three_prime_outfile) {
    print "\nERROR: Command line not parsed correctly. Check input.\n\n";
    usage();
    exit(1);
}

open my $in, '<', $infile or die "ERROR: Could not open file: $infile\n";
open my $fout, '>', $five_prime_outfile or die "ERROR: Could not open file: $five_prime_outfile\n";
open my $rout, '>', $three_prime_outfile or die "ERROR: Could not open file: $three_prime_outfile\n";

my (%rseqhash, %fseqhash);

my @aux = undef;
my ($name, $comm, $seq, $qual);
$threshold //= 100;
$merlength //= 25;

while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
    my $radapter = substr($seq, -$merlength, $merlength);
    my $fadapter = substr($seq, 0, $merlength);
    $rseqhash{$radapter}++;
    $fseqhash{$fadapter}++;        
}

for my $fkey (reverse sort { $fseqhash{$a} <=> $fseqhash{$b} } keys %fseqhash) {
    if ($fseqhash{$fkey} > $threshold) {
	say $fout join "\t", $fkey, $fseqhash{$fkey};
    }
}

for my $rkey (reverse sort { $rseqhash{$a} <=> $rseqhash{$b} } keys %rseqhash) {
    if ($rseqhash{$rkey} > $threshold) {
	say $rout join "\t", $rkey, $rseqhash{$rkey};
    }
}
close $fout;
close $rout;

exit;
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
USAGE: $script [-i] [-a] [-b] [-t] [-l] [-h] [-m]

Required:
    -i|infile               :       File of reads (Fasta for Fastq format).
    -a|five_prime_outfile   :       Name of file to write the five prime adapter counts.
    -b|three_prime_outfile  :       Name of file to write the three prime adapter counts.
    
Options:
    -t|threshold            :       Set the lower threshold of counts for matches to report (Default: 100).
    -l|merlength            :       The length to search on either end of a read for adapters (Default: 25).
    -h|help                 :       Print a usage statement.

EOF
}
