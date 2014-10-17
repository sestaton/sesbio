#!/usr/bin/env perl

## NB: faSomeRecords from Kent source is the fastest for this task, but 
##     this script works for Fastq, as well as Fasta.

use 5.010;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

# lexical vars
my $idlist;
my $infile;
my $outfile;
my $format;

GetOptions(
	   'id|idlist=s' => \$idlist,
	   'i|infile=s'  => \$infile,
	   'o|outfile=s' => \$outfile,
	   );

# check @ARGV
if (!$infile || !$outfile || !$idlist) {
    usage();
    exit(1);
}

my @aux = undef;
my ($name, $comm, $seq, $qual);

open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

my $fh  = get_fh($infile);
my $ids = id2hash($idlist);

while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
    $seq =~ s/.{60}\K/\n/g;
    if (exists $ids->{$name}) {
	if (defined $qual) {
	    say $out join "\n", "@".$name, $seq, "+", $qual;
	}
	else {
	    say $out join "\n", ">".$name, $seq;
	}
    }
}
close $fh;
close $out;

exit;
#
# methods
#
sub id2hash {
    my $idlist = shift;
    open my $fh, '<', $idlist or die "\nERROR: Could not open file: $!\n";

    my %ids;
    while (<$fh>) {
	chomp;
	#s/\s+/_/g;
	$ids{$_} = 1;
    }
    close $fh;
    return \%ids;
}

sub get_fh {
    my ($file) = @_;

    my $fh;
    if ($file =~ /\.gz$/) {
        open $fh, '-|', 'zcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /\.bz2$/) {
        open $fh, '-|', 'bzcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /^-$|STDIN/) {
	open $fh, '< -' or die "\nERROR: Could not open STDIN\n";
    }
    else {
        open $fh, '<', $file or die "\nERROR: Could not open file: $file\n";
    }

    return $fh;
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
    print STDERR<<EOF
USAGE: $script [-id] [-i] [-o]

Required:
    -id|idlist    :      An ID list of records (one per line).  
    -i|infile     :      A fasta/q file to pull sequences from.
    -o|outfile    :      A file to write the records to.

EOF
}
