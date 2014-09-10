#!/usr/bin/env perl

##TODO: add POD

use 5.010;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

#
# lexical vars
#
my $infile;
my $outfile;
my $numreads;

GetOptions(
           'i|infile=s'    => \$infile,
           'n|numreads=s'  => \$numreads,
          );

#
# check input
#
usage() and exit(1) if !$infile;

$numreads  //= '1000000';
my $num    = 0;
my $count  = 0;
my $fcount = 1;

my ($iname, $ipath, $isuffix) = fileparse($infile, qr/\.[^.]*/);
if ($iname =~ /\.fa.*/) {
    ($iname, $ipath, $isuffix) = fileparse($iname, qr/\.[^.]*/);
}

my $out;
my @split_files;
    
my $fname = $iname."_".$fcount.$isuffix;
open $out, '>', $fname or die "\nERROR: Could not open file: $fname\n";
    
push @split_files, $fname;
my $in  = get_fh($infile);
my @aux = undef;
my ($name, $comm, $seq, $qual);

while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
    if ($count % $numreads == 0 && $count > 0) {
	$fcount++;
	$fname = $iname."_".$fcount.$isuffix;
	open $out, '>', $fname or die "\nERROR: Could not open file: $fname\n";
	
	push @split_files, $fname;
    }
    if (defined $qual) {
	say $out join "\n", "@".$name, $seq, "+", $qual;
    }
    else {
	say $out join "\n", ">".$name, $seq;
    }
    $count++;
}
close $in; 
close $out;

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

sub get_fh {
    my ($file) = @_;

    my $fh;
    if ($file =~ /\.gz$/) {
        open $fh, '-|', 'zcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /\.bz2$/) {
        open $fh, '-|', 'bzcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    else {
        open $fh, '<', $file or die "\nERROR: Could not open file: $file\n";
    }

    return $fh;
}

sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-i] [-n] [-h] [-m] 

Required:
     -i|infile         :      A FASTA/Q file to split into smaller files.
                              (the input may be compressed with gzip or bzip2).

Options:
    -n|numreads       :       The number of reads to write to each file (Default: 1000000).
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}
