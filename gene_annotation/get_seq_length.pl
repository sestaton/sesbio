#!/usr/bin/env perl

# TODO: Add minimal POD
use 5.010;
use strict;
use warnings;
use Getopt::Long;
       
#
# Vars
#
my $usage = "\nUSAGE: perl $0 -i infile -o outfile

The outfile will contain a summary of all the contig lengths formatted as:

read_name\tlength(bp)

";

my $infile; 
my $outfile;

GetOptions(
          'i|infile=s'  => \$infile,
          'o|outfile=s' => \$outfile,
          );

if (!$infile || !$outfile) {
    die "\nERROR: Command line not parsed correctly. Exiting.\n",$usage;
}

my @aux = undef;
my ($name, $comm, $seq, $qual);

open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

my $fh = get_fh($infile);

while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
    say join "\t", $name, length($seq);
}

close $fh;
close $out;

exit;

#
# methods
#
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
