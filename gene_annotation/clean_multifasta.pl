#!/usr/bin/env perl

=head1 NAME 
                                                                       
 clean_multifasta.pl - remove characters that break EMBOSS and other programs 

=head1 SYNOPSIS    

 clean_multifasta.pl -i weirdseqs.fas -o cleanseqs.fas 

=head1 DESCRIPTION
                                                                   
 Spaces, semicolons, commas, parenthesis and other characters cause the
 fasta header to not be parsed correctly when using EMBOSS applications.
 This leads to sequences being skipped by the application and sequence
 IDs being modified without intent. This script also replaces non-ATGCN
 characters with an "N."

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
evan at evanstaton dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The multi-fasta file to be cleaned.

=item -o, --outfile

A file to place the cleaned sequences.

=back

=head1 OPTIONS

=over 2

=item -l, --idlength

The length the ID should be truncated to. This is useful for sequences from public
databases which may have many identifers, long descriptions, and odd characters.

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut      

#
# Includes
#
use 5.010;
use strict;
use warnings;
use File::Basename;
use Time::HiRes qw(gettimeofday);
use Getopt::Long;
use Pod::Usage;

# 
# Variables with scope
#
my $infile;
my $outfile;
my $length;
my $help;
my $man;

GetOptions(# Required
	   'i|infile=s'   => \$infile,
	   'o|outfile=s'  => \$outfile,
	   # Options
	   'l|idlength=i' => \$length,
	   'h|help'       => \$help,
	   'm|man'        => \$man,
	   );

pod2usage( -verbose => 2 ) if $man;

usage() and exit(0) if $help;

if (!$infile || !$outfile) {
    say "\nERROR: No input was given.\n";
    usage();
    exit(1);
}

open my $in, '<', $infile or die "ERROR: Could not open file: $infile\n";
open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

# counters
my $t0 = gettimeofday();
my $fasnum = 0;
my $headchar = 0;
my $non_atgcn = 0;

my @aux = undef;
my ($name, $comm, $seq, $qual);
my ($n, $slen, $qlen) = (0, 0, 0);

while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
    $fasnum++;
    if ($name =~ /\s+|\;|\:|\(|\)|\[|\]|\./ || (defined $comm && $comm =~ /\s+|\;|\:|\(|\)|\[|\]|\./)) {
	$headchar++;
	$name =~ s/\s+|\;|\:|\(|\)|\[|\]|\./_/g;
	$comm =~ s/\s+|\;|\:|\(|\)|\[|\]|\./_/g if $comm;
	$headchar++ if $comm;

	$name = join '_', $name, $comm if $comm;
	# TODO: shorten the header, optionally
	if ($length) {
	    $name = substr($name, 0, $length);
	}
    }

    if ($seq =~ /\r|\n/) { # Mac
	$seq =~ s/\r?//;   # PC
	$seq =~ s/\n//;    # Unix
    }
 
    my @nt = split //, $seq;
    for my $base (@nt) {
	if ($base !~ /[ACGTNacgtn]/) {
	    $non_atgcn++;
	    $base = "N";
	}
    }       
    my $dna = join '', @nt;
    $dna =~ s/\s//g;

    my $nucleic_bc = ($dna =~ tr/NACGTacgtn//);
    my $nonnucleic = (length($dna) - $nucleic_bc);
    #$dna =~ s/(.{60})/$1\n/gs;        
    $dna =~ s/.{60}\K/\n/g;

    say $out join "\n", ">".$name, $dna;
}
close $in;
close $out;

my $t1 = gettimeofday();
my $elapsed = $t1 - $t0;
my $time = sprintf("%.2f",$elapsed);

say "\n========== Done. $fasnum sequences read and cleaned in $time seconds."; 
say "========== $non_atgcn Non-ATGCN characters changed in sequence. $headchar characters changed in headers.\n";

exit;
#
# subs
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
    my ($name, $comm) = /^.(\S+)(?:\s+)(.*)/ ? ($1, $2) : 
	                /^.(\S+)/ ? ($1, '') : ('', '');
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
  print STDERR <<END
USAGE: $script -i seqsin.fas -o filterseqs.fas 

Required:
    -i|infile    :    Fasta file of reads or contigs to filter.
    -o|outfile   :    File to place the filtered reads or contigs.
    
Options:
    -l|idlength  :    The length to truncate the ID to (Default: print full ID).
    -h|help      :    Print usage statement.
    -m|man       :    Print full documentation.
END
}
