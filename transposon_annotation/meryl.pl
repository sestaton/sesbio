#!/usr/bin/env perl

=head1 NAME 
                                                                       
meryl.pl - Compute k-mer frequencies for a set of DNA sequences

=head1 SYNOPSIS    
 
 perl meryl.pl -i contig.fas -t target.fas -k 20 -o contig_target.gff

=head1 DESCRIPTION

 (...)

=head1 DEPENDENCIES

Non-core Perl modules used are IPC::System::Simple and Try::Tiny.

Tested with:

=over 2

=item *

Perl 5.18.0 (on Red Hat Enterprise Linux Server release 5.9 (Tikanga))

=item *

Perl 5.16.0 (on Mac OS X 10.6.8 (Snow Leopard))

=back

=head1 LICENSE

Copyright (C) 2013 S. Evan Staton

This program is distributed under the MIT (X11) License: http://www.opensource.org/licenses/mit-license.php

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, 
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies 
or substantial portions of the Software.

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

A Fasta file (contig or chromosome) to search.

=item -o, --outfile

The name a GFF3 file that will be created with the search results.

=back

=head1 OPTIONS

=over 2

=item -t, --target

A file of WGS reads to index and search against the input Fasta.

=item -k, --kmerlen

The k-mer length to use for building the index. Integer (Default: 20).

=item -s, --search

Search the input Fasta file against an existing index. The index must be
specified with this option.

=item -idx, --index

The name of the index to search against the input Fasta. Leave this option off if you want 
to build an index to search.

=item --log

Report the log number of counts instead of raw counts. This is often a good option with WGS
data because many regions have very, very high coverage.

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut 
##TODO: 1) Clean up POD
##      2) Add readfq sub and remove bioperl dep   DONE
##      3) add IPCSS and TryTiny                   DONE
##      4) test all of the changes to the code

use 5.010;
use strict;
use warnings;
use IPC::System::Simple qw(system);
use Try::Tiny;
use autodie qw(open);
use Getopt::Long;
use File::Basename;

# input sequence to search
# input sequence to index 
#
# output gff

my $infile;
my $outfile;
my $k;
my $db;
my $help;
my $man;

my $index;
my $search;
my $log;
my $quiet;

my $matches;

GetOptions(# Required
	   'i|infile=s'        => \$infile,
	   'o|outfile=s'       => \$outfile,
	   # Options
	   't|target=s'        => \$db,
	   'k|kmerlen=i'       => \$k,
	   'idx|index=s'       => \$index,
	   'search'            => \$search,
	   'log'               => \$log,
	   'quiet'             => \$quiet,
	   );

if (!$infile || !$outfile || !$index) {
    say "\nERROR: No input was given.";
    usage();
    exit(1);
}

$k //= 20;

my $meryl = findprog('meryl');
my $mapMers = findprog('mapMers-depth');

if ($search) {
    $matches = meryl_search($infile, $index, $k);   
} else {
    my ($idxname) = build_index($db, $index, $k);
    $matches = meryl_search($infile, $idxname, $k);
}

my ($seqid, $seqlen) = return_seq($infile);
 
open my $mers,'<',$matches or die "\nERROR: Could not open file: $matches\n";
open my $gff,'>',$outfile or die "\nERROR: Could not open file: $outfile\n";

say $gff "##gff-version 3";
say $gff "##sequence-region ",$seqid," 1 ",$seqlen;

my $merct = 0;

while(my $match = <$mers>) {
    chomp $match;
    my ($offset, $count) = split /\t/, $match;
    $offset =~ s/\s//g;
    $count =~ s/\s//g;
    eval { $count = sprintf("%.2f",log($count)) if $log; };
    $merct++;
    say $gff join "\t", $seqid, "meryl", "MDR", $offset, $offset, $count, ".", "+",
		     join ";", "Name=mapMers-depth_".$k."_mer","ID=mer:$merct","dbxref=SO:0000657"; 
}
close $mers;
close $gff;
unlink $matches;

exit;
#
# Subs
#
sub findprog {
    my $prog = shift;
    my $path = `which $prog 2> /dev/null`;
    chomp $path;
    if ( (! -e $path) && (! -x $path) ) {
	die "\nERROR: Cannot find $prog binary. Exiting.\n";
    } else {
	return $path;
    }
}

sub return_seq {
    my $infile = shift;

    my $in, '<', $infile;
    my %seq;  
    my $seqct = 0;

    my @aux = undef;
    my ($name, $comm, $seq, $qual);
    my ($n, $slen, $qlen) = (0, 0, 0);
    while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
	$seqct++;
	$seq{$name} = $seq;
	return ($name, length($seq));
    }
    close $in;

    if ($seqct > 1) {
	die "\nERROR: $seqct sequences present in $infile when only 1 sequence is expected. Exiting.\n";
    } 
    else {
	for (keys %seq) {
	    say "\n========> Running meryl on sequence: $_" unless $quiet;
	}
    }
}

sub build_index {
    my ($db, $indexname, $k) = @_;

    my @index = "$meryl ".
	        "-v -B ".
                "-m $k ".
		"-s $db ".
		"-o $indexname ";
		$index .= $index." 2>&1 > /dev/null" if $quiet;

    say "\n========> Creating meryl index for mersize $k for sequence: $db";
    my $exit_code;
    try {
	$exit_code = system([0..5], @index);
    }
    catch {
	say "ERROR: meryl failed with exit code $exit_code. Here is the exception: $_.\n";
    };
}

sub meryl_search {
    my ($infile, $indexname, $k) = @_;
    my ($seqfile,$seqdir,$seqext) = fileparse($infile, qr/\.[^.]*/);
    my ($indfile,$inddir,$indext) = fileparse($indexname, qr/\.[^.]*/);
    my $searchout = $seqfile."_".$indfile.".meryl_mapMers-depth.out";

    my @search = "$mapMers ".
	         "-m $k ".
		 "-mers $indexname ".
		 "-seq $infile ".
                 "> $searchout";
    say "\n========> Searching $infile with $indexname" unless $quiet;
    my $exit_code;
    try {
	$exit_code = system([0..5], @search);
    }
    catch {
	say "ERROR: meryl mapMers failed with exit code $exit_code. Here is the exception: $_.\n";
    }

    return $searchout;
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
    my ($name, $comm) = /^.(\S+)(?:\s+)(\S+)/ ? ($1, $2) : 
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

USAGE: $script -i contig.fas -t target.fas -k 20 -o contig_target.gff 

Required:
    -i|infile       :    Fasta file to search (contig or chromosome).
    -o|outfile      :    File name to write the gff to.

Options:
    -t|target       :    Fasta file of WGS reads to index.
    -k|kmerlen      :    Kmer length to use for building the index.
    -s|search       :    Just search the (--infile). Must specify an existing index.
    -idx|index      :    Name of the index (if used with --search option, otherwise leave ignore this option).
    --log           :    Return the log number of matches instead of the raw count.
    --quiet         :    Don't print progress or program output.
	
END
}
