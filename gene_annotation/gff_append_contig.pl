#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use Sort::Naturally;
use Bio::Tools::GFF;
use Data::Dump::Color;

my $usage  = "$0 orig.gff new.gff > updated.gff";
my $infile = shift or die $usage;
my $newgff = shift or die $usage;
open my $in, '<', $infile;
open my $gff, '<', $newgff;

my @aux = undef;
my ($id, $comm, $seq, $qual);
my (%chroms, %contigs);

while (my $line = <$in>) {
    chomp $line;
    if ($line =~ /contig/) {
        my @f = split /\t/, $line;
	$contigs{$f[0]} = join "||", @f;
    }

    if ($line =~ /^##FASTA$/) {
	while (($id, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
	    $chroms{$id} = $seq;
        }
    }
}
close $in;

while (my $line = <$gff>) {
    chomp $line;
    if ($line =~ /^##\w+|^###$/) {
	say $line;
    }
    else {
	if ($line =~ /\S/) {
	    my @f = split /\t/, $line;
	    if (exists $contigs{$f[0]}) {
		my @c = split /\|\|/, $contigs{$f[0]};
		say join "\t", @c;
		say "###";
		say $line unless $line =~ /^#/;
		delete $contigs{$f[0]};
	    }
	}
	else {
	    say "##FASTA";
	    for my $chr (nsort keys %chroms) {
		$chroms{$chr} =~ s/.{60}\K/\n/g;
		say join "\n", ">$chr", $chroms{$chr};
	    }
	}
    }
}
close $gff;

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
