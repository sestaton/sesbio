#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use File::Basename;
use Getopt::Long;
use Bio::Tools::GFF;
use Data::Dump;

my %opt;
my %tirs;
my $tirct = 0;

GetOptions(\%opt, 'infile|i=s', 'fasta|f=s');

usage() and exit(0) if !$opt{infile} or !$opt{fasta};

my $hash = seq_to_hash($opt{fasta});

my $gffio = Bio::Tools::GFF->new( -file => $opt{infile}, -gff_version => 3 );

my ($start, $end, $region, %feature);
while (my $feature = $gffio->next_feature()) {
    if ($feature->primary_tag eq 'terminal_inverted_repeat_element') {
        my @string = split /\t/, $feature->gff_string;
        ($region) = ($string[8] =~ /ID=?\s+?(terminal_inverted_repeat_element\d+)/);
        ($start, $end) = ($feature->start, $feature->end);
	$tirs{$region}{'full'} = join "-", $string[0], $feature->primary_tag, $string[3], $string[4];
    }
    next unless defined $start && defined $end;
    if ($feature->primary_tag eq 'terminal_inverted_repeat') {
	my @string = split /\t/, $feature->gff_string;
        if ($feature->start >= $start && $feature->end <= $end) {
            push @{$tirs{$region}{'tirs'}}, 
	        join "||", $string[0], $feature->primary_tag, $string[3], $string[4];
        }
    }
}

for my $tir (sort keys %tirs) {
    my ($source, $element, $start, $end) = split /\-/, $tirs{$tir}{'full'};
    my $length = ($end - $start) + 1;

    my $outfile = $tir.".fasta";
    open my $out, '>', $outfile;
    my $seq = substr $hash->{$source}, $start, $length;
    $seq =~ s/.{60}\K/\n/g;
    say $out join "\n", ">".$element."_".$start."-".$end, $seq;
    close $out;

    for my $tir_repeat (@{$tirs{$tir}{'tirs'}}) {
	my ($src, $tire, $s, $e) = split /\|\|/, $tir_repeat;
	my $tirlen = ($e - $s) + 1;
	my $tirseq = substr $hash->{$src}, $s, $tirlen;
	$tirseq =~ s/.{60}\K/\n/g;
	if ($tirct) { 
	    my $fiveprime_outfile = $tir."_5prime_tir.fasta";
	    open my $tout, '>', $fiveprime_outfile;
	    say $tout join "\n", ">".$tir."_5prime_tir_".$s."-".$e, $tirseq;
	    close $tout;
	}
	else {
	    my $threeprime_outfile = $tir."_3prime_tir.fasta";
            open my $tout, '>', $threeprime_outfile;
	    say $tout join "\n", ">".$tir."_3prime_tir_".$s."-".$e, $tirseq;
	    close $tout;
	    $tirct++;
	}
    }
    $tirct = 0;
}

sub seq_to_hash {
    my ($file) = @_;

    open my $in, '<', $file;
    my %hash;
    my @aux = undef;
    my ($name, $comm, $seq, $qual);

    while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
	$hash{$name} = $seq;
    }
    close $file;

    return \%hash;
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
  print STDERR <<END
USAGE: $script -i file.gff -f seqs.fas

Required:
 -i|infile    :    GFF file to extract TIR coordinates from
 -f|fasta     :    FASTA file to pull the TIR regions from.
    
Options:
 -h|help      :    Print usage statement (not implemented).
 -m|man       :    Print full documentation (not implemented).
END
}
