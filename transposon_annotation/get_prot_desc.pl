#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use Bio::DB::HTS::Kseq;
use Data::Dump::Color;

my $usage = "$0 blasttable fas match";
my $blast = shift or die $usage;
my $fasta = shift or die $usage;
my $match = shift or die $usage;

my $lens = get_lengths($fasta);
#dd $lens and exit;

my %type = ( 'DT'  => 'DNA transposon', 
	     'RLG' => 'Gypsy', 
	     'RLC' => 'Copia', 
	     'DHH' => 'Helitron', 
	     'RIL' => 'LINE', 
	     'RLX' => 'LTR-Unclassified',
             'RLT' => 'TRIM' );

say join "\t", 'Type', 'Element_with_gene', 'Protein', 'PID', 'AlignLength', 'QueryLength', 'AlignPercent', 'Description';
open my $b, '<', $blast;
while (my $line = <$b>) {
    chomp $line;
    my @f = split /\t/, $line;
    if ($f[1] =~ /^$match/) {
	my ($len, $desc) = split /\|\|/, $lens->{$f[0]};
	my $aln_perc = sprintf("%.2f",($f[3]/$len)*100);
	if ($f[2] >= 80 && $f[3] >= 80) {
	    my ($len, $desc) = split /\|\|/, $lens->{$f[0]};
	    say join "\t", $type{$match}, @f[1,0,2,3], $len, $aln_perc, $desc;
	}
    }
}

sub get_lengths {
    my ($fas) = @_;

    my %lens;
   
    open my $f, '<', $fas;
    {
	local $/ = "\n>";
	while (my $line = <$f>) {
	    chomp $line;
	    my ($id, $seq) = split /\n/, $line, 2;
	    my ($name) = ($id =~ />?(\S+)/);
	    my ($desc) = ($id =~ /def=(.*)/);
	    defined $seq && $seq =~ s/>//g;
	    $seq =~ s/\s+|\n//g;
	    my $length = length($seq);
	    $lens{$name} = join "||", $length, $desc;
	}
    }
    close $f;
 
    return \%lens;
}
