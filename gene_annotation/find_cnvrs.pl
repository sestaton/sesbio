#!/usr/bin/env perl

use 5.022;
use warnings;
use Set::IntervalTree;
use Sort::Naturally qw(nsort);
use List::UtilsBy   qw(nsort_by);
use Data::Dump::Color;
use experimental 'signatures';

my $usage = "$0 cnvs";
my $cnvs  = shift or die $usage;

my $lg4_tree  = Set::IntervalTree->new;
my $lg16_tree = Set::IntervalTree->new;

#Ha16 150857390 150861558 LG16_scf7180038271797
#Ha4  174602340 174608711 LG4_scaffold_500
$lg4_tree->insert(  'Ha4',  174602340, 174608711);
$lg16_tree->insert( 'Ha16', 150857390, 150861558);

my %gene_coords = (
    'Ha4'  => $lg4_tree,
    'Ha16' => $lg16_tree,    
);

my $copies = parse_calls($cnvs, \%gene_coords);
for my $cn (nsort_by { $copies->{$_} =~ m/CN(\d+)/ and $1 } keys %$copies) {
    say join "\t", $cn, $copies->{$cn};
}
#dd $copies;

sub parse_calls ($calls, $gene_coords) {
    open my $in, '<', $calls;
    my @samples;
    my %copies;

    while (my $line = <$in>) {
        chomp $line;
	if ($line =~ /^seqnames/) {
	    my ($sq, $s, $e, $w, $st, @t) = split /\s+/, $line;
	    @samples = @t;
	    #dd \@samples and exit;
	}
	else {
	    # seqnames start end width strand
	    my ($index, $chr, $start, $end, $width, $strand, @cnvs) = split /\s+/, $line;
	    my $res = $gene_coords->{$chr}->fetch($start, $end);
	    if (@$res) {
		say join "\t", $index, $chr, $start, $end, $width, $strand, $cnvs[0];
		@copies{@samples} = @cnvs;
	    }
	}
    }
    close $in;

    return \%copies;
}

