#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use Sort::Naturally;
use Bio::DB::HTS::Kseq;

my $usage = "$0 seq";
my $seqfile = shift or die $usage;

my $lengths = get_lengths($seqfile);

my $idx = 1;
for my $id (nsort keys %$lengths) {
    next if $id =~ /Chr00/; # unplaced scaffolds
    #chr - chr01 1 0 153905722 black
    say join q{ }, 'chr', '-', $id, $idx, '0', $lengths->{$id}, 'black'; 
    $idx++;
}

exit;
#
# methods
#
sub get_lengths {
    my ($seqfile) = @_;

    my %lengths;
    my $kseq = Bio::DB::HTS::Kseq->new($seqfile);
    my $iter = $kseq->iterator;
    
    while (my $seqobj = $iter->next_seq) {
	my $seq = $seqobj->{seq};
	my $id  = $seqobj->{name};
	my $length = length($seq);

	$lengths{$id} = $length;
    }

    return \%lengths;
}
