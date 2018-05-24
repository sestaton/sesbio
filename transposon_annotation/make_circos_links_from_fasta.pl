#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use File::Basename;
use File::Spec;
use Bio::DB::HTS::Kseq;
use Sort::Naturally qw(nsort);
use Scalar::Util    qw(openhandle);

my $usage = basename($0).' tephra_fragments.fasta';
my $file = shift or die $usage;
say STDERR "\nERROR: $file does not exist or is empty. Check input. Exiting.\n"
    unless -e $file && -s $file;

my ($chroms, $funder, $fover) = get_fragments($file);
my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);

my %handles;
for my $chr (keys %$chroms) {
    my $uoutfile = File::Spec->catfile($path, $name.'_'.$chr.'_under100bp.txt');
    open my $uout, '>', $uoutfile;
    $handles{under}{$chr} = $uout;

    my $ooutfile = File::Spec->catfile($path, $name.'_'.$chr.'_over5kb.txt');
    open my $oout, '>', $ooutfile;
    $handles{over}{$chr} = $oout;
}

for my $key (nsort keys %$funder) {
    my ($fchr, $fstart) = split /\|\|/, $key;
    my ($fend, $schr, $sstart, $send) = split /\|\|/, $funder->{$key};
    my $fh = $handles{under}{$fchr};
    say $fh join "\t", $fchr, $fstart, $fend, $schr, $sstart, $send;
}

for my $key (nsort keys %$fover) {
    my ($fchr, $fstart) = split /\|\|/, $key;
    my ($fend, $schr, $sstart, $send) = split /\|\|/, $fover->{$key};
    my $fh = $handles{over}{$fchr};
    say $fh join "\t", $fchr, $fstart, $fend, $schr, $sstart, $send;
}

for my $type (keys %handles) {
    for my $chr (keys %{$handles{$type}}) {
	close $handles{$type}{$chr};
    }
}

exit;
## methods
sub get_fragments {
    my ($file) = @_;

    my $kseq = Bio::DB::HTS::Kseq->new($file);
    my $iter = $kseq->iterator;

    my (%fover, %funder, %chroms);
    while (my $seqobj = $iter->next_seq) {
	my $id  = $seqobj->name;
	my $seq = $seqobj->seq;
	my $len = length($seq);
	#NB: ID format should be modified to be more abstract
	next if $id =~ /MT|CP|Chr0/;
	next unless $id =~ /fragment/;
	my ($schr, $sstart, $send, $fchr, $fstart, $fend) = 
	    ($id =~ /(MtrunA17(?:MT|CP)?(?:Chr\d+)?(?:\w+\d+)?)_(\d+)_(\d+)_fragment_.*(MtrunA17(?:MT|CP)?(?:Chr\d+)?(?:\w+\d+)?)_(\d+)_(\d+)$/);
	my $flen = $fend - $fstart + 1;
	my $key = join "||", $fchr, $fstart;
	my $val = join "||", $fend, $schr, $sstart, $send;
	$chroms{$fchr} = 1;

	if ($flen <= 100) {
	    $funder{$key} = $val;
	}
	elsif ($flen >= 5000) {
	    $fover{$key} = $val;
	}
    }

    return (\%chroms, \%funder, \%fover);
}
