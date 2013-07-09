#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::SeqIO;
use Benchmark qw(:all);
use autodie qw(open);

my $usage = "perl $0 infile";
my $infile = shift or die $usage;

#open my $in, '<', $infile;

my $count = 50;
my $seqct = 0;
my @aux = undef;
my ($id, $seq, $qual);
cmpthese($count, {
    'readseq' => sub {
	open my $in, '<', $infile;
	while (($id, $seq) = readseq(\*$in)) {
	    $seqct++ if defined $seq;
	}
	close $in;
	#say $seqct, " seqs from readseq";
    },
    'readfq' => sub {
	open my $in, '<', $infile;
	while (($id, $seq, $qual) = readfq(\*$in, \@aux)) {
	    $seqct++ if defined $seq;
        }
	close $in;
	#say $seqct, " seqs from readfq";
    },
    'bioseqio' => sub {
	open my $in, '<', $infile;
	my $seqio = Bio::SeqIO->new(-fh => $in, -format => 'fasta');
	while (my $seq = $seqio->next_seq) {
	    $seqct++ if defined $seq->seq;
	}
	#say $seqct, " seqs from bioseqio";
    },
	  });

#close $in;

#
# subroutines
#
sub readseq {
    my ($self) = @_;
    
    local $/ = "\n>";
    return unless my $entry = $self->getline;
    chomp $entry;

    my ($id, $seq) = split /\n/, $entry, 2;
    defined $seq && $seq =~ s/>//g;
    return ($id, $seq);
}

sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!defined(@$aux));
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
    my $name;
    if (/^.?(\S+\s(\d)\S+)/) {          # Illumina 1.8+
	$name = $1."/".$2;
    }
    elsif (/^.?(\S+)/) {            # Illumina 1.3+
	$name = $1;
    } else {
	$name = '';                 # ?
    }
    #my $name = /^.(\S+)/? $1 : ''; # Heng Li's original regex
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
    return ($name, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
}
