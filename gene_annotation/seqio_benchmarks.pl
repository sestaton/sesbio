#!/usr/bin/env perl

# Output of cmpthese() for a fasta file of ~500k sequences:
#
#            s/iter bioseqio  readseq   readfq
#bioseqio     14.5       --     -88%    -100%
#readseq      1.77     719%       --     -98%
#readfq   3.98e-02   36372%    4356%       --
#
# Output of timethese() on same file:
#
#Benchmark: timing 50 iterations of bioseqio, readfq, readseq...
#  bioseqio: 731 wallclock secs (726.93 usr +  1.18 sys = 728.11 CPU) @  0.07/s (n=50)
#    readfq:  2 wallclock secs ( 2.06 usr +  0.02 sys =  2.08 CPU) @ 24.04/s (n=50)
#   readseq: 91 wallclock secs (89.63 usr +  0.70 sys = 90.33 CPU) @  0.55/s (n=50)

use 5.010;
use strict;
use warnings;
use Bio::SeqIO;
use Benchmark qw(:all);
use autodie qw(open);

my $usage = "perl $0 infile";
my $infile = shift or die $usage;

my $count = 50;
my $seqct = 0;
my @aux = undef;
my ($id, $seq, $qual);

timethese($count, {
    'readseq' => sub {
	open my $in, '<', $infile;
	while (($id, $seq) = readseq(\*$in)) {
	    $seqct++ if defined $seq;
	}
	close $in;
    },
    'readfq' => sub {
	open my $in, '<', $infile;
	while (($id, $seq, $qual) = readfq(\*$in, \@aux)) {
	    $seqct++ if defined $seq;
        }
	close $in;
    },
    'bioseqio' => sub {
	open my $in, '<', $infile;
	my $seqio = Bio::SeqIO->new(-fh => $in, -format => 'fasta');
	while (my $seq = $seqio->next_seq) {
	    $seqct++ if defined $seq->seq;
	}
	close $in;
    },
	 });

#
# subs
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
