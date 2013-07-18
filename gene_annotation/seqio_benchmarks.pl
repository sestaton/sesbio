#!/usr/bin/env perl

# Output of cmpthese() for a fasta file of ~500k sequences:
#
#                     Rate    bioseqio transposome_seqio     readseq       readfq
#bioseqio          0.205/s          --              -82%        -88%        -100%
#transposome_seqio  1.15/s        461%                --        -32%         -98%
#readseq            1.70/s        732%               48%          --         -98%
#readfq             72.5/s      35246%             6196%       4151%           --
#
# Output of timethese() on same file:
#
#Benchmark: timing 50 iterations of bioseqio, readfq, readseq, transposome_seqio...
#  bioseqio: 245 wallclock secs (244.93 usr +  0.28 sys = 245.21 CPU) @  0.20/s (n=50)
#    readfq:  1 wallclock secs ( 0.69 usr +  0.00 sys =  0.69 CPU) @ 72.46/s (n=50)
#   readseq: 30 wallclock secs (29.27 usr +  0.18 sys = 29.45 CPU) @  1.70/s (n=50)
#transposome_seqio: 44 wallclock secs (43.88 usr +  0.20 sys = 44.08 CPU) @  1.13/s (n=50)

use 5.010;
use strict;
use warnings;
use Bio::SeqIO;
use SeqIO;       # my own SeqIO, based on readfq   
use Benchmark qw(:all);
use autodie qw(open);

my $usage = "perl $0 infile";
my $infile = shift or die $usage;

my $count = 50;
my $seqct = 0;
my @aux = undef;
my ($id, $seq, $qual);

cmpthese($count, {
    'transposome_seqio' => sub {
	 my $seqio = SeqIO->new( file => $infile );
	 my $fh = $seqio->get_fh;
	 while (my $seq = $seqio->next_seq($fh)) {
	     $seqct++ if $seq->has_seq;
	 }
    },
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
