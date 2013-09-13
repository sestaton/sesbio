#!/usr/bin/env perl

# Output of cmpthese() for a fasta file of 1m sequences:
#
#                  s/iter    bioseqio transposome_seqio      readseq       readfq
#bioseqio            28.1          --               -1%         -88%         -98%
#transposome_seqio   27.9          1%                --         -88%         -98%
#readseq             3.39        729%              724%           --         -85%
#readfq             0.499       5523%             5493%         579%           --
#
# Output of timethese() on same file:
#
#Benchmark: timing 10 iterations of bioseqio, readfq, readseq, transposome_seqio...
#  bioseqio: 278 wallclock secs (277.15 usr +  0.29 sys = 277.44 CPU) @  0.04/s (n=10)
#    readfq:  5 wallclock secs ( 4.93 usr +  0.01 sys =  4.94 CPU) @  2.02/s (n=10)
#   readseq: 34 wallclock secs (33.53 usr +  0.19 sys = 33.72 CPU) @  0.30/s (n=10)
#transposome_seqio: 277 wallclock secs (261.90 usr + 14.62 sys = 276.52 CPU) @  0.04/s (n=10)
#
# Output of cmpthese() for a fastq file of 1m sequences:
#
#                 s/iter          bioseqio transposome_seqio            readfq
#bioseqio             174                --              -68%             -100%
#transposome_seqio   56.5              209%                --              -99%
#readfq             0.673            25792%             8292%                --
#
# Output of timethese() for a fastq file of 1m sequences:
#
#Benchmark: timing 10 iterations of bioseqio, readfq, transposome_seqio...
#  bioseqio: 1762 wallclock secs (1757.15 usr +  0.61 sys = 1757.76 CPU) @  0.01/s (n=10)
#    readfq:  6 wallclock secs ( 6.69 usr +  0.03 sys =  6.72 CPU) @  1.49/s (n=10)
#transposome_seqio: 578 wallclock secs (560.79 usr + 16.02 sys = 576.81 CPU) @  0.02/s (n=10)

use 5.010;
use strict;
use warnings;
use Bio::SeqIO;
use SeqIO;
use Benchmark qw(:all);
use autodie qw(open);

my $usage = "perl $0 infile";
my $infile = shift or die $usage;

my $count = 10;
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
