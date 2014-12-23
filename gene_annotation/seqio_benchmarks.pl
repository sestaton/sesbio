#!/usr/bin/env perl

# Output of cmpthese() for a fasta file of 1m sequences:
#
#                  s/iter biome_seqio bioseqio transposome_seqio base_perl readfq
#biome_seqio         54.5          --      -4%              -49%      -86%   -98%
#bioseqio            52.2          4%       --              -47%      -86%   -98%
#transposome_seqio   27.8         96%      88%                --      -73%   -96%
#base_perl           7.48        629%     598%              272%        --   -86%
#readfq              1.06       5040%    4824%             2521%      605%     --
#
# Output of timethese() on same file:
#
#Benchmark: timing 10 iterations of base_perl, biome_seqio, bioseqio, readfq, transposome_seqio...
# base_perl: 73 wallclock secs (67.23 usr +  5.50 sys = 72.73 CPU) @  0.14/s (n=10)
# biome_seqio: 546 wallclock secs (539.17 usr +  7.35 sys = 546.52 CPU) @  0.02/s (n=10)
# bioseqio: 526 wallclock secs (518.96 usr +  6.13 sys = 525.09 CPU) @  0.02/s (n=10)
# readfq: 10 wallclock secs ( 9.94 usr +  0.60 sys = 10.54 CPU) @  0.95/s (n=10)
#transposome_seqio: 281 wallclock secs (184.67 usr + 95.54 sys = 280.21 CPU) @  0.04/s (n=10)
#
# Output of cmpthese() for a fastq file of 1m sequences:
#
#                  s/iter          bioseqio transposome_seqio            readfq
#bioseqio             465                --              -93%             -100%
#transposome_seqio   31.2             1390%                --              -96%
#readfq              1.38            33534%             2157%                --
#
# Output of timethese() for a fastq file of 1m sequences:
#
#Benchmark: timing 10 iterations of bioseqio, readfq, transposome_seqio...
#  bioseqio: 4540 wallclock secs (4528.11 usr + 11.55 sys = 4539.66 CPU) @  0.00/s (n=10)
#  readfq: 14 wallclock secs (13.10 usr +  0.87 sys = 13.97 CPU) @  0.72/s (n=10)
#transposome_seqio: 314 wallclock secs (209.94 usr + 104.05 sys = 313.99 CPU) @  0.03/s (n=10)

use 5.010;
use strict;
use warnings;
use Bio::SeqIO;
use Transposome::SeqIO;
use Benchmark qw(:all);
#use lib qw(/home/jmblab/statonse/apps/perlmod/biome/blib/lib);
use Biome::SeqIO;
use autodie qw(open);

my $usage = "perl $0 infile";
my $infile = shift or die $usage;

my $count = 10;
my $seqct = 0;
my @aux = undef;
my ($id, $comm, $seq, $qual);

my $results = timethese($count, {
    'biome_seqio' => sub {
        my $seqio = Biome::SeqIO->new( file => $infile, format => 'fasta' );
        while (my $seq = $seqio->next_Seq) {
            $seqct++ if defined $seq;
        }
    },
    'transposome_seqio' => sub {
	my $seqio = Transposome::SeqIO->new( file => $infile );
	my $fh = $seqio->get_fh;
	while (my $seq = $seqio->next_seq($fh)) {
	    $seqct++ if defined $seq;
	}
    },
    'base_perl' => sub {
	open my $in, '<', $infile;
	while (($id, $seq) = base_perl(\*$in)) {
	    $seqct++ if defined $seq;
	}
	close $in;
    },
    'readfq' => sub {
	open my $in, '<', $infile;
	while (($id, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
	    $seqct++ if defined $seq;
        }
	close $in;
    },
    'bioseqio' => sub {
	open my $in, '<', $infile;
	my $seqio = Bio::SeqIO->new(-fh => $in, -format => 'fastq');
	while (my $seq = $seqio->next_seq) {
	    $seqct++ if defined $seq;
	}
	close $in;
    },
});

cmpthese( $results );

#
# subs
#
sub base_perl {
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

