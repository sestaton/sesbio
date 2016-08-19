#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);

my $usage = "perl $0 seqfile.fastq\n";
my $infile = shift or die $usage;

open my $in, '<', $infile;
my @aux = undef;
my ($name, $comm, $seq, $qual);

my ($n_count, $total, $seq_ct, $n_ratio, $n_ratio_t, $n_perc, $n_perc_t, $gap_tot) = (0, 0, 0, 0, 0, 0, 0, 0);

my $gaplen = 73;
my $gapstr = 'N' x $gaplen;

say "Seq_ID\tLength\tfake-gap_count\tfake-gap_ratio\tfake-gap_perc";

while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
    my $seqlength = length($seq);
    if ($seqlength > 0) {
        $seq_ct++;
        $total += $seqlength;
	$n_count = () = $seq =~ /[catg]$gapstr[catg]/ig;
	$gap_tot = $n_count * $gaplen; # 73 is the artificial gap size
        $n_perc  = sprintf("%.2f",$gap_tot/$seqlength);
        $n_ratio = sprintf("%.2f",$gap_tot/($seqlength - $gap_tot));
        $n_ratio_t += $n_ratio;
        $n_perc_t  += $n_perc; 
        say join "\t", $name, $seqlength, $n_count, $n_ratio, $n_perc;
    }
}
close $in;

my $len_ave     = sprintf("%.0f",$total/$seq_ct);
my $n_ratio_ave = sprintf("%.2f",$n_ratio_t/$seq_ct);
my $n_perc_ave  = sprintf("%.2f",$n_perc_t/$seq_ct);

say "\nAverage length, N-ratio and N-percent for $infile:   $len_ave\t$n_ratio_ave\t$n_perc_ave";

#
# methods
#
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
