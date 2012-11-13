#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;

# lexical vars
my $idlist;
my $fas_in;
my $fas_out;

GetOptions(
	   'id|idlist=s'    => \$idlist,
	   'fi|fasta_in=s'  => \$fas_in,
	   'fo|fasta_out=s' => \$fas_out,
	   );

# check @ARGV
if (!$fas_in || !$fas_out || !$idlist) {
    usage();
    exit(1);
}

open(my $out, ">", $fas_out) or die "\nERROR: Could not open file: $fas_out\n";

my $idlist_hash = id2hash($idlist);
my $fa_hash = seq2hash($fas_in);

foreach my $gene (keys %$idlist_hash) {
    if (exists $fa_hash->{$gene}) {
	print $out join("\n",(">".$gene,$fa_hash->{$gene})),"\n";
    }
}
close($out);

exit;

#
# Subs
#
sub id2hash {
    my $idlist = shift;
    open(my $fh, '<', $idlist) or die "\nERROR: Could not open file: $!\n";

    my %hash;
    while(<$fh>) {
	chomp;
	$hash{$_} = 1;
    }
    close($fh);
    return(\%hash);
}

sub seq2hash {
    my $fas = shift;
    open(my $f, '<', $fas) or die "\nERROR: Could not open file: $fas\n";

    my ($name, $seq, $qual);
    my @aux = undef; 
    my %seqhash;
    my $seqct = 0;
    
    while (($name, $seq, $qual) = readfq(\*$f, \@raux)) {
	$seqct++;
	$seqhash{$name} = $seq;
    }
    
    print "$seqct sequences in $fas\n";
    return(\%seqhash);
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
    my $name = /^.(\S+.*)/? $1 : ''; 
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

sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-id] [-fi] [-fo]

Required:
    -id|idlist       :      An ID list of records (one per line).  
    -fi|fasta_in     :      A fasta/q file to pull sequences from.
    -fo|fasta_out    :      A file to write the records to.

EOF
}
