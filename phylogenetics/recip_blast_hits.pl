#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $query_blast;
my $subj_blast;
my $outfile;
my $help;

GetOptions(
	   'q|query=s'   => \$query_blast,
	   's|subject=s' => \$subj_blast,
	   'o|outfile=s' => \$outfile,
	   'h|help=s'    => \$help,
	   );

usage() and exit(0) if $help;

if (!$query_blast || !$subj_blast || !$outfile) {
    usage();
    exit(1);
}

my $recip_hit = 0;
my %qhash;

open my $query, '<', $query_blast or die "\nERROR: Could not open file: $query_blast\n";
open my $subj, '<', $subj_blast or die "\nERROR: Could not open file: $subj_blast\n";
open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

while (my $qline = <$query>) {
    chomp $qline;
    next if $qline =~ /^#/; # for mgblast or other variants
    my @query_fields = split /\t/, $qline;
    die "\nERROR: The input must be a tab-delimited BLAST report. Exiting.\n"
	unless @query_fields == 12;

    #$query_fields[0] =~ s/_.*//;
    my $hkey = join ",", @query_fields[0..1];
    $qhash{$hkey} = join "~~", @query_fields[2..3], @query_fields[10..11];
}

my %seen;
say $out "Query\tHit\tPID_query\tHSP_len_query\tEval_query\tBits_query\tPID_hit\tHSP_len_hit\tEval_hit\tBits_hit";

while (my $sline = <$subj>) {
    chomp $sline;
    next if $sline =~ /^#/;
    my @subj_fields = split /\t/, $sline;
    die"\nERROR: The input must be a tab-delimited BLAST report. Exiting.\n"
	unless @subj_fields == 12;
    #$subj_fields[1] =~ s/_.*//;
    next if exists $seen{$subj_fields[1]};
    my $key = join ",", @subj_fields[0,1];
    my $rkey = join ",",@subj_fields[1,0];

    if (exists $qhash{$rkey}) {
	my ($qq, $qh) = split /\,/, $rkey;
	my ($qpid, $qaln_len, $qeval, $qbits) = split /\~\~/, $qhash{$rkey};

	$recip_hit++;
	say $out join "\t", $qq,$qh,$qpid,$qaln_len,$qeval,$qbits,@subj_fields[2..3],@subj_fields[10..11];
	$seen{$subj_fields[1]} = 1;
    }
}

close $query;
close $subj;
close $out;

say "\nFound $recip_hit reciprocal hits in $query_blast and $subj_blast.\n";

exit;
#
# methods
#
sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -q blast_res_i_to_j.bln -s blast_res_j_to_i.bln -o blast_result [-h]

Required:
    -q|query         :    Tab-delimited BLAST file for query (generated with -outfmt 6, or -m 8 with legacy BLAST).
    -s|subject       :    Tab-delimited BLAST file for subject (generated with -outfmt 6, or -m 8 with legacy BLAST). 
                          (NB: It makes no difference which order you put subject and query, 
			   it only affects the order the results are printed.)
    -o|outfile       :    File name to write the reciprocal BLAST hit results to.

Options:
    -h|help          :    Print a usage statement.

END
}
