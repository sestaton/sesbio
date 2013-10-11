#!/usr/bin/env perl

use v5.12;
use strict;
use warnings;
use autodie qw(open);
use Data::Dump qw(dd dump);
use List::Util qw(sum);
use Getopt::Long;

my $usage = "$0 -i annot.tsv -c cls -o out -n seqct\n";

my $infile;
my $cls;
my $seqct;
my $outfile;
my %annot;
my %sfam_hitct;
my %fam_readct;
my $total_ct = 0;

GetOptions(
	   'i|infile=s'    => \$infile,
	   'o|outfile=s'   => \$outfile,
	   'n|seqnum=i'    => \$seqct,
	   'c|cls=s'       => \$cls,
	   );

die $usage if !$infile or !$outfile or !$cls;

my $cls_tot = get_cls_tot($cls);
#my ($seqs, $seqct) = fas2hash($fas_file);
my $repfrac = $cls_tot / $seqct;

open my $in, '<', $infile;
open my $out, '>', $outfile;

while (<$in>) {
    chomp;
    next if /^Cluster/;
    my @fields = split;
    $total_ct += $fields[1];
    if (scalar @fields == 4) {
	#CL95    866     pseudogene      SSU-rRNA_Ath
	$sfam_hitct{$fields[2]}++;
	$fam_readct{$fields[3]} += $fields[1];
	if (exists $annot{$fields[2]}{$fields[3]}) {
	    push @{$annot{$fields[2]}{$fields[3]}}, $fields[1];
	}
	else {
	    $annot{$fields[2]}{$fields[3]} = [$fields[1]];
	}
    }
    #CL109   673     integrated_virus        Caulimoviridae  Caulimovirus-4_STu
    elsif (scalar @fields == 5) {
	$sfam_hitct{$fields[3]}++;
        $fam_readct{$fields[4]} += $fields[1];
        if (exists $annot{$fields[3]}{$fields[4]}) {
            push @{$annot{$fields[3]}{$fields[4]}}, $fields[1];
        }
        else {
            $annot{$fields[3]}{$fields[4]} = [$fields[1]];
        }
    }
    #CL36    4732    transposable_element    ltr_retrotransposon     Copia   ATCOPIA24I
    else {
	my $fam = $fields[5];
	$fam =~ s/\-\d.*//;
	$sfam_hitct{$fields[4]}++;
	$fam_readct{$fam} += $fields[1];
	if (exists $annot{$fields[4]}{$fam}) {
	    push @{$annot{$fields[4]}{$fam}}, $fields[1];
	}
	else {
	    $annot{$fields[4]}{$fam} = [$fields[1]];
	}
    }
}
close $in;

#dd \%annot;
#dd \%sfam_hitct;
#dd \%fam_readct;
say $out "Superfamily\tFamily\tReadCt/TotalSeqs\tPercCov\tGPerc_corr";
for my $sf (reverse sort { $sfam_hitct{$a} <=> $sfam_hitct{$b} } keys %sfam_hitct) {
    for my $f (reverse sort { $fam_readct{$a} <=> $fam_readct{$b} } keys %fam_readct) {
	if (exists $annot{$sf}{$f}) {
	    my $read_ct = sum @{$annot{$sf}{$f}};
	    my $perc_cov = sprintf("%.12f",$read_ct/$total_ct);
	    my $corr_perc_cov = $perc_cov * $repfrac;
	    say $out join "\t", $sf, $f, $read_ct."/".$total_ct, $perc_cov, $corr_perc_cov;
	}
    }
}
close $out;

#
# subs
#
sub get_cls_tot {
    my $cls = shift;
    my $cls_tot = 0;

    open my $in, '<', $cls;
    while (<$in>) {
        chomp;
        if (/^>(\S+\s(\d+))/) {
            $cls_tot += $2;
        }
    }
    return $cls_tot;
}

