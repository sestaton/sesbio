#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Getopt::Long;
use autodie qw(open);
use Data::Dump;

my $usage = "$0 -i annotations.tsv -c cls_list -o outfile\n";
my $anno_file;
my $cls_file;
my $out_file;

GetOptions(
           'i|infile=s'   => \$anno_file,
           'c|cls_file=s' => \$cls_file,
           'o|out_file=s' => \$out_file,
          );

die $usage if !$anno_file or !$cls_file or !$out_file;
my %clsh;
my %annoh;
my $total_clustered = 0;
my $total_clusters = 0;
my $readtot = 500000;

my $colormap = { Gypsy      => "darkgreen",
		 Copia      => "aquamarine4",
		 R1         => "black",
		 Helitron   => "darkolivegreen3",
		 hAT        => "darkkhaki",
		 pseudogene => "azure",
		 SAT        => "azure",
		 EnSpm      => "azure2",
                 Kiri       => "azure2",
		 Penelope   => "cyan",
		 DIRS       => "lightgreen",
		 Polinton   => "chartreuse" };

open my $an, '<', $anno_file;
open my $cl, '<', $cls_file;
open my $out, '>', $out_file;

while (<$an>) {
    chomp;
    next if /^Cluster/;
    my @fields = split;
    #$total_ct += $fields[1];
    if (scalar @fields == 5) {
        #CL95    866     pseudogene      SSU-rRNA_Ath     0.87
	$annoh{$fields[0]} = $fields[2];
    }
    #CL109   673     integrated_virus        Caulimoviridae  Caulimovirus-4_STu     0.84
    elsif (scalar @fields == 6) {
	$annoh{$fields[0]} = $fields[3];
    }
    #G0	4669	transposable_element	ltr_retrotransposon	Gypsy	RLG-iketas	RLG-iketas-3_Contig173_BWAZ_227-1_66704_76056	0.10
    else {
	$annoh{$fields[0]} = $fields[4];
    }
}
close $an;

#dd \%annoh;
#exit;

while (<$cl>) {
    chomp;
    my ($clus, $readct) = split;
    $total_clusters++ if defined $clus;
    $total_clustered += $readct if defined $clus;
    $clsh{$clus} = $readct;
}
close $cl;

#dd \%clsh;
#exit;

say $out join "\t", "Cluster","ReadNum","GPerc","Superfamily","Color";
for my $cls (reverse sort { $clsh{$a} <=> $clsh{$b} } keys %clsh) {
    if (exists $annoh{$cls} && exists $colormap->{ $annoh{$cls} }) {
	my $cls_perc_clustered = $clsh{$cls} / $readtot * 100;
	say $out join "\t", $cls, $clsh{$cls}, $cls_perc_clustered, $annoh{$cls}, $colormap->{$annoh{$cls}};
    }
}
close $out;
my $perc_clustered = $total_clustered / $readtot * 100;
say "$total_clustered reads in $total_clusters clusters. $perc_clustered perc clustered.";
