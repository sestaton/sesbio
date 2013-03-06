#!/usr/bin/env perl

use strict;
use warnings;
use autodie qw(open);
use feature 'say';
use File::Spec qw(catfile rel2abs);
use File::Basename qw(fileparse);
use File::Path qw(make_path);
use Getopt::Long;
use Graph::UnionFind;
use Data::Dump qw(dd dump);

my $usage = "$0 -i cluster_fas_files_dir";
my $indir;
my $database;
my $outdir;
my $json;

GetOptions(
           'i|indir=s'  => \$indir,
           'o|outdir=s' => \$outdir,
           'd|db=s'     => \$database,
           'j|json=s'   => \$json,
           );

die $usage if !$indir;

## get input files                                                                                                                                           
opendir(DIR, $indir) || die "\nERROR: Could not open directory: $indir. Exiting.\n";
my @clus_fas_files = grep /\.fa.*$/, readdir(DIR);
closedir(DIR);

if (scalar(@clus_fas_files) < 1) {
    say "\nERROR: Could not find any fasta files in $indir. Exiting.\n";
    exit(1);
}
my $uf = Graph::UnionFind->new;
my %vertex;

chdir($indir);
my $read_pairs = find_pairs(\@clus_fas_files);

#dd $read_pairs;

my %mapped_pairs;
while (my ($cls, $reads) = each %$read_pairs) {
    for my $read (@$reads) {
	$read =~ s/\/\d$//;
	if (exists $mapped_pairs{$read}) {
	    push @{$mapped_pairs{$read}}, $cls;
	}
	else {
	    $mapped_pairs{$read} = [$cls];
	}
    }
}

my %seen;
my %cluster_merge_candidates;
my %cluster_connections;
for my $key (keys %mapped_pairs) {
    my $pairct = scalar @{$mapped_pairs{$key}};
    $seen{$_}++ for @{$mapped_pairs{$key}};
    my $r = (keys %seen)[0];            # just grab the first key
    if ($pairct > 1 && $seen{$r} < 2) { # this ensures the key is unique (i.e., we don't care about paired reads in the same cluster)
	my $pair = join "|", @{$mapped_pairs{$key}};
	$cluster_merge_candidates{$pair}++;
	#say "$key ==> $pairct ==> ",join ",", @{$mapped_pairs{$key}};
    }
    undef %seen;
}

my %linkage;
for my $clspair (reverse sort { $cluster_merge_candidates{$a} <=> $cluster_merge_candidates{$b} } keys %cluster_merge_candidates) {
    my ($i, $j) = split /\|/, $clspair;
    if ($cluster_merge_candidates{$clspair} > 1000) { # how many times were reads split between these clusters
	say join "\t", $i, $j, $cluster_merge_candidates{$clspair};
	if (exists $linkage{$i}) {
	    if (exists $linkage{$j}) {
		#warn "Merging $i and $j\n";
		$linkage{$_} = $linkage{$i} for grep $linkage{$_} eq $linkage{$j}, keys %linkage;
	    }
	    else {
		#warn "Adding $j to $i\n";
		$linkage{$j} = $linkage{$i};
	    }
	}
	elsif (exists $linkage{$j}) {
	    #warn "Adding $i to $j\n";
	    $linkage{$i} = $linkage{$j};
	}
	else {
	    #warn "New $i and $j to $i\n";
	    @linkage{$i, $j} = ($i) x 2;
	}	
    }
}
	

my %groups;
for my $key (keys %linkage) {
    push @{ $groups{ $linkage{$key} } }, $key;
}
 
say "======================================================================";
 
for my $key (sort keys %groups) {
    say join ',' => @{ $groups{$key} };
}

#
# subs
#
sub find_pairs {
    my $clus_fas_files = shift;

    my %read_pairs;
    local $/ = '>';
   
    for my $file (@$clus_fas_files) {
	my ($cluster, $readnum, @suffix) = split /\_/, $file;
	open(my $in, '<', $file);	
	while (my $line = <$in>) {
	    $line =~ s/>//g;
	    next if !length($line);
	    my ($seqid, @seqs) = split /\n/, $line;
	    my $seq = join '', @seqs;
	    push @{$read_pairs{$cluster}}, $seqid;  ## Not using the seq at this point
	}
	close($in);
    }

    return \%read_pairs;
}
