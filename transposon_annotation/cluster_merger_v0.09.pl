#!/usr/bin/env perl

## TODO: Need to write all clusters, those below threshold included, to new .cls file

use utf8;
use v5.12;
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
use POSIX qw(strftime);
use charnames qw(:full :short);

my $usage = "$0 -i cluster_file.cls -r cluster_grouping_report -f fasta";
my $infile;
my $fasta;
my $report;
my $cluster_size;

GetOptions(
           'i|infile=s'      => \$infile,
	   'f|fasta=s'       => \$fasta,
	   'r|report=s'      => \$report,
	   's|clustersize=i' => \$cluster_size,
           );

die $usage if !$infile or !$fasta or !$report; 

# set path for output
my $str = POSIX::strftime("%m_%d_%Y_%H_%M_%S", localtime); 
my ($iname, $ipath, $isuffix) = fileparse($infile, qr/\.[^.]*/);
my $cls_dir_base = $iname;
$cls_dir_base =~ s/\.[^.]+$//;
my $cls_with_merges = $cls_dir_base;
my $cls_dir = $cls_dir_base."_cls_fasta_files_$str";
$cls_with_merges .= "_merged_$str.cls";
my $cls_dir_path = $ipath.$cls_dir;
make_path($cls_dir_path, {verbose => 0, mode => 0711,}); # allows for recursively making paths
my $cls_with_merges_path = File::Spec->catfile($ipath, $cls_with_merges);
open(my $rep, '>', $report);
open(my $clsnew, '>', $cls_with_merges_path);

# find union in clusters
$cluster_size //= '500';
my ($seqs, $seqct) = fas2hash($fasta);
my ($read_pairs, $vertex, $uf) = find_pairs($infile, $rep);

my %cluster;
for my $v (keys %$vertex) {
    my $b = $$uf->find($v);
    die "$0: no block for $v" unless defined $b;
    push @{$cluster{$b}}, $v;
}

# generate groups based on cluster union
say $rep "=====> Cluster groupings (group_index\tclusters)";
my $group_index = 0;
for my $group (values %cluster) {
    my $groupseqnum; my @grpcp;
    for (@$group) { my $clsstrcp = $_; my ($id, $seqnum) = split /\_/, $clsstrcp, 2; $groupseqnum += $seqnum; push @grpcp, $id; }
    say $rep join "\t", $group_index, join ",", @grpcp;
    say $clsnew ">G$group_index $groupseqnum";
    my $group_file = "Cluster_grouping_".$group_index.".fas";
    my $group_file_path = File::Spec->catfile($cls_dir_path, $group_file);
    open(my $groupout, '>', $group_file_path);
    
    for my $clus (@$group) {
	if (exists $read_pairs->{$clus}) {
	    print $clsnew join " ",@{$read_pairs->{$clus}};
	    for my $read (@{$read_pairs->{$clus}}) {
		if (exists $seqs->{$read}) {
		    say $groupout join "\n", ">".$read, $seqs->{$read};
		}
		else {
		    say "WARNING: $read not found. This is a bug. Please report it.";
		}
            }
        }
	print $clsnew " ";
	delete $read_pairs->{$clus}
    }
    print $clsnew "\n";
    close($groupout);
    $group_index++;
}

# write out those clusters that weren't merged
say $rep "=====> Non-grouped clusters";
for my $non_paired_cls (keys %$read_pairs) {
    my ($non_paired_clsid, $non_paired_clsseqnum) = split /\_/, $non_paired_cls, 2;
    say $rep $non_paired_clsid;
    say $clsnew join "\n", ">$non_paired_clsid $non_paired_clsseqnum", join " ", @{$read_pairs->{$non_paired_cls}};

    if (scalar(@{$read_pairs->{$non_paired_cls}}) >= $cluster_size) {
	my $non_paired_clsfile .= $non_paired_cls.".fas";
	my $cls_file_path = File::Spec->catfile($cls_dir_path, $non_paired_clsfile);
	open(my $clsout, '>', $cls_file_path);

	for my $non_paired_read (@{$read_pairs->{$non_paired_cls}}) {
	    if (exists $seqs->{$non_paired_read}) {
		say $clsout join "\n", ">".$non_paired_read, $seqs->{$non_paired_read};
	    }
	    else {
		say "WARNING: $non_paired_read not found. This is a bug. Please report it.";
	    }
	}
	close($clsout);
    }
}
close($rep);
close($clsnew);

#
# subs
#
sub find_pairs {
    my ($cls_file, $rep) = @_;
    
    my $uf = Graph::UnionFind->new;

    say $rep "=====> Cluster connections above threshold";

    my %vertex;
    my %read_pairs;
    my %mapped_pairs;

    {
	local $/ = '>';
	
	open(my $in, '<', $cls_file);	
	while (my $line = <$in>) {
	    $line =~ s/>//g;
	    next if !length($line);
	    my ($clsid, $seqids) = split /\n/, $line;
	    $clsid =~ s/\s/\_/;
	    my @ids = split /\s+/, $seqids;
	    push @{$read_pairs{$clsid}}, $_ for @ids;
	}
	close($in);
    }

    while (my ($cls, $reads) = each %read_pairs) {
	for my $read (@$reads) {
	    my $readbase = $read;
	    $readbase =~ s/\/\d$//;
	    if (exists $mapped_pairs{$readbase}) {
		push @{$mapped_pairs{$readbase}}, {$read => $cls};
	    }
	    else {
		$mapped_pairs{$readbase} = [{$read => $cls}];
	    }
	}
    }
    
    my %cls_conn_ct;
    my ($cls_i, $cls_j);
    my @sep_reads;

    for my $allpairs (keys %mapped_pairs) {
	if (scalar(@{$mapped_pairs{$allpairs}}) < 2) {     # if no pair is found in another cluster, 
	    delete $mapped_pairs{$allpairs};               # remove this pair
	}
	else {
	    push @sep_reads, values %$_ for @{$mapped_pairs{$allpairs}};
	    ($cls_i, $cls_j) = sort @sep_reads;
	    if ($cls_i =~ /$cls_j/) {                      # remove reads that have pairs in the same cluster       
		delete $mapped_pairs{$allpairs};           # which is uninformative for merging clusters
	    }
	    else {
		#my $cls_merge_cand = join "|", $cls_i, $cls_j;
		#$cls_conn_ct{$cls_merge_cand}++;
		$cls_conn_ct{mk_key($cls_i,$cls_j)}++;
	    }
	}
	@sep_reads = ();
    }

    for my $p (reverse sort { $cls_conn_ct{$a} <=> $cls_conn_ct{$b} } keys %cls_conn_ct) {
	#my ($i, $j) = split /\|/, $p;
	my ($i, $j) = mk_vec($p);
	my $i_noct = $i; $i_noct =~ s/\_.*//;
	my $j_noct = $j; $j_noct =~ s/\_.*//;
        if ($cls_conn_ct{$p} >= 100) {    # threshold for merging clusters
	    say $rep join "\t", $i_noct, $j_noct, $cls_conn_ct{$p};
	    ++$vertex{$_} for $i, $j;
	    $uf->union($i, $j);
	}
    }
    return(\%read_pairs, \%vertex, \$uf);
}

sub fas2hash {
    my $fas_file = shift;
    open(my $fas, '<', $fas_file);

    my %seqhash;

    local $/ = '>';

    while (my $line = <$fas>) {
        $line =~ s/>//g;
        next if !length($line);
        my ($seqid, @seqs) = split /\n/, $line;
        my $seq = join '', @seqs;
        $seqhash{$seqid} = $seq;
        $seqct++ if defined $seq;
    }
    close($fas);
    
    return(\%seqhash, $seqct);
}

# the best way to make combined keys is not to concatenate them
# http://stackoverflow.com/a/15299397/1543853
#
# NB: in v5.16 charnames() will be loaded automatically
sub mk_key { join("\N{INVISIBLE SEPARATOR}", @_) }

sub mk_vec { map [split "\N{INVISIBLE SEPARATOR}"], @_ }
