#!/usr/bin/env perl

use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;
use Storable qw(freeze thaw);
use Capture::Tiny qw(:all);
use File::Basename;
use File::Path qw(make_path);
use File::Spec;
use File::Basename;
use List::Util qw(max);
use POSIX qw(strftime);
use feature 'say';

my $infile;
my $outfile;
my $percent_id;
my $percent_cov;
my $cluster_size;
my $fas_file;

my %match_pairs;
my %match_index;

#counters
my $total_hits = 0;
my $parsed_hits = 0;
my $index = 0;

GetOptions(
           'i|infile=s'         => \$infile,
	   'f|fas_file=s'       => \$fas_file,
           'o|outfile=s'        => \$outfile,
	   'cs|cluster_size=i'  => \$cluster_size,
	   'id|percent_id=f'    => \$percent_id,
	   'cov|percent_cov=f'  => \$percent_cov,
          );

# open the infile or die with a usage statement
if (!$infile || !$outfile || !$fas_file
    || !$percent_id || !$percent_cov) {
    print "\nERROR: Input not parsed correctly.\n";
    usage() and exit(1);
}

my ($match_pairs, $match_index) = parse_mgblast($infile, $percent_id, $percent_cov);

my $index_file = $outfile.".index";  # integer index for read IDs used in clustering
my $cluster_file = $outfile.".cls";  # cluster file in Repeat Explorer "cls" format
my $hitsort_int = $outfile.".int";   # integer index for clustering

open(my $hs_int, '>', $hitsort_int);
open(my $out, '>', $outfile);
open(my $indexout, '>', $index_file);

for my $id (keys %$match_index) {
    $match_index{$id} = $index; 
    say $indexout "$id $index";
    $index++; 
}
close($indexout);

for my $match (keys %$match_pairs) {
    my $match_score = max(@{$match_pairs{$match}});
    $parsed_hits++;
    next unless defined $match_score;
    my ($qry, $sbj) = split /\|/, $match;
    if (exists $match_index{$qry} && exists $match_index{$sbj}) {
	say $out join "\t", $qry, $sbj, $match_score;
	say $hs_int join "\t", $match_index{$qry}, $match_index{$sbj}, $match_score;
    }
}
close($hs_int);
close($out);

my $community = louvain_method($hitsort_int, $index_file);
my $cls_file = make_clusters($community, \%match_index, $cluster_file, $outfile, $index_file);

### generate plots from cluster file

my $seqhash = fas2hash($fas_file);
my $cls_dir_path = clusters2fasta($cls_file, $seqhash, $cluster_size);


### search each fasta file in $cls_dir_path against custom repeatdb and HMMER

#
# Subs
#
sub parse_mgblast {
    my ($infile, $percent_id, $percent_cov) = @_;

    open(my $in, '<', $infile);

    while (<$in>) { 
	chomp; 
	my ($q_name, $q_len, $q_start, $q_end, $s_name, $s_len, 
	    $s_start, $s_end, $pid, $score, $e_val, $strand) = split;
    
	my $pair = join "|", $q_name, $s_name;
	my $revpair = join "|", $q_name, $s_name;
	my $subj_hit_length = ($s_end - $s_start) + 1;
	my $subj_cov = $subj_hit_length/$s_len;
	
	if ($strand eq '-') {
	    $total_hits++;
	    my $neg_query_hit_length = ($q_start - $q_end) + 1;
	    my $neg_query_cov = $neg_query_hit_length/$q_len;
	    
	    if ( ($neg_query_cov >= $percent_cov) && ($subj_cov >= $percent_cov) && ($pid >= $percent_id) ) {
		if (exists $match_pairs{$pair}) {
		    push @{$match_pairs{$pair}}, $score;
		}
		else {
		    $match_pairs{$pair} = [$score];
		    $match_index{$q_name} = [] unless exists $match_index{$q_name}; # hash is much more memory efficient than array (1.6g v 2.3g)
		    $match_index{$s_name} = [] unless exists $match_index{$s_name}; # and faster (though only ~1s faster)
		}
	    }
	}
	else {
	    $total_hits++;
	    my $pos_query_hit_length = ($q_end - $q_start) + 1;
	    my $pos_query_cov = $pos_query_hit_length/$q_len;
	    
	    if ( ($pos_query_cov >= $percent_cov) && ($subj_cov >= $percent_cov) && ($pid >= $percent_id) ) {
		if (exists $match_pairs{$pair}) {
		push @{$match_pairs{$pair}}, $score;
		}
		else {
		    $match_pairs{$pair} = [$score];
		    $match_index{$q_name} = [] unless exists $match_index{$q_name};
		    $match_index{$s_name} = [] unless exists $match_index{$s_name};
		}
	    }
	}
    }
    close($in);
    
    return(\%match_pairs, \%match_index);
}

sub louvain_method {
    my ($hitsort_int, $index_file) = @_;
    my ($iname, $ipath, $isuffix) = fileparse($hitsort_int, qr/\.[^.]*/);
    my $cls_bin = $iname.".bin";                    # Community "bin" format
    my $cls_tree = $iname.".tree";                  # hierarchical tree of clustering results
    my $cls_tree_weights = $cls_tree.".weights";    # bit score, the weights applied to clustering
    my $cls_tree_log = $cls_tree.".log";            # the summary of clustering results at each level of refinement
    my $hierarchy_err = $cls_tree.".hierarchy.log"; # some other log
    system("louvain_convert -i $hitsort_int -o $cls_bin -w $cls_tree_weights");
    system("louvain_community $cls_bin -l -1 -w $cls_tree_weights -v >$cls_tree 2>$cls_tree_log");

    my $levels = `grep -c level $cls_tree_log`;
    chomp($levels);

    my @comm;
    for (my $i = 0; $i <= $levels-1; $i++) {
        my $cls_graph_comm = $cls_tree.".graph_node2comm_level_".$i;
	my $cls_graph_clusters = $cls_tree.".graph.clusters";
	my $cls_graph_membership = $cls_tree.".graph.membership";

        system("louvain_hierarchy $cls_tree -l $i > $cls_graph_comm");
  	push @comm, $cls_graph_comm;
    }
    return \@comm;
}

sub make_clusters {
    my ($graph_comm, $match_index, $cluster_file, $outfile, $index_file) = @_;

    my @graph_comm_sort = reverse sort { ($a =~ /(\d)$/) <=> ($b =~ /(\d)$/) } @$graph_comm;
    my $graph = shift @graph_comm_sort;
    my %rindex = reverse %$match_index; # this only works if you are *certain* there are no duplicate keys, which I am ;)
    my %clus;

    my $membership_file = $cluster_file.".membership.txt";
    open(my $mem,'>', $membership_file);
    open(my $in, '<', $graph);
    open(my $cls_out, '>', $cluster_file);

    while (my $line = <$in>) {
	chomp $line;
	my ($i, $j) = split /\s+/, $line;
	if (exists $clus{$j}) {
	    push @{$clus{$j}}, $i;
	}
	else {
	    $clus{$j} = [$i];
	}
    }
    close($in);

    my $cls_ct = 1;
    for my $cls (reverse sort { @{$clus{$a}} <=> @{$clus{$b}} } keys %clus) {
	my $clus_size = scalar @{$clus{$cls}};
	say $cls_out ">CL$cls_ct $clus_size";
	my @clus_members;
	for my $cls_member (@{$clus{$cls}}) {
	    say $mem "$cls_member $cls_ct";
	    if (exists $rindex{$cls_member}) {
		push @clus_members, $rindex{$cls_member};
	    }
	}
	say $cls_out join " ", @clus_members;
	$cls_ct++;
    }
    close($cls_out);
    close($mem);

    return $cluster_file;
    # output CLS for each level and clusterMembership.txt for each level
}

sub clusters2fasta {
    my ($cls_file, $seqhash, $cluster_size) = @_;

    local $/ = '>';

    my $cls_dir_path;

    $cluster_size = defined($cluster_size) ? $cluster_size : '500';

    open(my $cls, '<', $cls_file);
    local $/ = '>';

    my $str = POSIX::strftime("%m_%d_%Y_%H_%M_%S", localtime);
    my ($iname, $ipath, $isuffix) = fileparse($cls_file, qr/\.[^.]*/);
    my $cls_dir_base = $iname;
    $cls_dir_base =~ s/\.[^.]+$//;
    my $cls_dir = $cls_dir_base."_cls_fasta_files_$str";
    $cls_dir_path = $ipath.$cls_dir;
    make_path($cls_dir_path, {verbose => 0, mode => 0711,}); # allows for recursively making paths
    #mkdir($cls_dir_path, 0777) || die "ERROR: Could not make directory for clusters. Exiting.\n"; 

    while (my $line = <$cls>) {
        my ($clsid, $seqids) = split /\n/, $line;
        next unless defined $seqids;
        my @ids = split /\s+/, $seqids;
        my ($clid, $cls_seq_num) = $clsid =~ /(\S+)\s(\d+)/;
        if ($cls_seq_num >= $cluster_size) {
            my $indiv_clsfile = join "_", $clid, $cls_seq_num;
            $indiv_clsfile .= ".fas";  
            my $cls_file_path = File::Spec->catfile($cls_dir_path, $indiv_clsfile);
            open(my $clsout, '>', $cls_file_path);
            for my $seq (@ids) {
                if (exists $seqhash->{$seq}) {
                    say $clsout join "\n", ">".$seq, $seqhash->{$seq};
                }
            }
            close($clsout);
        }
    }

    return $cls_dir_path;
}

sub fas2hash {
    my $fas_file = shift;

    my %seqhash;
    open(my $fas, '<', $fas_file);

    local $/ = '>';

    while (my $line = <$fas>) {
        my ($seqid, $seq) = split /\n/, $line;;
        $seqhash{$seqid} = $seq;
    }
    close($fas);
    
    return \%seqhash;
}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -i inreport -o hitsort -id 90.00 -cov 0.55 -f fasta_file

Required:
    -i|infile        :    An mgblast report in tab format (-D 4).
    -f|fas_file      :    The (Fasta) file of sequences used in the all-vs-all blast.
    -o|outfile       :    File name to write the parsed results to in Repeat Explorer's hitsort format.
    -cs|cluster_size :    The minimum cluster size to convert to fasta (Default: 500).
    -id|percent_id   :    The percent identity threshold for matches.
    -cov|percent_cov :    The percent coverage for both the query and subject to be retained for clustering.

Options:
    -s|statsfile     :    A file to hold the cluster stats. [NOT IMPLEMENTED]

END
}
