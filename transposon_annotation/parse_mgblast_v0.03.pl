#!/usr/bin/env perl

use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;
use Storable qw(freeze thaw);
use feature 'say';
use lib qw(/home/jmblab/statonse/apps/perlmod/Capture-Tiny-0.19/blib/lib);
use Capture::Tiny qw(:all);
use File::Basename;

my $infile;
my $outfile;
my $percent_id;
my $percent_coverage;
my $usage = "USAGE: parse_mgblast -i inreport -o parsedreport -id 90.00 -cov 0.50\n";

my %match_pairs;
my %match_index;
my @matchIDs;

#counters
my $total_hits = 0;
my $parsed_hits = 0;
my $index = 0;

GetOptions(
           'i|infile=s'         => \$infile,
           'o|outfile=s'        => \$outfile,
	   'id|percent_id=f'    => \$percent_id,
	   'cov|percent_cov=f'  => \$percent_coverage,
          );

# open the infile or die with a usage statement
print $usage and exit(1) if !$infile or !$outfile;

open(my $in, '<', $infile);

while (<$in>) { 
  chomp; 
  my ($q_name, $q_len, $q_start, $q_end, $s_name, $s_len, 
      $s_start, $s_end, $pid, $score, $e_val, $strand) = split;

  my $pair = join "|", $q_name, $s_name;
  my $subj_hit_length = ($s_end - $s_start) + 1;
  my $subj_cov = $subj_hit_length/$s_len;

  if ($strand eq '-') {
      $total_hits++;
      my $neg_query_hit_length = ($q_start - $q_end) + 1;
      my $neg_query_cov = $neg_query_hit_length/$q_len;

      if ( ($neg_query_cov >= $percent_coverage) && ($subj_cov >= $percent_coverage) && ($pid >= $percent_id) ) {
          if (exists $match_pairs{$pair}) {
              push @{$match_pairs{$pair}}, $score;
          }
          else {
              $match_pairs{$pair} = [];
	      push @matchIDs, $q_name, $s_name;
          }
      }
  }
  else {
      $total_hits++;
      my $pos_query_hit_length = ($q_end - $q_start) + 1;
      my $pos_query_cov = $pos_query_hit_length/$q_len;

      if ( ($pos_query_cov >= $percent_coverage) && ($subj_cov >= $percent_coverage) && ($pid >= $percent_id) ) {
          if (exists $match_pairs{$pair}) {
              push @{$match_pairs{$pair}}, $score;
          }
          else {
              $match_pairs{$pair} = [];
	      push @matchIDs, $q_name, $s_name;
          }
      }
  }
}
close($in);

open(my $out, '>', $outfile);

my %seen = ();
my @unique_IDs = grep { ! $seen{$_} ++ } @matchIDs;
for my $id (@unique_IDs) { $match_index{$id} = $index; $index++; }

for my $match (reverse sort { @{$match_pairs{$a}} <=> @{$match_pairs{$b}} } keys %match_pairs) {
    $parsed_hits++;
    my $match_score = shift @{$match_pairs{$match}};
    next unless defined $match_score;
    my ($qry, $sbj) = split /\|/, $match;
    if (exists $match_index{$qry} && exists $match_index{$sbj}) {
	say $out join "\t", $match_index{$qry}, $match_index{$sbj}, $match_score;
    }
}
close($out);

louvain_method($outfile);
#my $sindex = freeze \%match_index;

#print "\n===== There are ",$total_hits," records in this mgblast report.\n";
#print "\n===== There are ",$parsed_hits, " parsed records in the output.\n\n";

#
# Subs
#
sub louvain_method {
    my $cls_int = shift;
    my ($iname, $ipath, $isuffix) = fileparse($cls_int, qr/\.[^.]*/);
    my $cls_bin = $iname.".bin";
    my $cls_tree = $iname.".tree";
    my $cls_tree_weights = $cls_tree.".weights";
    my $cls_tree_log = $cls_tree.".log";
    my $hierarchy_err = $cls_tree.".hierarchy.log";
    my (@comm_res, @hier_res);
    #my ($convert_out, $convert_err) = capture { system("louvain_convert -i $cls_int -o $cls_bin -w $cls_tree_weights"); };
    system("louvain_convert -i $cls_int -o $cls_bin -w $cls_tree_weights");
    #say "Convert out: $convert_out"; say "Convert err: $convert_err";
    #($cls_tree, $cls_tree_log) = capture { 
    system("louvain_community $cls_bin -l -1 -w $cls_tree_weights -v >$cls_tree 2>$cls_tree_log");
	#system("louvain_commmunity $cls_bin -l -1 -w $cls_tree_weights -v"); };
   
    # store the output and search it for "level", then create trees for each level
    #say "Community out: $community_out"; say "Community err: $community_err";

    my $levels = `grep -c level $cls_tree_log`;
    chomp($levels);
    
    #say "There are $levels levels in $cls_tree_log";
    for (my $i = 0; $i <= $levels-1; $i++) {
    #for my $i (0..$#levels) {
	my $cls_graph_comm = $cls_tree.".graph_node2comm_level_".$i;
	#($cls_graph_comm, $hierarchy_err) = capture { system("louvain_hierarchy $cls_tree -l $i"); };
	system("louvain_hierarchy $cls_tree -l $i > $cls_graph_comm");
	# take a look at the output and then convert to cluster IDs
	#makeCls.py graph_node2comm_leve_$i index cluster cluster_membership // This is where we map the ids back to the cluster indices
	#make_clusters($cls_graph_comm, \%match_index);
    }
}

sub make_clusters {
    my ($graph_comm, $match_index) = @_;
    
    # output CLS for each level and clusterMembership.txt for each level

}

