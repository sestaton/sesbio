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
use List::Util qw(max);

my $infile;
my $outfile;
my $percent_id;
my $percent_coverage;
my $usage = "USAGE: parse_mgblast -i inreport -o hitsort -id 90.00 -cov 0.50\n";

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
  my $revpair = join "|", $q_name, $s_name;
  my $subj_hit_length = ($s_end - $s_start) + 1;
  my $subj_cov = $subj_hit_length/$s_len;

  if ($strand eq '-') {
      $total_hits++;
      my $neg_query_hit_length = ($q_start - $q_end) + 1;
      my $neg_query_cov = $neg_query_hit_length/$q_len;

      if ( ($neg_query_cov >= $percent_coverage) && ($subj_cov >= $percent_coverage) && ($pid >= $percent_id) ) {
         if (exists $match_pairs{$pair} || exists $match_pairs{$revpair}) {
              push @{$match_pairs{$pair}}, $score;
          }
          else {
              $match_pairs{$pair} = [$score];
	      push @matchIDs, $q_name, $s_name;
          }
      }
  }
  else {
      $total_hits++;
      my $pos_query_hit_length = ($q_end - $q_start) + 1;
      my $pos_query_cov = $pos_query_hit_length/$q_len;

      if ( ($pos_query_cov >= $percent_coverage) && ($subj_cov >= $percent_coverage) && ($pid >= $percent_id) ) {
          if (exists $match_pairs{$pair} || exists $match_pairs{$revpair}) {
              push @{$match_pairs{$pair}}, $score;
          }
          else {
              $match_pairs{$pair} = [$score];
	      push @matchIDs, $q_name, $s_name;
          }
      }
  }
}
close($in);

my $index_file = $outfile.".index";
open(my $out, '>', $outfile);
open(my $indexout, '>', $index_file);

my %seen = ();
my @unique_IDs = grep { ! $seen{$_} ++ } @matchIDs;
for my $id (@unique_IDs) { 
    $match_index{$id} = $index; 
    say $indexout "$id $index";
    $index++; 
}
close($indexout);

#for my $match (reverse sort { @{$match_pairs{$a}} <=> @{$match_pairs{$b}} } keys %match_pairs) {
for my $match (keys %match_pairs) {
    my $match_score = max(reverse sort {$a <=> $b} @{$match_pairs{$match}});
    #say join "\t", $match, $match_score, join ",", @{$match_pairs{$match}};
    $parsed_hits++;
    #next unless defined $match_score;
    my ($qry, $sbj) = split /\|/, $match;
    #say join "\t", $qry, $sbj, $match_score, @{$match_pairs{$match}};
    if (exists $match_index{$qry} && exists $match_index{$sbj}) {
	say $out join "\t", $match_index{$qry}, $match_index{$sbj}, $match_score;
    }
}
close($out);



