#!/usr/bin/env perl

use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;
use Storable qw(freeze thaw);
use feature 'say';

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
           "i|infile=s"         => \$infile,
           "o|outfile=s"        => \$outfile,
	   "id|percent_id=f"      => \$percent_id,
	   "cov|percent_cov=f"    => \$percent_coverage,
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

my %seen = ();
my @unique_IDs = grep { ! $seen{$_} ++ } @matchIDs;
for my $id (@unique_IDs) { $match_index{$id} = $index; $index++; }

open(my $out, '>', $outfile);

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

my $sindex = freeze \%match_index;

#print "\n===== There are ",$total_hits," records in this mgblast report.\n";
#print "\n===== There are ",$parsed_hits, " parsed records in the output.\n\n";

