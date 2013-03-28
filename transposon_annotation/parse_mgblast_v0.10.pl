#!/usr/bin/env perl

use v5.12;
use utf8;
use strict;
use warnings;
use warnings FATAL => "utf8";
use autodie qw(open);
use Getopt::Long;
use File::Basename;
use List::Util qw(max);
use charnames qw(:full :short);
use Encode qw(encode decode);
use lib qw(/home/jmblab/statonse/apps/perlmod/DBM-Deep-2.0008/blib/lib);
use DBM::Deep;
BEGIN {
  @AnyDBM_File::ISA = qw( DB_File SQLite_File )
      unless @AnyDBM_File::ISA == 1; 
}
use AnyDBM_File;
use vars qw( $DB_BTREE &R_DUP );
use AnyDBM_File::Importer qw(:bdb);

# lexical vars
my $infile;
my $outfile;
my $percent_id;
my $percent_coverage;
my $usage = "USAGE: parse_mgblast -i inreport -o hitsort-int -id 90.00 -cov 0.50\n";

my $dbm = "matchpairs.dbm";
my $dbi = "matchindex.dbm";

if (-e $dbm) { unlink($dbm); }
if (-e $dbi) { unlink($dbi); }

my $db = DBM::Deep->new( file      => $dbm,
                         locking   => 1,
                         autoflush => 1,
                         type      => DBM::Deep::TYPE_HASH
                         );

my %match_index;
$DB_BTREE->{cachesize} = 100000;
$DB_BTREE->{flags} = R_DUP;

tie %match_index, 'AnyDBM_File', $dbi, O_RDWR|O_CREAT, 0666, $DB_BTREE
    or die "\nERROR: Could not open DBM file $dbi: $!\n";

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
    
    my $pair = mk_key($q_name, $s_name);
    my $revpair = mk_key($q_name, $s_name);
    my $enc_pair = encode("UTF-8", $pair, 1);
    my $enc_revpair = encode("UTF-8", $revpair, 1);
    my $subj_hit_length = ($s_end - $s_start) + 1;
    my $subj_cov = $subj_hit_length/$s_len;
    
    if ($q_start > $q_end) {
	$total_hits++;
	my $neg_query_hit_length = ($q_start - $q_end) + 1;
	my $neg_query_cov = $neg_query_hit_length/$q_len;
	
	if ( ($neg_query_cov >= $percent_coverage) && ($pid >= $percent_id) ) {
            if (exists $db->{$enc_pair}) {
                push @{$db->{$enc_pair}}, $score;
            }
            elsif (exists $db->{$enc_revpair}) {
                push @{$db->{$enc_revpair}}, $score;
            }
            else {
                $db->{$enc_pair} = [$score];
                $match_index{$q_name} = $index unless exists $match_index{$q_name}; # hash is much more memory efficient than array (1.6g v 2.3g)
                $index++;
                $match_index{$s_name} = $index unless exists $match_index{$s_name}; # and faster (though only ~1s faster)
                $index++;
            }
        }
    }
    else {
        $total_hits++;
        my $pos_query_hit_length = ($q_end - $q_start) + 1;
        my $pos_query_cov = $pos_query_hit_length/$q_len;
        
        if ( ($pos_query_cov >= $percent_coverage) && ($pid >= $percent_id) ) {
            if (exists $db->{$enc_pair}) {
                push @{$db->{$enc_pair}}, $score;
            }
            elsif (exists $db->{$enc_revpair}) {
                push @{$db->{$enc_revpair}}, $score;
            }
            else {
                $db->{$enc_pair} = [$score];
                $match_index{$q_name} = $index unless exists $match_index{$q_name};
                $index++;
                $match_index{$s_name} = $index unless exists $match_index{$s_name};
                $index++;
            }
        }
    }
}
close($in);

my $index_file = $outfile.".index";
open(my $out, '>', $outfile);
open(my $indexout, '>', $index_file);

for my $idx_mem (sort { $match_index{$a} <=> $match_index{$b} } keys %match_index) {
    say $indexout join " ", $idx_mem, $match_index{$idx_mem};
}
close($indexout);

#for my $match (keys %$db) {    # this uses a tremendous amount of memory; it stores all keys in an array
while (my ($match, $scores) = each %$db) {
    my $enc_match = encode("UTF-8", $match, 1);
    my $match_score = max(@$scores);
    #$parsed_hits++;
    my ($qry, $sbj) = mk_vec($enc_match);
    #my ($qry, $sbj) = split "|", $match;
    #my $revmatch = join "|", $sbj, $qry;
    my $revmatch = mk_key($sbj, $qry);
    my $enc_revmatch = encode("UTF-8", $revmatch, 1);
    if (exists $db->{$enc_revmatch}) { 
        my $rev_match_score = max(@{$db->{$enc_revmatch}});
        if ($rev_match_score > $match_score) {
            if (exists $match_index{$sbj} && exists $match_index{$qry}) {
                say $out join "\t", $match_index{$sbj}, $match_index{$qry}, $rev_match_score;
                delete $db->{$enc_match};
            }
        }
        else {
            delete $db->{$enc_revmatch};
        }   
    }
    else {
        if (exists $match_index{$qry} && exists $match_index{$sbj}) {
            say $out join "\t", $match_index{$qry}, $match_index{$sbj}, $match_score;
        }
    }
}
close($out);

untie %match_index;

# subs
sub mk_key { join "\N{INVISIBLE SEPARATOR}", @_ }

sub mk_vec { split "\N{INVISIBLE SEPARATOR}", shift }

