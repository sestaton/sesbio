#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use List::Util qw(max sum);
use Data::Dump qw(dd);
use File::Spec;
use File::Basename;
use Cwd;
use DBM::Deep;
use BerkeleyDB;

my $usage = "$0 blastfile\n";
my $blast = shift or die $usage;
my $inmemory = 0;
my $exists = 0;

my $cwd = getcwd();
my ($iname, $ipath, $isuffix) = fileparse($blast, qr/\.[^.]*/);
my $int_file = $iname;
my $idx_file = $iname;
my $hs_file  = $iname;
$int_file    .= "_louvain.int";
$idx_file    .= "_louvain.idx";
$hs_file     .= "_louvain.hs";
my $int_path = File::Spec->catfile($cwd, $int_file);
my $idx_path = File::Spec->catfile($cwd, $idx_file);
my $hs_path  = File::Spec->catfile($cwd, $hs_file);

# counters
my $total_hits  = 0;
my $parsed_hits = 0;
my $index       = 0;

open my $fh, '<', $blast;

if ($inmemory) {
    my %match_pairs;
    my %match_index;
    
    while (<$fh>) {
	chomp;
	my ($q_name, $s_name, $pid, $aln_len, $mistmatch, $gaps, $q_start, $q_end,
	    $s_start, $s_end, $e_val, $score) = split;
        
	my $pair            = mk_key($q_name, $s_name);
	my $revpair         = mk_key($s_name, $q_name);
	my $subj_hit_length = ($s_end - $s_start) + 1;

	if ($q_start > $q_end) {
	    $total_hits++;
	    my $neg_query_hit_length = ($q_start - $q_end) + 1;

	    if ( ($neg_query_hit_length >= 55) && ($pid >= 90) ) {
		if (exists $match_pairs{$pair}) {
		    push @{$match_pairs{$pair}}, $score;
		}
		elsif (exists $match_pairs{$revpair}) {
		    push @{$match_pairs{$revpair}}, $score;
		}
		else {
		    $match_pairs{$pair}   = [$score];
		    $match_index{$q_name} = $index unless exists $match_index{$q_name};
		    $index++;
		    $match_index{$s_name} = $index unless exists $match_index{$s_name};
		    $index++;
		}
	    }
	}
	else {
	    $total_hits++;
	    my $pos_query_hit_length = ($q_end - $q_start) + 1;
	    
	    if ( ($pos_query_hit_length>= 55) && ($pid >= 90) ) {
		if (exists $match_pairs{$pair}) {
		    push @{$match_pairs{$pair}}, $score;
		}
		elsif (exists $match_pairs{$revpair}) {
		    push @{$match_pairs{$revpair}}, $score;
		}
		else {
		    $match_pairs{$pair}   = [$score];
		    $match_index{$q_name} = $index unless exists $match_index{$q_name};
		    $index++;
		    $match_index{$s_name} = $index unless exists $match_index{$s_name};
		    $index++;
		}
	    }
	}
    }
    close $fh;

    open my $idx, '>', $idx_path or die "\n[ERROR]: Could not open file: $idx_path\n";
    
    for my $idx_mem (sort { $match_index{$a} <=> $match_index{$b} } keys %match_index) {
	say $idx join q{ }, $idx_mem, $match_index{$idx_mem};
    }
    close $idx;
    
    open my $int, '>', $int_path or die "\n[ERROR]: Could not open file: $int_path\n";
    open my $hs,  '>', $hs_path  or die "\n[ERROR]: Could not open file: $hs_path\n";

    while (my ($match, $scores) = each %match_pairs) {
	say $match;
	my $match_score = max(@$scores);
	my ($qry, $sbj) = mk_vec($match);
	my $revmatch = mk_key($sbj, $qry);
	if (exists $match_pairs{$revmatch}) {
	    my $rev_match_score = max(@{$match_pairs{$revmatch}});
	    if ($rev_match_score > $match_score) {
		if (exists $match_index{$sbj} && exists $match_index{$qry}) {
		    say $hs join "\t", $sbj, $qry, $rev_match_score;
		    say $int join "\t", $match_index{$sbj}, $match_index{$qry}, $rev_match_score;
		    delete $match_pairs{$match};
		}
	    }
	    else {
		delete $match_pairs{$revmatch};
	    }
	}
	else {
	    if (exists $match_index{$qry} && exists $match_index{$sbj}) {
		say $hs join "\t", $qry, $sbj, $match_score;
		say $int join "\t", $match_index{$qry}, $match_index{$sbj}, $match_score;
	    }
	}
    }
    close $int;
    close $hs;

    #say "$exists keys in match_pairs";
}
else {
    my $dbm = "mgblast_matchpairs.dbm";
    my $dbi = "mgblast_matchindex.dbm";
    
    unlink $dbm if -e $dbm;
    unlink $dbi if -e $dbi;
    
    #my $db = DBM::Deep->new( file      => $dbm,
    #			     locking   => 1,
    #			     autoflush => 1,
    #			     type      => DBM::Deep::TYPE_HASH );

    my %match_pairs;
    my %match_index;

    tie %match_pairs, 'DBM::Deep', { file => $dbm, locking => 1, autoflush => 0, type => DBM::Deep::TYPE_HASH };
    tie %match_index, 'BerkeleyDB::Btree', -Filename => $dbi, -Flags => DB_CREATE
	or die "\n[ERROR]: Could not open DBM file: $dbi: $! $BerkeleyDB::Error\n";
    
    while (<$fh>) {
        chomp;
        my ($q_name, $s_name, $pid, $aln_len, $mistmatch, $gaps, $q_start, $q_end,
            $s_start, $s_end, $e_val, $score) = split;
        
        my $pair            = mk_key($q_name, $s_name);
        my $revpair         = mk_key($s_name, $q_name);
        my $subj_hit_length = ($s_end - $s_start) + 1;

        if ($q_start > $q_end) {
            $total_hits++;
            my $neg_query_hit_length = ($q_start - $q_end) + 1;

            if ( ($neg_query_hit_length >= 55) && ($pid >= 90) ) {
                if (exists $match_pairs{$pair}) {
                    push @{$match_pairs{$pair}}, $score;
                }
                elsif (exists $match_pairs{$revpair}) {
                    push @{$match_pairs{$revpair}}, $score;
                }
                else {
                    $match_pairs{$pair}   = [$score];
                    $match_index{$q_name} = $index unless exists $match_index{$q_name};
                    $index++;
                    $match_index{$s_name} = $index unless exists $match_index{$s_name};
                    $index++;
                }
            }
	}
        else {
            $total_hits++;
            my $pos_query_hit_length = ($q_end - $q_start) + 1;
            
            if ( ($pos_query_hit_length>= 55) && ($pid >= 90) ) {
                if (exists $match_pairs{$pair}) {
                    push @{$match_pairs{$pair}}, $score;
                }
                elsif (exists $match_pairs{$revpair}) {
                    push @{$match_pairs{$revpair}}, $score;
                }
                else {
                    $match_pairs{$pair}   = [$score];
                    $match_index{$q_name} = $index unless exists $match_index{$q_name};
                    $index++;
                    $match_index{$s_name} = $index unless exists $match_index{$s_name};
                    $index++;
                }
            }
        }
    }
    close $fh;

    open my $idx, '>', $idx_path or die "\n[ERROR]: Could not open file: $idx_path\n";

    for my $idx_mem (sort { $match_index{$a} <=> $match_index{$b} } keys %match_index) {
        say $idx join q{ }, $idx_mem, $match_index{$idx_mem};
    }
    close $idx;

    #dd \%match_index;

    open my $int, '>', $int_path or die "\n[ERROR]: Could not open file: $int_path\n";
    open my $hs,  '>', $hs_path  or die "\n[ERROR]: Could not open file: $hs_path\n";

    #dd \%match_pairs;
    #exit;
    #while (my ($match, $scores) = each %match_pairs) {
	#say $match;
    #}
    #exit;


    #my $exists = 0;
    while (my ($match, $scores) = each %match_pairs) {
    #for my $match (keys %match_pairs) {
	say $match; # 101
        my $match_score = max(@$scores);
        my ($qry, $sbj) = mk_vec($match);
	#say $qry; 101
        my $revmatch = mk_key($sbj, $qry);
	#say $revmatch; 101
        if (exists $match_pairs{$revmatch}) {
	    #say $revmatch;
            my $rev_match_score = max(@{$match_pairs{$revmatch}});
            if ($rev_match_score > $match_score) {
                if (exists $match_index{$sbj} && exists $match_index{$qry}) {
		    #$exists++;
		    #say $sbj;
                    say $hs join "\t", $sbj, $qry, $rev_match_score;
                    say $int join "\t", $match_index{$sbj}, $match_index{$qry}, $rev_match_score;
                    #delete $match_pairs{$match};
                }
            }
	    #else {
                #delete $match_pairs{$revmatch};
            #}
        }
        else {
	    #say $qry;
            if (exists $match_index{$qry} && exists $match_index{$sbj}) {
		#$exists++;
		#say $sbj;
                say $hs join "\t", $qry, $sbj, $match_score;
                say $int join "\t", $match_index{$qry}, $match_index{$sbj}, $match_score;
            }
        }
    }
    #say "$exists keys in match_pairs";
    close $int;
    close $hs;

    untie %match_index;
    untie %match_pairs;
    unlink $dbi;
    unlink $dbm;
}
    
sub mk_key { return join "||", @_ }

sub mk_vec { return split /\|\|/, shift }
