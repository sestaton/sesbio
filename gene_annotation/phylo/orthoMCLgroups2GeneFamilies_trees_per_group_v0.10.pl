#!/usr/bin/env perl

#########

use v5.14;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use lib qw(/iob_home/jmblab/statonse/apps/perlmod/Data-Dump-1.21/blib/lib);
use Data::Dump qw(dd);
use lib qw(/iob_home/jmblab/statonse/apps/perlmod/cjfields-Bio-Kseq-dc7a71c/lib/site_perl/5.14.1/x86_64-linux-thread-multi/auto); #/Bio/Kseq/Kseq.bs
use lib qw(/iob_home/jmblab/statonse/apps/perlmod/cjfields-Bio-Kseq-dc7a71c/lib/site_perl/5.14.1/x86_64-linux-thread-multi);      #/Bio/Kseq.pm
use Bio::Kseq;
use lib qw(/iob_home/jmblab/statonse/apps/perlmod/Capture-Tiny-0.19/blib/lib);
use Capture::Tiny qw(:all);
BEGIN {
  @AnyDBM_File::ISA = qw( DB_File SQLite_File )
      unless @AnyDBM_File::ISA == 1; # if loaded already, AnyDBM_File::ISA has a length of one;
}
use AnyDBM_File;
use vars qw( $DB_BTREE &R_DUP );
use AnyDBM_File::Importer qw(:bdb);
use Cwd;

#
# lexical vars
#
my $groups;
my $groups_blast;
my $orthogroup_stats;
my $tree_stats;
my $pep_fas;
my $nt_fas;
my $og_pep_fas;
my $og_nt_fas;
my $help;
my $man;

GetOptions(
           'g|groups=s'            => \$groups,
           'gb|groups_blast=s'     => \$groups_blast,
           'pf|pep_fas=s'          => \$pep_fas,
           'nf|nt_fas=s'           => \$nt_fas,
           'onf|og_nt_fas=s'       => \$og_nt_fas,
           'opf|og_pep_fas=s'      => \$og_pep_fas,
           'os|orthogroup_stats=s' => \$orthogroup_stats,
	   'ts|tree_stats=s'       => \$tree_stats,
           'h|help'                => \$help,
           'm|man'                 => \$man,
	   );

#
# check input
#
if (!$groups || !$nt_fas || 
    !$groups_blast || !$pep_fas || 
    !$og_nt_fas || !$og_pep_fas ||
    !$orthogroup_stats || !$tree_stats) {
    usage();
    exit(1);
}

my %seqhash;
my %oghash;
my %statshash;
my @treestats;
my @clusterstats;

$DB_BTREE->{cachesize} = 100000;
$DB_BTREE->{flags} = R_DUP;
my $db_file = "orthoMCL_groups.bdb";
my $og_db_file = "orthoMCL_outgroup.bdb";
tie( %seqhash, 'AnyDBM_File', $db_file, 0666, $DB_BTREE);
tie( %oghash, 'AnyDBM_File', $og_db_file, 0666, $DB_BTREE);

# set PATH for programs we need
my $muscle = find_prog("muscle");
my $pal2nal = find_prog("pal2nal");
my $raxml = find_prog("raxml");

# create Kseq objects for reading
my $knseq = Bio::Kseq->new($nt_fas);
my $nt_it = $knseq->iterator;

my $kpseq = Bio::Kseq->new($pep_fas);
my $pep_it = $kpseq->iterator;

my $kognseq = Bio::Kseq->new($og_nt_fas);
my $ogn_it = $kognseq->iterator;

my $kogpseq = Bio::Kseq->new($og_pep_fas);
my $ogp_it = $kogpseq->iterator;

while (my $nseq = $nt_it->next_seq) {
    $seqhash{ $nseq->{name} } = [ $nseq->{seq} ];
}

while (my $pseq = $pep_it->next_seq) {
    if (exists $seqhash{ $pseq->{name} }) {
        push(@{$seqhash{ $pseq->{name} }}, $pseq->{seq});
    }
}

while (my $ognseq = $ogn_it->next_seq) {
    $oghash{ $ognseq->{name} } = [ $ognseq->{seq} ];
}

while (my $ogpseq = $ogp_it->next_seq) {
    if (exists $oghash{ $ogpseq->{name} }) {
	push(@{$oghash{ $ogpseq->{name} }}, $ogpseq->{seq});
    }
}

# parse orthoMCL groups file and blast report
my $clusters = parse_groups($groups, %seqhash);
my $clusters_best_match = parse_best_hits($groups_blast);

## look at data structures
#dd %$clusters;
#dd $clusters_best_match;
#dd %oghash;

my %clusters_with_recip_hit;
my %taxahash;
my %tree_members;
my $og_taxa = 'grap';

for my $cluster (sort keys %$clusters) {
    my $cluster_size = scalar @{$clusters->{$cluster}};
    push(@clusterstats, $cluster_size);
    $statshash{$cluster} = $cluster_size; 
    my $gene_file = "gene_cluster_".$cluster."_nt.fasta";
    my $pep_file = "gene_cluster_".$cluster."_pep.fasta";
    open(my $nt_group, ">>", $gene_file) or die "\nERROR: Could not open file: $gene_file\n";
    open(my $pep_group, ">>", $pep_file) or die "\nERROR: Could not open file: $pep_file\n";
    
    my $gene_ct = 0;
    
    for my $gene (@{$clusters->{$cluster}}) {
	my ($taxa, $gene_id) = split(/\|/,$gene);
	$taxahash{$taxa}++;   
	my $og_id = $clusters_best_match->{$cluster}{$gene};
	if (defined $og_id && exists $oghash{$og_id}) {
	    $clusters_with_recip_hit{$cluster} = $og_id;
	    $taxahash{$og_taxa} = 1; # if present, count of outgroup taxa in alignment is 1.
	}
	if (exists $seqhash{$gene}) {
	    $gene_ct++;
	    $tree_members{$gene} = 1;
	    print $nt_group join("\n",(">".$gene,${$seqhash{$gene}}[0])),"\n";
	    print $pep_group join("\n",(">".$gene,${$seqhash{$gene}}[1])),"\n";
	}
    }
    
    my $cluster_recip_hit_ct = keys %clusters_with_recip_hit;
    my $taxa_ct = keys %taxahash;
    
    if ($cluster_recip_hit_ct > 0) {
	while (my ($recip_hit_cluster, $recip_hit_id) =  each %clusters_with_recip_hit) {
	    ####### 
	    ## this is for RAxML
	    my $raxml_recip_hit_id = $recip_hit_id;
	    $raxml_recip_hit_id =~ s/\|.*//;
	    #$raxml_recip_hit_id = "$og_taxa|".$raxml_recip_hit_id;
	    #######
	    print $nt_group join("\n",(">".$raxml_recip_hit_id,${$oghash{ $recip_hit_id  }}[0])),"\n";
	    print $pep_group join("\n",(">".$raxml_recip_hit_id,${$oghash{ $recip_hit_id }}[1])),"\n";
	    close($nt_group);
	    close($pep_group);
	}
    }
    else {
	close($nt_group);
	close($pep_group);
	unlink($gene_file, $pep_file);
    }

    ###### Align $gene_file here
    if (-s $gene_file && -s $pep_file && $gene_ct > 3 && $gene_ct < 500 && $taxa_ct > 3) { # including the outgroup, we want 4 taxa

	my ($gene_aln, $pep_aln) = align($gene_file, $pep_file, $muscle, $gene_ct);
	
	###### refine alignment 
	my ($gene_phy, $gene_aln_fas) = pal2nal($pep_aln, $gene_file, $pal2nal);
	
	###### infer tree
	if (-s $gene_phy) {
	    my $gene_tree = infer_tree($gene_phy, $raxml, $clusters_with_recip_hit{$cluster});
	    next unless defined $gene_tree; ## try to catch errors from muscle or raxml that will cause Bio::Phylo to crash
	    my $tree_stats = summarize_trees($gene_tree, $gene_ct, \%tree_members);
	    push(@treestats,$tree_stats);
	}
    }
    else {
	unlink($gene_file, $pep_file);
    }
    undef %taxahash;
    undef %clusters_with_recip_hit;
    undef %tree_members;
    $gene_ct = 0;
}

# free some memory
undef %seqhash;
untie %seqhash;
undef %oghash;
untie %oghash;

#
# statistical calculations for each ortho-cluster
#
my $count = scalar @clusterstats;
my $mean = mean(@clusterstats);
my $median = median(@clusterstats);
my $min = min(@clusterstats);
my $max = max(@clusterstats);

open(my $cluster_report, ">", $orthogroup_stats) or die "\nERROR: Could not open file: $orthogroup_stats\n";
print $cluster_report "=-=" x 25, "\n";
print $cluster_report join("\t",("Cluster_count","Cluster_mean","Cluster_median","Cluster_min","Cluster_max")),"\n";
print $cluster_report join("\t",($count, $mean, $median, $min, $max)), "\n";
print $cluster_report "=-=" x 25, "\n";
print $cluster_report join("\t",("Cluster_number","Cluster_size")),"\n";

for my $key (reverse sort { $statshash{$a} <=> $statshash{$b} } keys %statshash) {
    print $cluster_report join("\t",($key,$statshash{$key})), "\n";
}
close($cluster_report);

my $gene_tree_ct = 0;
open(my $tree_report, ">", $tree_stats) or die "\nERROR: Could not open file: $tree_stats\n";
print $tree_report join("\t",("Gene_tree","Taxon_code","Ave_branch_length","Gene_tree_members")),"\n";
for my $tree_stats_ref (@treestats) {
    for my $gene_fam_tree (keys %$tree_stats_ref) {
	$gene_tree_ct++;
	for my $taxon (keys %{$tree_stats_ref->{$gene_fam_tree}}) {
	    my ($tr_members, $br_len_ave) = split(/\,/,($tree_stats_ref->{$gene_fam_tree}{$taxon}));
	    print $tree_report join("\t",($gene_fam_tree, $taxon, $tr_members, $br_len_ave)),"\n";
	}
    }
}
print $tree_report "\n$gene_tree_ct total trees with all three taxa and an outgroup match.\n";
close($tree_report);

exit;

#
# subs
#
sub parse_groups {
    # modified from:
    # http://cpansearch.perl.org/src/EASR/ONTO-PERL-1.19/lib/OBO/CCO/OrthoMCLParser.pm

    my $infile = shift;	
    open(my $fh, '<', $infile) or die "\nERROR: Could not open file: $infile\n.";
    
    my %clusters; 
    while(<$fh>){
	my ($cluster, $proteins) = split(/:\s+/xms);
	my $cluster_num;
	if ($cluster =~ /(\d+)/) { # work on this regex
	    $cluster_num = $1;
	}
	$cluster = $cluster_num;
	my @proteins = split(/\s/xms, $proteins);
	my $protein_ct = @proteins;
	foreach my $protein ( @proteins) {
	    if ($protein  =~ /((\w+)\|(\w+))/) {
		if (exists $clusters{$cluster}) {
		    push(@{$clusters{$cluster}}, $1);
		}
		else {
		    $clusters{$cluster} = [ $1 ];
		}
	    }		
	}
    }		
    close($fh);
    return \%clusters;
}

sub parse_best_hits {
    my $infile = shift;	
    open(my $fh, '<', $infile) or die "\nERROR: Could not open file: $infile\n.";
    
    my %clusters_best_match;
    while(<$fh>) {
	next if /^Query/ || /^\#/;
	my @hits = split(/\t/,$_);
	my ($clusternum, $read) = ($hits[0] =~ /(\d+)\|(\w+\|\S+)/);
	$clusters_best_match{$clusternum}{$read} = $hits[1];
    }
    close($fh);
    return \%clusters_best_match;
}

sub summarize_trees {
    my ($gene_tree, $gene_ct, $tree_members) = @_;
    eval { 
	use lib qw(/iob_home/jmblab/statonse/apps/perlmod/Bio-Phylo-0.50/blib/lib);
	require Bio::Phylo::IO;    
        };
    if ($@) {
        die "\nERROR: The Bio::Phylo Perl package is required to work with trees. Exiting.\n";
    }

    my %treehash;
    my %tree_stats;
    my $tree = Bio::Phylo::IO->parse(
        -format => 'newick',
        -file   => $gene_tree
        )->first;

    print "Tree file on line 315 is: ",$gene_tree,"\n";
    for my $node_name (keys %$tree_members) {
	my $node = $tree->get_by_name($node_name);
	my $length = $node->calc_path_to_root;

	my $taxa = $node_name;
	$taxa =~ s/\|.*//;
	if (exists $treehash{$taxa}) {
	    push(@{$treehash{$taxa}}, $length);
	}
	else {
	    $treehash{$taxa} = [ $length ];
	}
	print "Length to root for $node_name is: ",$length,"\n";
    }

    my @keys = keys %treehash;
    for my $key (@keys) {
	my @taxon_tree_lengths = @{$treehash{$key}};
	my $taxon_tree_members = scalar @taxon_tree_lengths;
	my $branch_len_sum = sum(@taxon_tree_lengths);
	my $branch_len_ave = $branch_len_sum / $taxon_tree_members;
	print "Total branch lengths for $key in $gene_tree is: ",$branch_len_sum,"\n";
	print "Average branch length for $key in $gene_tree is: ",$branch_len_ave,"\n";
	$tree_stats{$gene_tree}{$key} = $taxon_tree_members.",".$branch_len_ave; 
  }

    return \%tree_stats;
}

sub infer_tree {
    my ($gene_phy, $raxml, $og_seq_id) = @_;
    my $raxml_tree = $gene_phy;
    $raxml_tree =~ s/\.phy$/\_raxml\.tre/;
    my $raxml_out_tree = "RAxML_bestTree.".$raxml_tree; ## unfortunately, the filename changes with different options
    $og_seq_id =~ s/\|.*//;   # just create a simple id
    
    my ($raxml_out, $raxml_err, @raxml_res) = capture { 
        #system("$raxml -p 34623 -B 0.03 -x 24562 -N 100 -f a -K GTR -m GTRGAMMAI -T 8 -U -s $gene_phy -n $raxml_tree 2>&1 /dev/null");
	system("$raxml -p 34623 -f d -m GTRGAMMA -T 8 -U -o $og_seq_id -s $gene_phy -n $raxml_tree"); 
    };

    print "RAxML output: ",$raxml_out,"\n";
    print "RAxML error: ",$raxml_err,"\n";
    return $raxml_out_tree;
}

sub pal2nal {
    my ($pep_aln, $gene_file, $pal2nal) = @_;
    my $pal2nal_aln = $gene_file;
    $pal2nal_aln =~ s/\.fa.*/\_pal2nal\.aln/;
    my ($pal2nal_out, $pal2nal_err, @pal2nal_res) = capture { system("$pal2nal $pep_aln $gene_file -nogap > $pal2nal_aln"); };

    my $clustalw_to_fasta = sub {
        eval { require Bio::AlignIO;
        };
        if ($@) {
            die "BioPerl is needed to do alignment format conversion. Exiting.\n";
        }
        my $clustalw_aln = shift;
        my $fas_aln = $clustalw_aln;
        $fas_aln =~ s/\.aln$/\.fas/;
        my $aln_in = Bio::AlignIO->new(-file => $clustalw_aln, -format => 'clustalw');
        my $aln_out = Bio::AlignIO->new(-file => ">$fas_aln", -format => 'fasta');
        
        while(my $aln = $aln_in->next_aln ) {
            $aln_out->write_aln($aln);
        }
        
        return $fas_aln;
    };

    my $clustalw_to_phylip = sub {
        eval { require Bio::AlignIO; 
               require Bio::SimpleAlign;
        };
        if ($@) {
            die "BioPerl is needed to do alignment format conversion. Exiting.\n";
        }
        my $clustalw_aln = shift;
        my $phy_outtmp = $clustalw_aln;
        $phy_outtmp =~ s/\.aln$//;
        my $phy_out = $phy_outtmp;
        $phy_outtmp .= "_tmp.phy";
        $phy_out .= ".phy";
        my $aln_in = Bio::AlignIO->new(-file => $clustalw_aln, -format => 'clustalw');
        my $aln_out = Bio::AlignIO->new(-file => ">$phy_outtmp", -format => 'phylip', -longid => 1, -interleaved => 0);
        
        while(my $aln = $aln_in->next_aln ) {
            $aln_out->write_aln($aln);
        }
        
        # this is so RaxML won't complain about pal2nal's clustalw format (including quotes in the identifiers)
        open(my $phy, "<", $phy_outtmp) or die "\nERROR: Could not open file: $phy_outtmp\n";
        open(my $correct_phy, ">", $phy_out) or die "\nERROR: Could not open file: $phy_out\n";
        my ($seqnum, $length);

        my %phyhash;
        my %idhash;
        my $id;

	my $Nseqs = 0;
        while (my $line = <$phy>) {
            chomp $line;
            if ($line =~ /^\s(\d+)\s+(\d+)/){
                $seqnum = $1;
                $length = $2;
            }
            if ($line =~ /^(\'\w+.*\')\s+(\w+)/) { 
                my $seqname = $1;
                my $seq = $2;
                $seqname =~ s/\'//g;
                $seqname =~ s/\s//g;
                $seq =~ s/\s//g;
                if ($seq !~ /A|T|C|G/i) {
                    $Nseqs++;
                    next;
                }
                my $adjusted = $seqnum - $Nseqs;
                $phyhash{$seqname} = $seq;
                $idhash{$adjusted} = $length;
            }
        }
               
        my $ids = 0;
        for my $id (keys %idhash) {
            if ($ids < 1) {
                print $correct_phy " $id $idhash{$id}\n";
            }
            $ids++;
        }

        for my $seq (keys %phyhash) {
            print $correct_phy "$seq  $phyhash{$seq}\n";
        }
        close($phy);
        close($correct_phy);
        unlink ($phy_outtmp);

	return $phy_out;
    };

    my $gene_phy = $clustalw_to_phylip->($pal2nal_aln);
    my $gene_aln_fas = $clustalw_to_fasta->($pal2nal_aln);

    return($gene_phy, $gene_aln_fas); 
}

sub align {
    my ($gene_file, $pep_file, $muscle, $gene_ct) = @_;
    my $gene_aln = $gene_file;
    $gene_aln =~ s/\.fa.*//;
    $gene_aln .= ".aln";
    my $pep_aln = $pep_file;
    $pep_aln =~ s/\.fa.*//;
    $pep_aln .= ".aln";

    print "file for alignment is: ",$gene_file," and: ",$pep_file,"\n";
    #if ($gene_ct > 100) { ## muscle crashes with a core dump on large alignment so we just run the first two iterations
	my ($gene_aln_out, $gene_aln_err, @gene_aln_res) = capture { system("$muscle -in $gene_file -out $gene_aln -maxiters 2"); }; # -quiet
	my ($pep_aln_out, $pep_aln_err, @pep_aln_res) = capture { system("$muscle -in $pep_file -out $pep_aln -maxiters 2"); };
	print "muscle pep aln error on line 475 is: ",$pep_aln_err,"\n";
	print "muscle gene aln error on line 476 is: ",$gene_aln_err,"\n";
	print "muscle pep out on line 477 is: ",$pep_aln_out,"\n";
	print "muscle gene out on line 477 is: ",$gene_aln_out,"\n";
    #}
    #else {
	#my ($gene_aln_out, $gene_aln_err, @gene_aln_res) = capture { system("$muscle -in $gene_file -out $gene_aln"); };
        #my ($pep_aln_out, $pep_aln_err, @pep_aln_res) = capture { system("$muscle -in $pep_file -out $pep_aln"); };
	#print "muscle pep aln error on line 483 is: ",$pep_aln_err,"\n";
	#print "muscle gene aln error on line 484 is: ",$gene_aln_err,"\n";
	#print "muscle pep aln out on line 485 is: ",$pep_aln_out,"\n";
	#print "muscle gene aln out on line 486 is: ",$gene_aln_out,"\n";
    #}
    return($gene_aln, $pep_aln);
}

sub sum {
    my @array = @_;
    my $sum;
    foreach (@array) { $sum += $_; }
    return $sum;
}

sub min {
    my $min = shift;
    foreach ( @_ ) { $min = $_ if $_ < $min }
    return $min;
}

sub max {
    my $max = shift;
    foreach ( @_ ) { $max = $_ if $_ > $max }
    return $max;
}

sub mean { 
    my @array = @_; 
    my $sum = sum(@array); 
    my $count = scalar @array; 
    return sprintf("%.2f",$sum / $count); 
}

sub median {
    my @orig_array = @_;
    my @array = sort {$a <=> $b} @orig_array;
    if ($#array % 2 == 0) {
        my $median = $array[($#array / 2)];
        return $median;
    }
    else {
        my $median = $array[int($#array / 2)] + (($array[int($#array / 2) + 1] - $array[int($#array / 2)]) / 2);
        return $median;
    }
}
    
sub find_prog {
    my $prog = shift;
    my ($path, $err) = capture { system("which $prog"); };
    chomp($path);

    given ($path) {
	when (/$prog$/) { say "Using $prog located at $path." }
	when (/no $prog in/) { say "Could not find $prog. Try installing it or adding it to the PATH. Exiting."; exit(1); }
	when ('') { say "Could not find $prog. Try installing it or adding it to the PATH. Exiting."; exit(1); }
    }
    return $path;
}

sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-g] [-gb] [-nf] [-pf] [-onf] [-opf] [-os] [-ts]  

Required:
     -g|groups            :       The orthoMCL groups.txt file produced by the program mcl.
     -gb|groups_blast     :       A blast (tab-delimited) report of each cluster representative to an outgroup.
                                  (NB: format is: >clusternumber|taxoncode|sequenceid)
     -onf|og_nt_fas       :       Nucleotide file for the outgroup species in the groups_blast file (NB: sequences must be in same format as in blast file).
     -opf|og_pep_fas      :       Peptide file for the outgroup speces in the groups_blast file (NB: sequences must be in same format as in blast file).
     -nf|nt_fas           :       The goodGenes.fasta file (you must create this).
     -pf|pep_fas          :       The goodProteins.fasta file.
     -os|orthogroup_stats :       A report of the length statistics for each cluster.
     -ts|tree_stats       :       A report of the branch length statistics for each gene tree.

Options:
     -h|help              :       Print a usage statement.
     -m|man               :       Print the full documentation.

EOF
}

