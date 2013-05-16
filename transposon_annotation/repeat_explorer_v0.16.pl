#!/usr/bin/env perl

## TODO: 
##       
##
##      

use utf8;
use v5.12;
use strict;
use warnings;
use warnings FATAL => "utf8";
use open qw(:std :utf8);
use autodie qw(open);
use Getopt::Long;
use Capture::Tiny qw(:all);
use File::Basename;
use File::Path qw(make_path);
use File::Spec;
use File::Basename;
use POSIX qw(strftime);
use Graph::UnionFind;
use List::Util qw(sum max);
use JSON;
use Cwd;
BEGIN {
  @AnyDBM_File::ISA = qw( DB_File SQLite_File )
      unless @AnyDBM_File::ISA == 1;
}
use AnyDBM_File;                  
use vars qw( $DB_BTREE &R_DUP );  
use AnyDBM_File::Importer qw(:bdb);
use DBM::Deep;
use charnames qw(:full :short);
use Encode qw(encode decode);

# lexical vars
my $infile;
my $outdir;
my $percent_id;
my $percent_cov;
my $memory;
my $cluster_size;    # this influences what will get annotated, not plotted in the size-dist graph
my $merge_threshold;
my $fas_file;
my $cores;
my $report;
my $graph;
my $database;
my $json;
my $evalue;

GetOptions(
           'i|infile=s'          => \$infile,
	   'f|fas_file=s'        => \$fas_file,
           'o|outdir=s'          => \$outdir,
	   'im|in_memory'        => \$memory,
	   'p|cpus=i'            => \$cores,
	   'r|report=s'          => \$report,
	   'cs|cluster_size=i'   => \$cluster_size,
	   'id|percent_id=f'     => \$percent_id,
	   'cov|percent_cov=f'   => \$percent_cov,
	   'm|merge_threshold=i' => \$merge_threshold,
	   'g|graph'             => \$graph,
	   'd|database=s'        => \$database,
	   'j|repbase_json=s'    => \$json,
	   'e|evalue=f'          => \$evalue,
	   );

# open the infile or die with a usage statement
if (!$infile || !$outdir || !$fas_file || !$database
    || !$percent_id || !$percent_cov || !$report || !$json) {
    print "\nERROR: Input not parsed correctly.\n";
    usage() and exit(1);
}

### Set paths to be used
my $cwd = getcwd();
my $str = POSIX::strftime("%m_%d_%Y_%H_%M_%S", localtime);

### Parse blast into form for clustering
my ($idx_file, $int_file, $hs_file) = parse_blast($infile, $percent_id, $percent_cov, $outdir, $memory);
#my $hitsort_int = $outfile.".int";     ##// for debugging community
#my $index_file = $outfile.".index";    ##// for debugging community

### Clustering
my $community = louvain_method($int_file, $outdir);
my $cls_file  = make_clusters($community, $int_file, $idx_file, $outdir);
#my $cls_file = make_clusters($community, $int_file, $idx_file, $outdir);    ##// this is for debugging community

### Find union in clusters
my ($seqs, $seqct) = fas2hash($fas_file, $memory);
my ($read_pairs, $vertex, $uf) = find_pairs($cls_file, $report, $merge_threshold);
my ($cls_dir_path, $cls_with_merges_path, $cls_tot) = merge_clusters($infile, $str, $vertex, $seqs, $read_pairs, $report, $cluster_size, $outdir);
untie %$seqs if defined $memory;

### Annotate clusters, produce summary of merged and non-merged cluster size distribution
#annotate_clusters($cls_dir_path, $database, $json, $report, $outdir, $evalue);
my $anno_rp_path = annotate_clusters($cls_dir_path, $database, $json, $report, $outdir, $evalue, $seqct, $cls_tot);
#////////////// end here 3:19 - 4/3 - SES
clusters_annotation_to_summary($anno_rp_path, $anno_sum_rep_path, $total_readct,
			       $seqct, $rep_frac, \@blasts, \@superfams, $report);

exit; ## This is the end

#
# Subs
#
sub parse_blast {
    my ($infile, $percent_id, $percent_cov, $outdir, $memory) = @_;

    my ($iname, $ipath, $isuffix) = fileparse($infile, qr/\.[^.]*/);
    unless (-d $outdir) {
	make_path($outdir, {verbose => 0, mode => 0771,});
    }
    my $int_file = $iname;
    my $idx_file = $iname;
    my $hs_file = $iname;
    $int_file .= "_louvain.int";
    $idx_file .= "_louvain.idx";
    $hs_file .= "_louvain.hs";
    my $int_path = File::Spec->catfile($outdir, $int_file);
    my $idx_path = File::Spec->catfile($outdir, $idx_file);
    my $hs_path = File::Spec->catfile($outdir, $hs_file);

    # counters
    my $total_hits = 0;
    my $parsed_hits = 0;
    my $index = 0;

    open my $in, '<', $infile;
    open my $int, '>', $int_path;
    open my $idx, '>', $idx_path;
    open my $hs, '>', $hs_path;

    if (defined $memory) {
	my %match_pairs;
	my %match_index;

	while (<$in>) {
	    chomp;
	    my ($q_name, $q_len, $q_start, $q_end, $s_name, $s_len,
		$s_start, $s_end, $pid, $score, $e_val, $strand) = split;

	    my $pair = mk_key($q_name, $s_name);
	    my $revpair = mk_key($s_name, $q_name);
	    my $subj_hit_length = ($s_end - $s_start) + 1;
	    my $subj_cov = $subj_hit_length/$s_len;

	    if ($q_start > $q_end) {
		$total_hits++;
		my $neg_query_hit_length = ($q_start - $q_end) + 1;
		my $neg_query_cov = $neg_query_hit_length/$q_len;

		if ( ($neg_query_cov >= $percent_cov) && ($pid >= $percent_id) ) {
		    if (exists $match_pairs{$pair}) {
			push @{$match_pairs{$pair}}, $score;
		    }
		    elsif (exists $match_pairs{$revpair}) {
			push @{$match_pairs{$revpair}}, $score;
		    }
		    else {
			$match_pairs{$pair} = [$score];
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
                my $pos_query_cov = $pos_query_hit_length/$q_len;

                if ( ($pos_query_cov >= $percent_cov) && ($pid >= $percent_id) ) {
                    if (exists $match_pairs{$pair}) {
                        push @{$match_pairs{$pair}}, $score;
                    }
                    elsif (exists $match_pairs{$revpair}) {
                        push @{$match_pairs{$revpair}}, $score;
                    }
                    else {
                        $match_pairs{$pair} = [$score];
                        $match_index{$q_name} = $index unless exists $match_index{$q_name};
                        $index++;
                        $match_index{$s_name} = $index unless exists $match_index{$s_name};
                        $index++;
                    }
                }
            }
        }
        close $in;

        for my $idx_mem (sort { $match_index{$a} <=> $match_index{$b} } keys %match_index) {
            say $idx join " ", $idx_mem, $match_index{$idx_mem};
        }
        close $idx;

        while (my ($match, $scores) = each %match_pairs) {
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

	return($idx_path, $int_path, $hs_path);
    }
    else {
	my $dbm = "mgblast_matchpairs.dbm";
	my $dbi = "mgblast_matchindex.dbm";
	
	unlink $dbm if -e $dbm;
	unlink $dbi if -e $dbi;
	
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

	while (<$in>) {
	    chomp;
	    my ($q_name, $q_len, $q_start, $q_end, $s_name, $s_len,
		$s_start, $s_end, $pid, $score, $e_val, $strand) = split;

	    my $pair = mk_key($q_name, $s_name);
	    my $revpair = mk_key($s_name, $q_name);
	    my $enc_pair = encode("UTF-8", $pair, 1);
	    my $enc_revpair = encode("UTF-8", $revpair, 1);
	    my $subj_hit_length = ($s_end - $s_start) + 1;
	    my $subj_cov = $subj_hit_length/$s_len;

	    if ($q_start > $q_end) {
		$total_hits++;
		my $neg_query_hit_length = ($q_start - $q_end) + 1;
		my $neg_query_cov = $neg_query_hit_length/$q_len;

		if ( ($neg_query_cov >= $percent_cov) && ($pid >= $percent_id) ) {
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
	    else {
		$total_hits++;
		my $pos_query_hit_length = ($q_end - $q_start) + 1;
		my $pos_query_cov = $pos_query_hit_length/$q_len;

		if ( ($pos_query_cov >= $percent_cov) && ($pid >= $percent_id) ) {
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
	close $in;

	for my $idx_mem (sort { $match_index{$a} <=> $match_index{$b} } keys %match_index) {
	    say $idx join " ", $idx_mem, $match_index{$idx_mem};
	}
	close $idx;

	while (my ($match, $scores) = each %$db) {
	    my $enc_match = encode("UTF-8", $match, 1);
	    my $match_score = max(@$scores);
	    my ($qry, $sbj) = mk_vec($enc_match);
	    my $revmatch = mk_key($sbj, $qry);
	    my $enc_revmatch = encode("UTF-8", $revmatch, 1);
	    if (exists $db->{$enc_revmatch}) {
		my $rev_match_score = max(@{$db->{$enc_revmatch}});
		if ($rev_match_score > $match_score) {
		    if (exists $match_index{$sbj} && exists $match_index{$qry}) {
			say $hs join "\t", $sbj, $qry, $rev_match_score;
			say $int join "\t", $match_index{$sbj}, $match_index{$qry}, $rev_match_score;
			delete $db->{$enc_match};
		    }
		}
		else {
		    delete $db->{$enc_revmatch};
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

	untie %match_index;

	return($idx_path, $int_path, $hs_path);

    }
}

sub louvain_method {
    my ($int_file, $outdir) = @_;
    #chdir($outdir) || die "\nERROR: Could not change $outdir: $!\n";
    my ($iname, $ipath, $isuffix) = fileparse($int_file, qr/\.[^.]*/);
    my $cls_bin = $iname.".bin";                    # Community "bin" format
    my $cls_tree = $iname.".tree";                  # hierarchical tree of clustering results
    my $cls_tree_weights = $cls_tree.".weights";    # bit score, the weights applied to clustering
    my $cls_tree_log = $cls_tree.".log";            # the summary of clustering results at each level of refinement
    my $hierarchy_err = $cls_tree.".hierarchy.log"; # hierarchical tree building log (not actually used)

    my $cls_bin_path = File::Spec->catfile($outdir, $cls_bin);
    my $cls_tree_path = File::Spec->catfile($outdir, $cls_tree);
    my $cls_tree_weights_path = File::Spec->catfile($outdir, $cls_tree_weights);
    my $cls_tree_log_path = File::Spec->catfile($outdir, $cls_tree_log);
    my $hierarchy_err_path = File::Spec->catfile($outdir, $hierarchy_err);

    #try {
    system("louvain_convert -i $int_file -o $cls_bin_path -w $cls_tree_weights_path");
    #catch {
    #unless (defined $cls_bin && defined $cls_tree_weights) {
	#say "ERROR: louvain_convert failed. Exiting." and exit(1);
    #}
    system("louvain_community $cls_bin_path -l -1 -w $cls_tree_weights_path -v >$cls_tree_path 2>$cls_tree_log_path");
    unless (defined $cls_tree && defined $cls_tree_log) { # this is not a good way to catch errors
	say "ERROR: louvain_community failed. Exiting." and exit(1);
    }

    my $levels = `grep -c level $cls_tree_log_path`;
    chomp $levels;

    my @comm;
    for (my $i = 0; $i <= $levels-1; $i++) {
        my $cls_graph_comm = $cls_tree.".graph_node2comm_level_".$i;
	#my $cls_graph_clusters = $cls_tree_path.".graph.clusters";
	#my $cls_graph_membership = $cls_tree_path.".graph.membership";
	my $cls_graph_comm_path = File::Spec->catfile($outdir, $cls_graph_comm);

        system("louvain_hierarchy $cls_tree_path -l $i > $cls_graph_comm_path");
  	push @comm, $cls_graph_comm;
    }
    return \@comm;
}

sub make_clusters {
    my ($graph_comm, $int_file, $idx_file, $outdir) = @_;

    my ($iname, $ipath, $isuffix) = fileparse($int_file, qr/\.[^.]*/);
    my $cluster_file = $iname.".cls";
    my $membership_file = $cluster_file.".membership.txt";
    my $cluster_file_path = File::Spec->catfile($outdir, $cluster_file);
    my $membership_file_path = File::Spec->catfile($outdir, $membership_file);

    my @graph_comm_sort = reverse sort { ($a =~ /(\d)$/) <=> ($b =~ /(\d)$/) } @$graph_comm;
    my $graph = shift @graph_comm_sort;
    die "\nERROR: Community clustering failed. Exiting.\n" unless defined $graph;
    my $graph_path = File::Spec->catfile($outdir, $graph);
    my %clus;
    my %index;

    open my $idx, '<', $idx_file;
    while (my $idpair = <$idx>) {
        chomp $idpair;
        my ($readid, $readindex) = split /\s+/, $idpair;
        $index{$readindex} = $readid;
    }
    close $idx;

    open my $mem,'>', $membership_file_path;
    open my $in, '<', $graph_path;
    open my $cls_out, '>', $cluster_file_path;

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
    close $in;

    my $cls_ct = 1;
    for my $cls (reverse sort { @{$clus{$a}} <=> @{$clus{$b}} } keys %clus) {
        my $clus_size = scalar @{$clus{$cls}};
        say $cls_out ">CL$cls_ct $clus_size";
        my @clus_members;
        for my $cls_member (@{$clus{$cls}}) {
            say $mem "$cls_member $cls_ct";
            if (exists $index{$cls_member}) {
                push @clus_members, $index{$cls_member};
            }
        }
        say $cls_out join " ", @clus_members;
        $cls_ct++;
    }
    close $cls_out;
    close $mem;

    return $cluster_file;
}

sub merge_clusters {
    my ($infile, $str, $vertex, $seqs, $read_pairs, $report, $cluster_size, $outdir) = @_;

    my ($iname, $ipath, $isuffix) = fileparse($infile, qr/\.[^.]*/);
    my $cls_dir_base = $iname;
    #$cls_dir_base =~ s/\.[^.]+$//;     # this is messing up the filenames
    my $cls_with_merges = $cls_dir_base;
    my $cls_dir = $cls_dir_base."_cls_fasta_files_$str";
    $cls_with_merges .= "_merged_$str.cls";
    my $cls_dir_path = $ipath.$cls_dir;
    make_path($outdir."/".$cls_dir_path, {verbose => 0, mode => 0711,}); # allows for recursively making paths                                                                
    my $cls_with_merges_path = File::Spec->catfile($outdir, $cls_with_merges);
    open my $clsnew, '>', $cls_with_merges_path;

    my ($rpname, $rppath, $rpsuffix) = fileparse($report, qr/\.[^.]*/);
    my $rp_path = File::Spec->rel2abs($rppath.$rpname.$rpsuffix);
    open my $rep, '>>', $rp_path;

    $cluster_size //= 500;
    my $cls_tot = 0;

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
	$cls_tot += $groupseqnum;
	say $rep join "\t", $group_index, join ",", @grpcp;
	say $clsnew ">G$group_index $groupseqnum";
	my $group_file = "G$group_index"."_$groupseqnum".".fas";
	my $group_file_path = File::Spec->catfile($outdir."/".$cls_dir_path, $group_file);
	open my $groupout, '>', $group_file_path;
    
	for my $clus (@$group) {
	    if (exists $read_pairs->{$clus}) {
		print $clsnew join " ",@{$read_pairs->{$clus}};
		for my $read (@{$read_pairs->{$clus}}) {
		    if (exists $seqs->{$read}) {
			say $groupout join "\n", ">".$read, $seqs->{$read};
		    }
		    else {
			say "WARNING: $read not found. This is possibly a bug. Please report it.";
		    }
		}
	    }
	    print $clsnew " ";
	    delete $read_pairs->{$clus}
	}
	print $clsnew "\n";
	close $groupout;
	$group_index++;
    }

    # write out those clusters that weren't merged
    say $rep "=====> Non-grouped clusters";
    for my $non_paired_cls (keys %$read_pairs) {
	my ($non_paired_clsid, $non_paired_clsseqnum) = split /\_/, $non_paired_cls, 2;
	$cls_tot += $non_paired_clsseqnum;
	say $rep $non_paired_clsid;
	say $clsnew join "\n", ">$non_paired_clsid $non_paired_clsseqnum", join " ", @{$read_pairs->{$non_paired_cls}};

	if (scalar(@{$read_pairs->{$non_paired_cls}}) >= $cluster_size) {
	    my $non_paired_clsfile .= $non_paired_cls.".fas";
	    my $cls_file_path = File::Spec->catfile($outdir."/".$cls_dir_path, $non_paired_clsfile);
	    open my $clsout, '>', $cls_file_path;

	    for my $non_paired_read (@{$read_pairs->{$non_paired_cls}}) {
		if (exists $seqs->{$non_paired_read}) {
		    say $clsout join "\n", ">".$non_paired_read, $seqs->{$non_paired_read};
		}
		else {
		    say "WARNING: $non_paired_read not found. This is possibly a bug. Please report it.";
		}
	    }
	    close $clsout;
	}
    }
    close $rep;
    close $clsnew;

    return ($cls_dir_path, $cls_with_merges_path, $cls_tot);
}

sub fas2hash {
    my ($fas_file, $memory) = shift;
    open my $fas, '<', $fas_file;

    my %seqhash;
    if (defined $memory) {
	$DB_BTREE->{cachesize} = 100000;
	$DB_BTREE->{flags} = R_DUP;
	my $seq_dbm = "repeat_explorer_seqs.dbm";
	tie %seqhash, 'AnyDBM_File', $seq_dbm, O_RDWR|O_CREAT, 0666, $DB_BTREE
	     or die "\nERROR: Could not open DBM file $seq_dbm: $!\n";
    }
    my $seqct = 0;
    local $/ = '>';

    while (my $line = <$fas>) {
	$line =~ s/>//g;
        next if !length($line);
        my ($seqid, @seqs) = split /\n/, $line;
	my $seq = join '', @seqs;
        $seqhash{$seqid} = $seq;
        $seqct++ if defined $seq;
    }
    close $fas;
    
    return(\%seqhash, $seqct);
}

sub find_pairs {
    my ($cls_file, $report, $merge_threshold) = @_;
    
    my ($rpname, $rppath, $rpsuffix) = fileparse($report, qr/\.[^.]*/);
    my $rp_path = File::Spec->rel2abs($rppath.$rpname.$rpsuffix);
    my ($clname, $clpath, $clsuffix) = fileparse($cls_file, qr/\.[^.]*/);
    my $cls_file_path = File::Spec->rel2abs($clpath.$outdir."/".$clname.$clsuffix);
    open my $rep, '>', $rp_path;

    $merge_threshold //= 500;

    my $uf = Graph::UnionFind->new;

    say $rep "=====> Cluster connections above threshold";

    my %vertex;
    my %read_pairs;
    my %mapped_pairs;

    {
        local $/ = '>';
        
        open my $in, '<', $cls_file_path;   
        while (my $line = <$in>) {
            $line =~ s/>//g;
            next if !length($line);
            my ($clsid, $seqids) = split /\n/, $line;
            $clsid =~ s/\s/\_/;
            my @ids = split /\s+/, $seqids;
            #if (scalar(@ids) >= $cluster_size) {       # limit cluster size in .cls file here, if desired
            push @{$read_pairs{$clsid}}, $_ for @ids;
            #}
        }
        close $in;
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
		my $k = mk_key($cls_i, $cls_j);
                $cls_conn_ct{$k}++;
            }
        }
        @sep_reads = ();
    }

    for my $p (reverse sort { $cls_conn_ct{$a} <=> $cls_conn_ct{$b} } keys %cls_conn_ct) {
	my ($i, $j) = mk_vec($p);
        my $i_noct = $i; $i_noct =~ s/\_.*//;
        my $j_noct = $j; $j_noct =~ s/\_.*//;
        if ($cls_conn_ct{$p} >= $merge_threshold) {   
            say $rep join "\t", $i_noct, $j_noct, $cls_conn_ct{$p};
            ++$vertex{$_} for $i, $j;
            $uf->union($i, $j);
        }
    }
    close $rep;
    return(\%read_pairs, \%vertex, \$uf);
}

# the best way to make combined keys is not to concatenate them
# http://stackoverflow.com/a/15299397/1543853
#
# NB: in v5.16 charnames() will be loaded automatically
sub mk_key { join "\N{INVISIBLE SEPARATOR}", @_ }

sub mk_vec { split "\N{INVISIBLE SEPARATOR}", shift }

sub json2hash {
    my $json = shift;
   
    my $json_text;
    local $/;
    open my $in, '<', $json;
    $json_text = <$in>;
    close $in;

    my $repeats = JSON->new->utf8->space_after->decode($json_text);
    return $repeats;
}

sub parse_blast_to_top_hit {
    my ($blast_out, $blast_file_path) = @_;
    my %blhits;

    my $top_hit;
    my $top_hit_num = 0;
    my $hit_ct = 0;

    for my $hit (@$blast_out) {
        chomp $hit;
        my ($ct, $hittype) = split /\t/, $hit;
        next unless defined $ct;
        $blhits{$hittype} = $ct;
        $hit_ct++;
    }
    
    my $sum = sum values %blhits;
    if ($hit_ct > 0) {
        open my $out, '>', $blast_file_path;
        $top_hit = (reverse sort { $blhits{$a} <=> $blhits{$b} } keys %blhits)[0];
        keys %blhits; #reset iterator                                                                                                                         
        for my $hits (reverse sort { $blhits{$a} <=> $blhits{$b} } keys %blhits) {
            #say $out join "\t", $hits, $blhits{$hits};
	    my $perc = sprintf("%.2f", $blhits{$hits} / $sum);
	    say $out join "\t", $hits, $blhits{$hits}, $perc;
        }
        close $out;
        return \$hit_ct, \$top_hit, \%blhits;
    }
    else { ## if (!%blhits) {                                                                                                                                 
        unlink $blast_file_path;
        return undef, undef, undef;
    }
}

sub annotate_clusters {
    my ($cls_with_merges_dir, $database, $json, $report, $outdir, $evalue, $seqct, $cls_tot) = @_;

    my ($rpname, $rppath, $rpsuffix) = fileparse($report, qr/\.[^.]*/);
    my $rp_path = File::Spec->rel2abs($rppath.$rpname.$rpsuffix);
    open my $rep, '>>', $rp_path;
    my $anno_rep = $rpname."_annotations.tsv";
    my $anno_summary_rep = $rpname."_annotations_summary.tsv";
    my $anno_rp_path = File::Spec->rel2abs($rppath.$anno_rep);
    my $anno_sum_rep_path = File::Spec->rel2abs($rppath.$anno_summary_rep);
    my $total_readct = 0;
    $evalue //= 10;
    my $rep_frac = $cls_tot / $seqct;
    say $rep "======> Total seqs: ",$seqct;
    say $rep "======> Total clustered: ",$cls_tot;
    say $rep "======> Repeat fraction: ",$rep_frac;
    close $rep;
    my $top_hit_superfam = {};

    ## get input files                                                                                                                                        
    opendir my $dir, $outdir."/".$cls_with_merges_dir || die "\nERROR: Could not open directory: $cls_with_merges_dir. Exiting.\n";
    my @clus_fas_files = grep /\.fa.*$/, readdir $dir;
    closedir $dir;


    if (scalar @clus_fas_files < 1) {
        say "\nERROR: Could not find any fasta files in $cls_with_merges_dir. Exiting.\n";
        exit(1);
    }

     ## set path to output dir
    my $annodir = $outdir."/".$cls_with_merges_dir."_annotations";
    #my ($oname, $opath, $osuffix) = fileparse($annodir); # for working with directory names with periods
    #my $out_path = File::Spec->rel2abs($opath.$oname);
    my $out_path = File::Spec->rel2abs($annodir);
    make_path($annodir, {verbose => 0, mode => 0711,}); # allows for recursively making paths                                                                 
    my ($dname, $dpath, $dsuffix) = fileparse($database, qr/\.[^.]*/);
    my $db_path = File::Spec->rel2abs($dpath.$dname);
    my @blasts; # container for each report (hash) // need to rethink how duplicate entries will be handled                                                   
    my @superfams;

    open my $out, '>>', $anno_rp_path;
    chdir $cls_with_merges_dir;

    say $out join "\t", "Cluster", "Read_count", "Type", "Class", "Superfam", "(SINE_family; if present)","Top_hit";
    for my $file (@clus_fas_files) {
        my $query = $outdir."/".$cls_with_merges_dir."/".$file;
        my ($fname, $fpath, $fsuffix) = fileparse($query, qr/\.[^.]*/);
        my $blast_res = $fname;
        my ($filebase, $readct) = split /\_/, $fname, 2;
        $total_readct += $readct;
        $blast_res =~ s/\.[^.]+$//;
        $blast_res .= "_blast_$evalue.tsv";
        my $blast_file_path = File::Spec->catfile($out_path, $blast_res);

        my $blastcmd = "blastn -dust no -query $query -evalue $evalue -db $db_path -outfmt 6 | ".
                       "sort -k1,1 -u | ".                       # count each read in the report only once                                                 
                       "cut -f2 | ".                             # keep only the ssids        
                       "sort | ".                                # sort the list
                       "uniq -c | ".                             # reduce the list
                       "sort -bnr | ".                           # count unique items
                       "perl -lane 'print join(\"\\t\",\@F)'";   # create an easy to parse format
	my @blast_out = qx($blastcmd);

        my ($hit_ct, $top_hit, $blhits) = parse_blast_to_top_hit(\@blast_out, $blast_file_path);
        next unless defined $top_hit && defined $hit_ct;
        #say $file, " ==> ", $blast_res, " ==> ", dd $blhits;                                                                                                 
        push @blasts, $blhits;
        $top_hit_superfam = blast2annot($json, $filebase, $readct, $top_hit, $out);
        push @superfams, $top_hit_superfam unless !%$top_hit_superfam;
    }
    close $out;

    return $anno_rp_path;
}

sub clusters_annotation_to_summary  {
    my ($anno_rp_path, $anno_sum_rep_path, $total_readct, 
	$seqct, $rep_frac, $blasts, $superfams, $report) = @_;

    my ($rpname, $rppath, $rpsuffix) = fileparse($report, qr/\.[^.]*/);
    my $rp_path = File::Spec->rel2abs($rppath.$rpname.$rpsuffix);
    open my $rep, '>>', $rp_path;

    my %top_hit_superfam;
    @top_hit_superfam{keys %$_} = values %$_ for @$superfams;

    #dd \%top_hit_superfam;                                                                                                                                   
    for my $f (keys %top_hit_superfam) {
        if ($f =~ /^((RL.\-\w+)\-\d)/) {
            my $fam = $2;
            $top_hit_superfam{$fam} = $top_hit_superfam{$f};
            delete $top_hit_superfam{$f};
        }
    }
    #dd \%top_hit_superfam;                                                                                                                                   
    open my $outsum, '>', $anno_sum_rep_path;

    my %annot;
    my %fams;
    my $total_ct = 0;
    my $hashct = scalar @$blasts;
    my $hitct;
    for my $blast (@$blasts) {
        for my $fam (keys %$blast) {
            $total_ct += $blast->{$fam};
            if ($fam =~ /^((RL.\-|\_\w+)\-|\_\d)/) {
                my $famname = $2;
                if (exists $fams{$famname}) {
                    $fams{$famname} += $blast->{$fam};
                }
                else {
                    $fams{$famname} = $blast->{$fam};
                }
            }
            else {
                if (exists $fams{$fam}) {
                    $fams{$fam} += $blast->{$fam};
                }
                else {
                    $fams{$fam} = $blast->{$fam};
                }
            }
        }
    }
    my $total_gcov = 0;
    ### NEED TO SIMPLIFY THE STATS REPORTED, once I know they're correct 
    say $outsum join "\t", "ReadNum", "Superfamily", "Family", "ReadCt/ReadsWithHit", "HitPerc", "GPerc";
    for my $k (reverse sort { $fams{$a} <=> $fams{$b} } keys %fams) {
        if (exists $top_hit_superfam{$k}) {
	    my $hit_perc = sprintf("%.12f",$fams{$k}/$total_ct);
	    my $gperc_corr = $hit_perc * $rep_frac;
            $total_gcov += $gperc_corr;
            say $outsum join "\t", $seqct, $top_hit_superfam{$k}, $k, $fams{$k}."/".$total_ct, $hit_perc, $gperc_corr;
        }
    }
    close $outsum;
    say $rep "======> Total repeat fraction from annotations: ",$total_gcov;
    close $rep;
}

sub blast2annot {
    my ($json, $filebase, $readct, $top_hit, $out) = @_;

    my %top_hit_superfam;
    my $repeats = json2hash($json);

    for my $type (keys %$repeats) {
        if ($type eq 'pseudogene' || $type eq 'simple_repeat' || $type eq 'integrated_virus') {
            if ($type eq 'pseudogene' && $$top_hit =~ /rrna|trna|snrna/i) {
                say $out join "\t", $filebase, $readct, $type, $$top_hit;
            }
            elsif ($type eq 'simple_repeat' && $$top_hit =~ /msat/i) {
                say $out join "\t", $filebase, $readct, $type, "Satellite", "MSAT", $$top_hit;
            }
            elsif ($type eq 'simple_repeat' && $$top_hit =~ /sat/i) {
                say $out join "\t", $filebase, $readct, $type, "Satellite", "SAT", $$top_hit;
            }
            elsif ($type eq 'integrated_virus' && $$top_hit =~ /caul/i) {
                say $out join "\t", $filebase, $readct, $type, "Caulimoviridae", $$top_hit;
            }
            elsif ($type eq 'integrated_virus' && ($$top_hit eq 'PIVE' || $$top_hit eq 'DENSOV_HM')) {
                say $out join "\t", $filebase, $readct, $type, "DNA Virus", $$top_hit;
            }
            elsif ($type eq 'endogenous_retrovirus' && $$top_hit =~ /erv/i) {
                say $out join "\t", $filebase, $readct, $type, "Endogenous Retrovirus", $$top_hit;
            }
            next;
        }
	for my $class (keys %{$repeats->{$type}}) {
            while ( my ($superfam_index, $superfam) = each @{$repeats->{$type}{$class}} ) {
                for my $superfam_h (keys %$superfam) {
                    if ($superfam_h =~ /sine/i) {
                        while (my ($sine_fam_index, $sine_fam_h) = each @{$superfam->{$superfam_h}}) {
                            for my $sine_fam_mem (keys %$sine_fam_h) {
                                for my $sines (@{$repeats->{$type}{$class}[$superfam_index]{$superfam_h}[$sine_fam_index]{$sine_fam_mem}}) {
                                    for my $sine (@$sines) {
                                        if ($sine =~ /$$top_hit/) {
                                            $top_hit_superfam{$$top_hit} = $sine_fam_mem;
                                            say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $sine_fam_mem, $$top_hit;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    elsif ($superfam_h =~ /gypsy/i && $$top_hit =~ /^RLG/) {
                        $top_hit_superfam{$$top_hit} = $superfam_h;
                        say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $$top_hit;
                    }
                    elsif ($superfam_h =~ /copia/i && $$top_hit =~ /^RLC/) {
                        $top_hit_superfam{$$top_hit} = $superfam_h;
                        say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $$top_hit;
                    }
		    else {
                        for my $fam (@{$repeats->{$type}{$class}[$superfam_index]{$superfam_h}}) {
                            for my $mem (@$fam) {
                                if ($mem =~ /$$top_hit/i) {
                                    $top_hit_superfam{$$top_hit} = $superfam_h;
                                    say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $$top_hit;
                                }
                                #else { say "NOMATCH $mem => $$top_hit"; }                                                                                    
                            }
                        }
                    }
                }
            }
        }
    }
    return \%top_hit_superfam;
}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -i inreport -o hitsort -id 90.00 -cov 0.55 -f fasta_file 

Required:
    -i|infile           :    An mgblast report in tab format (-D 4).
    -f|fas_file         :    The (Fasta) file of sequences used in the all-vs-all blast.
    -o|outdir           :    A directory to place all clustering and annotation results.
    -p|cpus             :    The number of processors to use for plotting.
    -cs|cluster_size    :    The minimum cluster size to convert to fasta (Default: 500).
    -id|percent_id      :    The percent identity threshold for matches (ex. 90.00).
    -cov|percent_cov    :    The percent coverage for both the query and subject to be retained for clustering (ex. 0.55).
    -r|report           :    A file to hold the cluster stats.
    -d|database         :    Name a repeat database (made with makeblastdb) to use for annotation.
    -im|in_memory       :    Perform all calculations in memory (Default: No).
                             (NB: The is recommended as it will greatly increase the speed of the analysis.)
    -j|repbase_json     :    RepBase database in JSON format for annotation.

Options:
    -g|graph            :    Generate a graph file for each cluster, and a summary PDF for all the graphs.
                             (NB: This is very time-consuming and intensive.)
    -m|merge_threshold  :    Threshold for merging clusters based on split paired-end reads (Default: 500).
                             (NB: A logical threshold to pick is dependent on the number of reads.)
    -e|evalue           :    The BLASTN E-value threshold to be used for annotation (Default: 10).

END
}
