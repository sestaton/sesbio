#!/usr/bin/env perl

## TODO: add time stamp to GL and ncol directories DONE 
##       --------> Okay, now fix the printing (just printing "8") DONE
##       write method for annotating clusters
##       turn off printing in R code and do it with Perl
##       final R routine is still printing "NULL" to the screen (need to capture dev.off())
##       Don't store all the blast scores in an array, only keep the top hit DONE
##
##       //This version is an attempt to generate larger clusters by keeping all BLAST hits above a threshold\\
use utf8;
use v5.12;
use strict;
use warnings;
use warnings qw(FATAL utf8);
use open qw(:std :utf8);
use autodie qw(open);
use Getopt::Long;
use Capture::Tiny qw(:all);
use File::Basename;
use File::Path qw(make_path);
use File::Spec;
use File::Basename;
#use List::Util qw(max);       ## shouldn't need this routine if we only keep the top hit
use POSIX qw(strftime);
use Graph::UnionFind;
use Cwd;
BEGIN {
  @AnyDBM_File::ISA = qw( DB_File SQLite_File )
      unless @AnyDBM_File::ISA == 1; # if loaded already, AnyDBM_File::ISA has a length of one;
}
use AnyDBM_File;                  
use vars qw( $DB_BTREE &R_DUP );  
use AnyDBM_File::Importer qw(:bdb);
use charnames qw(:full :short);

# lexical vars
my $infile;
my $outfile;
my $percent_id;
my $percent_cov;
my $cluster_size;
my $fas_file;
my $cores;
my $report;

#counters
#my $total_hits = 0;
#my $parsed_hits = 0;
#my $index = 0;

GetOptions(
           'i|infile=s'         => \$infile,
	   'f|fas_file=s'       => \$fas_file,
           'o|outfile=s'        => \$outfile,
	   'p|cpus=i'           => \$cores,
	   'r|report=s'         => \$report,
	   'cs|cluster_size=i'  => \$cluster_size,
	   'id|percent_id=f'    => \$percent_id,
	   'cov|percent_cov=f'  => \$percent_cov,
          );

# open the infile or die with a usage statement
if (!$infile || !$outfile || !$fas_file
    || !$percent_id || !$percent_cov || !$report) {
    print "\nERROR: Input not parsed correctly.\n";
    usage() and exit(1);
}

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

## Parse blast into form for clustering
my ($hitsort, $index_file, $hitsort_int) = parse_blast($infile, $percent_id, $percent_cov, $outfile);

### Clustering
my $community = louvain_method($hitsort_int, $index_file);
my $cls_file = make_clusters($community, $hitsort, $index_file);

# find union in clusters
$cluster_size = defined($cluster_size) ? $cluster_size : '500';
my ($seqs, $seqct) = fas2hash($fas_file);
my ($read_pairs, $vertex, $uf) = find_pairs($infile, $rep);
merge_clusters($vertex, $seqs, $read_pairs, $clsnew, $cls_dir_path, $rep);

### generate plots from cluster file
#my $cls_dir_path = clusters2fasta($cls_file, $seqhash, $cluster_size, $str);
#untie %$seqhash;
#undef %$seqhash;

#my $gl_dir = clusters2graphs($cls_file, $outfile, $cores, $str);
#GL2summary($gl_dir, $cls_file);

### search each fasta file in $cls_dir_path against custom repeatdb and HMMER
### -- then, incorporate results into size dist summary plot (see old 454summary.R
### script use in plant journal paper)
my $cwd = getcwd();
write_cls_size_dist_summary($cls_with_merges, $seqct, $cwd);

#
# Subs
#
sub parse_blast {
    my ($infile, $percent_id, $percent_cov, $outfile) = @_;

    # counters
    my $total_hits = 0;
    my $parsed_hits = 0;
    my $index = 0;

    my $index_file = $outfile.".index";  # integer index for read IDs used in clustering
    #my $cluster_file = $outfile.".cls";  # cluster file in Repeat Explorer "cls" format
    my $hitsort_int = $outfile.".int";   # integer index for clustering           

    open(my $hs_int, '>', $hitsort_int);
    open(my $indexout, '>', $index_file);
    open(my $out, '>', $outfile);

    my %match_index;
    $DB_BTREE->{cachesize} = 100000;
    $DB_BTREE->{flags} = R_DUP;
    tie( %match_index, 'AnyDBM_File', ':memory:', 0666, $DB_BTREE);

    open(my $in, '<', $infile);

    while (<$in>) { 
	chomp; 
	my ($q_name, $q_len, $q_start, $q_end, $s_name, $s_len, 
	    $s_start, $s_end, $pid, $score, $e_val, $strand) = split;
    
	if ($strand eq '2') {
	    $total_hits++;
	    my $neg_query_hit_length = ($q_start - $q_end) + 1;
	    my $neg_query_cov = $neg_query_hit_length/$q_len;
	    
	    if ( ($neg_query_cov >= $percent_cov) && ($pid >= $percent_id) ) {
		say $out join "\t", $q_name, $s_name, $score;
		$match_index{$q_name} = $index unless exists $match_index{$q_name}; # hash is much more memory efficient than array (1.6g v 2.3g)
		$index++;
		$match_index{$s_name} = $index unless exists $match_index{$s_name}; # and faster (though only ~1s faster)
		$index++;
	    }
	}
	else {
	    $total_hits++;
	    my $pos_query_hit_length = ($q_end - $q_start) + 1;
	    my $pos_query_cov = $pos_query_hit_length/$q_len;
	    
	    if ( ($pos_query_cov >= $percent_cov) && ($pid >= $percent_id) ) {
		say $out join "\t", $q_name, $s_name, $score;
		$match_index{$q_name} = $index unless exists $match_index{$q_name};
		$index++;
		$match_index{$s_name} = $index unless exists $match_index{$s_name};
		$index++;
	    }
	}
    }
    close($out);

    open(my $hs, '<', $outfile);
    while (my $hitpair = <$hs>) {
	chomp $hitpair;
	my ($q, $s, $sc) = split /\t/, $hitpair;
	if (exists $match_index{$q} && exists $match_index{$s}) {
	    say $hs_int join "\t", $match_index{$q}, $match_index{$s}, $sc;
	}
    }
    close($hs);
    close($hs_int);
    
    for my $idx_mem (sort { $match_index{$a} <=> $match_index{$b} } keys %match_index) {
	say $indexout join " ", $idx_mem, $match_index{$idx_mem};
    }
    close($indexout);

    return($outfile, $index_file, $hitsort_int);
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
    unless (defined $cls_bin && defined $cls_tree_weights) {
	say "ERROR: louvain_convert failed. Exiting." and exit(1);
    }
    system("louvain_community $cls_bin -l -1 -w $cls_tree_weights -v >$cls_tree 2>$cls_tree_log");
    unless (defined $cls_tree && defined $cls_tree_log) {
	say "ERROR: louvain_community failed. Exiting." and exit(1);
    }

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
    my ($graph_comm, $hitsort, $index_file) = @_;

    my $cluster_file = $hitsort.".cls";

    my @graph_comm_sort = reverse sort { ($a =~ /(\d)$/) <=> ($b =~ /(\d)$/) } @$graph_comm;
    my $graph = shift @graph_comm_sort;
    my %clus;
    my %index;

    open(my $idx, '<', $index_file);
    while (my $idpair = <$idx>) {
        chomp $idpair;
        my ($readid, $readindex) = split /\s+/, $idpair;
        $index{$readindex} = $readid;
    }
    close($idx);

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
            if (exists $index{$cls_member}) {
                push @clus_members, $index{$cls_member};
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

sub merge_clusters {
    my ($vertex, $seqs, $read_pairs, $clsnew, $cls_dir_path, $rep) = @_;

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

    #return \%cluster, \$clsnew;
}

sub fas2hash {
    my $fas_file = shift;
    open(my $fas, '<', $fas_file);

    my %seqhash;
    $DB_BTREE->{cachesize} = 100000;
    $DB_BTREE->{flags} = R_DUP;
    #my $seq_dbfile = "repeat_explorer_seqs.bdb";
    tie( %seqhash, 'AnyDBM_File', ':memory:', 0666, $DB_BTREE);
   
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
    close($fas);
    
    return(\%seqhash, $seqct);
}

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
            #if (scalar(@ids) >= $cluster_size) {
            push @{$read_pairs{$clsid}}, $_ for @ids;
            #}
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
		my $k = mk_key($cls_i, $cls_j);
                $cls_conn_ct{$k}++;
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

# the best way to make combined keys is not to concatenate them
# http://stackoverflow.com/a/15299397/1543853
#
# NB: in v5.16 charnames() will be loaded automatically
sub mk_key { join "\N{INVISIBLE SEPARATOR}", @_ }

sub mk_vec { split "\N{INVISIBLE SEPARATOR}", shift }

sub clusters2graphs {
    my ($cls_file, $hitsort, $cores, $str) = @_;

    my $clusters2graph_rscript = $cls_file;
    $clusters2graph_rscript .= ".clusters2graph.rscript";
    my ($hname, $hpath, $hsuffix) = fileparse($hitsort, qr/\.[^.]*/);
    my $gl_dir = $hname."_GL_$str";
    my $ncol_dir = $hname."_ncol_$str";
    open(my $rscript, '>', $clusters2graph_rscript);

    say $rscript "suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(multicore))

## FUNCTIONS :
# this is modification of igraph function  - faster for large graphs
indexing=function(x,ind){
        y=as.numeric(factor(ind,levels=x))
        y
}

graph.data.frame2 <- function(d, directed=TRUE) {
        
        if (ncol(d) < 2) {
                stop(\"the data frame should contain at least two columns\")
        }

        # assign vertex ids
        if (nrow(d)>200000000){
                cat(\'getting unique...\')
                names <- uniqueM(c(d[,1], d[,2]))
                cat(\'done\\n\')

        }else{
                names <- unique( c(as.character(d[,1]), as.character(d[,2])) )

        }
        
        ids <- seq(along=names)

        # create graph
        g <- graph.empty(n=0, directed=directed)
        g <- add.vertices(g, length(ids), name=names)

        # create edge list
        from <- as.character(d[,1])

        to <- as.character(d[,2])
        #edges <- t(matrix(c(ids[from], ids[to]), nc=2))
        edges <- t(matrix(c(indexing(names,from),indexing(names,to)),nc=2))
        
        # edge attributes
        attrs <- list()
        if (ncol(d) > 2) {
                for (i in 3:ncol(d)) {
                        newval <- d[,i]
                        if (class(newval) == \"factor\") {
                                newval <- as.character(newval)
                        }
                        attrs[[ names(d)[i] ]] <- newval
                }
        }
        
        # add the edges
        g <- add.edges(g, edges, attr=attrs)
        g
}

uniqueM=function(x,n=100)
{ print(length(x)  )
        options(warn=1)
        x=unlist(mclapply(split(x,1:4),FUN=unique))
        options(warn=0)  
        print(length(x))  
        x=unique(x)
        print(length(x)  )
        x
    }

#######################
registerDoMC(cores=$cores)

vmax=50000

hitsort=\"$hitsort\"
#clusterlist=\"$cls_file\"
minClsSize=30

cat(\'\\n reading .cls file\\n\')
clusters=scan(file=\"$cls_file\",what=character(),sep=\"\\n\",comment.char=\">\",quiet=T)  # load cluster file
clusters=strsplit(clusters,split=\"[ \\t]\")
clusters=clusters[sapply(clusters,length)>=minClsSize]  # remove smaller clusters
    ncolDir=paste(dirname(hitsort),\"/$ncol_dir\",sep=\'\')

##get hitsort
cat(\"\\nreading hitsort file\\n\\n\")
edgelist=read.table(hitsort,sep=\'\\t\',header=F,as.is=T,colClasses=c(\"character\",\"character\",\"numeric\"))
# create output set directories ncol
dir.create(ncolDir)

#save individual clusters: ncol
    for (i in seq_along(clusters)){
     cat(paste(\'saving cluster no.\',i,\'\\n\'))
     subEdgelist=edgelist[(edgelist[,1] %in% clusters[[i]]+edgelist[,2] %in% clusters[[i]])==2,]
     write.table(subEdgelist,sep=\'\\t\',col.names=F,row.names=F,quote=F,file=paste(ncolDir,\"/CL\",i,\'.ncol\',sep=\"\"))
 }                                                                                                                                                           
rm(edgelist)           
 
#clusters
# create output set directories GL
    GLdir=paste(dirname(hitsort),\"/$gl_dir\",sep=\'\')
dir.create(GLdir)

# calculate layouts
# read ncol files again:
    out=foreach(i = seq_along(clusters),.inorder=FALSE) %dopar% {
        dir.create(ncolDir)
        
        # samples are created in advance!!
        ncolfile=paste(ncolDir,\"/CL\",i,\'.ncol\',sep=\"\")
        ncolSampledFile=paste(ncolDir,\"/CL\",i,\'_sample.ncol\',sep=\"\")
        
        if (file.exists(ncolSampledFile)){
                cat(paste(\'original cluster CL\',i, \'was above threshold!, sample of graph is used\\n\'))
                gd=read.table(file=ncolSampledFile,sep=\'\\t\',header=F,as.is=T,col.names=c(1,2,\'weight\'))
                cat(paste(\'ncol file for cluster CL\',i,\' loaded - \',Sys.time(),\"\\n\",sep=\"\"))
                g=graph.data.frame2(gd,directed=F)
                cat(paste(\'ncol file converted to graph for cluster CL\',i,\" \",Sys.time(),\"\\n\",sep=\"\"))
                rm(gd) # to free memory!
                cat(paste(\'graph of  cluster CL\',i,\' with \',vcount(g),\' nodes, \',ecount(g),\'edges created\',sep=\'\'),\"\\n\")
                CCs=clusters(g)
                #the largest cc
                maxCC=which.max(CCs\$csize)
                g=induced.subgraph(g,vids=which(CCs\$membership==maxCC))
                cat(paste(\'Largest connected component of graph sample of cluster CL\',i,\' - \',vcount(g),\' nodes, \',ecount(g),\'edges created\',sep=\'\'),\"\\n\")
	    }else{
                gd=read.table(file=ncolfile,sep=\'\\t\',header=F,as.is=T,col.names=c(1,2,\'weight\'))
                cat(paste(\'ncol file for cluster CL\',i,\' loaded - \',Sys.time(),\"\\n\",sep=\"\"))
                g=graph.data.frame2(gd,directed=F)
                cat(paste(\'ncol file converted to graph for cluster CL\',i,\' \',Sys.time(),\"\\n\",sep=\"\"))
                rm(gd) # to free memory!
                cat(paste(\'graph of  cluster CL\',i,\' with \',vcount(g),\' nodes, \',ecount(g),\'edges created\',sep=\'\'),\"\\n\")
                
	    }

            cat(paste(\'layout for cluster CL\',i,\' - calculation start at - \',Sys.time(),\"\\n\",sep=\"\"))
            set.seed(1)
            GL=list(G=g,L=layout.fruchterman.reingold(g,dim=3,verbose=F))
            cat(paste(\'layout for cluster CL\',i,\' - calculation finished at - \',Sys.time(),\"\\n\",sep=\"\"))
            save(GL,file=paste(GLdir,\"\/CL\",i,\'.GL\',sep=\"\"))
            cat(paste(\'layout for cluster CL\',i,\' with \',vcount(g),\' nodes, \',ecount(g),\' edges saved at \',Sys.time(),sep=\'\'),\"\\n\")

    }";

close($rscript);
         
system("/usr/local/R/2.15.2/lib64/R/bin/R --vanilla --slave --silent < $clusters2graph_rscript 2> /dev/null");
unlink($clusters2graph_rscript);

return $gl_dir;

}

sub GL2summary {
    my ($gl_dir, $cls_file) = @_;

    unless ($gl_dir =~ /\/$/) {
        $gl_dir .= "/";
    }
    
    ## edit these filenames to something better
    my $gl2summary_plot = $cls_file;
    $gl2summary_plot .= ".cluster_graph_summary.pdf";
    my $gl2summary_rscript = $cls_file;
    $gl2summary_rscript .= "gl2summary.rscript";
    open(my $rscript, '>', $gl2summary_rscript);

    say $rscript "suppressPackageStartupMessages(library(igraph))
plotg=function(GG,LL,wlim=NULL,...){
        
        e=get.edgelist(GG,names=F)
        w=E(GG)\$weight
        if (!is.null(wlim)) {e=e[w>wlim,]; w=w[w>wlim]}
        X0=LL[e[,1],1]
        Y0=LL[e[,1],2]
        X1=LL[e[,2],1]
        Y1=LL[e[,2],2]
        plot(range(LL[,1]),range(LL[,2]),xlab=\"\",ylab=\"\",axes=F,type=\"n\",...)
        brv=\'grey\'
        segments(X0,Y0,X1,Y1,lwd=.5,col=brv)
        points(LL,pch=18,cex=.4,...)
}
############################################################################
Diameter=''
tmp=capture.output(as.null(Diameter)) # we don't want to calculate this because it takes too long
# get .GL file names
setwd(\"$gl_dir\")
GLfiles=system(\"ls *.GL\",intern=T)
# assume it is numbered CLXX.GL - get number:
GLfilesIndex=as.numeric(gsub(\"[^0-9]\",\"\",GLfiles))
GLfiles=GLfiles[order(GLfilesIndex)]
                        
j=0
page=1
figPerPage=12
ll=c(4,4)
pngFiles=list()
newPage=seq(1,length(GLfiles),figPerPage)
pageLayout=matrix(c(matrix(1:8,ncol=2,byrow=T),matrix(9:16,ncol=2,byrow=T),matrix(17:24,ncol=2,byrow=T)),ncol=4,byrow=T)

for (i in GLfiles){
        j=j+1
        if (j %in% newPage){
                graphics.off()
                pngFiles[[j]]=paste(\'page\',sprintf(\"%04d\", page),\".png\",sep=\'\')
                png(pngFiles[[j]],width=2481,height=3507,pointsize=40)
                layout(pageLayout,heights=c(2,1,2,1,2,1,2,1))
                par(mar=c(0.5,0.5,0.5,0.5))
                page=page+1
        }
        cat (i,\'...\')
        load(i)
        plotg(GL\$G,GL\$L)
        par(mar=c(0.0,0.0,0.0,0.0))
        plot(0:10,0:10,type=\'n\',axes=F,xlab=\'\',ylab=\'\',main=paste(\'\\nCL\',j,sep=\'\'))
        text2plot=paste(\'Number of reads:   \',vcount(GL\$G),\'\\nNumber of pairs:   \',ecount(GL\$G),
                        \"\\nDensity:      \",signif(graph.density(GL\$G),4),\"\\nDiameter:      \",ifelse(is.null(Diameter),\"NA\",length(get.diameter(GL\$G,directed=F))),
                        \"\\nMean edge weigth:   \", signif(mean(E(GL\$G)\$weight),5),\"\\nMax. degree:     \",max(degree(GL\$G)),sep=\'\')
        text(1,4,labels=text2plot,pos=4)
        abline(h=0)
        par(mar=c(0.5,0.5,0.5,0.5))
        cat(\"done\\n\")
}
tmp2=capture.output(dev.off())
# make pdf file from pngs 
#warnings()
cmd=paste(\"convert page????.png $gl2summary_plot\")
system(cmd)
tmp=lapply(pngFiles,unlink)  # remove all png files";

close($rscript);

    system("/usr/local/R/2.15.2/lib64/R/bin/R --vanilla --slave --silent < $gl2summary_rscript 2> /dev/null");
    unlink($gl2summary_rscript);
    #return $gl2summary_plot
}

sub write_cls_size_dist_summary {
    my ($cls_file, $seqct, $cwd) = @_;

    my $cls_size_dist_plot = $cls_file;
    $cls_size_dist_plot .= ".png";
    my $cls_size_dist_rscript = $cls_file;
    $cls_size_dist_rscript .= ".rscript";
    open(my $rscript, '>', $cls_size_dist_rscript);

    say $rscript "cls=scan(file=\"$cls_file\",what=character(),sep=\"\\n\",comment.char=\">\",quiet=T)
                  cls=strsplit(cls,split=\"[ \\t]\")

png(filename=\"$cls_size_dist_plot\",width=800,height=500)

options(scipen=999) # this turns of printing in scientific notation
#n=${seqct}
#formatC(n, format = \"d\")
#format(n, scientific = FALSE)
#setwd(${cwd})     #################### if this is necessary set it to the cwd
# plot barplot:
NinClusters=length(unlist(cls))
#NinAll=length(seqs)
NinAll=$seqct
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))
barplot(clsLengthAll,width=clsLengthAll,space=0,ylim=c(0,max(clsLength)*1.2),ylab=\"number of reads [reads]\",xlab=\"number of reads [%]\",main=paste(NinAll, \"reads total\"))
rect(0,0,NinClusters,clsLength[[1]]*1.2,col=\"#FF000010\")
rect(NinClusters,0,NinAll,clsLength[[1]]*1.2,col=\"#00FF0010\")
text(NinClusters/2,clsLength[[1]]*1.05, labels=paste(NinClusters,\"reads in\\n\",length(cls),\"clusters\"))
text(NinClusters+NinSingles/2,clsLength[[1]]*1.05, labels=paste(NinSingles,\"singlets\"))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10))

tmp=capture.output(dev.off())";
    close($rscript);

    system("/usr/local/R/2.15.2/lib64/R/bin/R --vanilla --slave --silent < $cls_size_dist_rscript 2> /dev/null");
    unlink($cls_size_dist_rscript);
}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -i inreport -o hitsort -id 90.00 -cov 0.55 -f fasta_file -p 4

Required:
    -i|infile        :    An mgblast report in tab format (-D 4).
    -f|fas_file      :    The (Fasta) file of sequences used in the all-vs-all blast.
    -o|outfile       :    File name to write the parsed results to in Repeat Explorer\'s hitsort format.
    -p|cpus          :    The number of processors to use for plotting.
    -cs|cluster_size :    The minimum cluster size to convert to fasta (Default: 500).
    -id|percent_id   :    The percent identity threshold for matches.
    -cov|percent_cov :    The percent coverage for both the query and subject to be retained for clustering.

Options:
    -r|report        :    A file to hold the cluster stats. [NOT IMPLEMENTED]

END
}
