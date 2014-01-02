#!/usr/bin/env perl

##TODO:

use 5.014;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dump qw(dd);
use Bio::Kseq;
use Capture::Tiny qw(:all);
BEGIN {
  @AnyDBM_File::ISA = qw( DB_File SQLite_File )
      unless @AnyDBM_File::ISA == 1; # if loaded already, AnyDBM_File::ISA has a length of one;
}
use AnyDBM_File;
use vars qw( $DB_BTREE &R_DUP );
use AnyDBM_File::Importer qw(:bdb);

## add warnings for given/when

#
# lexical vars
#
my $infile;
my $outfile;
my $pep_fas;
my $nt_fas;
my $help;
my $man;

GetOptions(
           'i|infile=s'    => \$infile,
           'pf|pep_fas=s'  => \$pep_fas,
           'nf|nt_fas=s'   => \$nt_fas,    
           'o|outfile=s'   => \$outfile,
           'h|help'        => \$help,
           'm|man'         => \$man,
          );

#
# check input
#
if (!$infile || !$nt_fas || 
    !$pep_fas || !$outfile) {
    usage();
    exit(1);
}

my %seqhash;
my %statshash;
my @stats;
my @treestats;

$DB_BTREE->{cachesize} = 100000;
$DB_BTREE->{flags} = R_DUP;
my $db_file = "orthoMCL_groups.bdb";
#tie( %seqhash, 'AnyDBM_File', ':memory:', 0666, $DB_BTREE);
tie( %seqhash, 'AnyDBM_File', $db_file, 0666, $DB_BTREE);

# set PATH for programs we need
my $muscle = find_prog("muscle");
my $gblocks = find_prog("Gblocks");
my $pal2nal = find_prog("pal2nal");
my $raxml = find_prog("raxml");

my $knseq = Bio::Kseq->new($nt_fas); ## TODO: Check if file exists
my $nt_it = $knseq->iterator;

my $kpseq = Bio::Kseq->new($pep_fas);
my $pep_it = $kpseq->iterator;

while (my $nseq = $nt_it->next_seq) {
    $seqhash{ $nseq->{name} } = [ $nseq->{seq} ];
}

while (my $pseq = $pep_it->next_seq) {
    if (exists $seqhash{ $pseq->{name} }) {
	push @{$seqhash{ $pseq->{name} }}, $pseq->{seq};
    }
}

my $clusters = parse_groups($infile, %seqhash);

for my $cluster (sort keys %$clusters) {
    my $cluster_size = scalar @{$clusters->{$cluster}};
    push @stats, $cluster_size;
    if ($cluster_size  > 10 && $cluster_size < 12) { # for testing
	my @species = map { s/\|.*//r } @{$clusters->{$cluster}};  # need Perl 5.14+ to do non-destructive substitution
	my %ortho_group; 
	for my $gene (@{$clusters->{$cluster}}) {
	    my $species = $gene;
	    $species =~ s/\|.*//;
	    if (exists $ortho_group{$cluster}{$species}) {
		push @{$ortho_group{$cluster}{$species}}, $gene;
	    }
	    else {
		$ortho_group{$cluster}{$species} = [ $gene ];
	    }
	}
        	
        for my $group (keys %ortho_group) {
	    for my $taxa (keys %{$ortho_group{$group}}) {
		
		my $gene_file = "gene_cluster_".$cluster."_".$taxa."_nt.fasta";
		my $pep_file = "gene_cluster_".$cluster."_".$taxa."_pep.fasta";
		open my $nt_group, ">>", $gene_file or die "\nERROR: Could not open file: $gene_file\n";
		open my $pep_group, ">>", $pep_file or die "\nERROR: Could not open file: $pep_file\n";
	
		my $gene_ct = 0;
		for my $paralog (@{$ortho_group{$group}{$taxa}}) {
		    if (exists $seqhash{$paralog}) {
			$gene_ct++;
			say $nt_group join "\n", ">".$paralog, ${$seqhash{$paralog}}[0];
			say $pep_group join "\n", ">".$paralog, ${$seqhash{$paralog}}[1];
		    }
		}
		close $nt_group;
		close $pep_group;
    
		if (-s $gene_file && -s $pep_file && $gene_ct > 3) {
		    my ($gene_aln, $pep_aln) = align($gene_file, $pep_file, $muscle);
		    
                    ###### refine alignment 
		    my ($gene_phy, $gene_aln_fas) = pal2nal($pep_aln, $gene_file, $pal2nal);
		    #my $gblocks_aln = gblocks($gene_aln, $gblocks); # this works, but I'm not sure what to do the format yet
	
		    ###### infer tree
		    if (-s $gene_phy) {
			my $gene_tree = infer_tree($gene_phy, $raxml);
			#Print "Getting to here -> line 129.\n";
			#print "Gene tree on line 131 is: ", $gene_tree,"\n";
			my $treehash = summarize_trees($taxa, $gene_tree, $gene_ct);
			#print "Branch length for $gene_tree is: ", $total_length,"\n";
			push @treestats, $treehash;
			###### compute ka/ks, etc. on alignment
			#infer_substitution_patterns($gene_aln_fas);     # need to look at pal2nal aln
		    }
		}
		else {
		    unlink $gene_file, $pep_file;
		}
		$gene_ct = 0;
	    }
	}
    }
}

undef %seqhash;
untie %seqhash;

my $count = scalar @stats;
my $mean = mean(@stats);
my $median = median(@stats);
my $min = min(@stats);
my $max = max(@stats);

open my $report, ">", $outfile or die "\nERROR: Could not open file: $outfile\n";
say $report "=-=" x 25;
say $report join "\t", "Cluster_count","Cluster_mean","Cluster_median","Cluster_min","Cluster_max";
say $report join "\t", $count, $mean, $median, $min, $max;
say $report "=-=" x 25;
say $report join "\t", "Cluster_number","Cluster_size";

for my $key (reverse sort { $statshash{$a} <=> $statshash{$b} } keys %statshash) {
    say $report join "\t", $key,$statshash{$key};
}
close $report;

say join "\t","Species", "Tree_file", "Tree_stat";
for my $treehashref (@treestats) {
    for my $taxa (keys %$treehashref) {
	dd $taxa;
	#for my $tree (keys %{$treehashref{$taxa}}) {
	    #for my $stat (@{$treehashref{$taxa}{$tree}}) {
		#print join("\t",($taxa, $tree, $stat)), "\n";
	    #}
	#}
    }
}

exit;
#
# subs
#
sub parse_groups {
    # modified from:
    # http://cpansearch.perl.org/src/EASR/ONTO-PERL-1.19/lib/OBO/CCO/OrthoMCLParser.pm

    my $infile = shift;	
    open my $fh, '<', $infile or die "\nERROR: Could not open file: $infile\n.";
    
    my %clusters; 
    while (<$fh>){
	my ($cluster, $proteins) = split /:\s+/xms;
	my $cluster_num;
	if ($cluster =~ /(\d+)/) { # work on this regex
	    $cluster_num = $1;
	}
	$cluster = $cluster_num;
	my @proteins = split /\s/xms, $proteins;
	my $protein_ct = @proteins;
	for my $protein ( @proteins) {
	    if ($protein  =~ /((\w+)\|(\w+))/) {
		if (exists $clusters{$cluster}) {
		    push @{$clusters{$cluster}}, $1;
		}
		else {
		    $clusters{$cluster} = [ $1 ];
		}
	    }		
	}
    }		
    close $fh;
    return \%clusters;
}

sub infer_substitution_patterns {
    my $gene_aln_fas = shift;

    eval { require Bio::Tools::Run::Phylo::PAML::Codeml;
	   require Bio::Tools::Run::Alignment::Muscle;
	};
    if ($@) {
	die "\nERROR: The BioPerl-run package is required to work with alignments. Exiting.\n";
    }
    
    use Bio::Align::Utilities qw(aa_to_dna_aln);

    my $aln_factory = Bio::Tools::Run::Alignment::Muscle->new;
    my $kaks_factory = Bio::Tools::Run::Phylo::PAML::Codeml->new( -params => { 'runmode' => -2,
									       'seqtype' => 1,
								              } );
    
    my $seqio = Bio::SeqIO->new(-file => $gene_aln_fas, -format => 'fasta');

    my %seqs; 
    my @prots;
    while (my $seq = $seqio->next_seq) {
	$seqs{$seq->display_id} = $seq;
	# translate them into protein
	my $protein = $seq->translate();
	my $pseq = $protein->seq();
	if ( $pseq =~ /\*/ &&
	     $pseq !~ /\*$/ ) {
	    warn("provided a CDS sequence with a stop codon, PAML will choke!");
	    exit(0);
	}
	# Tcoffee can't handle '*' even if it is trailing
	$pseq =~ s/\*//g;
	$protein->seq($pseq);
	push @prots, $protein;
    }
 
    if ( @prots < 2 ) {
	warn("Need at least 2 CDS sequences to proceed");
	exit(0);
    }

    my $aa_aln = $aln_factory->align(\@prots);

    my $dna_aln = aa_to_dna_aln($aa_aln, \%seqs);

    # set the alignment object
    $kaks_factory->alignment($dna_aln);

    my @each = $dna_aln->each_seq();

    # run the KaKs analysis
    my ($rc,$parser) = $kaks_factory->run();
    my $result = $parser->next_result;
    my $MLmatrix = $result->get_MLmatrix();
    
    my @otus = $result->get_seqs();
    # this gives us a mapping from the PAML order of sequences back to
    # the input order (since names get truncated)
    my @pos = map {
	my $c = 1;
	for my $s ( @each ) {
	    last if $s->display_id eq $_->display_id;
	    $c++;
	}
	$c;
    } @otus;
 
    say join "\t", qw(SEQ1 SEQ2 Ka Ks Ka/Ks PROT_PERCENTID CDNA_PERCENTID);
    for my $i ( 0 .. $#otus -1 ) {
	for my $j ( $i+1 .. $#otus ) {
	    my $sub_aa_aln  = $aa_aln->select_noncont($pos[$i],$pos[$j]);
	    my $sub_dna_aln = $dna_aln->select_noncont($pos[$i],$pos[$j]);
	    say join "\t", $otus[$i]->display_id,
		       $otus[$j]->display_id,$MLmatrix->[$i]->[$j]->{'dN'},
		       $MLmatrix->[$i]->[$j]->{'dS'},
		       $MLmatrix->[$i]->[$j]->{'omega'},
		       sprintf("%.2f",$sub_aa_aln->percentage_identity),
	               sprintf("%.2f",$sub_dna_aln->percentage_identity);
	}
    }
    # return \%somehash;
}

sub summarize_trees {
    my ($taxa, $gene_tree, $gene_ct) = @_;
    eval { require Bio::Phylo::IO; };
    if ($@) {
	die "\nERROR: The Bio::Phylo Perl package is required to work with trees. Exiting.\n";
    }

    my %treehash;
    my $forest = Bio::Phylo::IO->parse(
	-format => 'newick',
	-file   => $gene_tree
	);

    my $tree = $forest->first;
    my $tree_length = $tree->calc_tree_length;
    my $path_length = $tree->calc_total_paths;
    #my $ltt = $tree->calc_ltt;  # tree must be ultrameric 
    my $node_num = $tree->calc_number_of_nodes;

    my $ave_branch_len = $tree_length / $node_num;
    my $branch_len_by_genes = $tree_length / $gene_ct;
    say "Total tree length for $gene_tree is: ",$tree_length;
    say "Total path length for $gene_tree is: ",$path_length;
    #say "LTT points are: ",$ltt;
    say "Total number of sequences for $gene_tree is: ",$gene_ct;
    say "Total node number for $gene_tree is: ",$node_num;
    say "Ave branch length for $gene_tree is: ",$ave_branch_len;
    say "Branch length divided by number of genes is: ",$branch_len_by_genes,"\n";

    for my $stat (qw($ave_branch_len $tree_length $path_length $node_num $branch_len_by_genes)) {
	push @{$treehash{$taxa}{$gene_tree}}, $stat;
    }

    return \%treehash;

}

sub infer_tree {
    my ($gene_phy, $raxml) = @_;
    my $raxml_tree = $gene_phy;
    $raxml_tree =~ s/\.phy$/\_raxml\.tre/;
    my $raxml_out_tree = "RAxML_bootstrap.".$raxml_tree;
    ## TODO: Add better method for executing commands
    my ($raxml_out, $raxml_err, @raxml_res) = capture { 
	system("$raxml -A S16 -b 4294382 -p 34623 -B 0.03 -f d -k -K GTR -m GTRGAMMAI -T 8 -U -N 100 -s $gene_phy -n $raxml_tree 2>&1 /dev/null"); 
    };

    print "RAxML error: ",$raxml_out,"\n";

    return $raxml_out_tree;
}

sub gblocks {
    my ($gene_aln, $gblocks) = @_;
    my $gblocks_aln .= $gene_aln."-gb";
    ## TODO: Add better method for executing commands
    my ($gblocks_out, $gblocks_err, @gblocks_res) = capture { system("$gblocks $gene_aln -t=d -s=y -p=y -b4=5 -e=-gb"); };
    return $gblocks_aln;
}

sub pal2nal {
    my ($pep_aln, $gene_file, $pal2nal) = @_;
    my $pal2nal_aln = $gene_file;
    $pal2nal_aln =~ s/\.fa.*/\_pal2nal\.aln/;
    my ($pal2nal_out, $pal2nal_err, @pal2nal_res) = capture { system("$pal2nal $pep_aln $gene_file -nogap > $pal2nal_aln"); };

    my $clustalw_to_fasta = sub {
	eval { require Bio::AlignIO; };
	if ($@) {
	    die "BioPerl is needed to do alignment format conversion. Exiting.\n";
	}
	my $clustalw_aln = shift;
	my $fas_aln = $clustalw_aln;
	$fas_aln =~ s/\.aln$/\.fas/;
	my $aln_in = Bio::AlignIO->new(-file => $clustalw_aln, -format => 'clustalw');
	my $aln_out = Bio::AlignIO->new(-file => ">$fas_aln", -format => 'fasta');
	
	while (my $aln = $aln_in->next_aln ) {
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
	
	while (my $aln = $aln_in->next_aln ) {
	    $aln_out->write_aln($aln);
	}
	
	# this is so RaxML won't complain about pal2nal's clustalw format (including quotes in the identifiers)
	open my $phy, "<", $phy_outtmp or die "\nERROR: Could not open file: $phy_outtmp\n";
	open my $correct_phy, ">", $phy_out or die "\nERROR: Could not open file: $phy_out\n";
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
		    #print "$seqname $seq\n";
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
	close $phy;
	close $correct_phy;
	unlink $phy_outtmp;
	
	return $phy_out;
    };

    my $gene_phy = $clustalw_to_phylip->($pal2nal_aln);
    my $gene_aln_fas = $clustalw_to_fasta->($pal2nal_aln);

    return($gene_phy, $gene_aln_fas); 
}

sub align {
    my ($gene_file, $pep_file, $muscle) = @_;
    my $gene_aln = $gene_file;
    $gene_aln =~ s/\.fa.*//;
    $gene_aln .= ".aln";
    my $pep_aln = $pep_file;
    $pep_aln =~ s/\.fa.*//;
    $pep_aln .= ".aln";

    ## TODO: Add better method for executing commands
    my ($gene_aln_out, $gene_aln_err, @gene_aln_res) = capture { system("$muscle -in $gene_file -out $gene_aln -quiet"); };
    my ($pep_aln_out, $pep_aln_err, @pep_aln_res) = capture { system("$muscle -in $pep_file -out $pep_aln -quiet"); };
    
    return($gene_aln, $pep_aln);
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
    my $sum; 
    my $count = scalar @array; 
    foreach (@array) { $sum += $_; } 
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
USAGE: $script [-i] [-nf] [-pf] [-o] 

Required:
     -i|infile         :       The orthoMCL groups.txt file produced by the program mcl.
     -nf|nt_fas        :       The goodGenes.fasta file (you must create this).
     -pf|pep_fas       :       The goodProteins.fasta file.
     -o|outfile        :       A report of the length statistics for each cluster.

Options:
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}

