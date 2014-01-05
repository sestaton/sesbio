#!/usr/bin/env perl

##TODO: Create a goodGenes.fasta file in the same format as goodProteins.fasta,
##      with the former being the nucleotide sequence of the genes. Use blast2orthologs.pl
##      as a guide for creating that file (or data structure). Use the nucleotide/peptide
##      alignments for calculating evolutionary rates. Need to trim alignments...
##      suppress warnings for given/when if v5.18+

use 5.010;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dump;
use Bio::SeqIO;
use Capture::Tiny qw(:all);
BEGIN {
  @AnyDBM_File::ISA = qw( DB_File SQLite_File )
      unless @AnyDBM_File::ISA == 1; # if loaded already, AnyDBM_File::ISA has a length of one;
}
use AnyDBM_File;
use vars qw( $DB_BTREE &R_DUP );
use AnyDBM_File::Importer qw(:bdb);

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

$DB_BTREE->{cachesize} = 100000;
$DB_BTREE->{flags} = R_DUP;
my $db_file = "orthoMCL_groups.bdb";
#tie( %seqhash, 'AnyDBM_File', ':memory:', 0666, $DB_BTREE);
tie( %seqhash, 'AnyDBM_File', $db_file, 0666, $DB_BTREE);

# set PATH for programs we need
my $muscle = find_prog("muscle");
my $pal2nal = find_prog("pal2nal");
my $raxml = find_prog("raxml");

my $nt_seq_in = Bio::SeqIO->new(-file => $nt_fas, -format => 'fasta');
my $pep_seq_in = Bio::SeqIO->new(-file => $pep_fas, -format => 'fasta');

while (my $nt_seq = $nt_seq_in->next_seq()) {
    my $seqname = $nt_seq->id;
    my $seq = $nt_seq->seq;
    $seqhash{$seqname} = [ $seq ];
}

while (my $pep_seq = $pep_seq_in->next_seq()) {
    my $pepname = $pep_seq->id;
    my $pepseq = $pep_seq->seq;
    if (exists $seqhash{$pepname}) {
	push @{$seqhash{$pepname}}, $pepseq;
    }
}

#dd %seqhash;

my $clusters = parse_groups($infile, %seqhash);

for my $cluster (sort keys %$clusters) {
    my $cluster_size = scalar @{$clusters->{$cluster}};
    push @stats, $cluster_size;
    if ($cluster_size  > 1 && $cluster_size < 20) { # for testing
	my $gene_file = "gene_cluster_".$cluster."_nt.fasta";
        my $pep_file = "gene_cluster_".$cluster."_pep.fasta";
        open my $nt_group, ">>", $gene_file or die "\nERROR: Could not open file: $gene_file\n";
        open my $pep_group, ">>", $pep_file or die "\nERROR: Could not open file: $pep_file\n";
        
        for my $gene (@{$clusters->{$cluster}}) {
            if (exists $seqhash{$gene} ) { 
                say $nt_group join "\n", ">".$gene,${$seqhash{$gene}}[0];
                say $pep_group join "\n", ">".$gene,${$seqhash{$gene}}[1];
             }
        }
        close $nt_group;
        close $pep_group;

	###### Align $gene_file here
	if (-s $gene_file && -s $pep_file) {
	    my ($gene_aln, $pep_aln) = align($gene_file, $pep_file, $muscle);
	    ###### construct tree and/or...
	    #my ($gene_tree, $pep_tree) = infer_tree($gene_aln, $pep_aln);
	    ###### compute ka/ks, etc. on alignment
	    my $pal2nal_aln = pal2nal($pep_aln, $gene_file, $pal2nal);
	}
	else {
	    unlink $gene_file, $pep_file;
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

sub infer_tree {
    my ($gene_aln, $raxml) = @_;
}

sub pal2nal {
    my ($pep_aln, $gene_file, $pal2nal) = @_;
    my $pal2nal_aln = $gene_file;
    $pal2nal_aln =~ s/\.fa.*/\_pal2nal\.aln/;
    my ($pal2nal_out, $pal2nal_err, @pal2nal_res) = capture { system("$pal2nal $pep_aln $gene_file -nogap > $pal2nal_aln"); };
    return $pal2nal_aln; 
}

sub align {
    my ($gene_file, $pep_file, $muscle) = @_;
    my $gene_aln = $gene_file;
    $gene_aln =~ s/\.fa.*//;
    $gene_aln .= ".aln";
    my $pep_aln = $pep_file;
    $pep_aln =~ s/\.fa.*//;
    $pep_aln .= ".aln";

    my ($gene_aln_out, $gene_aln_err, @gene_aln_res) = capture { system("$muscle -in $gene_file -out $gene_aln -quiet"); };
    my ($pep_aln_out, $pep_aln_err, @pep_aln_res) = capture { system("$muscle -in $pep_file -out $pep_aln -quiet"); };
    return ($gene_aln, $pep_aln);
}

sub min {
    my $min = shift;
    for ( @_ ) { $min = $_ if $_ < $min }
    return $min;
}

sub max {
    my $max = shift;
    for ( @_ ) { $max = $_ if $_ > $max }
    return $max;
}

sub mean { 
    my @array = @_; 
    my $sum; 
    my $count = scalar @array; 
    for (@array) { $sum += $_; } 
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
    chomp $path;

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

