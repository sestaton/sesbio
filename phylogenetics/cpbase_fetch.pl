#!/usr/bin/env perl

##TODO: Check if it is possible to get alignments. Since they are stored in the document path of /tmp
##      they may be generated on the fly (hence, not available on demand). It appears that are not available.
##      Need to confirm this an remove the alignments option.
##
##      Generate statistics from higher taxonomic levels

use 5.010;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Time::HiRes qw(gettimeofday);
use IPC::System::Simple qw(system);
use Try::Tiny;
use Pod::Usage;
use WWW::Mechanize;
use HTML::TableExtract;
use LWP::UserAgent;
use XML::LibXML;
use Data::Dump qw(dd);

# given/when emits warnings in v5.18+
no if $] >= 5.018, 'warnings', "experimental::smartmatch";

#
# lexical vars
#
my $db;
my $all;
my $type;
my $genus;
my $species;
my $outfile; ## log
my $statistics;
my $available;
my $help;
my $man;
my $alphabet;
my $assemblies;
my $alignments;
my $sequences;
my $gene_name;
my $gene_clusters;
my $rna_clusters;
my $lineage;
my $cpbase_response = "CpBase_database_response.html"; # HTML

GetOptions(
           'all'                     => \$all,
           'd|db=s'                  => \$db,
	   't|type=s'                => \$type,
           'g|genus=s'               => \$genus,
	   's|species=s'             => \$species,
	   'o|outfile=s'             => \$outfile,
	   'stats|statistics'        => \$statistics,
	   'available'               => \$available,
	   'mol|alphabet=s'          => \$alphabet,
	   'asm|assemblies'          => \$assemblies,
	   'aln|alignments'          => \$alignments,
	   'seq|sequences'           => \$sequences, 
	   'gn|gene_name=s'          => \$gene_name,
	   'gc|gene_clusters'        => \$gene_clusters,
	   'rc|rna_clusters'         => \$rna_clusters,
           'l|taxonomic_lineage'     => \$lineage,
	   'h|help'                  => \$help,
	   'm|man'                   => \$man,
    );

#pod2usage( -verbose => 1 ) if $help;
#pod2usage( -verbose => 2 ) if $man;

usage() and exit(0) if $help;

## set defaults for search
$type //= 'fasta';
$alphabet //= 'dna';

get_lineage_for_taxon($cpbase_response, $available) and exit(0) if $lineage;

if ($statistics) {
    if ($db || $gene_clusters || $rna_clusters || $assemblies || $alignments || $sequences) {
	say "\nERROR: The 'stats' option only works when given a genus and species name and no other options.";
	usage();
	exit(1);
    }
}

##This is a reminder for the TODO below; just say this option isn't implemented and exit
if ($alignments && !$gene_name || !$gene_clusters) {
    say "ERROR: Currently only alignments for a specific gene cluster or gene name ".
	"may be returned. Exiting.";
    exit(1);
}

##TODO add handling of option for getting all alignments
if ($gene_clusters && $gene_name && $alignments) {
    my $gene_stats = fetch_ortholog_sets($gene_name, $alignments, $alphabet, $type);
    say join "\t", "Gene","Genome","Locus","Product";
    for my $gene (keys %$gene_stats) {
	for my $genome (keys %{$gene_stats->{$gene}}) {
	    for my $locus (keys %{$gene_stats->{$gene}{$genome}}) {
		say join "\t", $gene, $genome, $locus, $gene_stats->{$gene}{$genome}{$locus};
	    }
	}
    }
    exit;
}

if ($rna_clusters && $gene_name) {
    fetch_rna_clusters($statistics, $sequences, $gene_name, $all);
    exit;
}

#
# check @ARGV
#
if (!$db) {
    say "\nERROR: A database to query must be given for getting stats or assemblies. Or, other arguments may be supplied. Exiting.";
    usage();
    exit(1);
}

if (!$genus && $species) {
    say "\nERROR: Can not query a species without a genus. Exiting.";
    usage();
    exit(1);
}

# make epithet
my $epithet;
$epithet = $genus."_".$species if $genus && $species;
my %stats;

#
# counters
#
my $t0 = gettimeofday();
my $records = 0;
my $genomes;

#
# Set the CpBase database to search and type
#
given ($db) {
    when (/algae/i) {             $db = "Algae"; }
    when (/red lineage/i) {       $db = "Red_Lineage"; }
    when (/rhodophyta/i) {        $db = "Rhodophyta"; }
    when (/stramenopiles/i) {     $db = "stramenopiles"; }
    when (/non_viridiplantae/i) { $db = "NOT_Viridiplantae"; }
    when (/viridiplantae/i) {     $db = "Viridiplantae"; }
    default {                     die "Invalid name for option db."; }
}

my $ua = LWP::UserAgent->new;
my $urlbase = "http://chloroplast.ocean.washington.edu/tools/cpbase/run&genome_taxonomy=$db";
my $response = $ua->get($urlbase);

# check for a response
unless ($response->is_success) {
    die "Can't get url $urlbase -- ", $response->status_line;
}

open my $out, '>', $cpbase_response or die "\nERROR: Could not open file: $!\n";
say $out $response->content;
close $out;

my $id_map = get_species_id($urlbase);

my $te = HTML::TableExtract->new( attribs => { border => 1 } );
$te->parse_file($cpbase_response);

for my $ts ($te->tables) {
    for my $row ($ts->rows) {
	my @elem = grep { defined } @$row;
	if ($elem[0] =~ /(\d+) Genomes/i) {
	    $genomes = $1;
	    unlink $cpbase_response;
	    say "$genomes $db genomes available in CpBase." and exit if $available;
	}
	else {
	    say "\nERROR: Be advised, this command will now attempt to download $genomes assemblies. ".
		"It would be nicer to specify one species. Exiting now.\n" and exit if $assemblies && $all;
	    my ($organism,$locus,$sequence_length,$assembled,$annotated,$added) = @elem;
	    $organism =~ s/\s+/_/g;
		    
	    if (exists $id_map->{$organism}) {
		if ($genus && $species && $organism =~ /\Q$epithet\E/i) {
		    my $id = $id_map->{$organism};
		    my $assem_stats = get_cp_data($id);
		    $stats{$organism} = $assem_stats;
		    my $file = $organism."_".$locus;
		    my $endpoint = "http://chloroplast.ocean.washington.edu/CpBase_data/$locus/files/$file";

		    if ($type eq 'genbank') {
			$file = $file.".gb";
			$endpoint = $endpoint.".gb";
		    }
		    elsif ($type eq 'fasta') {
			$file = $file.".fasta";
			$endpoint = $endpoint.".fasta";
		    }
		    fetch_file($file, $endpoint) if $assemblies;
		}
		elsif ($genus && $organism =~ /\Q$genus\E/i) {
		    my $id = $id_map->{$organism};
		    my $assem_stats = get_cp_data($id);
		    $stats{$organism} = $assem_stats;
		}
	    }
	}
    }
}

if ($statistics) {
    for my $genome (keys %stats) {
	say "====> Showing chloroplast genome statistics for: $genome";
	for my $stat (keys %{$stats{$genome}}) {
	    say join "\t", $stat, $stats{$genome}{$stat};
	}
    }
}
unlink $cpbase_response;

#
# methods
#
sub get_species_id {
    my ($urlbase) = @_;

    my %id_map;

    my $mech = WWW::Mechanize->new();
    $mech->get( $urlbase );
    my @links = $mech->links();
    for my $link ( @links ) {
	next unless defined $link->text;
	my ($g, $sp) = split /\s+/, $link->text;
	next unless defined $g && defined $sp;
	my $ep = $g."_".$sp;
	if ($link->url =~ /id=(\d+)/) {
	    $id_map{$ep} = $1;
	}
    }
    return \%id_map;
}

sub get_cp_data {
    my ($id) = @_;
    
    my %assem_stats;
    my $ua = LWP::UserAgent->new;
    my $cpbase_response = "CpBase_database_response_$id".".html";
    my $urlbase = "http://chloroplast.ocean.washington.edu/tools/cpbase/run?genome_id=$id&view=genome";
    my $response = $ua->get($urlbase);

    unless ($response->is_success) {
	die "Can't get url $urlbase -- ", $response->status_line;
    }

    open my $out, '>', $cpbase_response or die "\nERROR: Could not open file: $!\n";
    say $out $response->content;
    close $out;

    my $te = HTML::TableExtract->new( attribs => { border => 1 } );
    $te->parse_file($cpbase_response);
    
    for my $ts ($te->tables) {
	for my $row ($ts->rows) {
	    my @elem = grep { defined } @$row;
	    if ($elem[0] =~ /GC Content/) {
		$elem[1] =~ s/%$//;
		$assem_stats{gc_content} = $elem[1];
	    }
	    elsif ($elem[0] =~ /GC Skew/) {
		$assem_stats{gc_skew} = $elem[1];
	    }
	    elsif ($elem[0] =~ /Total Sequence Length/) {
		$elem[1] =~ s/\s.*//;
		$assem_stats{seq_len_total} = $elem[1];
	    }
	    elsif ($elem[0] =~ /Total CDS Bases/) {
		$elem[1] =~ s/\s.*//;
		$assem_stats{cds_bases_total} = $elem[1];
	    }
	    elsif ($elem[0] =~ /Average CDS Length/) {
		$elem[1] =~ s/\s.*//;
		$assem_stats{cds_len_ave} = $elem[1];
	    }
	    elsif ($elem[0] =~ /Total RNA Bases/) {
		$elem[1] =~ s/\s.*//;
		$assem_stats{rna_bases_total} = $elem[1];
	    }
	    elsif ($elem[0] =~ /Total RNA Bases/) {
		$elem[1] =~ s/\s.*//;
		$assem_stats{rna_bases_total} = $elem[1];
	    }
	    elsif ($elem[0] =~ /Average Repeat Length/) {
                $elem[1] =~ s/\s.*//;
                $assem_stats{repeat_bases_total} = $elem[1];
            }
	    elsif ($elem[0] =~ /Average Intergenic Distance/) {
                $elem[1] =~ s/\s.*//;
                $assem_stats{intergenic_dist_ave} = $elem[1];
            }
	}
    }
    unlink $cpbase_response;
    return \%assem_stats;
}

sub fetch_ortholog_sets {
    my ($gene_name, $alignments, $alphabet, $type) = @_;
    my %gene_stats;
    my $mech = WWW::Mechanize->new;
    my $urlbase = "http://chloroplast.ocean.washington.edu/tools/cpbase/run?view=u_feature_index"; 
    $mech->get($urlbase);
    my @links = $mech->links();

    for my $link ( @links ) {
        next unless defined $link && $link->url =~ /tools/;
	if ($link->url =~ /u_feature_id=(\d+)/) {
	    my $id = $1;
	    my $ua = LWP::UserAgent->new;
	    my $cpbase_response = "CpBase_database_response_gene_clusters_$id".".html";
	        
	    my $urlbase = "http://chloroplast.ocean.washington.edu/tools/cpbase/run?u_feature_id=$id&view=universal_feature"; 
	    my $response = $ua->get($urlbase);
	        
	    unless ($response->is_success) {
		die "Can't get url $urlbase -- ", $response->status_line;
	    }
	        
	    open my $out, '>', $cpbase_response or die "\nERROR: Could not open file: $!\n";
	    say $out $response->content;
	    close $out;
	        
	    my $te = HTML::TableExtract->new( attribs => { border => 1 } );
	    $te->parse_file($cpbase_response);
	        
	    for my $ts ($te->tables) {
		for my $row ($ts->rows) {
		    my @elem = grep { defined } @$row;
		    if (defined $elem[1] && $elem[1] eq $link->text) {
			next if $elem[0] =~ /Gene/i;
			my ($g, $sp) = split /\s+/, $elem[3] if defined $elem[3];

			if ($alignments && $all) {
			    $gene_stats{$elem[1]}{$elem[3]} = { $elem[0] => $elem[2]};
			    my ($file, $endpoint) = make_alignment_url_from_gene($link->text, $alphabet, $type);
			    fetch_file($file, $endpoint);
			    unlink $cpbase_response;
			}
			elsif ($alignments && 
			       defined $genus && 
			       $genus =~ /$g/ && 
			       defined $species && 
			       $species =~ /$sp/) {
			    $gene_stats{$elem[1]}{$elem[3]} = { $elem[0] => $elem[2]};
			    my ($file, $endpoint) = make_alignment_url_from_gene($link->text, $alphabet, $type);
			    fetch_file($file, $endpoint);
			    unlink $cpbase_response;
			}
			elsif ($alignments && 
			       defined $genus && 
			       $genus =~ /$g/ && 
			       defined $species && 
			       $species =~ /$sp/ && 
			       defined $gene_name && 
			       $gene_name =~ /$elem[0]/) {
			    $gene_stats{$elem[1]}{$elem[3]} = { $elem[0] => $elem[2]};
			    my ($file, $endpoint) = make_alignment_url_from_gene($link->text, $alphabet, $type);
			    fetch_file($file, $endpoint);
			    unlink $cpbase_response;
			}
			elsif ($alignments && 
			       !defined $genus && 
			       !defined $species && 
			       defined $gene_name && 
			       $gene_name =~ /$elem[1]/) {
			    $gene_stats{$elem[1]}{$elem[3]} = { $elem[0] => $elem[2]};
			    my ($file, $endpoint) = make_alignment_url_from_gene($link->text, $alphabet, $type);
			    fetch_file($file, $endpoint);
			    unlink $cpbase_response;
			}
			else { 
			    say join "\t", $genus, $g, $species, $sp; 
			}
		    }
		}
	    }
	    unlink $cpbase_response if -e $cpbase_response;
	}
    }
    return \%gene_stats;
}

sub make_alignment_url_from_gene {
    my ($gene, $alphabet, $type) = @_;

    my $file = $gene."_orthologs";
    my $endpoint = "http://chloroplast.ocean.washington.edu/CpBase_data/tmp/$gene";
    if ($alphabet =~ /dna/i && $type =~ /fasta/i) {
	$endpoint .= "_orthologs.nt.aln.fa";
	$file .= ".nt.aln.fa";
    }
    elsif ($alphabet =~ /protein/i && $type =~ /fasta/i) {
	$endpoint .= "_orthologs.aa.aln.fa";
	$file .= ".aa.aln.fa";
    }
    elsif ($alphabet =~ /dna/i && $type =~ /clustal/i) {
	$endpoint .= "_orthologs.nt.aln.clw";
	$file .= ".nt.aln.clw";
    }
    elsif ($alphabet =~ /protein/i && $type =~ /clustal/i) {
	$endpoint .= "_orthologs.aa.aln.clw";
	$file .= ".aa.aln.clw";
    }
    else {
	die "\nERROR: Could not determine parameter options for fetching ortholog clusters. alpha: $alphabet type: $type";
    }
    
    return ($file, $endpoint)
    }

sub fetch_rna_clusters {
    my ($statistics, $sequences, $gene_name, $all) = @_;
    my $rna_cluster_stats;
    my %rna_cluster_links;
    my $urlbase = "http://chloroplast.ocean.washington.edu/tools/cpbase/run?view=rna_cluster_index";
    my $mech = WWW::Mechanize->new();
    $mech->get( $urlbase );
    my @links = $mech->links();
    my $gene;
    for my $link ( @links ) {
        next unless defined $link->text;
	if ($link->url =~ /u_feature_id=(\d+)/) {
	    my $id = $1;
	    $gene = $link->text;
	    my $url = "http://chloroplast.ocean.washington/tools/cpbase/run?u_feature_id=$id&view=universal_feature";
	    $rna_cluster_links{$gene} = $url;
	}
    }
    
    if ($statistics && $all) {
	for my $gene (keys %rna_cluster_links) {
	    my $ua = LWP::UserAgent->new;
	    my $response = $ua->get($rna_cluster_links{$gene});
	        
	    unless ($response->is_success) {
		die "Can't get url $urlbase -- ", $response->status_line;
	    }
	        
	    open my $out, '>', $cpbase_response or die "\nERROR: Could not open file: $!\n";
	    say $out $response->content;
	    close $out;
	        
	    my $te = HTML::TableExtract->new( attribs => { border => 1 } );
	    $te->parse_file($cpbase_response);
	        
	    for my $ts ($te->tables) {
		for my $row ($ts->rows) {
		    my @elem = grep { defined } @$row;
		    say join q{, }, @elem;
		}
	    }
	}
    }
    elsif ($statistics && $gene_name) { 
	for my $gene (keys %rna_cluster_links) {
	    my $ua = LWP::UserAgent->new;
	    my $response = $ua->get($rna_cluster_links{$gene});
	        
	    unless ($response->is_success) {
		die "Can't get url $urlbase -- ", $response->status_line;
	    }
	        
	    open my $out, '>', $cpbase_response or die "\nERROR: Could not open file: $!\n";
	    say $out $response->content;
	    close $out;

	    my $te = HTML::TableExtract->new( attribs => { border => 1 } );
	    $te->parse_file($cpbase_response);
	        
	    for my $ts ($te->tables) {
		for my $row ($ts->rows) {
		    my @elem = grep { defined } @$row;
		    say join q{, }, @elem;
		}
	    }
	}
    }
    elsif ($sequences && $all) {
	my $file = $gene."_orthologs.nt.fasta";
	my $endpoint = "http://chloroplast.ocean.washington.edu/CpBase_data/tmp/$file";
	fetch_file($file, $endpoint);
    }
    elsif ($sequences && $gene_name) {
	if ($gene_name eq $gene) {
	    my $file = $gene."_orthologs.nt.fasta";
	    my $endpoint = "http://chloroplast.ocean.washington.edu/CpBase_data/tmp/$file";
	    fetch_file($file, $endpoint);
	}
    }
    elsif ($alignments && $all) {
	my $file = $gene."_orthologs.nt.aln.";
	my $suf; ##TODO
	$suf = "clw" if $type =~ /cl/i;
	$suf = "fa" if $type =~ /fa/i;
	$file .= $suf;
	my $endpoint = "http://chloroplast.ocean.washington.edu/CpBase_data/tmp/$file";
	fetch_file($file, $endpoint);
    }
    elsif ($alignments && $gene_name) {
	if ($gene_name eq $gene) {
	    my $suf; ##TODO
	    my $file = $gene."_orthologs.nt.aln.";
	    $suf = "clw" if $type =~ /cl/i;
	    $suf = "fa" if $type =~ /fa/i;
	    $file .= $suf;
	    my $endpoint = "http://chloroplast.ocean.washington.edu/CpBase_data/tmp/$file";
	    fetch_file($file, $endpoint);
	}
    }
    unlink $cpbase_response;
    ## reorder control flow to if gene_name
    ##                             if sequences
    ##                             elsif alignments
    #return \%rna_cluster_stats;
}

sub get_lineage_for_taxon {
    my ($cpbase_response, $available) = @_;

    my %taxa;
    my $ua = LWP::UserAgent->new;
    my $urlbase = "http://chloroplast.ocean.washington.edu/tools/cpbase/run";
    my $response = $ua->get($urlbase);

    unless ($response->is_success) {
	die "Can't get url $urlbase -- ", $response->status_line;
    }

    open my $out, '>', $cpbase_response or die "\nERROR: Could not open file: $!\n";
    say $out $response->content;
    close $out;

    my $te = HTML::TableExtract->new( attribs => { border => 1 } );
    $te->parse_file($cpbase_response);

    for my $ts ($te->tables) {
	for my $row ($ts->rows) {
	    my @elem = grep { defined } @$row;
	    if ($available) {
		if ($elem[0] =~ /(\d+) Genomes/i) {
		    $genomes = $1;
		    unlink $cpbase_response;
		    say "$genomes genomes available in CpBase." and exit(0);
		}
	    }
	    my ($organism, $locus, $sequence, $length, $assembled, $annotated, $added) = @elem;
	    my ($genus, $species) = split /\s+/, $organism;
	    next unless defined $species && length($species) > 3;
	    next if $species eq 'hybrid';
	    my $taxonid = fetch_taxonid($genus, $species);
	    if (defined $taxonid) {
		my ($lineage, $order, $family) = get_lineage_from_taxonid($taxonid);
		say join "\t", $order, $family, $genus, $species if defined $order && defined $family;
	    } 
	}
    }

    exit;
}

sub get_lineage_from_taxonid {
    my ($id) = @_;
    my $esumm = "esumm_$id.xml"; 
 
    my $ua = LWP::UserAgent->new;
    my $urlbase  = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$id";
    my $response = $ua->get($urlbase);

    unless ($response->is_success) {
	die "Can't get url $urlbase -- ", $response->status_line;
    }

    open my $out, '>', $esumm or die "\nERROR: Could not open file: $!\n";
    say $out $response->content;
    close $out;

    my $parser = XML::LibXML->new;
    my $doc    = $parser->parse_file($esumm);
    my ($order, $family, $lineage);

    for my $node ( $doc->findnodes('//TaxaSet/Taxon') ) {
	($lineage) = $node->findvalue('Lineage/text()');
	if ($lineage =~ /viridiplantae/i) {
	    ($family) = map  { s/\;$//; $_; }
	                grep { /(\w+aceae)/ } 
	                map  { split /\s+/  } $lineage;
	
	    ($order)  = map  { s/\;$//; $_; }
	                grep { /(\w+ales)/ }
	                map  { split /\s+/  } $lineage;

	}
	else {
	    ## need method to get order/family from non-viridiplantae
	    say $lineage;
	}
    }
    unlink $esumm;

    return ($lineage, $order, $family);
}

sub fetch_taxonid {
    my ($genus, $species) = @_;

    my $esearch = "esearch_$genus"."_"."$species.xml";
    my $ua = LWP::UserAgent->new;
    my $urlbase  = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=$genus%20$species";
    my $response = $ua->get($urlbase);

    unless ($response->is_success) {
	die "Can't get url $urlbase -- ", $response->status_line;
    }

    open my $out, '>', $esearch or die "\nERROR: Could not open file: $!\n";
    say $out $response->content;
    close $out;

    my $id;
    my $parser = XML::LibXML->new;
    my $doc    = $parser->parse_file($esearch);
    
    for my $node ( $doc->findnodes('//eSearchResult/IdList') ) {
	($id) = $node->findvalue('Id/text()');
    }
    
    unlink $esearch;
    
    return $id;
}

sub fetch_file {
    my ($file, $endpoint) = @_;

    unless (-e $file) {
	my $exit_code;
	try {
	    $exit_code = system([0..5], "wget -O $file $endpoint");
	}
	catch {
	    say "\nERROR: wget exited abnormally with exit code: $exit_code. Here is the exception: $_\n";
	};
    }
}

sub usage { 
    my $script = basename( $0, () );
    print STDERR <<END

USAGE: perl $script [-l] [-g] [-s] [-d] [-asm] [-aln] [-gn] [-gc] [-rc] [-t] [-mol] [-stats] [--available] [-h] [-m]

Options for metadata from CpBase:
  l|taxonomic_lineage :      Return the order, family, genus and species (tab separated) of all entries.
                             NB: For some taxa, the order and family returned are actually different taxonomic ranks,
                                 and this will have to be resolved on a per species basis. This has to do with the
                                 Entrez data returned from NCBI, though there may be a general solution.
  available           :      Print the number of species available in CpBase and exit.

Options for chlorplast assemblies:
  d|db                :      The database to search.
                             Must be one of: viridiplantae, non_viridiplanate, 'red lineage', rhodophyta, or stramenopiles.
  all                 :      Download files of the specified type for all species in the database.
  g|genus             :      The name of a genus query.
  s|species           :      The name of a species to query.
  asm|assemblies      :      Specifies that the chlorplast genome assemblies should be fetched.
  t|type              :      Format of the sequence file to fetch. Options are: genbank or fasta (Default: fasta).
  stats|statistics    :      Get statistics for the specified species.
  available           :      Print the number of species available in the database and exit.

Options for chloroplast orthologs:
  all                 :      Download files of the specified type for all species in the database.
  aln|alignments      :      Download ortholog alignments for a gene, or all genes.
  t|type              :      Format of the alignment file to fetch. Options are: clustalw or fasta (Default: fasta).
  gn|gene_name        :      The name of a specific gene to fetch ortholog cluster stats or alignments for.
  gc|gene_clusters    :      Fetch gene cluster information.
  mol|alphabet        :      The type of alignments to return. Options are: DNA or protein (Default: DNA).
  stats|statistics    :      Get statistics for the specified species.

Options for RNA orthologs:
  all                 :      Download files of the specified type for all species in the database.
  rc|rna_clusters     :      Download RNA clusters for the specified genes.
  seq|sequences       :      Download RNA cluster ortholog sequences for each gene (if --all) or specific genes (if --gene_name).
  t|type              :      Format of the alignment file to fetch. Options are: clustalw or fasta (Default: fasta).
  mol|alphabet        :      The type of alignments to return. Options are: DNA or protein (Default: DNA).
  stats|statistics    :      Get statistics for the specified species.
  h|help              :      Print a help statement.
  m|man               :      Print the full manual. 

END
}
