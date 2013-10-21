#!/usr/bin/env perl

#TODO: add POD

#
# library imports
#
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
use Data::Dump qw(dd);

# given/when emits warnings in v5.18+
no if $] >= 5.018, 'warnings', "experimental::switch";

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
my $assemblies;
my $cpbase_response = "CpBase_database_response.html"; # HTML

#
# set opts
#
GetOptions(
           'all'              => \$all,
           'd|db=s'           => \$db,
           't|type=s'         => \$type,
           'g|genus=s'        => \$genus,
	   's|species=s'      => \$species,
	   'o|outfile=s'      => \$outfile,
           'stats|statistics' => \$statistics,
           'available'        => \$available,
           'asm|assemblies'   => \$assemblies,
	   'h|help'           => \$help,
	   'm|man'            => \$man,
	  );

#pod2usage( -verbose => 1 ) if $help;
#pod2usage( -verbose => 2 ) if $man;

#
# check @ARGV
#
if (!$db) {
   say "\nERROR: A database to query must be given. Exiting.";
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
    when (/viridiplantae/i) {     $db = "Viridiplantae"; }
    when (/non viridiplantae/i) { $db = "NOT_Viridiplantae"; }
    default {                     die "Invalid name for option db."; }
}

$type //= 'fasta';

#
# create the UserAgent
# 
my $ua = LWP::UserAgent->new;

#
# perform the request
#
my $urlbase = "http://chloroplast.ocean.washington.edu/tools/cpbase/run&genome_taxonomy=$db";
my $response = $ua->get($urlbase);

#
# check for a response
unless ($response->is_success) {
    die "Can't get url $urlbase -- ", $response->status_line;
}

#
# open and parse the results
#
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
	    say "$genomes $db genomes available in CpBase." and exit if $available;
	}
	else {
	    my ($organism,$locus,$sequence_length,$assembled,$annotated,$added) = @elem;
	    $organism =~ s/\s+/_/g;
	    
	    if (exists $id_map->{$organism}) {
		if ($genus && $species && $organism =~ /\Q$epithet\E/i) {
		    my $id = $id_map->{$organism};
		    my $assem_stats = get_cp_data($id);
		    $stats{$organism} = $assem_stats;
		    #dd \%stats;
		    # annotations
		    # http://chloroplast.ocean.washington.edu/tools/cpbase/run?genome_id=750&view=genome
		    ## write code to fetch assemblies here
		    #say join "\t", $type, $locus, $organism;
		    fetch_files($type, $locus, $organism) if $assemblies;
		}
		elsif ($genus && $organism =~ /\Q$genus\E/i) {
		    #say $organism;
		    my $id = $id_map->{$organism};
		    my $assem_stats = get_cp_data($id);
		    $stats{$organism} = $assem_stats;
		}
	    }
	}
    }
}
dd \%stats;
unlink $cpbase_response;

#
# subroutines
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
	    #printf "%s %s, %d\n", $g, $sp, $1;
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
	    #say join q{,}, @elem;
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
    return \%assem_stats;
}

sub fetch_files {
    my ($type, $locus, $organism) = @_;

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

    my $exit_code;
    try {
	$exit_code = system([0..5], "wget -O $file $endpoint");
    }
    catch {
	say "\nERROR: wget exited abnormally with exit code $exit_code. Here is the exception: $_\n";
    };
}

sub usage { 
    my $script = basename( $0, () );
    print STDERR <<END

USAGE: perl $script [-g] [-s] [-d]

Required Arguments:
  d|db              :      The database to search. Must be one of: Viridiplantae, non_viridiplanate, 
  g|genus           :      The name of a genus query.
  s|species         :      The name of a species to query.    
  
Options:
  seq|sequences     :      Specifies that the raw EST sequences should be fetched.
  aln|alignments    :      Specifies that the assemblies aligned to Arabidopsis should be fetched.
  asm|assemblies    :      Specifies that the EST assemblies should be fetched.
  all               :      Download files of the specified type for all species in the database.
  h|help            :      Print a help statement.
  m|man             :      Print the full manual. 

END
}
