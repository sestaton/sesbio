#!/usr/bin/env perl

##TODO: Use new Phytomine API via Webservice::InterMine (https://metacpan.org/pod/Webservice::InterMine)

use strict;
use warnings;
use autodie;
use File::Basename;
use XML::LibXML;
use Getopt::Long;

my %opts;
GetOptions(\%opts, 
	   'u|user=s'     => \$opts{user}, 
	   'p|password=s' => \$opts{password},
	   's|species=s'  => \$opts{species},
	   'f|xmlfile=s'  => \$opts{xmlfile},
    );

if (!$opts{user} || !$opts{password} || !$opts{xmlfile}) {
    print "\nERROR: Command line not parsed correctly. Required arguments missing.\n\n";
    print "USAGE: ".basename($0)." -u USER -p PASSWORD -f phytozome_dirlisting.xml\n\n". 
        "\tOptionally, you may select a single species to download:\n\n".
	basename($0)." -u USER -p PASSWORD -f phytozome_dirlisting.xml -s 'genus species'\n\n";
    exit(1);
}

my $cookie = generate_cookie($opts{user}, $opts{password});
my $xmldoc = XML::LibXML->load_xml(location => $opts{xmlfile}, no_blanks => 1);
my ($gs_map, $abrv_map) = species_listing();

unless ($opts{species}) {
    print "WARNING: No species name was specified so all genomes will be downloaded.....";
    sleep 5;
    print "Okay, downloading all genomes.\n\n";
}

for my $sample ($xmldoc->findnodes('/organismDownloads/folder')) {
    next if $sample->{name} eq 'global_analysis' || $sample->{name} eq 'early_release';
    for my $property ($sample->findnodes('./folder')) {
	if ($property->{name} eq 'assembly') {
	    for my $file (grep { $_->{filename} !~ /masked/ } $property->findnodes('./file')) {
		if ($opts{species}) {
		    my $sp = $abrv_map->{ $sample->{name} };
		    if (exists $gs_map->{ $opts{species} } && defined $sp && $sp eq $opts{species}) {
			fetch_file($sample->{name}, $file->{url}, $file->{filename}, $cookie, $abrv_map);
			unlink $cookie;
			exit(0);
		    }
		}
		else {
		    fetch_file($sample->{name}, $file->{url}, $file->{filename}, $cookie, $abrv_map);
		}
	    }
	}
    }
}
unlink $cookie;

sub fetch_file {
    my ($sample, $url, $file, $cookie, $abrv_map) = @_;
    
    my $urlbase  = 'http://genome.jgi.doe.gov';
    my ($dir)    = ($url =~ /url=(\S+)$/);
    my $endpoint = $urlbase.$dir;
    my $genera   = $abrv_map->{$sample};
    print "=====> Fetching $file for $genera..."; # at: $endpoint";
    system("curl -s $endpoint -b $cookie > $file 2> /dev/null") == 0
	or die $!;
    print "Done.\n";
}

sub generate_cookie {
    my ($user, $pass) = @_;
    ##TODO: check if response was successful
    my $signon = 'https://signon.jgi.doe.gov/signon/create';
    my $cookie = 'jgi_cookie';
    system("curl -s $signon --data-urlencode 'login=$user' --data-urlencode 'password=$pass' -c $cookie 1> /dev/null") == 0
	or die $!;

    return $cookie;
}

sub species_listing {
    my %species = (
	Zmays               => 'Zea mays',
	Vvinifera           => 'Vitis vinifera',
	Vcarteri            => 'Volvox carteri',
	Tcacao              => 'Theobroma cacao',
	Sviridis            => 'Setaria viridis',
	Stuberosum          => 'Solanum tuberosum',
	Spurpurea           => 'Salix purpurea',
	Spolyrhiza          => 'Spirodela polyrhiza',
	Smoellendorffii     => 'Selaginella moellendorffii',
	Slycopersicum       => 'Solanum lycopersicum',
	Sitalica            => 'Setaria italica',
	Sbicolor            => 'Sorghum bicolor',
	Rcommunis           => 'Ricinus communis',
	Pvulgaris           => 'Phaseolus vulgaris',
	Pvirgatum           => 'Panicum virgatum',
	Ptrichocarpa        => 'Populus trichocarpa',
	Ppersica            => 'Prunus persica',
	Ppatens             => 'Physcomitrella patens',
	Phallii             => 'Panicum hallii ',
	Osativa             => 'Oryza sativa',
	Olucimarinus        => 'Ostreococcus lucimarinus',
	Mtruncatula         => 'Medicago truncatula',
	#MspRCC299           => 'Micromonas sp. RCC299',
	MpusillaCCMP1545    => 'Micromonas pusilla', # CCMP1545',
	Mguttatus           => 'Mimulus guttatus',
	Mesculenta          => 'Manihot esculenta',
	Mdomestica          => 'Malus domestica',
	Macuminata          => 'Musa acuminata',
	Lusitatissimum      => 'Linum usitatissimum',
	Kmarnieriana        => 'Kalanchoe marnieriana',
	Graimondii          => 'Gossypium raimondii',
	Gmax                => 'Glycine max',
	Fvesca              => 'Fragaria vesca',
	Esalsugineum        => 'Eutrema salsugineum',
	Egrandis            => 'Eucalyptus grandis',
	CsubellipsoideaC169 => 'Coccomyxa subellipsoidea', # C-169',
	Csinensis           => 'Citrus sinensis',
	Csativus            => 'Cucumis sativus',
	Crubella            => 'Capsella rubella',
	Creinhardtii        => 'Chlamydomonas reinhardtii',
	Cpapaya             => 'Carica papaya',
	Cgrandiflora        => 'Capsella grandiflora',
	Cclementina         => 'Citrus clementina',
	Bstricta            => 'Boechera strict',
	Bstacei             => 'Brachypodium stacei',
	BrapaFPsc           => 'Brassica rapa FPsc',
	Bdistachyon         => 'Brachypodium distachyon',
	Atrichopoda         => 'Amborella trichopoda',
	Athaliana           => 'Arabidopsis thaliana', # columbia',
	Alyrata             => 'Arabidopsis lyrata',
	Acoerulea           => 'Aquilegia coerulea' ,# Goldsmith',
	);

    my %sample_map = reverse %species;
    return (\%sample_map, \%species);
}
