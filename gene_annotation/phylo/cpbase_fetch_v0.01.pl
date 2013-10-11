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
use LWP::UserAgent;
use Time::HiRes qw(gettimeofday);
use HTML::TreeBuilder;
use Data::Dump qw(dd);
use IPC::System::Simple qw(system);
use Try::Tiny;
use Pod::Usage;

# given/when emits warnings in v5.18+
no if $] >= 5.018, 'warnings', "experimental::switch";

#
# lexical vars
#
my $db;
my $all;
my $genus;
my $species;
my $outfile; ## log
my $help;
my $man;
my $sequences;
my $alignments;
my $assemblies;
my $cpbase_response = "CpBase_database_response.html"; # HTML

#
# set opts
#
GetOptions(
           'all'              => \$all,
           'd|db=s'           => \$db,
	   'g|genus=s'        => \$genus,
	   's|species=s'      => \$species,
	   'o|outfile=s'      => \$outfile,
	   'seq|sequences'    => \$sequences,
           'aln|alignments'   => \$alignments,
           'asm|assemblies'   => \$assemblies,
	   'h|help'           => \$help,
	   'm|man'            => \$man,
	  );

#pod2usage( -verbose => 1 ) if $help;
#pod2usage( -verbose => 2 ) if $man;

#
# check @ARGV
#
#if (!$assemblies && !$sequences && !$alignments) {
#   say "\nERROR: Command line not parsed correctly. Exiting.";
#   usage();
#   exit(1);
#}

die "USAGE: perl $0 dbname\n" if !$db;


#
# counters
#
my $t0 = gettimeofday();
my $records = 0;

#
# set terms for search
#
my ($gen, $sp);
if ($genus) {
    $gen = substr($genus, 0, 4);
}
if ($species) {
    $sp = substr($species, 0, 4);
}

#
# Set the CpBase database to search
#
given ($db) {
    when (/algae/i) {             $db = "Algae"; }
    when (/red lineage/i) {       $db = "Red_Lineage"; }
    when (/rhodophyta/i) {        $db = "Rhodophyta"; }
    when (/stramenopiles/i) {     $db = "stramenopiles"; }
    when (/viridiplantae/i) {      $db = "Viridiplantae"; }
    when (/non viridiplantae/i) { $db = "NOT_Viridiplantae"; }
    default {                     die "Invalid name for option db."; }
}

#
# create the UserAgent
# 
my $ua = LWP::UserAgent->new;
my $tree = HTML::TreeBuilder->new;

#
# perform the request
#
my $urlbase = "http://chloroplast.ocean.washington.edu/tools/cpbase/run&genome_taxonomy=$db";
my $response = $ua->get($urlbase);

#
# check for a response
#
unless ($response->is_success) {
    die "Can't get url $urlbase -- ", $response->status_line;
}

#
# cpen and parse the results
#
open my $out, '>', $cpbase_response or die "\nERROR: Could not open file: $!\n";
say $out $response->content;
close $out;
$tree->parse_file($cpbase_response);

for my $tag ($tree->look_down(_tag => 'td')) {
    if ($tag->attr('href')) {
	say $tag->as_text;
    }
}

#unlink $cpbase_response;

#
# subroutines
#
sub usage { 
    my $script = basename( $0, () );
    print STDERR <<END

USAGE: perl $script [-seq] [-aln] [-asm] [-g] [-s] [--all]

Required Arguments:
  o|outfile         :      File to place the results (NOT IMPLEMENTED).
  seq|sequences     :      Specifies that the raw EST sequences should be fetched.
  aln|alignments    :      Specifies that the assemblies aligned to Arabidopsis should be fetched.
  asm|assemblies    :      Specifies that the EST assemblies should be fetched.

Options:
  all               :      Download files of the specified type for all species in the database.
  g|genus           :      The name of a genus query.
  s|species         :      The name of a species to query.
  h|help            :      Print a help statement.
  m|man             :      Print the full manual. 

END
}

sub set_url_for_db { 
    my ($db) = @_;

# Viridiplantae epithet list

#Acidosasa purpurea
#Acorus americanus
#Acorus calamus
#Acutodesmus obliquus
#Adiantum capillus-veneris
#Aethionema cordifolium
#Aethionema grandiflorum
#Ageratina adenophora
#Agrostis stolonifera
#Alsophila spinulosa
#Amborella trichopoda
#Angiopteris evecta
#Anomochloa marantoidea
#Anthoceros formosae
#Anthriscus cerefolium
#Arabidopsis thaliana
#Arabis hirsuta
#Ardisia polysticta
#Artemisia frigida
#Arundinaria gigantea
#Atropa belladonna
#Bambusa emeiensis
#Bambusa oldhamii
#Barbarea verna
#Bismarckia nobilis
#Boea hygrometrica
#Brachypodium distachyon
#Brassica napus
#Bryopsis hypnoides
#Buxus microphylla
#Calamus caryotoides
#Calycanthus floridus var. glaucus
#Camellia sinensis
#Capsella bursa-pastoris
#Capsicum annuum
#Carica papaya
#Castanea mollissima
#Catharanthus roseus
#Cathaya argyrophylla
#Cedrus deodara
#Cephalotaxus oliveri
#Cephalotaxus wilsoniana
#Ceratophyllum demersum
#Chaetosphaeridium globosum
#Chara vulgaris
#Cheilanthes lindheimeri
#Chlamydomonas reinhardtii
#Chloranthus spicatus
#Chlorella variabilis
#Chlorella vulgaris
#Chlorokybus atmophyticus
#Chrysanthemum indicum
#Chrysanthemum x morifolium
#Cicer arietinum
#Cistanche deserticola
#Citrus sinensis
#Coccomyxa subellipsoidea C-169
#Coffea arabica
#Coix lacryma-jobi
#Colocasia esculenta
#Corynocarpus laevigata
#Crucihimalaya wallichii
#Cryptomeria japonica
#Cucumis melo subsp. melo
#Cucumis sativus
#Cunninghamia lanceolata
#Cycas revoluta
#Cycas taitungensis
#Cymbidium aloifolium
#Cymbidium mannii
#Cymbidium sinense
#Cymbidium tortisepalum
#Cymbidium tracyanum
#Dasypogon bromeliifolius
#Datura stramonium
#Daucus carota
#Dendrocalamus latiflorus
#Dioscorea elephantipes
#Draba nemorosa
#Drimys granadensis
#Dunaliella salina
#Elaeis guineensis
#Eleutherococcus senticosus
#Elodea canadensis
#Ephedra equisetina
#Equisetum arvense
#Equisetum hyemale
#Erodium carvifolium
#Erodium guttatum
#Erodium texanum
#Eucalyptus globulus subsp. globulus
#Eucalyptus grandis
#Fagopyrum esculentum subsp. ancestrale
#Ferrocalamus rimosivaginus
#Festuca altissima
#Festuca arundinacea
#Festuca ovina
#Festuca pratensis
#Floydiella terrestris
#Fragaria chiloensis
#Fragaria mandshurica
#Fragaria vesca subsp. bracteata
#Fragaria vesca subsp. vesca
#Fragaria virginiana
#Francoa sonchifolia
#Geranium palmatum
#Ginkgo biloba
#Glycine max
#Glycine tomentella
#Gnetum montanum
#Gnetum parvifolium
#Gonium pectorale
#Gossypium areysianum
#Gossypium barbadense
#Gossypium capitis-viridis
#Gossypium darwinii
#Gossypium gossypioides
#Gossypium herbaceum subsp. africanum
#Gossypium hirsutum
#Gossypium incanum
#Gossypium mustelinum
#Gossypium raimondii
#Gossypium robinsonii
#Gossypium somalense
#Gossypium thurberi
#Gossypium tomentosum
#Guizotia abyssinica
#Helianthus annuus
#Heliconia collinsiana
#Hevea brasiliensis
#Hordeum vulgare subsp. vulgare
#Huperzia lucidula
#Illicium oligandrum
#Indocalamus longiauritus
#Ipomoea purpurea
#Isoetes flaccida
#Jacobaea vulgaris
#Jasminum nudiflorum
#Jatropha curcas
#Keteleeria davidiana
#Lactuca sativa
#Larix decidua
#Lathyrus sativus
#Leersia tisserantii
#Lemna minor
#Lepidium virginicum
#Leptosira terrestris
#Liriodendron tulipifera
#Lobularia maritima
#Lolium multiflorum
#Lolium perenne
#Lotus japonicus
#Magnolia denudata
#Magnolia grandiflora
#Magnolia kwangsiensis
#Magnolia officinalis subsp. biloba
#Manihot esculenta
#Mankyua chejuensis
#Marchantia polymorpha
#Medicago truncatula
#Megaleranthis saniculifolia
#Mesostigma viride
#Micromonas sp. RCC299
#Millettia pinnata
#Monomastix sp. OKE-1
#Monsonia speciosa
#Morus indica
#Nandina domestica
#Nasturtium officinale
#Nelumbo lutea
#Nelumbo nucifera
#Neottia nidus-avis
#Nephroselmis olivacea
#Nicotiana sylvestris
#Nicotiana tabacum
#Nicotiana tomentosiformis
#Nicotiana undulata
#Nothoceros aenigmaticus
#Nuphar advena
#Nymphaea alba
#Oedogonium cardiacum
#Oenothera argillicola
#Oenothera biennis
#Oenothera elata subsp. hookeri
#Oenothera glazioviana
#Oenothera parviflora
#Olea europaea
#Olea europaea subsp. cuspidata
#Olea europaea subsp. europaea
#Olea europaea subsp. maroccana
#Olea woodiana subsp. woodiana
#Olimarabidopsis pumila
#Oltmannsiellopsis viridis
#Oncidium hybrid cultivar
#Ophioglossum californicum
#Oryza meridionalis
#Oryza nivara
#Oryza rufipogon
#Oryza sativa Indica Group
#Oryza sativa Japonica Group
#Ostreococcus tauri
#Pachycladon cheesemanii
#Pachycladon enysii
#Panax ginseng
#Panicum virgatum
#Parachlorella kessleri
#Pedinomonas minor
#Pelargonium x hortorum
#Pellia endiviifolia
#Pentactina rupicola
#Phalaenopsis aphrodite subsp. formosana
#Phalaenopsis equestris
#Pharus latifolius
#Phaseolus vulgaris
#Phoenix dactylifera
#Phyllostachys edulis
#Phyllostachys nigra var. henonis
#Phyllostachys propinqua
#Physcomitrella patens subsp. patens
#Picea abies
#Picea morrisonicola
#Picea sitchensis
#Pinus contorta
#Pinus gerardiana
#Pinus koraiensis
#Pinus krempfii
#Pinus lambertiana
#Pinus massoniana
#Pinus monophylla
#Pinus nelsonii
#Pinus taeda
#Pinus thunbergii
#Piper cenocladum
#Pisum sativum
#Platanus occidentalis
#Pleodorina starrii
#Podocarpus totara
#Populus alba
#Populus trichocarpa
#Prinsepia utilis
#Prunus persica
#Pseudendoclonium akinetum
#Pseudophoenix vinifera
#Pseudotsuga sinensis var. wilsoniana
#Psilotum nudum
#Pteridium aquilinum subsp. aquilinum
#Ptilidium pulcherrimum
#Pycnococcus provasolii
#Pyramimonas parkeae
#Pyrus pyrifolia
#Quercus rubra
#Ranunculus macranthus
#Rhizanthella gardneri
#Rhynchoryza subulata
#Ricinus communis
#Saccharum hybrid cultivar NCo 310
#Saccharum hybrid cultivar SP80-3280
#Salvia miltiorrhiza
#Schizomeris leibleinii
#Selaginella moellendorffii
#Sesamum indicum
#Silene conica
#Silene latifolia
#Silene noctiflora
#Solanum bulbocastanum
#Solanum lycopersicum
#Solanum tuberosum
#Sorghum bicolor
#Spinacia oleracea
#Spirodela polyrhiza
#Staurastrum punctulatum
#Stigeoclonium helveticum
#Syntrichia ruralis
#Taiwania cryptomerioides
#Taiwania flousiana
#Taxus mairei
#Tectona grandis
#Tetracentron sinense
#Theobroma cacao
#Trachelium caeruleum
#Trebouxiophyceae sp. MX-AZ01
#Trifolium subterraneum
#Trithuria inconspicua
#Triticum aestivum
#Trochodendron aralioides
#Typha latifolia
#Utricularia gibba
#Vaccinium macrocarpon
#Vigna angularis
#Vigna radiata
#Vigna unguiculata
#Vitis vinifera
#Welwitschia mirabilis
#Wolffia australiana
#Wolffiella lingulata
#Zea mays
#Zingiber spectabile
#Zygnema circumcarinatum

}
