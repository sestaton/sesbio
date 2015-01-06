#!/usr/bin/env perl

=head1 NAME 
                                                                       
kew_fetch_plantCvalues.pl - Fetch C-values for a plant family  

=head1 SYNOPSIS    
 
kew_fetch_plantCvalues.pl -f familyname -e email -o resultsfile

=head1 DESCRIPTION
                                                                   
This is a  web client to fetch C-values for a plant family from the Kew Royal Botanic Gardens database (http://data.kew.org/cvalues/) that 
returns a file sorted by ascending genome size for each species in the database as in the example shown below. The data in each column are family, subfamily, tribe, genus, species, chromosome number, ploidy level, and genome size (listed as 1C (Mbp) estimates). Note that some species have no ploidy level or chromosome number estimate as in the first two species in the example below. This is because these values are not listed in the Kew database. Also note that the subfamily and tribe for some species is not listed, as with the third and fourth species below, because these species are not in the NCBI Entrez Taxonomy database. Please note that neither of these features are because of an error with this script, but rather a reflection of what information is present in each database. In addition,the subfamily and tribe may actually be incorrect because these values are missing from NCBI so the taxonomic rank that is expected to be subfamily and tribe is printed.

ASTERACEAE      Cichorioideae   Hypochaeridinae Leontodon       longirostris                    391
ASTERACEAE      Asteroideae     Senecioninae    Pericallis      appendiculata   60      6       533
ASTERACEAE                                      Taraxacum       linearisquameum 16      2       866
ASTERACEAE                                      Centaurea       cuneifolia      18      2       870


(...)

=head1 DEPENDENCIES

This client uses URI to format data for a request, and LWP::UserAgent 
and HTTP GET to perform a request.

Tested with:

(to be added later)

=head1 LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -d, --db

The database to search. Can be one of Angiosperm, Gymnosperm,
Pteridophyte, Bryophyte, Algae. Case is not important but the 
database MUST be spelled correctly.

=item -f, --family

The name of the plant family to search for.

=item -e, --email

An email must be used to fill out the query form online.

=back

=head1 OPTIONS

=over 2

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut      

##TODO: get return type (2C Mbp, 1C Mbp, etc...) Currently, 1C Mbp is returned 
##      add method for donating or supporting Kew gardens, at least in the documentation
##
##      remove use of bio-db-taxonomy

use 5.010;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Bio::DB::Taxonomy;
use URI;
use LWP::UserAgent;
use HTTP::Request::Common qw(GET);
use Pod::Usage;
use Time::HiRes qw(gettimeofday);

# given/when emits warnings in v5.18+
no if $] >= 5.018, 'warnings', "experimental::smartmatch";

my $family;
my $email;
my $outfile;
my $help;
my $man;
my $db;
my $Kew_response = "LWP_Client_Kew_Royal_Botanic_Gardens_Plant-Cvalues_Database.html"; # XHTML

#
# Set Opts
#
GetOptions(
	   'f|family=s'     => \$family,
	   'e|email=s'      => \$email,
	   'o|outfile=s'    => \$outfile,
	   'd|db:s'         => \$db,
	   'h|help'         => \$help,
	   'm|man'          => \$man,
	  );

pod2usage( -verbose => 1 ) if $help;
pod2usage( -verbose => 2 ) if $man;

#
# Check @ARGV
#
if (!$family || !$email || !$outfile) {
   say "\nERROR: Command line not parsed correctly. Exiting.";
   usage();
   exit(0);
}

#
# Set the Kew database to search
#
given ($db) {
    when (/angiosperm/i) {   $db = "Angiosperm"; }
    when (/gymnosperm/i) {   $db = "gymnosperm"; }
    when (/pteridophyte/i) { $db = "pteridophyte"; }
    when (/bryophyte/i) {    $db = "bryophyte"; }
    when (/algae/i) {        $db = "algae"; }
    default {                die "Invalid name for option db."; }
}

#
# Counters
#
my $t0 = gettimeofday();
my $records = 0;

#
# Set the entrez database to search
#
my $entrezdb = Bio::DB::Taxonomy->new(-source => 'entrez');

#
# Set the base URL
#
my $URL = geturlfordb($db,$email,$family);

#
# Create the UserAgent
# 
my $ua = LWP::UserAgent->new;

#
# Define the request
#
my $request = HTTP::Request->new(GET => $URL);

#
# Perform the request
#
my $response = $ua->request($request,$Kew_response);

#
# Check for a response
#
unless ($response->is_success) {
    die "Can't get $URL -- ", $response->status_line;
}

#
# Open and parse the results
#
open my $xhtml, '<', $Kew_response or die "\nERROR: Could not open file: $Kew_response";
open my $Kew_results, '>', $outfile or die "\nERROR: Could not open file: $outfile";

while (my $cvalues = <$xhtml>) {
    chomp $cvalues;
    my $FAM = uc($family);
    while ($cvalues =~ m/($FAM)<.*?><.*?>(\w+)<.*?><.*?>(\w+.*?)<.*?><.*?>(\d+|)<.*?><.*?>(\d+|)<.*?><.*?>(\d+)/g) {
	$records++;
	my $kewfam = $1;
	my $genus = $2;
	my $species = $3;
	my $chrnum = $4;
	my $ploidy = $5;
	my $cval = $6;
	$species =~ s/\(.*//;
	if ($family =~ /Asteraceae/i) {
	    my $taxonid = $entrezdb->get_taxonid("$genus $species");
	    if ( defined($taxonid) ) {
		my $node    = $entrezdb->get_Taxonomy_Node(-taxonid => $taxonid);
                #print "\n$node->rank\n";
                #print "\n$node->get_all_Descendants";
		my $sf = $node;
		my $tr = $node;  
		for ( 1..4 ) { 
		    $sf = $entrezdb->get_Taxonomy_Node(-taxonid => $sf->parent_id);
		}
		for ( 1..2 ) { 
		    $tr = $entrezdb->get_Taxonomy_Node(-taxonid => $tr->parent_id);
		}
		say $Kew_results join "\t", $kewfam, $sf->scientific_name, $tr->scientific_name, $genus, $species, $chrnum, $ploidy, $cval; 
	    } else {
		say $Kew_results join "\t", $kewfam, "", "", $genus, $species, $chrnum, $ploidy, $cval;
	    }
	} else {
	    say $Kew_results join "\t", $kewfam, $genus, $species, $chrnum, $ploidy, $cval;
	}
    }
}

close $xhtml;
close $Kew_results;
unlink $Kew_response;

my $t1 = gettimeofday();
my $elapsed = $t1 - $t0;
my $time = sprintf("%.2f",$elapsed/60);

say "Fetched $records records in $time m:s.";

exit;
#
# Methods
#
sub usage {
    my $script = basename( $0, () );
    print STDERR <<END

$script --db Angiosperm --family Asteraceae --email name\@domain.com -o myresults.txt

Required Arguments:
  -d|db            :       The database to query. Type $script --man for more details.
  -f|family        :       Name of the plant family to search.
  -e|email         :       Valid email to attach to the request.
  -o|outfile       :       File to place the results.

Options:
  -h|help          :       Print a help statement.
  -m|man           :       Print the full manual. 

END

}

sub geturlfordb {
    my ($database, $em, $fam) = @_;
    my $url = URI->new('http://data.kew.org/cvalues/CvalServlet');
    
    if ($database =~ m/angiosperm/i) {
	$url->query_form(
			 'generatedby'         => $database,  
			 'querytype'           => "-1",    
			 'txtEmail'            => $em,       # uri_escape() from URI::Escape is called automatically. Nice!
			 'chkFamily'           => "on",
			 'selectfamilytype'    => "phyloFam",
			 'chkGenus'            => "on",
			 'chkSpecies'          => "on",
			 'chkChromosome'       => "on",
			 'chkPloidyLevel'      => "on",
			 'cmbFourcSelect'      => "1",   # this is the return type, currently 1C Mbp
			 'chkOrigReference'    => "on",
			 'cmbPrimeOrAllData'   => "0",
			 'txtFamily'           => $fam,
			 'qryfamilytype'       => "phyloFam",
			 'txtGenus'            => "",
			 'txtSpecies'          => "",
			 'txtAuthority'        => "",
			 'txtFromFourc'        => "",
			 'txtToFourc'          => "",
			 'txtFromChromosome'   => "",
			 'txtToChromosome'     => "",
			 'cmbEstimationMethod' => "0",
			 'txtFromPloidy'       => "",
			 'txtToPloidy'         => "",
			 'rdoVoucher'          => "2",
			 'cmbLifeCycle'        => "0",    # different for algae
			 'rdoAngioGroup'       => "3",
			 'cmbSort'             => "0",
			);
	return $url;
    }
    if ($database =~ m/gymnosperm/i) {
	$url->query_form(
			 'generatedby'         => $database,
			 'querytype'           => "-1",
			 'txtEmail'            => $em,
			 'chkFamily'           => "on",
			 'selectfamilytype'    => "pubFam",
			 'chkGenus'            => "on",
			 'chkSpecies'          => "on",
			 'chkChromosome'       => "on",
			 'chkPloidyLevel'      => "on",
			 'cmbFourcSelect'      => "1",
			 'chkOrigReference'    => "on",
			 'cmbPrimeOrAllData'   => "0",
			 'txtFamily'           => $fam,
			 'qryfamilytype'       => "pubFam",
			 'txtGenus'            => "",
			 'txtSpecies'          => "",
			 'txtAuthority'        => "",
			 'txtFromFourc'        => "",
			 'txtToFourc'          => "",
			 'txtFromChromosome'   => "",
			 'txtToChromosome'     => "",
			 'cmbEstimationMethod' => "0",
			 'txtFromPloidy'       => "",
			 'txtToPloidy'         => "",
			 'rdoVoucher'          => "2",
			 'cmdGymnoGroup'       => "0",
			 'cmbFlagellaNumber'   => "0",
			 'cmbSort'             => "0",
			);
	return $url;
    }
    if ($database =~ m/pteridophyte/i) {
	$url->query_form(
			 'generatedby'         => $database,
			 'querytype'           => "-1",
			 'txtEmail'            => $em,
			 'chkFamily'           => "on",
			 'selectfamilytype'    => "pubFam",
			 'chkGenus'            => "on",
			 'chkSpecies'          => "on",
			 'chkChromosome'       => "on",
			 'chkPloidyLevel'      => "on",
			 'cmbFourcSelect'      => "1",
			 'chkOrigReference'    => "on",
			 'cmbPrimeOrAllData'   => "0",
			 'txtFamily'           => $fam, #+Ophioglossaceae
			 'qryfamilytype'       => 'pubFam',
			 'txtGenus'            => "",
			 'txtSpecies'          => "",
			 'txtAuthority'        => "",
			 'txtFromFourc'        => "",
			 'txtToFourc'          => "",
			 'txtFromChromosome'   => "",
			 'txtToChromosome'     => "",
			 'cmbEstimationMethod' => "0",
			 'txtFromPloidy'       => "",
			 'txtToPloidy'         => "",
			 'rdoVoucher'          => "2",
			 'cmbSporeType'        => "0",
			 'cmbPteridoGroup'     => "0",
			 'cmbSporangium'       => "0",
			 'cmbSpermType'        => "0",
			 'cmbSort'             => "0",
			 );
	return $url;
    }
    if ($database =~ m/bryophyte/i) {
	$url->query_form(
			 'generatedby'         => $database,
			 'querytype'           => "-1",
			 'txtEmail'            => $em,
			 'chkFamily'           => "on",
			 'selectfamilytype'    => "pubFam",
			 'chkGenus'            => "on",
			 'chkSpecies'          => "on",
			 'chkChromosome'       => "on",
			 'chkPloidyLevel'      => "on",
			 'cmbFourcSelect'      => "1",
			 'chkOrigReference'    => "on",
			 'cmbPrimeOrAllData'   => "0",
			 'txtFamily'           => $fam,
			 'qryfamilytype'       => "pubFam",
			 'txtGenus'            => "",
			 'txtSpecies'          => "",
			 'txtAuthority'        => "",
			 'txtFromFourc'        => "",
			 'txtToFourc'          => "",
			 'txtFromChromosome'   => "",
			 'txtToChromosome'     => "",
			 'cmbEstimationMethod' => "0",
			 'txtFromPloidy'       => "",
			 'txtToPloidy'         => "",
			 'rdoVoucher'          => "2",
			 'cmbBryoGroup'        => "0",
			 'cmbSort'             => "0",
			 );
	return $url;
    }
    if ($database =~ m/algae/i) {
	$url->query_form(
			 'generatedby'         => $database,
			 'querytype'           => "-1",
			 'txtEmail'            => $em,
			 'chkFamily'           => "on",
			 'selectfamilytype'    => "pubFam",
			 'chkGenus'            => "on",
			 'chkSpecies'          => "on",
			 'chkChromosome'       => "on",
			 'chkPloidyLevel'      => "on",
			 'cmbFourcSelect'      => "1",
			 'chkOrigReference'    => "on",
			 'cmbPrimeOrAllData'   => "0",
			 'txtFamily'           => $fam,
			 'qryfamilytype'       => "pubFam",
			 'txtGenus'            => "",
			 'txtSpecies'          => "",
			 'txtAuthority'        => "",
			 'txtFromFourc'        => "",
			 'txtToFourc'          => "",
			 'txtFromChromosome'   => "",
			 'txtToChromosome'     => "",
			 'cmbEstimationMethod' => "0",
			 'txtFromPloidy'       => "",
			 'txtToPloidy'         => "",
			 'rdoVoucher'          => "2",
			 'cmdAlgaeGroup'       => "0",
			 'txtAlgaeOrder'       => "",
			 'cmbSort'             => "0",
			 );
	return $url;
    }
}
