use 5.010;
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
    say "\nERROR: Command line not parsed correctly. Required arguments missing.\n";
    say "USAGE: ".basename($0)." -u USER -p PASSWORD -f phytozome_dirlisting.xml\n";
    exit(1);
}

my $cookie = generate_cookie($opts{user}, $opts{password});
my $xmldoc = XML::LibXML->load_xml(location => $opts{xmlfile}, no_blanks => 1);

for my $sample ($xmldoc->findnodes('/organismDownloads/folder')) {
    next if $sample->{name} eq 'global_analysis' || $sample->{name} eq 'early_release';
    #say $sample->{name}; ## species level
    for my $property ($sample->findnodes('./folder')) {
	#say $property->{name}; ## annotation/assembly
	if ($property->{name} eq 'assembly') {
	    for my $file (grep { $_->{filename} !~ /masked/ } $property->findnodes('./file')) {
		## apply species-level filter
		fetch_file($sample->{name}, $file->{url}, $file->{filename}, $cookie);
	    }
	}
    }
}
unlink $cookie;

sub fetch_file {
    my ($sample, $url, $file, $cookie) = @_;
    
    my $urlbase = 'http://genome.jgi.doe.gov';
    my ($dir) = ($url =~ /url=(\S+)$/);
    my $endpoint = $urlbase.$dir;
    say "=====> Fetching $file for sample $sample at: $endpoint";
    system("curl $endpoint -b $cookie > $file") == 0
	or die $!;
}

sub generate_cookie {
    my ($user, $pass) = @_;
    my $signon = 'https://signon.jgi.doe.gov/signon/create';
    my $cookie = 'jgi_cookie';
    system("curl $signon --data-urlencode 'login=$user' --data-urlencode 'password=$pass' -c $cookie > /dev/null") == 0
	or die $!;

    return $cookie;
}
