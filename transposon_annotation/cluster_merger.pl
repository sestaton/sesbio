#!/usr/bin/env perl

use strict;
use warnings;
use autodie qw(open);
use feature 'say';
use File::Spec qw(catfile rel2abs);
use File::Basename qw(fileparse);
use File::Path qw(make_path);
use Getopt::Long;
use JSON;
use Data::Dump qw(dd dump);

my $usage = "$0 -i cluster_fas_files_dir";
my $indir;
my $database;
my $outdir;
my $json;

GetOptions(
           'i|indir=s'  => \$indir,
           'o|outdir=s' => \$outdir,
           'd|db=s'     => \$database,
           'j|json=s'   => \$json,
           );

die $usage if !$indir;

## get input files                                                                                                                                           
opendir(DIR, $indir) || die "\nERROR: Could not open directory: $indir. Exiting.\n";
my @clus_fas_files = grep /\.fa.*$/, readdir(DIR);
closedir(DIR);

if (scalar(@clus_fas_files) < 1) {
    say "\nERROR: Could not find any fasta files in $indir. Exiting.\n";
    exit(1);
}

## set path to output dir
#my ($oname, $opath, $osuffix) = fileparse($outdir, qr/\.[^.]*/);
#my $db = File::Spec->rel2abs($repeats);
#my $out_path = File::Spec->rel2abs($opath.$oname);
chdir($indir);
my $read_pairs = find_pairs(\@clus_fas_files);

#dd $read_pairs;

my %mapped_pairs;
while (my ($key, $vals) = each %$read_pairs) {
    for my $val (@$vals) {
	$val =~ s/\/\d$//;
	if (exists $mapped_pairs{$val}) {
	    push @{$mapped_pairs{$val}}, $key;
	}
	else {
	    $mapped_pairs{$val} = [$key];
	}
    }
}

my %seen;
for my $key (keys %mapped_pairs) {
    my $pairct = scalar @{$mapped_pairs{$key}};
    $seen{$_}++ for @{$mapped_pairs{$key}};
    my $r = (keys %seen)[0];
    if ($pairct > 1 && $seen{$r} < 2) {
	say "$key ==> $pairct ==> ",join ",", @{$mapped_pairs{$key}};
    }
    undef %seen;
}

#
# subs
#
sub find_pairs {
    my $clus_fas_files = shift;

    my %read_pairs;
    local $/ = '>';
   
    for my $file (@$clus_fas_files) {
	my ($cluster, $readnum, @suffix) = split /\_/, $file;
	open(my $in, '<', $file);	
	while (my $line = <$in>) {
	    $line =~ s/>//g;
	    next if !length($line);
	    my ($seqid, @seqs) = split /\n/, $line;
	    my $seq = join '', @seqs;
	    push @{$read_pairs{$cluster}}, $seqid;  ## Not using the seq at this point
	}
	close($in);
    }

    return \%read_pairs;
}
