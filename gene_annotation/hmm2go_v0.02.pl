#!/usr/bin/perl -w

use strict; 
use Getopt::Long;
use Data::Dumper;
#use File::Copy;

my $usage = "\n$0 -i annot -p pfam2go -o outfile <--map>\n";
my $infile;
my $pfam2go; 
my $outfile;
my $mapping;
my $mapfile;
my $map_fh;

GetOptions(
	   'i|infile=s'  => \$infile,
	   'p|pfam2go=s' => \$pfam2go,
	   'o|outfile=s' => \$outfile,
	   'map'         => \$mapping,
	   );

if (!$infile) {
    die "\nERROR: no infile found.\n",$usage;
}
if (!$pfam2go) {
    die "\nERROR: No pfam2go db found.\n",$usage;
}
if (!$outfile) {
    die "\nERROR: No outfile found.\n",$usage;
}
open(my $in, '<', $infile) or die "\nERROR: Could not open file: $infile\n";
open(my $pfams, '<', $pfam2go) or die "\nERROR: Could not open file: $pfam2go\n";

open(my $out, '>', $outfile) or die "\nERROR: Could not open file: $outfile\n";

if ($mapping) {
    $mapfile = $outfile;
    $mapfile =~ s/\..*//g;
    $mapfile .= "_GOterm_mapping.txt";
    open($map_fh, '>', $mapfile) or die "\nERROR: Could not open file: $mapfile\n";
}

my %pfamids;
while(<$in>) {
    chomp;
    next if /^\#/;
    my ($target_name, $accession, $query_name, $accession_q, $E_value_full, $score_full, $bias_full, $E_value_best, $score_best, $bias_best, $exp, $reg, $clu, $ov, $env, $dom, $rev, $inc, $description_of_target) = split;
    my $query_eval = join(",",$query_name,$E_value_full,$description_of_target);
    $accession =~ s/\..*//;
    $pfamids{$query_eval} = $accession;
}
close($in);

#print Dumper %pfamids;

#my %mappedterms;
my %goterms;
my $go_ct = 0;
my $map_ct = 0;

while(my $mapping = <$pfams>) {
    chomp $mapping;
    next if $mapping =~ /^!/;
    if ($mapping =~ /Pfam:(\S+) (\S+ \> )(GO\:\S+.*\;) (GO\:\d+)/) {
	my $pf = $1;
	my $pf_name = $2;
	my $pf_desc = $3;
	my $go_term = $4;
	$pf_name =~ s/\s.*//;
	$pf_desc =~ s/\s\;//;
	foreach my $key (keys %pfamids) { 
	    my ($query, $eval, $desc) = split(/\,/,$key);
	    if ($pfamids{$key} eq $pf) {
		#print "$key => $pfamids{$key}\n";
		print $out join("\t",($query,$pf,$pf_name,$pf_desc,$go_term,$desc)),"\n";
		#my $definition = join("\t",($pf,$pf_name,$pf_desc,$go_term,$desc));
		#$mappedterms{$query} = $definition;
		if ($mapping) {
		    if (exists $goterms{$query}) {
			$go_ct++ if defined($go_term);
			$goterms{$query} .= ",".$go_term;
		    } else {
			$goterms{$query} = $go_term;
		    }
		}
		last;
	    }
	}
    }
}
close($pfams);
close($out);

if ($mapping) {
    while(my ($key, $value) = each(%goterms)) {
	$map_ct++;
	print $map_fh join("\t",$key,$value),"\n";
    }
    print "\n$map_ct query sequences with $go_ct GO terms mapped in file $mapfile.\n\n";
}

#my $map_ct = map_go_terms($out, $map_fh, %goterms);

#print "\n$map_ct query sequences with $go_ct GO terms mapped in file $map_fh.\n\n";

exit;

#
#
#
#sub map_go_terms {
#    my ($out, $map_fh, %goterms) = @_;

#    my %goterms;
#    my $map_ct = 0;
#     my $go_ct = 0;

#    while(my ($mapkey, $mapval) = each(%mappedterms)) {
#	my @hmm_fields = split(/\t/,$mapval);
#	print $out join("\t",$mapkey,$mapval),"\n";
#	if ($mapping) {
#	    if (exists $goterms{$mapkey}) {
#		$go_ct++ if $hmm_fields[3];
#		$goterms{$mapkey} .= ",".$hmm_fields[3];
#	    } else {
#		$goterms{$mapkey} = $hmm_fields[3];
#	    }
#	}
#    }

#    while(<$out>) {
#	chomp;
#	my @fields = split(/\t/,$_);
#	if (exists $goterms{$fields[0]}) {
#	    $goterms{$fields[0]} .= ",".$fields[4];
#	} else {
#	    $goterms{$fields[0]} = $fields[4];
#	}
#    }

#    while(my ($key, $value) = each(%goterms)) {
#	$map_ct++;
#	print $map_fh join("\t",$key,$value),"\n";
#    }
#    return ($map_ct, $go_ct);
#}



