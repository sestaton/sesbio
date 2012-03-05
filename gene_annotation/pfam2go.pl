#!/usr/bin/perl 

use strict; 
use warnings;
use Getopt::Long;


my $usage = "$0 -i annot -p pfam2go -o outfile <--revigo>\n";
my $infile;
my $pfam2go; 
my $outfile;
my $revigo;

GetOptions(
	   'i|infile=s'  => \$infile,
	   'p|pfam2go=s' => \$pfam2go,
	   'o|outfile=s' => \$outfile,
	   'revigo'      => \$revigo,
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

my %pfamids;
while(my $line = <$in>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    my @ids = split("\t",$line);
    if (defined($ids[1])) {    # some lines just have the contig id (first field), which causes problems with split
	my $pfid = $ids[1];
	my $eval = pop(@ids) || "";
	#print $eval,"\n";
	if ($pfid =~ /\,/) {
	    my @many = split(/\,/,$pfid);
	    foreach my $id (@many) {
		$pfamids{$id} = $eval;
	    }
	} else {
	    $pfamids{$pfid} = $eval;
	}
    }
}
close($in);

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
	for my $key (sort(keys %pfamids)) { 
	    if ($key eq $pf) {
		if ($revigo) {        
		    print $out join("\t",($go_term, $pfamids{$key})),"\n";
		} else {
		    print $out join("\t",($pf,$pf_name,$pf_desc,$go_term)),"\n";
		}
		last;
	    }
	}
    }
}

close($pfams);
close($out);

exit;
