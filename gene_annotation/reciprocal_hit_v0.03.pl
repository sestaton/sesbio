#!/usr/bin/perl -w

use strict;
use File::Basename;
use Getopt::Long;

my $query_blast;
my $subj_blast;
my $outfile;
my $help;

GetOptions(
	   'q|query=s'   => \$query_blast,
	   's|subject=s' => \$subj_blast,
	   'o|outfile=s' => \$outfile,
	   'h|help=s'    => \$help,
	   );

usage() and exit(0) if $help;

if (!$query_blast || !$subj_blast || !$outfile) {
    usage();
    exit(1);
}

my $recip_hit = 0;
my %qhash;

open(my $query, '<', $query_blast) or die "\nERROR: Could not open file: $query_blast\n";
open(my $subj, '<', $subj_blast) or die "\nERROR: Could not open file: $subj_blast\n";
open(my $out, '>', $outfile) or die "\nERROR: Could not open file: $outfile\n";

while (my $qline = <$query>) {
    chomp $qline;
    next if $qline =~ /^#/;
    my @query_fields = split(/\t/,$qline);
    my $hkey = join(",",($query_fields[0],$query_fields[2]));
    $qhash{$hkey} = join("\t",($query_fields[3],$query_fields[4],$query_fields[5]));
}

print $out "Hit\tHit\tHSP_len\tPID\tPCov\tHSP_len\tPID\n";

while (my $sline = <$subj>) {
    chomp $sline;
    next if $sline =~ /^#/;
    my @subj_fields = split(/\t/,$sline);
    while( my ($qid, $qhit) = each(%qhash)) {
	my ($qq, $qh) = split(",",$qid);
	if ($qq =~ /$subj_fields[1]/ && $qh =~ /$subj_fields[0]/) {
	    $recip_hit++;
	    print $out "HIT $recip_hit : ", join("\t",$qq,$subj_fields[0],$qhit,$subj_fields[3],$subj_fields[2]),"\n";
	    #print $out "HIT $recip_hit : ", join("\t",@subj_fields),"\n";
	}
    }
}

close($query);
close($subj);
close($out);

print "\nWe found $recip_hit !\n\n";


exit;

#
# subs
#
sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -q blast_res_i_to_j.bln -s blast_res_j_to_i.bln -o blast_result [-h]

Required:
    -q|query         :    Tab-delimited BLAST file for query.
    -s|subject       :    Tab-delimited BLAST file for subject. 
                          (NB: It makes no difference which order you put subject and query, 
			   it only affects the order the results are printed.)
    -o|outfile       :    File name to write the reciprocal BLAST hit results to.

Options:
    -h|help          :    Print a usage statement.

END
}
