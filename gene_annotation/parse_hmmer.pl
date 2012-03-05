#!/usr/bin/perl -w

#TODO: Take multiple hmm names at the command line
#
use strict;
use lib qw(/home/jmblab/statonse/apps/perlmod/bioperl-live-latest/); # the bioperl-hmmer3 location has to be specified
use Bio::SearchIO;                                                   # because we want hmmer3.pm, not hmmer.pm
use Getopt::Long;                                                    # ///Actually, bioperl-hmmer3 won't work alone, must get latest bioperl-live

#--------------- +
# VARIABLE SCOPE |
#----------------+
my $infile; 
my $outfile;
my $seqfile;
my $hmmname;
my $hmmer_version;

#----------+
# OPTIONS  |
#----------+
GetOptions(#required arguments at this time
	   'i|infile=s'      => \$infile,
	   'o|outfile=s'     => \$outfile,
	   's|seq=s'         => \$seqfile,
	   'h|hmmname=s'     => \$hmmname,
	   'hv|hmmer_ver=i'  => \$hmmer_version,
	  );

# open the infile or die with a usage statement
if (!$infile && !$outfile) {
    print "\nERROR: Command line not parsed correctly. Exiting.\n";
    &usage();
    exit(1);
}


if ($seqfile) {
    #my $seqfile = $outfile . "_SEQ";

    open (SEQ, ">$seqfile");
}
if ($hmmname) {
    	$hmmname = $hmmname.".fasta";
	open (HMMSEQ, ">$hmmname");
}
#if ($hmmname =~ m/\,/) {
#    my @hmms = split(",", $hmmname);
#    foreach $hmmname (@hmms) {
#	$hmmname = $hmmname.".fasta";
#	open (HMMSEQ, ">$hmmname");
#    } 
#} else {
#	$hmmname = $hmmname.".fasta";
#	open (HMMSEQ, ">$hmmname");
#}
    
# get the hmmer3 report
my $hmmer3in = Bio::SearchIO->new(-file => "$infile", -format => 'hmmer3');

if ($outfile) {
    print OUTFILE join("\t", ("#query", "query_length", "number_of_hits", "hit_name", "hit_score", "hit_significance", "hsp_length", "hsp_query_start", "hsp_query_end", "hsp_hit_start", "hsp_hit_end")), "\n";
}
while( my $result = $hmmer3in->next_result() ) {
    
    my $query      = $result->query_name();
    my $qlen       = $result->query_length();
    my $num_hits   = $result->num_hits();
       
    while( my $hit = $result->next_hit() ) {
	
	my $hitid    = $hit->name();
	my $score    = $hit->raw_score();
	my $signif   = $hit->significance();
	
	while( my $hsp = $hit->next_hsp ) {
	   
	    my $hsplen    = $hsp->length('total');
	    my $hstart    = $hsp->start('hit');
	    my $hstop     = $hsp->end('hit');
	    my $qstart    = $hsp->start('query');
	    my $qstop     = $hsp->end('query');
	    my $qstring   = $hsp->query_string;
	    my $hstring   = $hsp->hit_string;

	    my $seqid = ">".$query."|".$hitid."_".$qstart."-".$qstop;
	    # 
	    if ($outfile) {
		print OUTFILE join("\t", ($query, $qlen, $num_hits, $hitid, $score, $signif, $hsplen, $qstart, $qstop, $hstart, $hstop)), "\n";
		
	    }#print SEQ join("\n", ($seqid, $qstring));
	    if ($seqfile) {
		print SEQ $seqid."\n".$qstring."\n";
		
	    }
	    if ($hmmname) {
		if ($hmmname =~ m/$hitid/) {
		    print HMMSEQ $seqid."\n".$qstring."\n";
		}
	    }
	}
    }
}
close(OUTFILE);
close(SEQ);
close(HMMSEQ);
exit;

#
# subs
#
sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script -i infile <-o outfile> <-s seqfile> <-h hmmname>                                                                                        
                                                                                                                                                                           Return a tab-delimited report to the given outfile in the form of:                                                                                                        
query query_length number_of_hits hit_name hit_score hit_significance hsp_length hsp_query_start hsp_query_end hsp_hit_start hsp_hit_end                                 

Required:
    -i|infile     :     HMMER report to parse.
    -o|outfile    :     File to hold the parsed results.
    -hv|hmmer_ver :     HMMER version used to create report (Can be either 2 or 3).

Options:
    -s|seqfile    :     Write all query sequence matches to domains to a file.
    -h|hmmname    :     Print all matches to a particular domain name.


END
}

#-----------+
# CHANGELOG |
#-----------+
# 11/23/10
# - added option to get the sequence matching
# a specific hmm name
