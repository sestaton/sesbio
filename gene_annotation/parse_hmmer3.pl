#!/usr/bin/perl -w
#_________________________________________________________________+
#                                                                 |
# parse_hmmer3.pl - Get report or sequences from hmmer3
#_________________________________________________________________+
#                                                                 |
# Description: Parse hmmscan (formerly hmmpfam in hmmer2) and     |
# hmmersearch. To get the sequences, alignments must be present   |
# in the report.                                                  |
#                                                                 |
#                                                                 |
# Author: S. Evan Staton                                          |
# Contact: statonse<at>uga.edu                                    |
# Started: 11/22/10                                               |
# Updated: 11/23/10      
#                                                                 |
# Suggested usage: perl parse_hmmer3.pl -i infile <-o outfile> \  |
#                                 <-s seqfile> <-h hmmname>       | 
#                                                                 | 
# The options <-s seqfile>,<-o outfile>,<-h hmmname> can be used  |
# all together on separately.                                     |                             
#_________________________________________________________________+
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
my $usage = "USAGE: parse_hmmer3.pl -i infile <-o outfile> <-s seqfile> <-h hmmname>

\tReturn a tab-delimited report to the given outfile in the form of: 
\tquery query_length number_of_hits hit_name hit_score hit_significance hsp_length hsp_query_start hsp_query_end hsp_hit_start hsp_hit_end

\tIf the -s option is given the HSP query sequences are written to <seqfile> in fasta format.

\tIf the -h option is given the HSP query sequences matching the specified hmm name are written to hmmname.fasta.
\tCurrently only one hmm name can be evaluated at a time, this will change.\n";

#----------+
# OPTIONS  |
#----------+
GetOptions(#required arguments at this time
	   "i|infile=s"      => \$infile,
	   "o|outfile=s"     => \$outfile,
	   "s|seq=s"         => \$seqfile,
	   "h|hmmname=s"     => \$hmmname,
	  );

# open the infile or die with a usage statement
if ($infile && $outfile) {
    open (OUTFILE, ">$outfile");
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
    
#esle {
if (!$infile){
        die "\n","ERROR: No infile was given at the command line.\n\n",$usage,"\n"; 
    }
    # Might just want to get the sequence and not write the report everytime
    #if (!$outfile){
    #    die "\n","ERROR: No outfile was given at the command line.\n\n",$usage,"\n";
    #}
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

#-----------+
# CHANGELOG |
#-----------+
# 11/23/10
# - added option to get the sequence matching
# a specific hmm name
