#!/usr/bin/env perl

##TODO: 1) Filter matches by repeat ratio 
           - This could generate artifacts since there will be repetitive regions of most targets

##NB: this method with vmatch returns the best match, often resulting in only short, identical matches.

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;
use File::Basename;
use File::Temp;
use Try::Tiny;
use File::Copy;
use POSIX qw(strftime);
use IPC::System::Simple qw(system capture);
use Cwd;
#use Data::Dump qw(dd);

# lexical vars
my $infile;
my $subject;
my $outfile;
my $merlen;
my $matchlen;
my $identity;
my $keep;
my $filter;
my $clean;
my $help;
my $seqmap;

GetOptions(
	   'i|infile=s'             => \$infile,
	   's|subject=s'            => \$subject,
	   'o|outfile=s'            => \$outfile,
	   'm|merlen=i'             => \$merlen,
	   'l|match_len=i'          => \$matchlen,
	   'pid|percent_identity=i' => \$identity,
           'k|keep'                 => \$keep,
           'f|filter_hitrange'      => \$filter,
	   'c|clean'                => \$clean,
	   'h|help'                 => \$help,
	   );

# help?
usage() and exit(0) if $help;

if ($filter && !$keep) {
    say "\nERROR: Filtering the match parts of each sequence is only possible when keeping matches. Exiting.\n";
    usage();
    exit(1);
}

# check @ARGV
if (!$infile || !$subject || !$outfile) {
    usage();
    exit(1);
}

# set defaults
$merlen //= 20;
$matchlen //= 50;

my $str = POSIX::strftime("%m_%d_%Y_%H_%M_%S", localtime);
my ($sfile, $sdir, $sext) = fileparse($subject, qr/\.[^.]*/);
my ($ifile, $idir, $iext) = fileparse($infile, qr/\.[^.]*/);
my $vsearch_out = $ifile."_".$sfile."_$str".".vmatch";
my $db = $sfile."_mkvtreedb";
$db .= "_$str";

my $cwd = getcwd();
my $tmpiname = $ifile."_XXXX";
my $o_tmp;
if ($keep) {
    $o_tmp = File::Temp->new( TEMPLATE => $tmpiname,
			      DIR      => $cwd,
			      SUFFIX   => $iext,
			      UNLINK   => 0);
}

# set paths to programs used
my $mkvtree = find_prog("mkvtree");
my $vmatch  = find_prog("vmatch");
my $records = find_prog("faSomeRecords");

#
# Create the index
#
try {
    my @mkvout = capture([0..5], "$mkvtree -db $subject -indexname $db -dna -allout -v -pl");
    for my $mko (@mkvout) {
	if ($mko =~ /^mkvtree\:/) {
	    say "\nERROR: $mko. Exiting.";
	    exit(1);
	}
    }
}
catch {
    say "\nERROR: mkvtree appears to have exited abnormally. Here is the exception: $_\n";
    exit(1);
};
 
# 
# Run Vmatch for the query
#
my $vmatch_cmd; 
my $vmerSearchSeqnum = $ifile."_vmermatches_$str";
if ($keep) {
    $vmatch_cmd = "$vmatch -p -d -showdesc 0 -q $infile -l $merlen";
    $vmatch_cmd .= " -identity $identity" if $identity;
    $vmatch_cmd .= " $db";

    try {
	my @vmatchout = capture([0..5], "$vmatch_cmd > $vsearch_out");
	for my $vmline (@vmatchout) {
	    if ($vmline =~ /^vmatch\:/) {
		say "\nERROR: $vmline. Exiting."; exit(1);
	    }
	} 
	system([0..5], "grep -v \"^#\" $vsearch_out | sed -e 's\/^\[ \\t\]*\/\/g' | perl -lane 'print \$F\[5\]' | sort -u > $vmerSearchSeqnum");
    }
    catch {
	say "\nERROR: Vmatch appears to have exited abnormally. Here is the exception: $_\n";
	exit(1);
    };
}
else {
    $vmatch_cmd = "$vmatch -p -d -s -showdesc 0 -qnomatch $matchlen -q $infile -l $merlen";
    $vmatch_cmd .= " -identity $identity" if $identity;
    $vmatch_cmd .= " $db";

    try {
	my @vmatchout = capture([0..5], "$vmatch_cmd | grep -v \"^#\" | sed 's\/\\s\.*\/\/' > $outfile");
	for my $vmline (@vmatchout) {
	    if ($vmline =~ /^vmatch\:/) {
		say "\nERROR: $vmline. Exiting.";
		exit(1);
	    }
	}
    }
    catch {
	say "\nERROR: Vmatch appears to have exited abnormally. Here is the exception: $_\n";
	exit(1);
    };
}

#
# Return the match list
#
my %match_range;

if ($keep) {
    # capture the match range
    open my $vsearch, '<', $vsearch_out;

    while (<$vsearch>) {
        chomp;
        next if /^#/;                                                                                                                                           
        s/^\s+//g;
        my @f = split;
        unless (exists $match_range{$f[5]}) {
            $match_range{$f[5]} = join "|", $f[4], $f[6];
        }
    }
    close $vsearch;

    # return the sequences matching the subject index
    try {
	system([0..5], "$records $infile $vmerSearchSeqnum $o_tmp");
    }
    catch {
	say "\nERROR: faSomeRecords appears to have exited abnormally. Here is the exception: $_\n";
	exit(1);
    };
}

#dd \%match_range; exit;

#
# Convert back to uppercase, select match range if keeping matches
#
my ($scrSeqCt, $validscrSeqCt) = (0, 0);
if ($keep && $filter) {
    open my $fas, '<', $o_tmp;
    open my $out, '>', $outfile;
    
    {
	local $/ = '>';
    
	while (my $line = <$fas>) {
	    $line =~ s/>//g;
	    next if !length($line);
	    my ($seqid, @seqs) = split /\n/, $line; 
	    my $seq = join '', @seqs;
	    my $useq = uc($seq);
	    $scrSeqCt++ if defined $seq;
	    $seqid =~ s/\s.*//;
	    if (exists $match_range{$seqid}) {
		my ($match_len,$match_offset) = split /\|/, $match_range{ $seqid };
		if ($match_len >= $matchlen) {
		    $validscrSeqCt++;
		    my $seq_match = substr $useq, $match_offset, $match_len;
		    say $out join "\n", ">".$seqid, $seq_match;
		}
	    }
	}
    }
    close $fas;
    close $out;
}

#
# Compute screening results
#
my $qrySeqCt_o = capture([0..5], "grep -c '>' $infile"); 
chomp $qrySeqCt_o;

my ($scrSeqCt_o, $scrSeqCt_e);
unless ($keep && $filter) {
    $scrSeqCt_o = capture([0..5], "grep -c '>' $o_tmp");
    chomp $scrSeqCt_o;
}

if ($keep && $filter) {
    my $totSeqScrPerc = sprintf("%.2f",$validscrSeqCt/$qrySeqCt_o * 100);
    
    say "\n$totSeqScrPerc % ($validscrSeqCt","/","$qrySeqCt_o) of reads were screened from ", basename($subject),
	    " in ",basename($infile), ". $validscrSeqCt reads written to $outfile.\n";
}
elsif ($keep && !$filter) {
    my $totSeqScrPerc = sprintf("%.2f",$scrSeqCt_o/$qrySeqCt_o * 100);

    say "\n$totSeqScrPerc % ($scrSeqCt_o","/","$qrySeqCt_o) of reads were screened from ", basename($subject),
    " in ",basename($infile), ". $scrSeqCt_o reads written to $outfile.\n";
}
else {
    my $totSeqScr = $qrySeqCt_o - $scrSeqCt_o;
    my $totSeqScrPerc = sprintf("%.2f",$totSeqScr/$qrySeqCt_o * 100);
    
	say "\n$totSeqScrPerc % ($totSeqScr","/","$qrySeqCt_o) of reads were screened in ", basename($infile), 
    ". $scrSeqCt_o reads written to $outfile.\n";
}

# clean up
if ($clean) {
    unlink glob "${db}*";
    unlink $vsearch_out if -e $vsearch_out;
    unlink $o_tmp if -e $o_tmp;
    unlink $vmerSearchSeqnum if -e $vmerSearchSeqnum;
}

#
# methods
#
sub find_prog {
    my $prog = shift;
    my $path = capture([0..5], "which $prog");
    chomp $path;
    
    ## given/when moved to experimental in v5.18
    #given ($path) {
	#when ($path =~ /$prog$/) { return $path; }
	#default { say "\nERROR: Could not find $prog PATH. Exiting."; exit(1); }    
    #}
    if ($path =~ /$prog$/) {
	return $path;
    }
    else {
	say "\nERROR: Could not find $prog in PATH. Exiting.";
	exit(1);
    }
}

sub fas2idmap {
    my ($infile) = @_;

    my %seqmap;
    open my $in, '<', $infile;
    
    my $seqct = 0;

    local $/ = '>';
    
    while (my $line = <$in>) {
	$line =~ s/>//g;
	next if !length($line);
	my ($seqid, @seqs) = split /\n/, $line;
	my $seq = join '', @seqs;
	$seqid =~ s/\s.*//;
	$seqmap{$seqid} = $seqct;
	$seqct++ if defined $seq;
    }
    close $in;

    return \%seqmap;
}

sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script -i infile -o outfile -s subject -m merlen -l matchlen [-h] [-u] [-pid] [-k] [f] [-c]

Required:
 -i|infile             :       A multifasta to screen for contamination.
 -s|subject            :       A subject file to use as the target.
 -o|outfile            :       A file to put the screened sequences.

Options:
 -k|keep               :       Keep the sequences matching the subject intead of screening them (Default: no).
 -f|filter_hitrange    :       Keep only the match parts of each sequence; can only be used with [--keep](Default: no).
 -m|merlen             :       Length to use for matching the index (Default: 20). 
 -l|matchlen           :       Minimum length of the query to keep (Default: 50).
 -pid|percent_identity :       Set the minimum percent identity (integer) threshold (Default: exact match).
 -c|clean              :       Remove the index files created my mkvtree (Default: no).
 -h|help               :       Print a usage statement.

END
}
