#!/usr/bin/env perl

##TODO: 1) Filter matches by length,       DONE 
##      2) Filter matches by repeat ratio  (This could generate artifacts since there will be repetitive regions of most targets)

use utf8;
use v5.12;
use strict;
use warnings;
use warnings FATAL => "utf8";
use open qw(:std :utf8);
use autodie qw(open);
use Getopt::Long;
use File::Basename;
use File::Temp;
use Try::Tiny;
use Capture::Tiny qw(:all);
use File::Copy;
use POSIX qw(strftime);
use IPC::System::Simple qw(system);
use Cwd;
use Data::Dump qw(dd);

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
my $qdb = $ifile."_mkvtreedb";
my $db_exists = file_exists($db);
say "$db_exists ", if $db_exists;
$db .= "_$str" if $db_exists;

#if ($keep) {
#    my $qdb_exists = file_exists($db);
#    say "$qdb_exists ", if $qdb_exists;
#    $qdb .= "_$str" if $qdb_exists;
#    $seqmap = fas2idmap($infile);
#    #dd $seqmap; exit;
#}

my $cwd = getcwd();
my $tmpiname = $ifile."_XXXX";
my $o_tmp = File::Temp->new( TEMPLATE => $tmpiname,
			     DIR      => $cwd,
			     SUFFIX   => $iext,
			     UNLINK   => 0);

# set paths to programs used
my $mkvtree = find_prog("mkvtree");
my $vmatch = find_prog("vmatch");
my $samtools = find_prog("samtools");

#
# Create the index
#
my ($mkvtree_o, $mkvtree_e);
try {
    ($mkvtree_o, $mkvtree_e) = capture { system("$mkvtree -db $subject -indexname $db -dna -allout -v -pl"); };
    if ($mkvtree_e =~ /^mkvtree\:/) {
	say "\nERROR: $mkvtree_e. Exiting."; exit(1);
    }
}
catch {
    say "\nERROR: mkvtree appears to have exited abnormally. Here is the exception: $_\n" and exit(1);
    say "mkvtree output: $mkvtree_o";
    say "mkvtree error:  $mkvtree_e";
};
 
# 
# Run Vmatch for the query
# 
my ($vmatch_o, $vmatch_e, $vmatch_cmd);
my $vmerSearchSeqnum = $ifile."_vmermatches_$str";
if ($keep) {
    $vmatch_cmd = "$vmatch -showdesc 0 -q $infile -l $merlen";
    $vmatch_cmd .= " -identity $identity" if $identity;
    $vmatch_cmd .= " $db";

    try {
	($vmatch_o, $vmatch_e) = capture { system("$vmatch_cmd > $vsearch_out"); };
	system("grep -v \"^#\" $vsearch_out | sed -e 's\/^\[ \\t\]*\/\/g' | perl -lane 'print \$F\[5\]' | sort -u > $vmerSearchSeqnum");
	if ($vmatch_e =~ /^vmatch\:/) {
	    say "\nERROR: $vmatch_e. Exiting."; exit(1);
	}
    }
    catch {
	say "\nERROR: Vmatch appears to have exited abnormally. Here is the exception: $_\n" and exit(1);
	say "vmatch output: $vmatch_o";
	say "vmatch error:  $vmatch_e";
    };
}
else {
    $vmatch_cmd = "$vmatch -s -showdesc 0 -qnomatch $matchlen -q $infile -l $merlen";
    $vmatch_cmd .= " -identity $identity" if $identity;
    $vmatch_cmd .= " $db";

    try {
	($vmatch_o, $vmatch_e) = capture { system("$vmatch_cmd | grep -v \"^#\" | sed 's\/\\s\.*\/\/' > $o_tmp"); };
	if ($vmatch_e =~ /^vmatch\:/) {
	    say "\nERROR: $vmatch_e. Exiting."; exit(1);
	}
    }
    catch {
	say "\nERROR: Vmatch appears to have exited abnormally. Here is the exception: $_\n" and exit(1);
	say "vmatch output: $vmatch_o";
	say "vmatch error: $vmatch_e";
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
    my ($fa_o, $fa_e);
    try {
	($fa_o, $fa_e) = capture { system("/usr/local/kent/latest/bin/faSomeRecords $infile $vmerSearchSeqnum $o_tmp") };
    }
    catch {
	say "\nERROR: faSomeRecords appears to have exited abnormally. Here is the exception: $_\n" and exit(1);
	say "faSomeRecords output: $fa_o";
	say "faSomeRecords error:  $fa_e";
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
my ($qrySeqCt_o, $qrySetCt_e) = capture { system("grep -c '>' $infile"); }; 
chomp $qrySeqCt_o;

my ($scrSeqCt_o, $scrSeqCt_e);
unless ($keep && $filter) {
    ($scrSeqCt_o, $scrSeqCt_e) = capture { system("grep -c '>' $o_tmp"); };
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
    my ($clean_o, $clean_e) = capture { system("rm ${db}*"); }; 
    unlink $vsearch_out;
    unlink $o_tmp;
}

#
# Subs
#
sub find_prog {
    my $prog = shift;
    my ($path, $err) = capture { system("which $prog"); };
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

sub file_exists {
    # http://stackoverflow.com/a/8584761
    my ($qfn) = @_;
    my $rv = -e $qfn;
    die "Unable to determine if file exists: $!"
	if !defined($rv) && !$!{ENOENT};
    return $rv;
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
