#!/usr/bin/env perl

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

# lexical vars
my $infile;
my $subject;
my $outfile;
my $merlen;
my $matchlen;
my $toupper;
my $identity;
my $keep;
my $filter;
my $clean;
my $help;

GetOptions(
	   'i|infile=s'             => \$infile,
	   's|subject=s'            => \$subject,
	   'o|outfile=s'            => \$outfile,
	   'm|merlen=i'             => \$merlen,
	   'l|match_len=i'          => \$matchlen,
	   'u|toupper'              => \$toupper,
	   'pid|percent_identity=i' => \$identity,
           'k|keep'                 => \$keep,
           'f|filter'               => \$filter,
	   'c|clean'                => \$clean,
	   'h|help'                 => \$help,
	   );

# help?
usage() and exit(0) if $help;

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
my $db = $sfile."_mkvtreedb";
my $qdb = $ifile."_mkvtreedb";
my $db_exists = file_exists($db);
say "$db_exists ", if $db_exists;
$db .= "_$str" if $db_exists;
if ($keep) {
    my $qdb_exists = file_exists($db);
    say "$qdb_exists ", if $qdb_exists;
    $qdb .= "_$str" if $qdb_exists;
}

my $cwd = getcwd();
my $tmpiname = $ifile."_XXXX";
#my $vmatch_e;
my $o_tmp = File::Temp->new( TEMPLATE => $tmpiname,
			     DIR => $cwd,
			     SUFFIX => $iext,
			     UNLINK => 0);

# set paths to programs used
my $mkvtree = find_prog("mkvtree");
my $vmatch = find_prog("vmatch");
my $vseqselect = find_prog("vseqselect");

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
    say "mkvtree error: $mkvtree_e";
};
 
# 
# Run Vmatch for the query
# 
my ($vmatch_o, $vmatch_e, $vmatch_cmd);
my $vmerSearchSeqnum = $ifile."_vmermatches2";
if ($keep) {
    $vmatch_cmd = "$vmatch -q $infile -l $merlen";
    $vmatch_cmd .= " -identity $identity" if $identity;
    $vmatch_cmd .= " $db";

    try {
	($vmatch_o, $vmatch_e) =
	    capture { system("$vmatch_cmd | grep -v \"^#\" | sed -e 's\/^\[ \\t\]*\/\/g' | perl -lane 'print \$F\[5\]' | sort -u > $vmerSearchSeqnum"); };
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

#find_matches($vmatch_o, $filter) if $keep;
##grep -v "^#" $vmerSearch | sed -e 's/^[ \t]*//g' | perl -lane 'print $F[1]' > $vmerSearchSeqnum

# select the matching subject sequences from the index
#vseqselect -seqnum $vmerSearchSeqnum $db > $vmerSeq # this requires the sequence numbers not the IDs

#
# Return the match list
#
if ($keep) {
    # create the query index
    my ($mkvtree_o, $mkvtree_e);
    try {
	($mkvtree_o, $mkvtree_e) = capture { system("$mkvtree -db $infile -indexname $qdb -dna -allout -v -pl"); };
	if ($mkvtree_e =~ /^mkvtree\:/) {
	    say "\nERROR: $mkvtree_e. Exiting."; exit(1);
	}
    }
    catch {
	say "\nERROR: mkvtree appears to have exited abnormally. Here is the exception: $_\n" and exit(1);
	say "mkvtree output: $mkvtree_o";
	say "mkvtree error: $mkvtree_e";
    };

    # return the sequences matching the subject index
    my ($vseqselect_o, $vseqselect_e);
    try {
	($vseqselect_o, $vseqselect_e) = capture { system("$vseqselect -seqnum $vmerSearchSeqnum $qdb > $o_tmp") };
    }
    catch {
	say "\nERROR: vseqselect appears to have exited abnormally. Here is the exception: $_\n" and exit(1);
	say "vmatch output: $vseqselect_o";
	say "vmatch error: $vseqselect_e";
    };
}

#
# Convert back to Uppercase
#
my $scrSeqCt = 0;
if ($toupper) {
    open my $fas, '<:utf8', $o_tmp;
    open my $out, '>:utf8', $outfile;

    {
	local $/ = '>';
	
	while (my $line = <$fas>) {
	    $line =~ s/>//g;
	    next if !length($line);
	    my ($seqid, @seqs) = split /\n/, $line; 
	    my $seq = join '', @seqs;
	    my $useq = uc($seq);
	    $scrSeqCt++ if defined $seq;
	    say $out join "\n", ">".$seqid, $useq;
	}
    }
    close $fas;
    close $out;

    unlink $o_tmp;
}
else {
    move($o_tmp, $outfile);
    unlink $o_tmp;
}

if ($toupper) {
    my ($qrySeqCt_o, $qrySetCt_e) = capture { system("grep -c '>' $infile"); }; 
    chomp $qrySeqCt_o;

    if ($keep) {
	#my $totSeqScr = $qrySeqCt_o - $scrSeqCt;
        my $totSeqScrPerc = sprintf("%.2f",$scrSeqCt/$qrySeqCt_o * 100);

        say "\n$totSeqScrPerc % ($scrSeqCt","/","$qrySeqCt_o) of reads were screened from ", basename($subject),
	    " in ",basename($infile), ". $scrSeqCt reads written to $outfile.\n";
    }
    else {
	my $totSeqScr = $qrySeqCt_o - $scrSeqCt;
	my $totSeqScrPerc = sprintf("%.2f",$totSeqScr/$qrySeqCt_o * 100);

	say "\n$totSeqScrPerc % ($totSeqScr","/","$qrySeqCt_o) of reads were screened in ", basename($infile), 
	    ". $scrSeqCt reads written to $outfile.\n";
    }
}

# clean up
if ($clean) {
    my ($clean_o, $clean_e) = capture { system("rm ${db}*"); }; 
    if ($keep) {
	my ($clean_o, $clean_e) = capture { system("rm ${qdb}*"); };
    }
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

sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script -i infile -o outfile -s subject -m merlen -l matchlen [-h] [-u] [-pid]

Required:
 -i|infile             :       A multifasta to screen for contamination.
 -s|subject            :       A subject file to use as the target.
 -o|outfile            :       A file to put the screened sequences.

Options:
 -k|keep               :       Keep the sequences matching the subject intead of screening them (Default: no).
 -m|merlen             :       Length to use for matching the index (Default: 20). 
 -l|matchlen           :       Minimum length of the query to keep (Default: 50).
 -u|toupper            :       Print all the sequences as uppercase (Defalut: no).
                               (NB: Vmatch prints sequences lowercase by default,
			        so use this option if uppercase is desired.)
 -pid|percent_identity :       Set the minimum percent identity (integer) threshold (Default: exact match).
 -c|clean              :       Remove the index files created my mkvtree (Default: no).
 -h|help               :       Print a usage statement.

END
}
