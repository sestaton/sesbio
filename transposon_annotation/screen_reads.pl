#!/usr/bin/env perl

use v5.12;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Temp;
use Capture::Tiny qw(:all);
use File::Copy;
use Cwd;

# lexical vars
my $infile;
my $subject;
my $outfile;
my $merlen;
my $matchlen;
my $toupper;
my $help;

GetOptions(
	   'i|infile=s'    => \$infile,
	   's|subject=s'   => \$subject,
	   'o|outfile=s'   => \$outfile,
	   'm|merlen=i'    => \$merlen,
	   'l|match_len=i' => \$matchlen,
	   'u|toupper'     => \$toupper,
	   'h|help'        => \$help,
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

my ($ifile, $idir, $iext) = fileparse($infile, qr/\.[^.]*/);
my $db = $ifile."_mkvtreedb";
my $cwd = getcwd();
my $tmpiname = $ifile."_XXXX";
#my $vmatch_e;
my $o_tmp = File::Temp->new( TEMPLATE => $tmpiname,
			     DIR => $cwd,
			     SUFFIX => $iext,
			     UNLINK => 0);

# set paths to programs used
my $seqret = find_prog("seqret");
my $mkvtree = find_prog("mkvtree");
my $vmatch = find_prog("vmatch");

#
# Create the index
#
my ($mkvtree_o, $mkvtree_e) = capture { system("$mkvtree -db $subject -indexname $db -dna -allout -v -pl"); };

# 
# Run Vmatch for the query
# 
my ($vmatch_o, $vmatch_e) = capture { system("$vmatch -s -showdesc 0 -qnomatch $matchlen -q $infile -l $merlen $db | grep -v \"^#\" | sed 's\/\\s\.*\/\/' > $o_tmp"); };

#
# Convert back to Uppercase
#
if (defined $toupper) {
    my ($seqret_o, $seqret_e) = capture { system("$seqret -sequence $o_tmp -supper1 -outseq $outfile -auto"); };
    unlink($o_tmp);
}
move($o_tmp, $outfile);
unlink($o_tmp);

my $qrySeqCt = `grep -c '>' $infile`; chomp($qrySeqCt);
my $scrSeqCt = `grep -c '>' $outfile`; chomp($scrSeqCt);
my $totSeqScr = $qrySeqCt - $scrSeqCt;
my $totSeqScrPerc = sprintf("%.2f",$totSeqScr/$qrySeqCt * 100);

say "\n$totSeqScrPerc % ($scrSeqCt","/","$qrySeqCt) of reads were screened in $infile. $scrSeqCt reads written to $outfile.";

# clean up
system("rm ${db}*");

#
# Subs
#
sub find_prog {
    my $prog = shift;
    my ($path, $err) = capture { system("which $prog"); };
    chomp($path);
    
    if ($path =~ /$prog$/) {
	return $path;
    }

    if ($path !~ /$prog$/) {
	say "Couldn\'t find $prog in PATH. Will keep looking.";
	$path = "/usr/local/emboss/latest/bin/$prog";           # path at zcluster
    }

    # Instead of just testing if the program exists and is executable 
    # we want to make sure we have permissions, so we try to 
    # invoke programs and examine the output. 
    my ($prog_path, $prog_err) = capture { system("$path --help"); };

    given ($prog_err) {
	when (/Version\: EMBOSS/) { say "Using $prog located at: $path"; }
	when (/^-bash: \/usr\/local\/emboss\/bin\/$prog\: No such file or directory$/) { die "Could not find $prog. Exiting.\n"; }
	when ('') { die "Could not find getorf. Exiting.\n"; }
	default { die "Could not find $prog. Trying installing EMBOSS or adding it's location to your PATH. Exiting.\n"; }
    }
    return $path;
}
sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script -i infile -o outfile -s subject -m len -l len [-h] [-u]

Required:
 -i|infile     :       A multifasta to screen for contamination.
 -s|subject    :       A subject file to use as the target.
 -o|outfile    :       A file to put the screened sequences.

Options:
 -l|merlen     :       Length to use for matching the index (Default: 20). 
 -m|matchlen   :       Minimum length of the query to keep (Default: 50).
 -u|toupper    :       Print all the sequences as uppercase (Defalut: no).
                       (NB: Vmatch prints sequences lowercase by default,
			    so use this option if uppercase is desired.)
 -h|help       :       Print a usage statement.

END
}
