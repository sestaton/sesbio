#!/usr/bin/env perl

use utf8;
use v5.12;
use strict;
use warnings;
use warnings FATAL => "utf8";
#use charnames qw(:full :short); # not actually using this right now
use open qw(:std :utf8);
use autodie qw(open);
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

my ($sfile, $sdir, $sext) = fileparse($subject, qr/\.[^.]*/);
my ($ifile, $idir, $iext) = fileparse($infile, qr/\.[^.]*/);
my $db = $sfile."_mkvtreedb";
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
my $scrSeqCt = 0;
if (defined $toupper) {
    open(my $fas, '<:utf8', $o_tmp);
    open(my $out, '>:utf8', $outfile);

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
    close($fas);
    close($out);

    unlink($o_tmp);
}
else {
    move($o_tmp, $outfile);
    unlink($o_tmp);
}

my $qrySeqCt = `grep -c '>' $infile`; chomp($qrySeqCt);
my $totSeqScr = $qrySeqCt - $scrSeqCt;
my $totSeqScrPerc = sprintf("%.2f",$totSeqScr/$qrySeqCt * 100);

say "\n$totSeqScrPerc % ($totSeqScr","/","$qrySeqCt) of reads were screened in $infile. $scrSeqCt reads written to $outfile.\n";

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
    else {
	say "\nERROR: Could not find $prog PATH. Exiting."; exit(1);
    }    
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
