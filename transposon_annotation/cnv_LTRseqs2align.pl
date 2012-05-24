#!/usr/bin/perl -w 
#_______________________________________________________________________+
#                                                                       |
# cnv_LTRseqs2align.pl               
#_______________________________________________________________________+
#                                                                       |
# Description: 
#                                                                       |
# Author: S. Evan Staton                                                |
# Contact: statonse<at>uga.edu                                          |
# Started: 1.6.11                                                      |                                                
# Updated:                                                      
#                                                                       |
# Suggested usage:                                                      |

#_______________________________________________________________________+
# TODO: add counting and timing when printing messages to screen

use strict;
use Getopt::Long;

# define vars
my $usage = "USAGE: cnv_LTRseqs2align.pl -i /path/to/dir/of/LTRseq_files <options>

\tRequired:
\t -i|indir      :     The input directory of combined LTR fasta files.

\tOptions:
\t --quiet        :     Be quiet, don't print all of the clustal progress to screen.
\t --phy          :     Create a Phylip formatted file (for use with PAML baseml).
\t --compress     :     Compress the output directory. (All alignment and phylip files
\t                      will be placed in a single uncompressed directory regardless).";

my $in_dir; # = shift or die "\n$usage\n\n";
my $quiet;
my $phylip;
my $compress;

GetOptions(# Required arguments
	   "i|indir=s"         => \$in_dir,                    # indir (for now)
#           "o|outfile=s"        => \$outfile,                 # make this the alignment report       

 	   "quiet"             => \$quiet,
           "phy"               => \$phylip,
	   "compress"          => \$compress,
	  );


# NOT SURE THIS IS THE RIGHT SYNTAX FOR open() in this case
if (!$in_dir) {
    die "\n","ERROR: No input directory was given at the command line\n\n",$usage,"\n\n";
} else {
    if ($in_dir) {
	opendir( DIR, $in_dir ); # || print "Can't open directory:\n$in_dir\n"; 
	my @LTR_fasta_files = grep /\.fasta$/, readdir DIR ;                            # we are working with the .fasta ending
	closedir( DIR );

	if (scalar(@LTR_fasta_files) >= 1) {

	    print "\n\n============= Starting alignments ...\n";
	    chdir($in_dir);

	} else {
	    die "\nERROR: No fasta files were found in the input directory! The must end with .fasta.\n\n",$usage,"\n\n";
	}

	foreach my $LTR_fas (@LTR_fasta_files) {
	    
	    my $aligndir = "LTR_alignment_files";
	    unless  ( -e $aligndir ) {
		mkdir($aligndir);
	    }

	    #if ($phylip) {
		my $phylipdir = "LTR_alignment_phylip_files";
		unless  ( -e  $phylipdir ) {
		    mkdir($phylipdir);
		}
	    #}

	    # need to make separate dirs for each infile
	    my $LTR_outdir = $LTR_fas;
	    $LTR_outdir =~ s/_LTRseqs\.fasta/_aligned/;
	    mkdir($LTR_outdir);
	    chdir($LTR_outdir);
	    
	    system("cp ../$LTR_fas .");

	    # create outfile or else the alignment will be destroyed : NOT NECESSARY
	    my $aligned = $LTR_fas;
	    $aligned =~ s/\.fasta$/\.aln/;

	    
	    if ($quiet) { # be quiet; don't tell me what clustalw is doing
	    	
		#  
		open( OUT , ">$aligned" ) || die "Could not open file: $aligned\n";
		my $clustal_aln_cmd = `clustalw -infile=$LTR_fas 2>&1 > /dev/null`;
	    	print OUT $clustal_aln_cmd;
		close(OUT);
		unlink($clustal_aln_cmd);
		
		system("cp *aln ../$aligndir");

	    } else { # print clustalw output
		if (!$quiet) {
		    my $clustal_aln_cmd = "clustalw -infile=$LTR_fas -outfile=./$aligned";
		    system($clustal_aln_cmd);
		    system("cp *aln ../$aligndir");
		}
	    }
	    
	    # create a file to use with PAML baseml
	    if ($phylip) {
		
		#my $clustal_phy_cmd = "clustalw -infile=$LTR_fas -output=phylip";
		#system($clustal_phy_cmd);
		#
		# The clusalw phylip format is not correct for PAML so, I'll use bioperl 
		# to get the spacing right and then add an I to tell PAML it's an interleaved file.
		format_phylip($aligned);
		
		system("cp *dnd *phy ../$phylipdir");
		
	    }
	    unlink($LTR_fas);       # remove the copy; it would be safer to remove a link, but it works!
	    chdir("..");

	    if ($compress) {

		my $compressdir = $LTR_outdir.".tar.gz";
		system("tar -czf $compressdir $LTR_outdir");
		system("rm -rf $LTR_outdir");

	    }
	}
    }
}

sub format_phylip {

    use Bio::AlignIO;
    my $clustalw_aligned = shift;
    
    my $inter = $clustalw_aligned;
    $inter =~ s/\.aln$/\.phy/;
    my $tmp = $inter;
    $tmp .= "-tmp";

    #print "\nchecking file creation...",$tmp,"\n";
    
    my $aln_in = Bio::AlignIO->new(-file   => $clustalw_aligned,
				   -format => 'clustalw');
    
    my $aln_out = Bio::AlignIO->new(-file   => ">$inter",
				    -format => 'phylip'); 
				    #-interleave => 0);
    
    while (my $aln = $aln_in->next_aln()) {
	$aln_out->write_aln($aln);
    }
    
    open(PHY,  $inter) || die "ERROR: Could not open file: $inter\n";
    open(TMP, ">$tmp") || die "ERROR: Could not open file: $tmp\n";
    while (<PHY>) {
	chomp;
	if (/^\s\d\s\d+/) {
	    print TMP $_." "."I"."\n";
	} else {
	    print TMP $_."\n";
	}
    }
    close(TMP);
    close(PHY);
    system("mv $tmp $inter");
    unlink($tmp);

}
