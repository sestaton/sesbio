#!/usr/bin/env perl
#--------------------------------------------------------------+
# cnv_LTRdigest_gff3_to_valid.pl - Convert to standard gff3    |
#--------------------------------------------------------------+
# DESCRIPTION:  The LTRdigest GFF3 may work alone but will     |
# likely cause problems when used with annotations from        |
# other sources so, just convert it.                           |
#                                                              |
# AUTHOR: S. Evan Staton                                       |
# CONTACT: statonse<at>uga.edu                                 |
# STARTED: 12.2.10                                             |
# UPDATED: 12.3.10                                             |
#                                                              |
# USAGE: cnv_LTRdigest_gff3_to_valid.pl -i in.gff3 -o out.gff3 |       
#                                                              |
#                                                              |
#--------------------------------------------------------------+
# TODO: Change some of the terms to a standard such as
# that used by the cnv_ltrfinder2gff.pl (DAWGPAWS)
#
# Don't need options right now but will leave them in case I 
# decide to work on changing some of the terms, which could be
# an option. Specifically, change feature "protein_match" to
# "mature_protein_region" and "pfamname" to "Name." These chages
# are necessary to get the ltr_simple glyph to display elements
# correctly. Also add (perhaps optionally) SO terms. 
#
# Convert directory of files instead of one at a time.

use 5.010;
use strict;
use warnings;
use Getopt::Long;

my $infile;
my $outfile;
my $usage = "USAGE: $0 -i in.gff3 -o out.gff3\n";

GetOptions(
           'i|infile=s'           => \$infile,
           'o|outfile=s'          => \$outfile,
	  );

# open the infile or die with a usage statement
if (!$infile || !$outfile) {
    print $usage and exit(1);
}

open my $in, '<', $infile or die "\nERROR: Can't open file: $infile\n";
open my $seq, '<', $infile or die "\nERROR: Can't open file: $infile\n";
open my $out, '>', $outfile or die "\nERROR: Can't open file: $outfile\n";

my @contig = grep {/# (\w)/} <$seq>;
close $seq;
my $contigID = get_contig(@contig);

my @gff = <$in>;
for my $line (@gff) {
    chomp $line;
    if ($line =~ m/^##gff-version /) {
	say $out $line;
    }
    if ($line =~ m/^##sequence-region /) {
	my @seq_region = split(/\s+/, $line);
	say $out join q{ }, $seq_region[0], $contigID, $seq_region[2], $seq_region[3];
    }
    if ($line =~ m/^seq/) {
	my @gff_fields = split(/\s+/,$line);
	my $correctID = $gff_fields[0];
	$correctID =~ s/(\w*)/$contigID/;

	#my @corrected_fields = $correctID."\t".        # Column 1: "seqid"
        #                       $gff_fields[1]."\t".    # Column 2: "source"
	#			$gff_fields[2]."\t".    # Column 3: "type"        ==> repeat_region is not a correct SO term.
	#			$gff_fields[3]."\t".    # Column 4: "start" 
	#			$gff_fields[4]."\t".    # Column 5: "end"
	#			$gff_fields[5]."\t".    # Column 6: "score"
	#			$gff_fields[6]."\t".    # Column 7: "strand"
	#			$gff_fields[7]."\t".    # Column 8: "phase"
	#			$gff_fields[8]."\n";    # Column 9: "attributes"  ==> Need to fix here too. (Parent=repeat_region2)
	
	say $out join "\t", $correctID, @gff_fields[1..8];
    }
}

sub get_contig {
    my @name = @_;
    for (@name) {
	my ($comm, $contig) = split /\s+/, $_;
	return $contig;
    }
}

close $in;
close $out;

#----------------+
# CHANGELOG      |
#----------------+
# 12/2/10
# - Started working on sequence-region line. Not printing
# regions in correct order
#
# 12/2/10
# - Basic conversion is working correctly. 
