#!/usr/bin/perl -w
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

use strict;
use Getopt::Long;

my $infile;
my $outfile;
my $usage = "USAGE: cnv_LTRdigest_gff3_to_valid.pl -i in.gff3 -o out.gff3";

GetOptions(
           "i|infile=s"           => \$infile,
           "o|outfile=s"          => \$outfile,
	  );

# open the infile or die with a usage statement
if ($infile && $outfile) {
    open (INFILE, "<$infile") || print "ERROR: Can't open $infile\n";
    open (SEQNAME, "<$infile");
    open (OUTFILE, ">$outfile");
    
}
else {
    if (!$infile){
        die "\n","ERROR: No infile was given at the command line\n\n",$usage,"\n\n"; 
    }
    if (!$outfile){
        die "\n","ERROR: No outfile was given at the command line\n\n",$usage,"\n\n";
    }
}

my @contig = grep {/# (\w)/} <SEQNAME>;
close(SEQNAME);
my $contigID = get_contig(@contig);

my @gff = <INFILE>;
foreach my $line (@gff) {
   
    chomp $line;
    if ($line =~ m/^##gff-version /) {

	print OUTFILE $line,"\n";
    }
    if ($line =~ m/^##sequence-region /) {
	
	my @seq_region = split(/\s+/, $line);
	
	print OUTFILE $seq_region[0]," ",$contigID," ",$seq_region[2]," ",$seq_region[3],"\n";
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
	
				
	print OUTFILE join("\t",($correctID,$gff_fields[1],$gff_fields[2],$gff_fields[3],$gff_fields[4],$gff_fields[5],
				 $gff_fields[6],$gff_fields[7],$gff_fields[8])), "\n";
	
    }
  
}

sub get_contig {
  
    my @name = @_;
    for(@name)  {
    
    my ($comm, $contig) = split(/\s+/,$_);

    return $contig;
    
    }
    
}

close(INFILE);
close(OUTFILE);

exit;

#----------------+
# CHANGELOG      |
#----------------+
# 12/2/10
# - Started working on sequence-region line. Not printing
# regions in correct order
#
# 12/2/10
# - Basic conversion is working correctly. 
