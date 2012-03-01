#!/usr/bin/perl -w

=head1 NAME 
                                                                       
 clean_multifasta.pl - remove characters that break EMBOSS and other programs 

=head1 SYNOPSIS    

 clean_multifasta.pl -i weirdseqs.fas -o cleanseqs.fas 

=head1 DESCRIPTION
                                                                   
 Spaces, semicolons, commas, parenthesis and other characters cause the
 fasta header to not be parsed correctly when using EMBOSS applications.
 This leads to sequences being skipped by the application and sequence
 IDs being modified without intent. This script also replaces non-ATGCN
 characters with an "N."

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The multi-fasta file to be cleaned.

=item -o, --outfile

A file to place the cleaned sequences.

=back

=head1 OPTIONS

=over 2

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut      

#
# Includes
#
use strict;
use File::Basename;
use Time::HiRes qw(gettimeofday);
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;

# 
# Variables with scope
#
my $infile;
my $outfile;
my $help;
my $man;

GetOptions(# Required
	   'i|infile=s'   => \$infile,
	   'o|outfile=s'  => \$outfile,
	   # Options
	   'h|help'       => \$help,
	   'm|man'        => \$man,
	   );

pod2usage( -verbose => 2 ) if $man;

if ($help) {
    &usage();
    exit(0);

if (!$infile || !$outfile) {
    print "\nERROR: No input was given.\n";
    &usage();
    exit(1);
}

my $seq_in  = Bio::SeqIO->new( -format => 'fasta', 
                               -file => $infile); 

open( my $SEQSOUT , '>', $outfile ) 
    or die "\nERROR: Could not open file: $outfile\n";

# counters
my $t0 = gettimeofday();
my $fasnum = 0;
my $headchar = 0;
my $non_atgcn = 0;

my %seqhash;
my %cleanseqhash;

while(my $seqs = $seq_in->next_seq()) {
    $fasnum++;
    my $id = $seqs->id."_".$seqs->desc;
    $seqhash{$id} = $seqs->seq;   
}

if ( $fasnum >= 1 ) {
    print "\n========== Cleaning up $fasnum fasta files ...\n";
} else {
    die "\nERROR: No sequences found!\n";
}

while (my ($seqname, $seq) = each(%seqhash)) {

    if ($seqname =~ m/\s+|\;|\:|\(|\)|\./g) {
	$headchar++;
	$seqname =~ s/\s/\_/g;
	$seqname =~ s/\;/\_/g;
	$seqname =~ s/\:/\_/g;
	$seqname =~ s/\(/\_/g;
	$seqname =~ s/\)/\_/g;
	$seqname =~ s/\./\_/g;
	$seqname =~ s/\|/\_/g;
	if ($seqname =~ m/\_+/g) {
	    $seqname =~ s/\_+/\_/g;
	}
	# TODO: shorten the header, optionally
    }

    if ($seq =~ m/\r|\n/) { # Mac
	$seq =~ s/\r?//;    # PC
	$seq =~ s/\n//;     # Unix
    }
 
    my @nt = split(//,$seq);
    foreach my $base (@nt) {
	if ($base !~ m/A|C|G|T|N|a|c|g|t|n/) {
	    $non_atgcn++;
	    $base = "N";
	}
    }       
    my $dna = join('',@nt);
    $dna =~ s/\s//g;

    my $nucleic_bc = ($dna =~ tr/NACGTacgtn//);
    my $nonnucleic = (length($dna) - $nucleic_bc);
    $dna =~ s/(.{60})/$1\n/gs;        

    print $SEQSOUT ">"."$seqname\n"."$dna\n";
    $cleanseqhash{$seqname} = $dna;
}

close($SEQSOUT);
my $cleanfasnum = scalar(keys %cleanseqhash);

my $t1 = gettimeofday();
my $elapsed = $t1 - $t0;
my $time = sprintf("%.2f",$elapsed);

print "\n========== Done. $fasnum sequences read. $cleanfasnum sequences cleaned in $time seconds.\n"; 
print "========== $non_atgcn Non-ATGCN characters changed in sequence. $headchar characters changed in headers.\n\n";

exit;

#
# subs
#
sub usage {
  my $script = basename($0);
  print STDERR <<END
USAGE: $script -i seqsin.fas -o filterseqs.fas 

Required:
    -i|infile    :    Fasta file of reads or contigs to filter.
    -o|outfile   :    File to place the filtered reads or contigs.
    
Options:
    -h|help      :    Print usage statement.
    -m|man       :    Print full documentation.
END
}
