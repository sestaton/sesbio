#!/usr/bin/perl -w


use strict;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;

# input sequence to search
# input sequence to index 
#
# output gff

my $infile;
my $outfile;
my $k;
my $db;
my $help;
my $man;

my $index;
my $search;
my $log;
my $quiet;

my $matches;

GetOptions(# Required
	   'i|infile=s'        => \$infile,
	   'o|outfile=s'       => \$outfile,
	   # Options
	   't|target=s'        => \$db,
	   'k|kmerlen=i'       => \$k,
	   'idx|index=s'       => \$index,
	   'search'            => \$search,
	   'log'               => \$log,
	   'quiet'             => \$quiet,
	   );

if (!$infile || !$outfile || !$index) {
    print "\nERROR: No input was given.\n";
    &usage;
    exit(1);
}

my $meryl = findprog('meryl');
my $mapMers = findprog('mapMers-depth');

if ($search) {
    $matches = meryl_search($infile, $index, $k);   
} else {
    my ($idxname) = build_index($db, $index, $k);
    $matches = meryl_search($infile, $idxname, $k);
}

my ($seqid,$seqlen) = return_seq($infile);
 
open(my $mers,'<',$matches) or die "\nERROR: Could not open file: $matches\n";
open(my $gff,'>',$outfile) or die "\nERROR: Could not open file: $outfile\n";

print $gff "##gff-version 3\n";
print $gff "##sequence-region ",$seqid," 1 ",$seqlen,"\n";

my $merct = 0;

while(my $match = <$mers>) {
    chomp $match;
    my ($offset, $count) = split(/\t/,$match);
    $offset =~ s/\s//g;
    $count =~ s/\s//g;
    eval { $count = sprintf("%.2f",log($count)) if $log; };
    $merct++;
    print $gff join("\t",
		    ($seqid, "meryl", "MDR", $offset, $offset, $count, ".", "+",
		     join(";",
			  ("Name=mapMers-depth_".$k."_mer","ID=mer:$merct","dbxref=SO:0000657")
			  )
		     )
		    ),"\n"; 
}
close($mers);
close($gff);
unlink($matches);

exit;
#
# Subs
#
sub findprog {

    my $prog = shift;
    my $path = `which $prog 2> /dev/null`;
    chomp $path;
    if ( (! -e $path) && (! -x $path) ) {
	die "\nERROR: Cannot find $prog binary. Exiting.\n";
    } else {
	return ($path);
    }

}

sub return_seq {

    my $infile = shift;
    my $seq_in  = Bio::SeqIO->new( -format => 'fasta', 
				   -file => $infile);

    my %seq;  
    my $seqct = 0;
    while(my $fas = $seq_in->next_seq()) {
	$seqct++;
	$seq{$fas->id} = $fas;
	return ($fas->id, $fas->length);
    }

    if ($seqct > 1) {
	die "\nERROR: $seqct sequences present in $infile when only 1 sequence is expected. Exiting.\n";
    } else {
	foreach (keys(%seq)) {
	    print "\n========> Running meryl on sequence: $_\n" unless $quiet;
	}
    }

}

sub build_index {

    my ($db, $indexname, $k) = @_;

    my $index = "$meryl ".
	        "-v -B ".
                "-m $k ".
		"-s $db ".
		"-o $indexname ";
		$index .= $index." 2>&1 > /dev/null" if $quiet;

    print "\n========> Creating meryl index for mersize $k for sequence: $db\n";
    system($index);

}

sub meryl_search {

    my ($infile, $indexname, $k) = @_;
    my ($seqfile,$seqdir,$seqext) = fileparse($infile, qr/\.[^.]*/);
    my ($indfile,$inddir,$indext) = fileparse($indexname, qr/\.[^.]*/);
    my $searchout = $seqfile."_".$indfile.".meryl_mapMers-depth.out";

    my $search = "$mapMers ".
	         "-m $k ".
		 "-mers $indexname ".
		 "-seq $infile ".
                 "> $searchout";
    print "\n========> Searching $infile with $indexname\n" unless $quiet;
    system($search);
    return($searchout);

}

sub usage {
    my $script = basename($0);
  print STDERR <<END

USAGE: $script -i contig.fas -t target.fas -k 20 -o contig_target.gff 

Required:
    -i|infile       :    Fasta file to search (contig or chromosome).
    -o|outfile      :    File name to write the gff to.

Options:
    -t|target       :    Fasta file of WGS reads to index.
    -k|kmerlen      :    Kmer length to use for building the index.
    -s|search       :    Just search the (--infile). Must specify an existing index.
    -idx|index      :    Name of the index (if used with --search option, otherwise leave ignore this option).
    --log           :    Return the log number of matches instead of the raw count.
    --quiet         :    Don't print progress or program output.
	
END
}
