#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
#use File::Temp;

# TODO: Filter simple, dinuc, trinuc, etc. repeats           DONE
#       Filter by threshhold (raw counts or log counts)      DONE
#       SEE "filter_oligocounts.pl"
#       SEE "cnv_oligocounts2log.pl" to convert existing
#       gff file
 
#       (Makes more sense to use separate script because 
#        there is no way to know the range of counts,
#        though it may be a useful option to just do the 
#        search and print the statistical properties
#        of the k-mers.) 
#
#       Add option --gff.                                    DONE
#       Remove option --outfile and create
#       outfile names with File::Spec.

#       Print information about the outfile names.
#       Print time or message that program is completed.

#       Use File::Temp for temp fasta file creation, which 
#       would be much safer because in the event that a 
#       multifasta and one of the sequences has the same 
#       name, the original would also be deleted with the 
#       --clean option.

my $infile;
my $outfile;
my $k;
my $db;
my $help;
my $man;

my $esa;
my $index;
my $search;
my $log;
my $gff;
my $quiet;
my $filter;
my $matches;
my $ratio;
my $clean;
my $debug;

GetOptions(# Required
	   'i|infile=s'        => \$infile,
	   'o|outfile=s'       => \$outfile,
	   'e|esa'             => \$esa,
	   # Options
	   't|target=s'        => \$db,
	   'k|kmerlen=i'       => \$k,
	   'idx|index=s'       => \$index,
	   's|search'          => \$search,
	   'r|repeat-ratio=f'  => \$ratio,
	   'filter'            => \$filter,
	   'log'               => \$log,
	   'gff'               => \$gff,
	   'quiet'             => \$quiet,
	   'clean'             => \$clean,
	   'debug'             => \$debug,
	   'h|help'            => \$help,
	   'm|man'             => \$man,
	   );

if (!$infile || !$outfile || !$index) {
    print "\nERROR: No input was given.\n";
    usage();
    exit(1);
}

my $gt = findprog('gt');

if ($search && $filter && !$ratio) {
    warn "\nWARNING: Using a simple repeat ratio of 0.80 for filtering since one was not specified.\n";
}

# return reference of seq hash here and do tallymer search for each fasta in file
my ($seqhash,$seqreg,$seqct) = split_mfasta($infile);

if ($search) {
    for my $key (sort keys %$seqhash) {
	say "\n========> Running Tallymer Search on sequence: $key" unless $quiet;
	my $oneseq = getFh($key);
	$matches = tallymer_search($oneseq, $index);
	if ($gff) {
	    for my $seqregion (sort keys %$seqreg) {
		if ($key eq $seqregion) {
		    tallymersearch2gff($seqct,$matches,$outfile,$seqregion,$seqreg->{$seqregion});
		    delete $seqreg->{$seqregion};
		}
	    }
	}
	unlink $oneseq if $clean && $seqct > 1; # what we really want is: && $infile ne $oneseq
	delete $seqhash->{$key};
    }
} else {
    build_suffixarray($db);
    my ($idxname) = build_index($db, $index);
    for my $key (sort keys %$seqhash) {
        say "\n========> Running Tallymer Search on sequence: $key" unless $quiet;
	my $oneseq = getFh($key);
	$matches = tallymer_search($oneseq, $idxname);
	if ($gff) {
	    for my $seqregion (sort keys %$seqreg) {
		if ($key eq $seqregion) {
		    tallymersearch2gff($seqct,$matches,$outfile,$seqregion,$seqreg->{$seqregion});
		    delete $seqreg->{$seqregion};
		    # clean all the files from the suffix array creation here
		    # my $clean_cmd = "rm *.al1 *.des *.esq *.lcp *.llv *.prj *.sds *.ssp *.suf";
		}
	    }
	}
	unlink $oneseq if $clean && $seqct > 1; # what we really want is: && $infile ne $oneseq
	delete $seqhash->{$key};
    }
}

exit;
#
# Subs
#
sub findprog {
    my $prog = shift;
    my $path = `which $prog 2> /dev/null`;
    chomp $path;
    if ( (! -e $path) && (! -x $path) ) {
        die "\nERROR: Cannot find $prog binary. Exiting.\n\n";
    } else {
        return $path;
    }
}

sub split_mfasta {
    my $seq = shift;
    my $seq_in  = Bio::SeqIO->new( -format => 'fasta', -file => $seq);

    my %seqregion;
    my %seq;
    my $seqct = 0;

    while(my $fas = $seq_in->next_seq()) {
	$seqct++;
	$seq{$fas->id} = $fas->seq;
	$seqregion{$fas->id} = $fas->length;
    }

    if ($seqct > 1) {
	say "\n========> Running Tallymer Search on $seqct sequences." unless $quiet;
    } 
   
    return(\%seq,\%seqregion,$seqct);
	
}

sub getFh {
    my ($key) = shift;
    my $singleseq = $key.".fasta";           # fixed bug adding extra underscore 2/10/12
    $seqhash->{$key} =~ s/(.{60})/$1\n/gs;   # may speed things up marginally to not format the sequence
    open my $tmpseq, '>', $singleseq or die "\nERROR: Could not open file: $singleseq\n";
    say $tmpseq ">".$key;
    say $tmpseq $seqhash->{$key};
    close $tmpseq;
    return $singleseq;
    
}

sub build_suffixarray {
    my $db = shift;
    my $suffix = "$gt suffixerator ".
	         "-dna ".
                 "-pl ".
                 "-tis ".
                 "-suf ".
                 "-lcp ".
                 "-v ".
                 "-parts 4 ".
                 "-db $db ".
                 "-indexname $db";
    $suffix .= $suffix." 2>&1 > /dev/null" if $quiet;
    system($suffix);
}

sub build_index {
    my $db = shift;
    my $indexname = shift;
    my $index = "$gt tallymer ".
	        "mkindex ".
                "-mersize $k ".
		"-minocc 10 ".
		"-indexname $indexname ".
		"-counts ".
		"-pl ".
		"-esa $db";
    $index .= $index." 2>&1 > /dev/null" if $quiet;

    say "\n========> Creating Tallymer index for mersize $k for sequence: $db";
    system($index);
}

sub tallymer_search {
    my ($infile, $indexname) = @_;
    my ($seqfile,$seqdir,$seqext) = fileparse($infile, qr/\.[^.]*/);
    my ($indfile,$inddir,$indext) = fileparse($indexname, qr/\.[^.]*/);
    #say "========> $seqfile";
    #say "========> $indfile";

    my $searchout = $seqfile."_".$indfile.".tallymer-search.out";

    my $search = "$gt tallymer ".
	         "search ".
		 "-output qseqnum qpos counts sequence ".
                 "-tyr $indexname ".
                 "-q $infile ".
                 "> $searchout";
    say "\n========> Searching $infile with $indexname" unless $quiet;
    #say "\n========> Outfile is $searchout" unless $quiet;    # The Tallymer search output. 
    system($search);
    return $searchout;

}

sub tallymersearch2gff {
    my ($seqct, $matches, $outfile, $seqid, $end) = @_;
    my ($file,$dir,$ext) = fileparse($outfile, qr/\.[^.]*/);
    my $out;
    $out = $file."_".$seqid.".gff3" if $seqct > 1;
    $out = $outfile if $seqct == 1;

    open my $mers,'<',$matches or die "\nERROR: Could not open file: $matches\n";
    open my $gff,'>',$out or die "\nERROR: Could not open file: $out\n";
    
    say $gff "##gff-version 3";
    say $gff "##sequence-region ",$seqid," 1 ",$end;
    
    while(my $match = <$mers>) {
	chomp $match;
	my ($seqnum, $offset, $count, $seq) = split /\t/, $match;
	$offset =~ s/^\+//;
	$seq =~ s/\s//;
	my $merlen = length($seq);

	if ($filter) {
	    my $repeatseq = filter_simple($seq, $merlen);
	    unless (exists $repeatseq->{$seq} ) {
		printgff($count, $seqid, $offset, $merlen, $seq, $gff);
	    }
	} else {
	    printgff($count, $seqid, $offset, $merlen, $seq, $gff);
	}
    }
    close $gff;
    close $mers;
    unlink $matches if $clean;
}

sub filter_simple {
    my ($seq, $len) = @_;
    my %di = ('AA' => 0, 'AC' => 0, 
	      'AG' => 0, 'AT' => 0, 
	      'CA' => 0, 'CC' => 0, 
	      'CG' => 0, 'CT' => 0, 
	      'GA' => 0, 'GC' => 0, 
	      'GG' => 0, 'GT' => 0, 
	      'TA' => 0, 'TC' => 0, 
	      'TG' => 0,'TT' => 0);
    my %mono = ('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0);
    my $dict = 0;
    my $monoct = 0;
    my $diratio;
    my $monoratio;

    my %simpleseqs;
    my $repeat_ratio = defined($ratio) ? $ratio : "0.80";

    for my $mononuc (keys %mono) {
        #my $monoct = (uc($seq) =~ tr/$mononuc//);
	while ($seq =~ /$mononuc/ig) { $monoct++ };
        $monoratio = sprintf("%.2f",$monoct/$len);
        #print "Mer: $mononuc\tMer: $seq\nMononuc count: $monoct\n"; # for debug, if these simple repeats are of interest they
	#print "Mer: $mononuc\tMer: $seq\nMononuc ratio: $monoratio\n"; # are stored, along with their repeat ratio, in the hash below
	if ($monoratio >= $repeat_ratio) {
	    $simpleseqs{$seq} = $monoratio;
	}
	$monoct = 0;
    }

    for my $dinuc (keys %di) {
	while ($seq =~ /$dinuc/ig) { $dict++ };
	$diratio = sprintf("%.2f",$dict*2/$len);
	#print "Mer: $dinuc\tMer: $seq\nDinuc count: $dict\n"; # for debug, if these simple repeats are of interest they
	#print "Mer: $dinuc\tMer: $seq\nDinuc ratio: $diratio\n"; # are stored, along with their repeat ratio, in the hash below
	if ($diratio >= $repeat_ratio) {
	    $simpleseqs{$seq} = $diratio;
	}
	$dict = 0;
    }

    return(\%simpleseqs);
	    
}

sub printgff {
    my ($count, $seqid, $offset, $merlen, $seq, $gff) = @_;

    if ($log) {
	# may want to consider a higher level of resolution that 2 sig digs
	eval { $count = sprintf("%.2f",log($count)) }; warn $@ if $@; 
    }
    
    say $gff join "\t", $seqid, "Tallymer", "MDR", $offset, $offset, $count, ".", "+",
                        join ";", "Name=Tallymer_".$merlen."_mer","ID=$seq","dbxref=SO:0000657";
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
    -e|esa          :    Build the suffix array from the WGS reads (--target).
    -s|search       :    Just search the (--infile). Must specify an existing index.
    -r|ratio        :    Repeat ratio to use for filtering simple repeats (must be used with --filter).
    -idx|index      :    Name of the index (if used with --search option, otherwise leave ignore this option).
    --filter        :    Filter out simple repeats including di- and mononucleotides. (In testing phase)
    --log           :    Return the log number of matches instead of the raw count.
    --gff           :    Create a GBrowse-compatible GFF3 file. 
    --clean         :    Remove all the files generated by this script. This does not currently touch any of
                         the Tallymer suffix or index files. 			 
    --quiet         :    Don't print progress or program output.
	
END
}