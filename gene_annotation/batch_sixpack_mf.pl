#!/usr/bin/perl -w

=head1 NAME 
                                                                       
batch_sixpack.pl - Translate many fasta files and keep the longest ORF from each

=head1 SYNOPSIS    

batch_sixpack.pl -i seqs.fas -o seqs_trans.faa

=head1 DESCRIPTION
                                                                   
Translate a nucleotide multi-fasta file in all 6 frames and select 
the longest translated ORF (for each ORF). The minimum ORF length to report translations for
can be given as an option.

=head1 DEPENDENCIES

This script uses EMBOSS and BioPerl so, both must be installed. 
EMBOSS v6.2+ must be installed (the latest is v6.4 as of this writing), 
as well as BioPerl v1.6.1+ (the latest is v1.6.9 as of this writing).

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The fasta files to be translated.

=item -o, --outfile

A file to place the translated sequences.

=back

=head1 OPTIONS

=over 2

=item -orflen, --orflength

The minimum ORF length for which to show a translation (default is 80).
Lowering this value will not likely result in any significant hits 
(though there may be a reason to do so).

=item --clean

Remove the files created by EMBOSS during the conversion.

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut    

#
# Includes
#
use strict;

use Cwd;
use Getopt::Long;
use File::Basename;
use File::Copy;
use File::Temp;
use Bio::SeqIO;
use Pod::Usage;

#
# Vars
#
my $infile; 
my $outfile;
my $orflen;
my $clean;
my $help;
my $man;

#
# Counters
#
my $seqstot = 0;
my $transseqstot = 0;
my $pepct = 0;
my $fcount = 0;

GetOptions(#Required
	   "i|infile=s"            => \$infile,
	   "o|outfile=s"           => \$outfile,
	   #Options
	   "orflen|orflength=i"   => \$orflen,
	   "clean"                => \$clean,
	   "h|help"               => \$help,
	   "m|man"                => \$man,
	  );

pod2usage( -verbose => 2 ) if $man;

#
# Check @ARGVs
#  
usage() and exit(0) if $help;

if (!$infile || !$outfile) {
    print "\nERROR: No input was given.\n";
    usage();
    exit(1);
}

my $fasnum = `grep -c ">" $infile`;
if ($fasnum >= 1) {
    print "\n========== Translating sequences with minimum ORF length of $orf.\n";
} else {
    die "\nERROR: No fasta files were found! Must end with .fasta, .fas, or .fna.\n";
}

my $orf = defined($orflen) ? $orflen : "80";
my $seq_in = Bio::SeqIO->new(-file => $infile, -format => 'fasta');
my ($iname, $ipath, $isuffix) = fileparse($infile, qr/\.[^.]*/);

while(my $seq = $seq_in->next_seq) {

    $fcount++;

    my $tmpiname = $iname."_".$fcount."_XXXX";

    my $cwd = getcwd();
    my $fname = File::Temp->new( TEMPLATE => $tmpiname,
				 DIR => $cwd,
				 SUFFIX => ".faa",
				 UNLINK => 0);

    open(my $fh, '>', $fname) or die "\nERROR: Could not open file: $!\n";

    print $fh ">".$seq->id."\n".$seq->seq."\n";
  
    close($fh);

    #my $pepfile = $iname.".faa";                            
    my $pepfiletmp = $fname."_tmp";                      
    my $transout = $iname."_".$fcount.".sixpack";

    my $transcmd = "/usr/local/emboss/latest/bin/sixpack ".
	           "-sequence $fname ".
		   "-outseq $pepfiletmp ".
		   "-outfile $transout ".
		   "--orfminsize $orf ". 
		   "-auto ";     
    system($transcmd);

    if (-s $pepfiletmp) {
	$transseqstot++;

	#my $pepfilesort = $pepfiletmp."_sorted";
	my $seqin = Bio::SeqIO->new('-file' => $pepfiletmp, '-format' => 'fasta');
        my $seqout = Bio::SeqIO->new('-file' => ">>$outfile", '-format' => 'fasta');

	my @seq_array;
	while(my $seq = $seqin->next_seq) {
	    push(@seq_array,$seq);
	}

	@seq_array = sort { $b->length <=> $a->length } @seq_array;
	foreach my $seqs (@seq_array) {
	    $seqout->write_seq($seqs) if $pepct < 1;
	    $pepct++;
	}
	$pepct = 0;
	@seq_array = ();
	
	#move("$pepfile","../$outdir") || die "Copy failed: $!";
	unlink($pepfiletmp);
	unlink($fname);
	#unlink($pepfilesort);
    }
    if ($clean) {
	unlink($transout);
    }
}

print "\n========== $fcount sequences processed.\n";
print "\n========== $transseqstot sequences translated with ORFs above $orf.\n";

exit;

#
# Subs
#
sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script -i infile -o outfile 

Required:
 -i|infile            :       A multifasta file. Each individual sequence will be translated.
 -o|outfile           :       A file to put the translated sequences.

Options:
 -orflen|orflength   :       An interger that will serve as the lower threshold
                             length for ORFs to consider prior to translating.
 --clean             :       Remove all of the alignment files created by EMBOSS
                             if only the translated sequences are desired (Recommeded).

END
}
