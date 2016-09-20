#!/usr/bin/env perl

## NB: This takes a FASTA of LTR retrotransposons and a reference FASTA file,
##     and the result will be be a file of only the internal coding regions,
##     not the LTRs.

use 5.010;
use strict;
use warnings;
use File::Basename;
use autodie qw(open);
use Bio::DB::HTS::Kseq;
use Bio::DB::HTS::Faidx;
use Getopt::Long;
use Data::Dump::Color;

my $ltr_fasta;
my $ref_fasta;
my $te_fasta;
my $out_file;
my $help;

GetOptions(
	   'l|ltr_fasta=s'  => \$ltr_fasta,
	   'r|ref_fasta=s'  => \$ref_fasta,
	   't|te_fasta=s'   => \$te_fasta,
	   'o|outfile=s'    => \$out_file,
	   'h|help'         => \$help,
	   );

usage() and exit(0) if $help;

if (!$ltr_fasta || !$ref_fasta || !$te_fasta || !$out_file) {
    say "\nERROR: Command line not parsed correctly. Check input.";
    usage();
    exit(1);
}

open my $out, '>', $out_file;

my $faidx = Bio::DB::HTS::Faidx->new($ref_fasta);

my $kseq = Bio::DB::HTS::Kseq->new($ltr_fasta);
my $iter = $kseq->iterator();

my %coords;
while (my $seqobj = $iter->next_seq) {
    my $name = $seqobj->name;
    #>3prime_RLG_family0_LTR_retrotransposon41304_HanXRQChr07_25192872_25193586
    if ($name =~ /^([35]prime)_(RL[CGX]_family\d+_LTR_retrotransposon\d+)_(\w+\d+)_(\d+)_(\d+)/) {
	my ($ltr, $elem, $chr, $start, $end) = ($1, $2, $3, $4, $5);
	$coords{$elem}{$ltr} = join "||", $chr, $start, $end;
    }
}

my $ckseq = Bio::DB::HTS::Kseq->new($te_fasta);
my $citer = $ckseq->iterator();

while (my $seqobj = $citer->next_seq) {
    my $name = $seqobj->name;
    if ($name =~ /(RL[CGX]_family\d+_LTR_retrotransposon\d+)_(\w+\d+)_(\d+)_(\d+)/) {
	my ($elem, $chr, $start, $end) = ($1, $2, $3, $4);
	my $element_length = $end - $start + 1;
	if (exists $coords{$elem}) {
	    my ($lchr, $lstrt, $lend) = split /\|\|/, $coords{$elem}{'3prime'};
	    my ($rchr, $rstrt, $rend) = split /\|\|/, $coords{$elem}{'5prime'};

	    my $lltr_len  = $lend - $lstrt + 1;
            my $rltr_len  = $rend - $rstrt + 1;

	    ## adjust for strand
	    if ($lstrt > $rstrt) {
		$rstrt = $lstrt;
		$lend = $rend;
	    }

	    my $int_start = $lend + 1;
	    my $int_end   = $rstrt - 1;
	    my $intlen    = $int_end - $int_start + 1;

	    my ($intseq, $int_len) = $faidx->get_sequence("$chr:$int_start-$int_end");
	    # for debugging
	    #say "SEQ:      $lchr";
	    #say "ELEM:     $elem";
	    #say "SEQLEN:   $element_length";
	    #say "LLTRLEN:  $lltr_len";
	    #say "INTSTRT:  $int_start";
	    #say "INTEND:   $int_end";
	    #say "RLTRLEN:  $rltr_len";
	    #say "INTLEN:   $int_len";
	    #say "EXPLEN:   $intlen";
	    #say "-" x 20; 
		
	    $intseq =~ s/.{60}\K/\n/g;
	    my $intid = join "_", $elem, $lchr, $int_start, $int_end, "Int";
	    say $out join "\n", ">".$intid, $intseq;
	}
    }
}
close $out;

sub get_fh {
    my ($file) = @_;

    my $fh;
    if ($file =~ /\.gz$/) {
        open $fh, '-|', 'zcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /\.bz2$/) {
        open $fh, '-|', 'bzcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /^-$|STDIN/) {
	open $fh, '< -' or die "\nERROR: Could not open STDIN\n";
    }
    else {
        open $fh, '<', $file or die "\nERROR: Could not open file: $file\n";
    }

    return $fh;
}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -l ltr_fasta -r ref_fasta -t te_fasta -o outfile

Required:
     -l|ltr_file     :   File of exemplar LTR sequences.
     -r|ref_fasta    :   File of reference sequences used for identifying LTR sequences.
     -t|te_fasta     :   File of full-length LTR sequences identified by LTRdigest.
     -o|outfile      :   File to write internal LTR element regions to.

Options:
    -h|help          :   Print a usage statement.

END
}
