#!/usr/bin/env perl

## NB: This takes a GFF of LTR retrotransposons and a reference FASTA file,
##     and the result will be be a file of only the internal coding regions,
##     not the LTRs.

use 5.010;
use strict;
use warnings;
use File::Basename;
use autodie qw(open);
use Getopt::Long;

my $annot_file;
my $ref_fasta;
my $te_fasta;
my $out_file;
my $help;

GetOptions(
	   'a|annot_file=s'  => \$annot_file,
	   'r|ref_fasta=s'   => \$ref_fasta,
	   't|te_fasta=s'    => \$te_fasta,
	   'o|outfile=s'     => \$out_file,
	   'h|help'          => \$help,
	   );

usage() and exit(0) if $help;

if (!$annot_file || !$ref_fasta || !$te_fasta || !$out_file) {
    say "\nERROR: Command line not parsed correctly. Check input.";
    usage();
    exit(1);
}
							  
open my $in, '<', $annot_file;
open my $out, '>>', $out_file;
my $tefas = get_fh($te_fasta);

my $seqstore = store_seq($ref_fasta, $annot_file);
my $seqidmap = seqid_map();

my @aux = undef;
my ($name, $comm, $seq, $qual);

my %ltr_map;
my ($mapped_ct, $notmapped_ct) = (0, 0);

my $header = <$in>;
while (<$in>) {
    chomp;
    my @f = split;
    my ($element_start, $element_end, $element_length, $sequence, 
	$lLTR_start, $lLTR_end, $lLTR_length, $rLTR_start, $rLTR_end) = @f[0..8];

    $sequence = correct_ids($sequence, $seqidmap);

    my $idkey = mk_key($element_start, $sequence);
    $ltr_map{$idkey} = mk_key($lLTR_start, $lLTR_end, $rLTR_start, $rLTR_end, $element_length);
}
close $in;


while (($name, $comm, $seq, $qual) = readfq(\*$tefas, \@aux)) {
    if ($name =~ /^(RL.[-_][a-zA-Z]{1,7}[-_]\d+)[-_](.*)/) {
	my ($famname, $contigname) = ($1, $2);
	my ($start, $end) = $contigname =~ /_(\d+)_(\d+)$/;
	my $pat = "_".$start."_".$end;
	$contigname =~ s/$pat//;
	$contigname = correct_ids($contigname, $seqidmap); 
	my $idkey = mk_key($start, $contigname);

	if (exists $ltr_map{$idkey}) {
	    $mapped_ct++;
	    my ($lLTR_start, $lLTR_end, $rLTR_start, $rLTR_end, $element_length) = mk_vec($ltr_map{$idkey});

	    my $int_start = $lLTR_end + 1;
	    my $int_end   = $rLTR_start - 1;
	    my $int_len   = $int_end - $int_start;

	    my $intseq = substr $seqstore->{$contigname}, $lLTR_end, $int_len;
	    $intseq =~ s/.{60}\K/\n/g;
	    my $intid = join "_", $famname, $contigname, $int_start, $int_end, "Int";
	    say $out join "\n", ">".$intid, $intseq;
	    # for debugging
	    say "SEQ:    $contigname";
	    say "SEQLEN: $element_length";
	    say "INTLEN: $int_len";
	    say "-" x 20; 
	}
	else {
	    $notmapped_ct++;
	}
    }
}
close $tefas;
close $out;

# for debugging
#say "NOT MAPPED: $notmapped_ct";
#say "Number of elements mapped: $mapped_ct;";

#
# methods
#
sub mk_key { return join "||", map { $_ // " " } @_ }

sub mk_vec { return split /\|\|/, shift }

sub correct_ids {
    my ($contigname, $seqidmap) = @_;

    if (exists $seqidmap->{$contigname}) {
	return $seqidmap->{$contigname};
    }
    else {
	return $contigname;
    }
}
    
sub seqid_map {
    # this is necessary because LTRdigest truncates the ref names in the GFF
    my %idmap  = (
		  'Contig112_HLAB-P189P24'   => 'Contig112_HLAB-P189P',
		  'Contig169_HLAE-P339N08'   => 'Contig169_HLAE-P339N',
		  'Contig23_HLAE-P245O15'    => 'Contig23_HLAE-P245O1',
		  'Contig25_HLAE-P245O15'    => 'Contig25_HLAE-P245O1',
		  'Contig283_HLAB-P94O19'    => 'Contig283_HLAB-P94O1',
		  'Contig29_HLAB-P347L09'    => 'Contig29_HLAB-P347L0',
		  'Contig36_HLAE-P245O15'    => 'Contig36_HLAE-P245O1', 
		  'Contig40_HLAB-P347K03'    => 'Contig40_HLAB-P347K0',
		  'Contig500_HLAE-P102A12'   => 'Contig500_HLAE-P102A',
		  'Contig84_HLAE-P408L01'    => 'Contig84_HLAE-P408L0',
		  'Contig86_HLAB-P392C18'    => 'Contig86_HLAB-P392C1',
		  'RL11_sub1_30x_c2-a_c2588' => 'RL11_sub1_30x_c2-a_c',
		  'Contig173_BWAZ_227-17'    => 'Contig173_BWAZ_227-1',
		  'Contig37_HLAE-P245O15'    => 'Contig37_HLAE-P245O1',
		  'RL11_sub1_30x_c3_c2588'   => 'RL11_sub1_30x_c3_c25',
		  'RL11_sub1_30x_c1_c2588'   => 'RL11_sub1_30x_c1_c25',
		  );

    my %reversed_map = reverse %idmap;

    return \%reversed_map;
}

sub store_seq {
    my ($ref_fasta) = @_;
    
    my %seqstore;
    my $reffh = get_fh($ref_fasta);
    
    my @aux = undef;
    my ($name, $comm, $seq, $qual);

    while (($name, $comm, $seq, $qual) = readfq(\*$reffh, \@aux)) {
	$seqstore{$name} = $seq;
    }
    close $reffh;

    return \%seqstore;
}

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

sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!@$aux);
    return if ($aux->[1]);
    if (!defined($aux->[0])) {
        while (<$fh>) {
            chomp;
            if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                $aux->[0] = $_;
                last;
            }
        }
        if (!defined($aux->[0])) {
            $aux->[1] = 1;
            return;
        }
    }
    my ($name, $comm);
    defined $_ && do {
        ($name, $comm) = /^.(\S+)(?:\s+)(\S+)/ ? ($1, $2) : 
	                 /^.(\S+)/ ? ($1, '') : ('', '');
    };
    my $seq = '';
    my $c;
    $aux->[0] = undef;
    while (<$fh>) {
        chomp;
        $c = substr($_, 0, 1);
        last if ($c eq '>' || $c eq '@' || $c eq '+');
        $seq .= $_;
    }
    $aux->[0] = $_;
    $aux->[1] = 1 if (!defined($aux->[0]));
    return ($name, $comm, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $comm, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -i -a annotation_file -r ref_fasta -t te_fasta -o outfile

Required:
     '-a|annot_file'   :   Tab-delimited annotation file ("*tabout.csv") produced by LTRdigest.
     '-r|ref_fasta'    :   File of reference sequences used for identifying LTR sequences.
     '-t|te_fasta'     :   File of full-length LTR sequences identified by LTRdigest.
     '-o|outfile'      :   File to write internal LTR element regions to.

Options:
    -h|help            :   Print a usage statement.

END
}
