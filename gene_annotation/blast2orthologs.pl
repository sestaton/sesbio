#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;
use Getopt::Long;
#use Data::Dumper;
use File::Basename;

my $blast_i; # blast for species i against target x
my $blast_j; # blast for species j against target x
my $blast_f; # etc. 

my $fas_i;   # fasta for species i 
my $fas_j;   # fasta for species j
my $fas_f;   # etc. 

my $species_i; # need something less abstract to write to the 
my $species_j; # sequence files for each hit. i,j,f don't
my $species_f; # mean very much in an alignment or tree.

GetOptions(
	   'ib|blast_i=s'   => \$blast_i,
	   'jb|blast_j=s'   => \$blast_j,
	   'fb|blast_f=s'   => \$blast_f,
	   'if|fas_i=s'     => \$fas_i,
	   'jf|fas_j=s'     => \$fas_j,
	   'ff|fas_f=s'     => \$fas_f,
	   'is|species_i=s' => \$species_i,
	   'js|species_j=s' => \$species_j,
	   'fs|species_f=s' => \$species_f,
	   );

if (!$blast_i || !$blast_j || !$blast_f || 
    !$fas_i || !$fas_j || !$fas_f || 
    !$species_i || !$species_j || !$species_f) {
    print "\nERROR: Command line not parsed correctly. Exiting.\n";
    usage();
    exit(1);
}

my $bl_hash_i = blast2hash($blast_i);
my $bl_hash_j = blast2hash($blast_j);
my $bl_hash_f = blast2hash($blast_f);

my $fa_hash_i = seq2hash($fas_i);
my $fa_hash_j = seq2hash($fas_j);
my $fa_hash_f = seq2hash($fas_f);

my $bl_match_hash = {};

foreach my $gene (keys %$bl_hash_f) {
    if (exists $bl_hash_i->{$gene} && exists $bl_hash_j->{$gene}) {

	$bl_match_hash->{$gene}->{'speciesi'} = $bl_hash_i->{$gene};
	$bl_match_hash->{$gene}->{'speciesj'} = $bl_hash_j->{$gene};
	$bl_match_hash->{$gene}->{'speciesf'} = $bl_hash_f->{$gene};

    }
}

#print Dumper $saff_fa_hash;
#print Dumper $lett_fa_hash;

foreach my $match (keys %$bl_match_hash) {

    my $gene_file = $match.".fasta";
    open(my $out, '>>', $gene_file) or die "\nERROR: Could not open file: $!\n";

    my $geneID_i = $bl_match_hash->{$match}->{'speciesi'};
    my $geneID_j = $bl_match_hash->{$match}->{'speciesj'};
    my $geneID_f = $bl_match_hash->{$match}->{'speciesf'};

    if (exists $fa_hash_i->{ $geneID_i }) {
	print $out ">".$species_i."_".$geneID_i,"\n";
        print $out $fa_hash_i->{ $geneID_i },"\n"; 
    } 
    else {
	print "$species_i\t$match\t$geneID_i\n";
    }

    if (exists $fa_hash_j->{ $geneID_j }) {
	print $out ">".$species_j."_".$geneID_j,"\n";
	print $out $fa_hash_j->{ $geneID_j },"\n"; 
    } 
    else {
	print "$species_j\t$match\t$geneID_j\n";
    }

    if (exists $fa_hash_f->{ $geneID_f }) {
	print $out ">".$species_f."_".$geneID_f,"\n";
        print $out $fa_hash_f->{ $geneID_f },"\n";
    } 
    else {
	print "$species_f\t$match\t$geneID_f\n";
    }

    close($out);

}

exit;

#
# Subs
#
sub blast2hash {
    my $bl = shift;

    open(my $fh, '<', $bl) or die "\nERROR: Could not open file: $!\n";

    my %hash;

    while(my $line = <$fh>) {
	chomp $line;
	next if $line =~ /^Query/ || $line =~ /^#/;
	my @fields = split(/\t/,$line);
	my $contigID = $fields[0];
	#$contigID =~ s/\_\d\_ORF\d//;
	my $geneID = $fields[1];
	$geneID =~ s/\|.*//;
	$hash{$geneID} = $contigID;
    }
    close($fh);

    return(\%hash);
}

sub seq2hash {
    my $fas = shift;

    my %seqhash;
    my $seqct = 0;

    my $seq_in = Bio::SeqIO->new(-file => $fas, -format => 'fasta');

    while(my $seq = $seq_in->next_seq) {
	$seqct++;
	$seqhash{$seq->id} = $seq->seq;
    }
    
    print "$seqct sequences in $fas\n";
    return(\%seqhash);
}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script [-ib] [-jb] [-fb] [-if] [-jf] [-ff] [-is] [-js] [-fs]

Required:
    -ib|blast_i     :    BLAST table (e.g., -m 8) for species i against target x.
    -jb|blast_j     :    BLAST table (e.g., -m 8) for species j against target x.
    -fb|blast_f     :    BLAST table (e.g., -m 8) for species f against target x.
    -if|fas_i       :    Fasta file for species i (NB: IDs must match those in BLAST report).
    -jf|fas_j       :    Fasta file for species j (NB: IDs must match those in BLAST report).
    -ff|fas_f       :    Fasta file for species f (NB: IDs must match those in BLAST report).
    -is|species_i   :    Common name for species i.
    -js|species_j   :    Common name for species j.
    -fs|species_s   :    Common name for species f.

Options:
   

END
}
