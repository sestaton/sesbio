#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

my $blast_i;   # blast for species i against target x
my $blast_j;   # blast for species j against target x
my $blast_f;   # blast for species f against target x

my $fas_i;     # fasta for species i 
my $fas_j;     # fasta for species j
my $fas_f;     # fasta for species f

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
    say "\nERROR: Command line not parsed correctly. Exiting.";
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

for my $gene (keys %$bl_hash_f) {
    if (exists $bl_hash_i->{$gene} && exists $bl_hash_j->{$gene}) {
	$bl_match_hash->{$gene}->{$species_i} = $bl_hash_i->{$gene};
        $bl_match_hash->{$gene}->{$species_j} = $bl_hash_j->{$gene};
	$bl_match_hash->{$gene}->{$species_f} = $bl_hash_f->{$gene};
    }
}

for my $match (keys %$bl_match_hash) {

    my $gene_file = $match.".fasta";
    open my $out, ">>", $gene_file or die "\nERROR: Could not open file: $!\n";

    for my $gene_i_copy (@{$bl_match_hash->{$match}->{$species_i}}) {
	if (exists $fa_hash_i->{ $gene_i_copy }) {
	    say $out ">".$species_i."_".$gene_i_copy;
	    say $out $fa_hash_i->{ $gene_i_copy }; 
	} 
	else { # for debug
	    say "$species_i\t$match\t$gene_i_copy";
	}
    }

    for my $gene_j_copy (@{$bl_match_hash->{$match}->{$species_j}}) {
	if (exists $fa_hash_j->{ $gene_j_copy }) {
	    say $out ">".$species_j."_".$gene_j_copy;
	    say $out $fa_hash_j->{ $gene_j_copy }; 
	} 
	else { # for debug
	    say "$species_j\t$match\t$gene_j_copy";
	}
    }

    for my $gene_f_copy (@{$bl_match_hash->{$match}->{$species_f}}) {
	if (exists $fa_hash_f->{ $gene_f_copy }) {
	    say $out ">".$species_f."_".$gene_f_copy;
	    say $out $fa_hash_f->{ $gene_f_copy };
	} 
	else {
	    say "$species_f\t$match\t$gene_f_copy";
	}
    }
    
    close $out;

}

exit;
#
# Subs
#
sub blast2hash {
    my $bl = shift;

    open my $fh, '<', $bl or die "\nERROR: Could not open file: $!\n";

    my %hash;

    while(my $line = <$fh>) {
	chomp $line;
	next if $line =~ /^Query/ || $line =~ /^#/;
	my @fields = split /\t/, $line;
	my $contigID = $fields[0];
	#$contigID =~ s/\_\d\_ORF\d//;           # this is for cleaning up ids from sixpack translation to match velvet contig IDs
	my $geneID = $fields[1];
	$geneID =~ s/\|.*//;
	# This is where multiple hits are kept for each gene
	if (exists $hash{$geneID}) {
	    push(@{ $hash{$geneID} }, $contigID); # I am using the old push @{} syntax for backwards compatability (anything prior to Perl 5.14)
	}
	else {
	    $hash{$geneID} = [ $contigID ];
	}
    }
    close $fh;

    return(\%hash);
}

sub seq2hash {
    my $fas = shift;

    my %seqhash;
    my $seqct = 0;
    
    open my $seq_in, "<", $fas or die "\nERROR: Could not open file: $fas\n";  

    my ($name, $seq, $qual);
    my @aux = undef;

    while (($name, $seq, $qual) = readfq(\*$seq_in, \@aux)) {
	$seqct++;
	$seqhash{$name} = $seq;
    }
    close $seq_in;

    say "$seqct sequences in $fas";
    return(\%seqhash);
}

sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!defined(@$aux));
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
    #my $name;
    #if (/^.?(\S+\s\S+.*)/) {          # Illumina 1.8+, now more greedy 8/1 SES
	#$name = $1;
    #}
    #elsif (/^.?(\S+)/) {              # Illumina 1.3+
    #    $name = $1;
    #} 
    #else {
        #$name = '';                   # ?
    #}
    my $name = /^.(\S+)/? $1 : '';   # Heng Li's original regex
    #}
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
    return ($name, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
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
