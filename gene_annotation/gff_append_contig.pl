#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie;
use Sort::Naturally;
use List::UtilsBy qw(nsort_by);
use Bio::Tools::GFF;

my $usage  = "$0 orig.gff new.gff > updated.gff";
my $infile = shift or die $usage;
my $newgff = shift or die $usage;
open my $in, '<', $infile;

my @aux = undef;
my ($id, $comm, $seq, $qual);
my (%chroms, %contigs);

while (my $line = <$in>) {
    chomp $line;
    if ($line =~ /contig/) {
        my @f = split /\t/, $line;
	$contigs{$f[0]} = join "||", @f;
    }
    if ($line =~ /^##FASTA$/) {
	while (($id, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
	    $chroms{$id} = $seq;
        }
    }
}
close $in;

my %seen;
my ($header, $features) = collect_gff_features($newgff);
say $header;

for my $ref (nsort keys %$features) {
    for my $id (nsort_by { m/\w+\.(\d+)\.\d+/ and $1 } keys %{$features->{$ref}}) {
	my ($parentid, $start, $stop) = split /\./, $id;
	for my $parent (keys %{$features->{$ref}{$id}}) {
	    my @parent_feats = split /\|\|/, $parent;
	    $parent_feats[8] = _format_parent_attribute($parent_feats[8]);
	    unless (exists $seen{$ref}) {                   
		my @c = split /\|\|/, $contigs{$ref};
		say join "\t", @c;                                                                         
		say "###";
		$seen{$ref} = 1;
	    }
	    say join "\t", @parent_feats;
	    
	    for my $feat (@{$features->{$ref}{$id}{$parent}}) {
		my @part_feats = split /\|\|/, $feat;
		$part_feats[8] = _format_part_attribute($part_feats[8]);
		say join "\t", @part_feats;
	    }
	}
	say "###";
    }
}

say "##FASTA";                                                              
for my $chr (nsort keys %chroms) {
    $chroms{$chr} =~ s/.{60}\K/\n/g;
    say join "\n", ">$chr", $chroms{$chr};                                                                    
}

exit;
#
# methods
#
sub collect_gff_features {
    my ($gff) = @_;

    my $header;
    open my $in, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    while (<$in>) {
        chomp;
        next if /^###$/;
        if (/^##?\w+/) {
            $header .= $_."\n";
        }
        else {
            last;
        }
    }
    close $in;
    chomp $header;

    my $gffio = Bio::Tools::GFF->new( -file => $gff, -gff_version => 3 );
    $gffio->ignore_sequence(1); # faster, more memory-efficient to not parse seqs in bioperl

    my ($start, $end, $region, $parent, $ref, %features);
  FEATURE:
    while (my $feature = $gffio->next_feature()) {
        if ($feature->primary_tag =~ /protein_match|expressed_sequence_match|gene/) {
            my @string = split /\t/, $feature->gff_string;
            ($region) = ($string[8] =~ /ID=?\s+?(protein_match\d+|expressed_sequence_match\d+|gene\d+)/);
            ($ref, $start, $end) = ($string[0], $feature->start, $feature->end);
            $parent = join "||", @string;
        }
        next FEATURE unless defined $start && defined $end && defined $ref;
        if ($feature->primary_tag !~ /protein_match|expressed_sequence_match|gene/) {
            if ($feature->start >= $start && $feature->end <= $end) {
                push @{$features{$ref}{$region.".".$start.".".$end}{$parent}}, 
		join "||", split /\t/, $feature->gff_string;
            }
        }
    }

    return ($header, \%features);
}

sub _format_parent_attribute {
    my ($str) = @_;

    $str =~ s/\s\;\s/\;/g;
    $str =~ s/\s+/=/g;
    $str =~ s/\s+$//;
    $str =~ s/=$//;
    $str =~ s/=\;/;/g;
    $str =~ s/\"//g;
    
    return $str;
}

sub _format_part_attribute {
    my ($str) = @_;

    $str =~ s/\s\;\s/\;/g;
    $str =~ s/\s\;/\;/g;
    $str =~ s/\s+$//;
    $str =~ s/=$//;
    $str =~ s/=\;/;/g;
    $str =~ s/\"//g;
    $str =~ s/^(\w+)\s/$1=/;
    $str =~ s/\;(\w+)\s/\;$1=/g;
    
    return $str;
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
