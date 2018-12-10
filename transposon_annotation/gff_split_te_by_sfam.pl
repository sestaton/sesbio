#!/usr/bin/env perl

#TODO: - Handle compressed input.
#      - Consider removing the text output since it is redundant with the JSON output now; this 
#        would simplify the usage.

use 5.010;
use strict;
use warnings;
use autodie;
use Sort::Naturally;
use Bio::DB::HTS::Faidx;
use Bio::GFF3::LowLevel qw(gff3_parse_feature gff3_format_feature);
use List::Util          qw(uniq);
use Tephra::Annotation::Util;
use Statistics::Descriptive;
use JSON ();
use Getopt::Long;
#use Data::Dump::Color;

my %opts;
GetOptions(\%opts, 'gfffile|i=s', 'genome|g=s', 'json|j=s', 'split|s', 'compactjson');

my $usage = "\nUSAGE: $0 -i tephra_transposons.gff3 -g genome.fas -j stats.json <--split>

Notes on usage:

* The optional '--split' argument will split each superfamily into separate GFF3 and FASTA files.

* The default JSON output is pretty-printed for reading, but the '--compactjson' argument will
print a compact, space-free file if this is not for human consumption.

* A human readable structured representation of the stats will be printed to STDOUT, 
which is a good thing to save with the GFF3 like so:

$0 -i tephra_transposons.gff3 -g genome.fas \
  -j tephra_transposons.gff3.stats.json --split > tephra_transposons.gff3.stats\n\n";

unless ($opts{gfffile} && -e $opts{gfffile}) {
    say STDERR "\nERROR: --gfffile argument is missing or the file does not exist. Check input.\n";
    say STDERR $usage;
    exit(1);
}

if ($opts{split} && !$opts{genome}) { 
    say STDERR "\nERROR: The --genome argument is required with the --split option. Check input.\n";
    say STDERR $usage;
    exit(1);
}

my ($header, $features) = collect_gff_features($opts{gfffile});
my ($sfams, $coords) = collate_features_by_superfamily($features);
my ($stats, $total) = compute_feature_stats($sfams);
write_te_stats($stats, $total);

if ($opts{json}) {
    open my $jfh, '>', $opts{json};
    my $utf8_encjson = $opts{compactjson} ? JSON->new->encode($stats) : JSON->new->pretty->encode($stats);
    say $jfh $utf8_encjson;
    close $jfh;
}

if ($opts{split}) { 
    split_gff3_by_superfamily($sfams, $coords, $opts{genome});
}

exit;
#
# methods
#
sub compute_feature_stats {
    my ($sfams) = @_;

    my $util = Tephra::Annotation::Util->new;
    my $repeat_map = $util->build_repeat_map;

    my %stats;
    for my $sfam (nsort keys %$sfams) {
	my $sf_code = $util->map_superfamily_name_to_code($sfam);
	my $sf_lineage = $repeat_map->{$sf_code};

	$stats{ $sf_lineage->{class} }{ $sf_lineage->{repeat_name} }{order} = $sf_lineage->{order};

	my %seen;
	for my $chr (nsort keys %{$sfams->{$sfam}}) {
	    for my $elements (@{$sfams->{$sfam}{$chr}}) {
		if (ref($elements) ne 'ARRAY') {
		    my $gff_feats = gff3_format_feature($elements);
		    chomp $gff_feats;
		    my @feats = split /\t/, $gff_feats;
		    if ($feats[2] =~ /helitron|non_LTR_retrotransposon/i) {
			my ($fam_id) = ($feats[8] =~ /family=(^[A-Z]{3}_\w+\d+?)\;/);
			my $length = $feats[4] - $feats[3] + 1;
			push @{$stats{ $sf_lineage->{class} }{ $sf_lineage->{repeat_name} }{lengths}}, $length;
			push @{$stats{ $sf_lineage->{class} }{ $sf_lineage->{repeat_name} }{families}}, $fam_id;
		    }
		}
		else {
		    for my $feat (@$elements) {
			if ($feat->{type} =~ /^(?:LTR|LARD|TRIM)_retrotransposon|terminal_inverted_repeat_element|MITE/) {
			    my $fam_id = @{$feat->{attributes}{family}}[0];
			    my $length = $feat->{end} - $feat->{start} + 1;
			    push @{$stats{ $sf_lineage->{class} }{ $sf_lineage->{repeat_name} }{lengths}}, $length;
			    push @{$stats{ $sf_lineage->{class} }{ $sf_lineage->{repeat_name} }{families}}, $fam_id;
			}
			elsif ($feat->{type} =~ /protein_match/) {
			    $stats{ $sf_lineage->{class} }{ $sf_lineage->{repeat_name} }{ 'protein_matches' }++;
			}
		    }
		}
	    }
	}
    }

    my %reduced;
    my $total = 0;
    my $order;
    for my $class (nsort keys %stats) {
	for my $name (sort keys %{$stats{$class}}) {
	    my @nosing_fams = grep { defined && ! /singleton/ } @{$stats{$class}{$name}{families}};
	    my @fams = uniq(@nosing_fams);

	    my $stat = Statistics::Descriptive::Full->new;
	    $stat->add_data(@{$stats{$class}{$name}{lengths}});

	    my $ave = $stat->mean;
	    my $min = $stat->min;
	    my $max = $stat->max;
	    my $ct  = $stat->count;
	    my $sd  = $stat->standard_deviation;
	    $total += $ct;

	    $reduced{$class}{ $stats{$class}{$name}{order} }{$name}{families} = @fams;
            $reduced{$class}{ $stats{$class}{$name}{order} }{$name}{count} = $ct;
            $reduced{$class}{ $stats{$class}{$name}{order} }{$name}{length}{mean} = $ave;
            $reduced{$class}{ $stats{$class}{$name}{order} }{$name}{length}{min} = $min;
            $reduced{$class}{ $stats{$class}{$name}{order} }{$name}{length}{max} = $max;
	    $reduced{$class}{ $stats{$class}{$name}{order} }{$name}{length}{stddev} = $sd;

	    if (defined $stats{$class}{$name}{protein_matches}) {
		$reduced{$class}{ $stats{$class}{$name}{order}  }{$name}{protein_matches} = $stats{$class}{$name}{protein_matches};
	    }
	}
    }

    return (\%reduced, $total);
}

sub collate_features_by_superfamily {
    my ($features) = @_;

    my $util = Tephra::Annotation::Util->new;

    my (%sfams, %coords);
    my ($seq_id, $source, $fstart, $fend, $sfam);
    for my $chr (nsort keys %$features) {
	for my $rreg (
	    map  { $_->[0] }
	    sort { $a->[1] <=> $b->[1] }
	    map  { [ $_, /\w+(?:\d+)?\|\|(\d+)\|\|\d+/ ] } 
	    keys %{$features->{$chr}} ) {
            
	    if ($rreg =~ /helitron|non_LTR_retrotransposon/i) {
		my $feat = @{$features->{$chr}{$rreg}}[0];
		my $fam_id = @{$feat->{attributes}{family}}[0];		

		my $sf;
		if (defined $fam_id && $fam_id !~ /^0$/) { 
		    $sf = $util->map_superfamily_name($fam_id);
		}
		elsif ($rreg =~ /helitron/) {
		    $sf = 'Helitron';
		}
		else {
		    $sf = 'LINE';
		}
		
		push @{$sfams{$sf}{$chr}}, $feat;
		$coords{$sf}{ @{$feat->{attributes}{ID}}[0] } = join "||", $feat->{seq_id}, $fam_id, $feat->{start}, $feat->{end};
	    }
	    elsif ($rreg =~ /repeat_region/) {
		for my $feature (@{$features->{$chr}{$rreg}}) {
		    if ($feature->{type} =~ /^(?:LTR|LARD|TRIM)_retrotransposon|terminal_inverted_repeat_element|MITE/) { 
			my $id = defined @{$feature->{attributes}{family}}[0] ? @{$feature->{attributes}{family}}[0] 
			    : @{$feature->{attributes}{ID}}[0];
			$sfam = $util->map_superfamily_name($id);
			#say STDERR join q{ }, $id, $sfam;
			push @{$sfams{$sfam}{$chr}}, $features->{$chr}{$rreg};
			$coords{$sfam}{ @{$feature->{attributes}{ID}}[0] } = 
			    join "||", $feature->{seq_id}, $id, $feature->{start}, $feature->{end};
		    }
		}
	    }
	}   
    }

    return (\%sfams, \%coords);
}

sub split_gff3_by_superfamily {    
    my ($sfams, $coords, $genome) = @_;

    my $faidx = Bio::DB::HTS::Faidx->new($genome);

    for my $sfam_id (nsort keys %$coords) {
	my $outfile;
	if ($sfam_id =~ /\//) {
	    my $tmp_id = $sfam_id;
	    $tmp_id =~ s/\//-/g;
	    $outfile = $tmp_id.'_tephra_features.fasta';
	}
	else {
	    $outfile = $sfam_id.'_tephra_features.fasta';
	}
	open my $outf, '>', $outfile or die $!;

	for my $id (nsort keys %{$coords->{$sfam_id}}) {
	    my ($chr, $fam_id, $start, $end) = split /\|\|/, $coords->{$sfam_id}{$id};
	    my $idx = "$chr:$start-$end";
	    my ($seq, $len) = $faidx->get_sequence($idx);
	    my $seq_id = join "_", $fam_id, $id, $chr, $start, $end;
	    $seq =~ s/.{60}\K/\n/g;

	    say $outf join "\n", ">$seq_id", $seq;
	}
	close $outf;
    }

    for my $sfam (nsort keys %$sfams) {
	my $sfam_id = $sfam;
	$sfam_id =~ s/\//-/;
	my $outfile = $sfam_id.'_tephra_features.gff3';
	open my $out, '>', $outfile or die $!;
	say $out $header;
	
	for my $chr (nsort keys %{$sfams->{$sfam}}) {
	    for my $elements (@{$sfams->{$sfam}{$chr}}) {
		if (ref($elements) ne 'ARRAY') {
		    my $gff_feats = gff3_format_feature($elements);
		    chomp $gff_feats;
		    say $out $gff_feats;
		}
		else {
		    for my $feat (@$elements) {
			my $gff_feats = gff3_format_feature($feat);
			chomp $gff_feats;
			say $out $gff_feats;
		    }
		}
	    }
	}
	close $out;
    }

    return;
}

sub collect_gff_features {
    my ($gff) = @_;

    my ($in, $gffio, $header);
    if ($gff =~ /\.gz/) { 
	open $in, '-|', 'zcat', $gff or die "\nERROR: Could not open file: $gff\n";
    }
    else {
	open $in, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    }

    while (my $line = <$in>) {
	chomp $line;
	next if $line =~ /^###$/;
	if ($line =~ /^##?\w+/) {
	    $header .= $line."\n";
	}
	else {
	    last;
	}
    }
    #close $in;
    #close $in or $? != 0 or die "close: $!";
    chomp $header;


    if ($gff =~ /\.gz/) {
	open $gffio, '-|', 'zcat', $gff or die "\nERROR: Could not open file: $gff\n";
    }
    else { 
	open $gffio, '<', $gff or die "\nERROR: Could not open file: $gff\n";
    }

    my ($start, $end, $region, $key, %features);
    while (my $line = <$gffio>) {
        chomp $line;
        next if $line =~ /^#/;
	my $feature = gff3_parse_feature( $line );
	next if $feature->{type} =~ /similarity|gene/;
	if ($feature->{type} =~ /helitron|non_LTR_retrotransposon/) { 
	    $region = @{$feature->{attributes}{ID}}[0];
	    $key = join "||", $region, $feature->{start}, $feature->{end};
	    push @{$features{$feature->{seq_id}}{$key}}, $feature;
	    next;
	}
        if ($feature->{type} eq 'repeat_region') {
            $region = @{$feature->{attributes}{ID}}[0];
            ($start, $end) = @{$feature}{qw(start end)};
	    $key = join "||", $region, $start, $end;
        }
	if ($feature->{type} !~ /repeat_region|gene|exon|intron|_utr|cds|rna|similarity|helitron|non_LTR_retrotransposon|solo_LTR/i) {
            if ($feature->{start} >= $start && $feature->{end} <= $end) {
		push @{$features{$feature->{seq_id}}{$key}}, $feature;
            }
        }
    }
    #close $gffio;

    return ($header, \%features);
}

sub write_te_stats {
    my ($stats, $total) = @_;

    my $util = Tephra::Annotation::Util->new;

    my $pad = "\t" x 3;
    my $s = "Total transposon number: $total\n\n";

    my %seen;
    for my $class (nsort keys %$stats) {

	$s .= "- $class\n";
	for my $order (nsort keys %{$stats->{$class}}) { 
	    for my $name (nsort keys %{$stats->{$class}{$order}}) {
		my $sf_name = $name;
		$sf_name =~ s/\s+/_/g;
		my $sf_code = $util->map_superfamily_name_to_code($sf_name);
		if ($class =~ /Class I$/) {
		    if ($order =~ /LTR/i) { 
			$s .= "\t - Long Terminal Repeat (LTR) retrotransposons\n"
			unless $seen{$order};
			$s .= "\t\t$name ($sf_code):\n";
		    }
		    elsif ($order =~ /L1|LINE/i) {
			$s .= "\t - non-Long Terminal Repeat (non-LTR) retrotransposons\n"
			    unless $seen{$order};
			$s .= "\t\t$name ($sf_code):\n";
			
		    }
		    $seen{$order} = 1;
		}
		else {
		    if ($order =~ /helitron/i) {
			unless ($seen{$order}) {
			    $s .= "\t - Subclass II\n";
			    $s .= "\t\t - Helitron transposons\n";
			}
			$s .= "$pad - $name ($sf_code)\n";
			
		    }
		    elsif ($order =~ /TIR/i) {
			unless ($seen{$order}) {
			    $s .= "\t - Subclass I\n";
			    $s .= "\t\t - Terminal Inverted Repeat (TIR) transposons\n";
			}
			$s .= "$pad - $name ($sf_code)\n";
		    }
		    $seen {$order} = 1;
		}
		
		$pad .= "\t" if $class =~ /Class II/;;
		$s .= "${pad}Total number: $stats->{$class}{$order}{$name}{count}\n";
		if (defined $stats->{$class}{$order}{$name}{protein_matches}) {
		    $s .= "${pad}Elements with protein matches: $stats->{$class}{$order}{$name}{protein_matches}\n";
		}
		else {
		    $s .= "${pad}Elements with protein matches: Undetermined\n";
		}
		$s .= "${pad}Number of families: $stats->{$class}{$order}{$name}{families}\n\n";
		$s .= "${pad}Length distribution:\n";
		$s .= "$pad\tMean: $stats->{$class}{$order}{$name}{length}{mean}\n";
		$s .= "$pad\tMinimum: $stats->{$class}{$order}{$name}{length}{min}\n";
		$s .= "$pad\tMaximum: $stats->{$class}{$order}{$name}{length}{max}\n";
		$s .= "$pad\tStandard deviation: $stats->{$class}{$order}{$name}{length}{stddev}\n";
		$pad =~ s/\t$// if $class =~ /Class II/;
	    }
	}
    }

    chomp $s;
    say $s;

    return;
}
