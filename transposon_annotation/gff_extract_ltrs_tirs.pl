#!/usr/bin/env perl

use 5.014;
use strict;
use warnings;
use autodie;
use File::Spec;
use File::Basename;
use File::Path          qw(make_path);
use Bio::GFF3::LowLevel qw(gff3_parse_feature);
use File::Temp     qw(tempfile);
use Scalar::Util   qw(looks_like_number);
use Bio::DB::HTS::Faidx;
use Sort::Naturally;
use Carp 'croak';
use Tephra::Alignment::Utils;

my $usage = basename($0)." genome.fas genome.gff3\n";
my $genome = shift or die $usage;
my $gff = shift or die $usage;

my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
my $ltrs_out = File::Spec->catfile($path, $name.'_ltrs.fasta');
my $tirs_out = File::Spec->catfile($path, $name.'_tirs.fasta');
open my $lfh, '>', $ltrs_out;
open my $tfh, '>', $tirs_out;

my $index = index_ref($genome);
my ($header, $features) = collect_gff_features($gff);
my ($ltr_files, $ltr_dir) = extract_ltr_sequences($index, $gff, $header, $features, $lfh);
my ($tir_files, $tir_dir) = extract_tir_sequences($index, $gff, $header, $features, $tfh);

exit;
#
# methods
#
sub collate {
    my ($file_in, $fh_out) = @_;

    my $lines = do { 
        local $/ = undef; 
        open my $fh_in, '<', $file_in or die "\n[ERROR]: Could not open file: $file_in\n";
        <$fh_in>;
    };
    print $fh_out $lines;

    return;
}

sub index_ref {
    my ($fasta) = @_;

    my $index = Bio::DB::HTS::Faidx->new($fasta);

    return $index;
}

sub write_element_parts {
    my ($index, $chromosome, $start, $end, $out, $id) = @_;

    my ($seq, $length) = get_full_seq($index, $chromosome, $start, $end);

    $seq =~ s/.{60}\K/\n/g;
    say $out join "\n", ">".$id, $seq;

    return;
}

sub get_full_seq {
    my ($index, $chromosome, $start, $end) = @_;

    my $location;
    if (defined $start && defined $end && looks_like_number($start) && looks_like_number($end)) {
	$location = "$chromosome:$start-$end";
    }
    elsif (defined $index && defined $chromosome) {
	$location = "$chromosome";
    }
    else {
	croak "\n[ERROR]: Index or Chromosome passed to Tephra::Role::Util::get_full_seq are undefined.\n";
    }

    my ($seq, $length) = $index->get_sequence($location);
    croak "\n[ERROR]: Something went wrong fetching sequence for '$location'. Got zero length.\n"
	unless $length;

    return ($seq, $length);
}

sub collect_gff_features {
    my ($gff) = @_;

    my $header;
    open my $in, '<', $gff or die "\n[ERROR]: Could not open file: $gff\n";
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

    open my $gffio, '<', $gff or die "\n[ERROR]: Could not open file: $gff\n";

    my ($seq_id, $start, $end, $region, $key, %features);
    while (my $line = <$gffio>) {
        chomp $line;
        next if $line =~ /^#/;
	    my $feature = gff3_parse_feature( $line );
        if ($feature->{type} eq 'repeat_region') {
            $region = @{$feature->{attributes}{ID}}[0];
            ($seq_id, $start, $end) = @{$feature}{qw(seq_id start end)};
	    $key = join "||", $region, $start, $end;

        }
	if ($feature->{type} ne 'repeat_region') {
            if ($feature->{start} >= $start && $feature->{end} <= $end) {
		push @{$features{$seq_id}{$key}}, $feature;
            }
        }
    }
    close $gffio;

    return ($header, \%features);
}

sub get_parent_coords {
    my ($parent, $coord_map) = @_;

    my ($seq_id, $start, $end) = split /\|\|/, $coord_map->{$parent};
    my $pkey = join "||", $parent, $seq_id, $start, $end;

    return ($seq_id, $pkey);
}

sub extract_ltr_sequences {
    my ($index, $gff, $header, $features, $lfh) = @_;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $dir = File::Spec->catdir($path, $name.'_ltrages');
    unless ( -d $dir ) {
        make_path( $dir, {verbose => 0, mode => 0771,} );
    }

    my ($family, %ltrs, %seen, %coord_map);
    for my $seq_id (keys %$features) {
	for my $rep_region (keys %{$features->{$seq_id}}) {
	    for my $ltr_feature (@{$features->{$seq_id}{$rep_region}}) {
		if ($ltr_feature->{type} =~ /(?:LTR|TRIM|LARD)_retrotransposon/) {
		    my $elem_id = @{$ltr_feature->{attributes}{ID}}[0];
		    $family  = @{$ltr_feature->{attributes}{family}}[0];
		    my ($start, $end) = @{$ltr_feature}{qw(start end)};
		    my $key = join "||", $family, $elem_id, $seq_id, $start, $end;
		    $ltrs{$key}{'full'} = join "||", @{$ltr_feature}{qw(seq_id type start end)};
		    $coord_map{$elem_id} = join "||", @{$ltr_feature}{qw(seq_id start end)};
		}
		elsif ($ltr_feature->{type} eq 'long_terminal_repeat') { # &&
		    my $parent = @{$ltr_feature->{attributes}{Parent}}[0];
		    my ($chr_id, $pkey) = get_parent_coords($parent, \%coord_map);
		        my ($ltr_type, $ltr_start, $ltr_end, $ltr_strand) = 
			    @{$ltr_feature}{qw(type start end strand)};
		    $ltr_strand //= '?';
		    my $ltrkey = join "||", $chr_id, $ltr_type, $ltr_start, $ltr_end, $ltr_strand;
		    my $parent_key = join "||", $family, $pkey;
		    push @{$ltrs{$parent_key}{'ltrs'}}, $ltrkey; # unless exists $seen{$ltrkey};
		    #$seen{$ltrkey} = 1;
		}
            }
        }
    }
    #dd \%ltrs and exit;

    my @files;
    my $ltrct = 0;
    my $orientation;
    for my $ltr (nsort keys %ltrs) {
        my ($family, $element, $chr, $rstart, $rend) = split /\|\|/, $ltr;
        my ($seq_id, $type, $start, $end) = split /\|\|/, $ltrs{$ltr}{'full'};
        my $ltr_file = join "_", $family, $element, $seq_id, $start, $end, 'ltrs.fasta';
        my $ltrs_out = File::Spec->catfile($dir, $ltr_file);
        die "\n[ERROR]: $ltrs_out exists. This will cause problems downstream. Please remove the previous ".
            "results and try again. Exiting.\n" if -e $ltrs_out;
        push @files, $ltrs_out;
        open my $ltrs_outfh, '>>', $ltrs_out or die "\n[ERROR]: Could not open file: $ltrs_out\n";

	if (@{$ltrs{$ltr}{'ltrs'}} != 2) {
            say STDERR "[ERROR]: $ltr contains ",scalar(@{$ltrs{$ltr}{'ltrs'}})," sequences";
            exit;
	}

        for my $ltr_repeat (@{$ltrs{$ltr}{'ltrs'}}) {
            #ltr: Contig57_HLAC-254L24||long_terminal_repeat||60101||61950||+
            my ($src, $ltrtag, $s, $e, $strand) = split /\|\|/, $ltr_repeat;
            my $ltrid;
            if ($ltrct) {
                $orientation = '5prime' if $strand eq '-';
                $orientation = '3prime'  if $strand eq '+';
                $orientation = 'unk-prime-r' if $strand eq '?';
                $ltrid = join "_", $orientation, $family, $element, $src, $s, $e;
                write_element_parts($index, $src, $s, $e, $ltrs_outfh, $ltrid);
                $ltrct = 0;
            }
            else {
                $orientation = '5prime' if $strand eq '+';
                $orientation = '3prime' if $strand eq '-';
                $orientation = 'unk-prime-f' if $strand eq '?';
                $ltrid = join "_", $orientation, $family, $element, $src, $s, $e;
                write_element_parts($index, $src, $s, $e, $ltrs_outfh, $ltrid);
                $ltrct++;
            }
        }
	close $ltrs_outfh;

	collate($dir.'/'.$ltr_file, $lfh);
	unlink $dir.'/'.$ltr_file;
    }
    return (\@files, $dir);
}

sub extract_tir_sequences {
    my ($index, $gff, $header, $features, $tfh) = @_;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $dir = File::Spec->catdir($path, $name.'_tirages');
    unless ( -d $dir ) {
        make_path( $dir, {verbose => 0, mode => 0771,} );
    }

    #dd $features;
    my ($family, %tirs, %seen, %coord_map);
    for my $seq_id (keys %$features) {
	for my $rep_region (keys %{$features->{$seq_id}}) {
	    for my $tir_feature (@{$features->{$seq_id}{$rep_region}}) {
		if ($tir_feature->{type} =~ /terminal_inverted_repeat_element|MITE/i) {
		    my $elem_id = @{$tir_feature->{attributes}{ID}}[0];
		    next unless defined $elem_id;
		    $family = @{$tir_feature->{attributes}{family}}[0];
		    my ($start, $end) = @{$tir_feature}{qw(start end)};
		        my $key = defined $family ? join "||", $family, $elem_id, $seq_id, $start, $end 
			    : join "||", 'DTX', $elem_id, $seq_id, $start, $end;
		    $tirs{$key}{'full'} = join "||", @{$tir_feature}{qw(seq_id type start end)};
		    $coord_map{$elem_id} = join "||", @{$tir_feature}{qw(seq_id start end)};
		}
		elsif ($tir_feature->{type} eq 'terminal_inverted_repeat') {
		    my $parent = @{$tir_feature->{attributes}{Parent}}[0];
		    my ($chr_id, $pkey) = get_parent_coords($parent, \%coord_map);
		        my ($tir_type, $tir_start, $tir_end, $tir_strand) = 
			    @{$tir_feature}{qw(type start end strand)};
		    $tir_strand //= '?';
		    my $tirkey = join "||", $chr_id, $tir_type, $tir_start, $tir_end, $tir_strand;
		    $pkey = defined $family ? join "||", $family, $pkey : join "||", 'DTX', $pkey;
		    push @{$tirs{$pkey}{'tirs'}}, $tirkey; # unless exists $seen{$tirkey};
		    #$seen{$tirkey} = 1;
		}
	    }
        }
    }

    my @files;
    my $tirct = 0;
    my $orientation;
    for my $tir (nsort keys %tirs) {
        my ($family, $element, $seq_id, $rstart, $rend) = split /\|\|/, $tir;
        my ($chr_id, $type, $start, $end) = split /\|\|/, $tirs{$tir}{'full'};
        my $tir_file = join "_", $family, $element, $seq_id, $start, $end, 'tirs.fasta';
        my $tirs_out = File::Spec->catfile($dir, $tir_file);
        die "\n[ERROR]: $tirs_out exists. This will cause problems downstream. Please remove the previous ".
            "results and try again. Exiting.\n" if -e $tirs_out;
        push @files, $tirs_out;
        open my $tirs_outfh, '>>', $tirs_out or die "\n[ERROR]: Could not open file: $tirs_out\n";
	
	if (@{$tirs{$tir}{'tirs'}} != 2) {
	    say STDERR "[ERROR]: $tir contains ",scalar(@{$tirs{$tir}{'tirs'}})," sequences";
	    exit;
	}

        for my $tir_repeat (@{$tirs{$tir}{'tirs'}}) {
            #Contig57_HLAC-254L24||terminal_inverted_repeat||60101||61950||+
            my ($src, $tirtag, $s, $e, $strand) = split /\|\|/, $tir_repeat;
            my $tirid;
            if ($tirct) {
                $orientation = '5prime' if $strand eq '-';
                $orientation = '3prime'  if $strand eq '+';
                $orientation = 'unk-prime-r' if $strand eq '?';
                write_tir_parts($index, $src, $element, $s, $e, $tirs_outfh, $orientation, $family);
                $tirct = 0;
            }
            else {
                $orientation = '5prime' if $strand eq '+';
                $orientation = '3prime' if $strand eq '-';
                $orientation = 'unk-prime-f' if $strand eq '?';
                write_tir_parts($index, $src, $element, $s, $e, $tirs_outfh, $orientation, $family);
                $tirct++;
            }
        }
	close $tirs_outfh;

	collate($dir.'/'.$tir_file, $tfh);
	unlink $dir.'/'.$tir_file;
    }
    return (\@files, $dir);
}

sub write_tir_parts {
    my ($index, $loc, $elem, $start, $end, $out, $orient, $family) = @_;
    my ($seq, $length) = get_full_seq($index, $loc, $start, $end);
    # need to reverse-complement the inverted seq
    my $utils = Tephra::Alignment::Utils->new;
    $seq = $utils->revcom($seq) if $orient =~ /3prime|prime-r/;

    my $id;
    $id = join "_", $family, $elem, $loc, $start, $end if !$orient;
    $id = join "_", $orient, $family, $elem, $loc, $start, $end if $orient; # for unique IDs with clustalw

    say $out join "\n", ">".$id, $seq;

    return;
}
