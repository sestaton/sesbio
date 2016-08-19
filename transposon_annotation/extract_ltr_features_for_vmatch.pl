#!/usr/bin/env perl

use 5.010;
use Sort::Naturally;
use Number::Range;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::DB::HTS::Kseq;
use Bio::DB::HTS::Faidx;
use Bio::GFF3::LowLevel qw(gff3_parse_feature);
use File::Path          qw(make_path);
use IPC::System::Simple qw(system);
use Parallel::ForkManager;
use Cwd;
use Try::Tiny;
use Data::Dump::Color;

my $usage  = "$0 genome outdir sf_gff\n";
my $fasta  = shift or die $usage;
my $dir    = shift or die $usage;
my $infile = shift or die $usage;

extract_features($fasta, $dir, $infile);

sub extract_features {
    my ($fasta, $dir, $infile) = @_;

    my $index = _index_seq($fasta);

    my ($name, $path, $suffix) = fileparse($infile, qr/\.[^.]*/);
    my $type = ($name =~ /(?:gypsy|copia|unclassified)$/i);
    die "\nERROR: Unexpected input. Should match /gypsy|copia|unclassified$/i. Exiting."
	unless defined $type;

    my $resdir = File::Spec->catdir($dir, $name);
    unless ( -d $resdir ) {
	make_path( $resdir, {verbose => 0, mode => 0771,} );
    }
    
    my $comp = File::Spec->catfile($resdir, $name."_complete.fasta");
    my $ppts = File::Spec->catfile($resdir, $name."_ppt.fasta");
    my $pbs  = File::Spec->catfile($resdir, $name."_pbs.fasta");
    my $five_pr_ltrs  = File::Spec->catfile($resdir, $name."_5prime-ltrs.fasta");
    my $three_pr_ltrs = File::Spec->catfile($resdir, $name."_3prime-ltrs.fasta");

    open my $allfh, '>>', $comp or die "\nERROR: Could not open file: $comp\n";
    open my $pptfh, '>>', $ppts or die "\nERROR: Could not open file: $ppts\n";
    open my $pbsfh, '>>', $pbs or die "\nERROR: Could not open file: $pbs\n";
    open my $fivefh, '>>', $five_pr_ltrs or die "\nERROR: Could not open file: $five_pr_ltrs\n";
    open my $threfh, '>>', $three_pr_ltrs or die "\nERROR: Could not open file: $three_pr_ltrs\n";

    open my $gffio, '<', $infile or die $!;

    my (%feature, %ltrs, %coord_map);
    while (my $line = <$gffio>) {
	chomp $line;
	next if $line =~ /^#/;
	my $feature = gff3_parse_feature( $line );
	if ($feature->{type} eq 'LTR_retrotransposon') {
	    my $elem_id = @{$feature->{attributes}{ID}}[0];
	    my $key = join "||", $elem_id, $feature->{start}, $feature->{end};
	    $ltrs{$key}{'full'} = join "||", @{$feature}{qw(seq_id type start end)};
	    $coord_map{$elem_id} = join "||", @{$feature}{qw(seq_id start end)};
	}
	if ($feature->{type} eq 'long_terminal_repeat') {
	    my $parent = @{$feature->{attributes}{Parent}}[0];
	    my ($seq_id, $pkey) = _get_parent_coords($parent, \%coord_map);
	    if ($seq_id eq $feature->{seq_id}) {
		my $ltrkey = join "||", @{$feature}{qw(seq_id type start end strand)};
		push @{$ltrs{$pkey}{'ltrs'}}, $ltrkey;
	    }
	}
	elsif ($feature->{type} eq 'primer_binding_site') {
	    my $name = $feature->{attributes}{trna};
	    my $parent = @{$feature->{attributes}{Parent}}[0];
	    my ($seq_id, $pkey) = _get_parent_coords($parent, \%coord_map);
            if ($seq_id eq $feature->{seq_id}) {
		$ltrs{$pkey}{'pbs'} =
		    join "||", $feature->{seq_id}, $feature->{type}, $name, $feature->{start}, $feature->{end};
	    }
	}
	elsif ($feature->{type} eq 'protein_match') {
	    my $name = @{$feature->{attributes}{name}}[0];
	    my $parent = @{$feature->{attributes}{Parent}}[0];
	    my ($seq_id, $pkey) = _get_parent_coords($parent, \%coord_map);
            if ($seq_id eq $feature->{seq_id}) {
		push @{$ltrs{$pkey}{'pdoms'}{$name}},
		    join "||", @{$feature}{qw(seq_id type start end strand)};
	    }
	}
	elsif ($feature->{type} eq 'RR_tract') {
	    my $parent = @{$feature->{attributes}{Parent}}[0];
	    my ($seq_id, $pkey) = _get_parent_coords($parent, \%coord_map);
            if ($seq_id eq $feature->{seq_id}) {
		$ltrs{$pkey}{'ppt'} =
		    join "||", @{$feature}{qw(seq_id type start end)};
	    }
	}
    }

    my (%pdoms, %seen_pdoms);
    my $ltrct = 0;
    for my $ltr (sort keys %ltrs) {
	my ($element, $rstart, $rend) = split /\|\|/, $ltr;
	# full element
	my ($source, $prim_tag, $fstart, $fend) = split /\|\|/, $ltrs{$ltr}{'full'};
	subseq($index, $source, $element, $fstart, $fend, $allfh);

	# pbs
	if ($ltrs{$ltr}{'pbs'}) {
	    my ($pbssource, $pbstag, $trna, $pbsstart, $pbsend) = split /\|\|/, $ltrs{$ltr}{'pbs'};
	    subseq($index, $pbssource, $element, $pbsstart, $pbsend, $pbsfh);
	}

	# ppt
	if ($ltrs{$ltr}{'ppt'}) {
	    my ($pptsource, $ppttag, $pptstart, $pptend) = split /\|\|/, $ltrs{$ltr}{'ppt'};
	    subseq($index, $source, $element, $pptstart, $pptend, $pptfh);
	}

	for my $ltr_repeat (@{$ltrs{$ltr}{'ltrs'}}) {
	    my ($src, $ltrtag, $s, $e, $strand) = split /\|\|/, $ltr_repeat;
	    my $lfname = $ltr;
	    if ($ltrct) {
		$lfname .= "_5prime-ltr.fasta" if $strand eq '+';
		$lfname .= "_3prime-ltr.fasta" if $strand eq '-';
		subseq($index, $src, $element, $s, $e, $fivefh);
		$ltrct = 0;
	    }
	    else {
		$lfname .= "_3prime-ltr.fasta" if $strand eq '+';
		$lfname .= "_5prime-ltr.fasta" if $strand eq '-';
		subseq($index, $src, $element, $s, $e, $threfh);
		$ltrct++;
	    }
	}

	if ($ltrs{$ltr}{'pdoms'}) {
	    for my $model_name (keys %{$ltrs{$ltr}{'pdoms'}}) {
		for my $ltr_repeat (@{$ltrs{$ltr}{'pdoms'}{$model_name}}) {
		    my ($src, $pdomtag, $s, $e, $str) = split /\|\|/, $ltr_repeat;
		    #"Ha10||protein_match||UBN2||132013916||132014240",
		    next if $model_name =~ /transpos(?:ase)?|mule|(?:dbd|dde)?_tnp_(?:hat)?|duf4216/i; 
		    # Do not classify elements based spurious matches 
		    push @{$pdoms{$src}{$element}{$model_name}}, join "||", $s, $e, $str;
		}
	    }
	}
    }
    close $allfh;
    close $pptfh;
    close $pbsfh;
    close $fivefh;
    close $threfh;

    ## This is where we merge overlapping hits in a chain and concatenate non-overlapping hits
    ## to create a single domain sequence for each element
    for my $src (keys %pdoms) {
	for my $element (keys %{$pdoms{$src}}) {
	    my ($pdom_s, $pdom_e, $str);
	    for my $pdom_type (keys %{$pdoms{$src}{$element}}) {
		my (%lrange, %seqs, $union);
		my $pdom_file = File::Spec->catfile($resdir, $pdom_type."_pdom.fasta");
		open my $fh, '>>', $pdom_file or die "\nERROR: Could not open file: $pdom_file\n";
		for my $split_dom (@{$pdoms{$src}{$element}{$pdom_type}}) {
		    ($pdom_s, $pdom_e, $str) = split /\|\|/, $split_dom;
		    push @{$lrange{$src}{$element}{$pdom_type}}, "$pdom_s..$pdom_e";
		}
		
		if (@{$lrange{$src}{$element}{$pdom_type}} > 1) {
		    {
			no warnings; # Number::Range warns on EVERY single interger that overlaps
			my $range = Number::Range->new(@{$lrange{$src}{$element}{$pdom_type}});
			$union    = $range->range;
		    }
		        
		    for my $r (split /\,/, $union) {
			my ($ustart, $uend) = split /\.\./, $r;
			my $seq = $self->subseq_pdoms($fasta, $src, $element, $ustart, $uend);
			my $k = join "_", $ustart, $uend;
			$seqs{$k} = $seq;
		    }
		        
		    concat_pdoms($src, $element, \%seqs, $fh);
		}
		else {
		    my ($nustart, $nuend, $str) = split /\|\|/, @{$pdoms{$src}{$element}{$pdom_type}}[0];
		    subseq($index, $src, $element, $nustart, $nuend, $fh);
		}
		close $fh;
		%seqs   = ();
		%lrange = ();
		unlink $pdom_file if ! -s $pdom_file;
	    }
	}
    }

    for my $file ($comp, $ppts, $pbs, $five_pr_ltrs, $three_pr_ltrs) {
	unlink $file if ! -s $file;
    }

    return $resdir
}

sub subseq_pdoms {
    my ($index, $loc, $elem, $start, $end) = @_;

    my $location = "$loc:$start-end";
    my ($seq, $length) = $index->get_sequence($location);

    return $seq;
}

sub concat_pdoms {
    my ($src, $elem, $seqs, $fh_out) = @_;
    my @ranges = map { split /\_/, $_ } keys %$seqs;
    my $start  = min(@ranges);
    my $end    = max(@ranges);
    my $id     = join "_", $elem, $src, $start, $end;

    my $concat_seq;
    for my $seq (values %$seqs) {
	$concat_seq .= $seq;
    }

    $concat_seq =~ s/.{60}\K/\n/g;
    say $fh_out join "\n", ">$id", $concat_seq;
}

sub subseq {
    my ($index, $loc, $elem, $start, $end, $out) = @_;

    my $location = "$loc:$start-end";
    my ($seq, $length) = $index->get_sequence($location);

    my $id = join "_", $elem, $loc, $start, $end;

    if ($seq) {
	$seq =~ s/.{60}\K/\n/g;
	say $out join "\n", ">$id", $seq;
    }
}

sub _get_parent_coords {
    my ($parent, $coord_map) = @_;

    my ($seq_id, $start, $end) = split /\|\|/, $coord_map->{$parent};
    my $pkey = join "||", $parent, $start, $end;

   return ($seq_id, $pkey);
}

sub _index_seq {
    my ($fasta) = @_;

    my $index = Bio::DB::HTS::Faidx->new($fasta);
    return $index;
}
