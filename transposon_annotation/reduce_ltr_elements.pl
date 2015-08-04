#!/usr/bin/env perl

use 5.020;
use warnings;
use autodie;
use File::Basename;
use Statistics::Descriptive;
use Bio::Tools::GFF;
use Sort::Naturally qw(nsort);
use List::Util      qw(sum max);
use List::MoreUtils qw(first_index);
use List::UtilsBy   qw(nsort_by);
use Set::IntervalTree;
use Data::Dump;
use experimental 'signatures';

my $usage = "$0 ltrdigest85.gff3 ltrdigest99.gff3";
my $fullgff = shift or die $usage;
my $partgff = shift or die $usage;

my ($all_feats, $intervals) = collect_features($fullgff);
my ($part_feats, $part_int) = collect_features($partgff);

my $best_elements = get_overlaps($all_feats, $part_feats, $intervals);
sort_features($fullgff, $best_elements);

#
# methods
#
sub sort_features ($gff, $best_elements) {
    my ($header, %features);
    open my $in, '<', $gff;
    while (<$in>) {
	chomp;
	if (/^#/) {
	    $header .= $_."\n";
	}
	else {
	    last;
	}
    }
    close $in;
    chomp $header;

    say $header;
    for my $str (@$best_elements) {
	for my $chr (keys %$str) {
	    push @{$features{$chr}}, @{$str->{$chr}};
	}
    }

    for my $chromosome (nsort keys %features) {
	for my $elements (nsort @{$features{$chromosome}}) {
	    for my $ltr (sort { $a =~ /repeat_region(\d+)/ <=> $b =~ /repeat_region(\d+)/ } keys %$elements) {
		for my $entry (@{$elements->{$ltr}}) {
		    my @feats = split /\|\|/, $entry;
		    $feats[8] =~ s/\s\;\s/\;/g;
		    $feats[8] =~ s/\s+$//;
		    $feats[8] =~ s/\"//g;
		    $feats[8] =~ s/(\;\w+)\s/$1=/g;
		    $feats[8] =~ s/\s;/;/;
		    $feats[8] =~ s/^(\w+)\s/$1=/;
		    say join "\t", @feats;
		}
	    }
	}
    }
    
}
    
sub collect_features ($gff) {
    my %intervals;

    my $gffio = Bio::Tools::GFF->new( -file => $gff, -gff_version => 3 );

    my ($source, $start, $end, $length, $region, %features);
    while (my $feature = $gffio->next_feature()) {
	if ($feature->primary_tag eq 'repeat_region') {
	    my @string = split /\t/, $feature->gff_string;
	    $source = $string[0];
	    ($region) = ($string[8] =~ /ID=?\s+?(repeat_region\d+)/);
	    ($start, $end) = ($feature->start, $feature->end);
	    $length = $end - $start + 1;
	}
	next $feature unless defined $start && defined $end;
	if ($feature->primary_tag ne 'repeat_region') {
	    if ($feature->start >= $start && $feature->end <= $end) {
		$intervals{$region} = join ".", $start, $end, $length;
		my $region_key = join ".", $region, $start, $end, $length;
		push @{$features{$source}{$region_key}}, join "||", split /\t/, $feature->gff_string;
	    }
	}
    }

    my $filtered = filter_compound_elements(\%features);
    return $filtered, \%intervals;
}

sub get_overlaps ($allfeatures, $partfeatures, $intervals) {
    my (@best_elements, %chr_intervals);

    for my $source (keys %$allfeatures) {
	my $tree = Set::IntervalTree->new;
	
	for my $rregion (keys %{$allfeatures->{$source}}) {
	    my ($reg, $start, $end, $length) = split /\./, $rregion;
	    $tree->insert($reg, $start, $end);
	    $chr_intervals{$source} = $tree;
	}
    }

    for my $source (keys %$partfeatures) {
	for my $rregion (keys %{$partfeatures->{$source}}) {
	    my (%scores, %sims);

	    my $tree;
	    if (exists $chr_intervals{$source}) {
		my ($reg, $start, $end, $length) = split /\./, $rregion;
		my $res = $chr_intervals{$source}->fetch($start, $end);
		#say "Found: ",scalar(@$res)," overlaps.";
		
		my ($score99, $sim99) = summarize_features($partfeatures->{$source}{$rregion});
		#say join q{ }, "LTR99 element: ", $reg, $start, $end, $length, $score99, $sim99;
		
		$scores{$source}{$rregion} = $score99;
		$sims{$source}{$rregion} = $sim99;
		
		## collect all features, then compare scores/sim...
		for my $over (@$res) {
		    my ($s, $e, $l) = split /\./, $intervals->{$over};
		    my $region_key = join ".", $over, $s, $e, $l;
		    my ($score85, $sim85) = summarize_features($allfeatures->{$source}{$region_key});
		    #say join q{ }, "LTR85 element: ", $over, $s, $e, $l, $score85, $sim85;
		    
		    $scores{$source}{$region_key} = $score85;
		    $sims{$source}{$region_key} = $sim85;
		}

		my $best_element = get_ltr_score_dups(\%scores, \%sims, $allfeatures, $partfeatures);
		push @best_elements, $best_element;
	    }
	}
    }
    
    return \@best_elements;
}

sub get_ltr_score_dups ($scores, $sims, $allfeatures, $partfeatures) {
    my %sccounts;
    my %sicounts;
    my %best_element;

    my $max_score = max(values %$scores);
    my $max_sims  = max(values %$sims);

    for my $source (sort keys %$scores) {
	for my $score_key (sort keys %{$scores->{$source}}) {
	    my $score_value = $scores->{$source}{$score_key};
	    push @{$sccounts{$score_value}}, $score_key;
	}
    }

    for my $source (sort keys %$sims) {
	for my $sim_key (sort keys %{$sims->{$source}}) {
	    my $sim_value = $sims->{$source}{$sim_key};
	    push @{$sicounts{$sim_value}}, $sim_key;
	}
    }

    my ($best_score_key, $best_sim_key);
    for my $src (keys %$scores) {
	$best_score_key = (reverse sort { $scores->{$src}{$a} <=> $scores->{$src}{$b} } keys %{$scores->{$src}})[0];
    }

    for my $src (keys %$sims) {
	$best_sim_key = (reverse sort { $sims->{$src}{$a} <=> $sims->{$src}{$b} } keys %{$sims->{$src}})[0];
    }

    #say "DEBUG: $best_score_key";
    #say "DEBUG: $best_sim_key";
    #dd $scores and exit;
    #dd \%sccounts and exit;
    #dd $sims and exit;
    for my $source (keys %$scores) {
	if (@{$sccounts{ $scores->{$source}{$best_score_key} }} == 1 &&
	    @{$sicounts{ $sims->{$source}{$best_sim_key} }} == 1  &&
	    $best_score_key eq $best_sim_key) {
	    my $bscore = $scores->{$source}{$best_score_key};
	    my $bsim   = $sims->{$source}{$best_score_key};
	    if (exists $partfeatures->{$source}{$best_score_key}) {
		#say join q{ }, "BEST: LTR99 element: ", $best_score_key, $bscore, $bsim;
		push @{$best_element{$source}}, { $best_score_key => $partfeatures->{$source}{$best_score_key} };
	    }
	    elsif (exists $allfeatures->{$source}{$best_score_key}) {
		#say join q{ }, "BEST: LTR85 element: ", $best_score_key, $bscore, $bsim;
		push @{$best_element{$source}}, { $best_score_key => $allfeatures->{$source}{$best_score_key} };
	    }
	    else {
		say "\nERROR: Something went wrong....'$best_score_key' not found in hash. This is a bug. Exiting.";
		exit(1);
	    }
	}
	elsif (@{$sccounts{ $scores->{$source}{$best_score_key} }} >= 1 && 
	       @{$sicounts{ $sims->{$source}{$best_sim_key} }} >=  1) {
	    my $best = @{$sicounts{ $sims->{$source}{$best_sim_key} }}[0];
	    my $bscore = $scores->{$source}{$best_score_key};
	    my $bsim   = $sims->{$source}{$best_score_key};
	    if (exists $partfeatures->{$source}{$best}) {
		#say join q{ }, "BEST: LTR99 element: ", $best, $bscore, $bsim;
		push @{$best_element{$source}},  { $best_score_key => $partfeatures->{$source}{$best} };
	    }
	    elsif (exists $allfeatures->{$source}{$best}) {
		#say join q{ }, "BEST: LTR85 element: ", $best, $bscore, $bsim;
		push @{$best_element{$source}}, { $best_score_key => $allfeatures->{$source}{$best} };
	    }
	    else {
		say "\nERROR: Something went wrong....'$best' not found in hash. This is a bug. Exiting.";
		exit(1);
	    }
	    
	}
	else {
	    say "DEBUG: $best_score_key"; dd $scores;
	    say "DEBUG 2: $best_sim_key"; dd $sims;
	    say "DEBUG 3: "; dd \%sccounts;
	    say "DEBUG 4: "; dd \%sicounts;
	}
    }
    
    return \%best_element;
}

sub summarize_features ($feature) {
    my ($three_pr_tsd, $five_pr_tsd, $ltr_sim);
    my ($has_pbs, $has_ppt, $has_pdoms, $has_ir, $tsd_eq, $tsd_ct) = (0, 0, 0, 0, 0, 0);
    for my $feat (@$feature) {
	my @part = split /\|\|/, $feat;
	if ($part[2] eq 'target_site_duplication') {
	    if ($tsd_ct > 0) {
		$five_pr_tsd = $part[4] - $part[3] + 1;
	    }
	    else {
		$three_pr_tsd = $part[4] - $part[3] + 1;
		$tsd_ct++;
	    }
	}
	elsif ($part[2] eq 'LTR_retrotransposon') {
	    ($ltr_sim) = ($part[8] =~ /ltr_similarity \"(\d+.\d+)\"/); 
	}
	$has_pbs = 1 if $part[2] eq 'primer_binding_site';
	$has_ppt = 1 if $part[2] eq 'RR_tract';
	$has_ir  = 1 if $part[2] eq 'inverted_repeat';
	$has_pdoms++ if $part[2] eq 'protein_match';
    }
    dd $feature and exit(1) unless defined $ltr_sim;
    $tsd_eq = 1 if $five_pr_tsd == $three_pr_tsd;

    my $ltr_score = sum($has_pbs, $has_ppt, $has_pdoms, $has_ir, $tsd_eq);
    return ($ltr_score, $ltr_sim);
}

sub filter_compound_elements ($features) {
    my @pdoms;
    my $is_gypsy   = 0;
    my $is_copia   = 0;
    my $has_pdoms  = 0;
    my $len_thresh = 25000; # are elements > 25kb real? probably not

    for my $source (keys %$features) {
	for my $ltr (nsort_by { m/repeat_region(\d+)/ and $1 } keys %{$features->{$source}}) {
	    my ($rreg, $s, $e, $l) = split /\./, $ltr;
	    my $region = @{$features->{$source}{$ltr}}[0];
	    
	    for my $feat (@{$features->{$source}{$ltr}}) {
		my @feats = split /\|\|/, $feat;
		$feats[8] =~ s/\s\;\s/\;/g;
		$feats[8] =~ s/\s+/=/g;
		$feats[8] =~ s/\s+$//;
		$feats[8] =~ s/=$//;
		$feats[8] =~ s/=\;/;/g;
		$feats[8] =~ s/\"//g;
		if ($feats[2] =~ /protein_match/) {
		    $has_pdoms = 1;
		    @pdoms = ($feats[8] =~ /name=(\S+)/);
		    if ($feats[8] =~ /name=RVT_1|name=Chromo/i) {
			$is_gypsy  = 1;
		    }
		    if ($feats[8] =~ /name=RVT_2/i) {
			$is_copia  = 1;
		    }
		}
	    }
	    
	    if ($is_gypsy && $is_copia) {
		delete $features->{$ltr};
	    }
	    
	    my %uniq;
	    if ($has_pdoms) {
		for my $element (@pdoms) {
		    $element =~ s/\;.*//;
		    next if $element =~ /chromo/i; # we expect these elements to be duplicated
		    delete $features->{$ltr} if $uniq{$element}++;
		}
	    }
	    
	    if ($l >= $len_thresh) {
		delete $features->{$ltr};
	    }
	
	    @pdoms = ();
	    $is_gypsy  = 0;
	    $is_copia  = 0;
	    $has_pdoms = 0;
	}
    }

    return $features;
}

sub get_source {
    my ($ref) = @_;

    #dd $ref and exit;
    for my $feat (@$ref) {
	my @feats = split /\|\|/, $feat;
	return ($feats[0]);
    }
}
