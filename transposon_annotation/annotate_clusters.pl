#!/usr/bin/env perl

## Load JSON file to get hits for each type DONE
## push matches into array, should stay ordered

##TODO deprecate in favor of transposome methods

use 5.012;
use strict;
use warnings;
use autodie qw(open);
use File::Spec qw(catfile rel2abs);
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long;
use JSON;
use Data::Dump qw(dd dump);
use List::Util qw(sum);
use Cwd;

my $usage = "$0 -i cluster_fas_files_dir -d repeat_db -j repbase.json -o annot_clusters_dir -r report.tsv\n";
my $indir; 
my $database;
my $outdir;
my $json;
my $report;

GetOptions(
	   'i|indir=s'  => \$indir,
	   'o|outdir=s' => \$outdir,
	   'd|db=s'     => \$database,
	   'j|json=s'   => \$json,
	   'r|report=s' => \$report,
	   );

#die $usage if !$indir or !$outdir or !$database or !$json or !$report;
if (!$indir || !$database || !$outdir ||
    !$json || !$report) {
    print $usage;
    exit(1);
}

## get input files
opendir(DIR, $indir) || die "\nERROR: Could not open directory: $indir. Exiting.\n";
my @clus_fas_files = grep /\.fa.*$/, readdir(DIR);
closedir(DIR);


if (scalar(@clus_fas_files) < 1) {
    say "\nERROR: Could not find any fasta files in $indir. Exiting.\n";
    exit(1);
}

## set path to output dir
my ($oname, $opath, $osuffix) = fileparse($outdir, qr/\.[^.]*/);
my $out_path = File::Spec->rel2abs($opath.$oname);
make_path($outdir, {verbose => 0, mode => 0711,}); # allows for recursively making paths

my ($dname, $dpath, $dsuffix) = fileparse($database, qr/\.[^.]*/);
my $db_path = File::Spec->rel2abs($dpath.$dname);

my $cwd = getcwd();
my ($rpname, $rppath, $rpsuffix) = fileparse($report, qr/\.[^.]*/);
my $rp_path = File::Spec->rel2abs($rppath.$rpname.$rpsuffix);
open my $out, '>>', $rp_path;

chdir $indir;

say $out join "\t", "Cluster", "Read_count", "Type", "Class", "Superfam", "SINE_family","Top_hit";  ## this is not accurate and needs to be updated

for my $file (sort { $a =~ /CL(\d+)/ <=> $b =~ /CL(\d+)/ } @clus_fas_files) {
    my ($fname, $fpath, $fsuffix) = fileparse($file, qr/\.[^.]*/);
    my $blast_res = $fname;
    my ($filebase, $readct, @blname) = split /\_/, $fname;
    $blast_res =~ s/\.[^.]+$//;
    $blast_res .= "_blast.tsv";
    my $blast_file_path = File::Spec->catfile($out_path, $blast_res);

    my $blastcmd = "blastn -dust no -query $file -evalue 1e-5 -db $db_path -outfmt 6 | ".
	           "cut -f2 | ".                           # keep only the ssids
		   "sort | ".                              # sort the list
		   "uniq -c | ".                           # reduce the list
		   "sort -bnr | ".                         # count unique items
		   "perl -lane 'print join(\"\\t\",\@F)'"; # create an easy to parse format
    my @blast_out = qx($blastcmd);

    my ($hit_ct, $top_hit) = parse_blast(\@blast_out, $blast_file_path);
    next unless defined $top_hit;
    blast2annot($json, $filebase, $readct, $top_hit, $out);
    if (not defined $hit_ct) {
	#say $$hit_ct," sequences have a BLAST hit ",$$top_hit," is the top hit ($$color) in ",$file;
	say "No hits for $file";
    }
}
close $out;

clusters_annotation_to_summary($rp_path);

#
# Subs
#
sub clusters_annotation_to_summary {
    my $report = shift;

    my ($rname, $rpath, $rsuffix) = fileparse($report, qr/\.[^.]*/);
    my $report_sum = $rname."_summary.tsv";
    my $report_sum_path = File::Spec->rel2abs($rpath.$report_sum);
    open my $outsum, '>', $report_sum_path;

    my %annot;
    my %sfam_hitct;
    my %fam_readct;
    my $total_ct = 0;

    open my $in, '<', $report;

    while (<$in>) {
	chomp;
	next if /^Cluster/;
	my @fields = split;
	$total_ct += $fields[1];
	if (scalar @fields == 4) {
	    $sfam_hitct{$fields[2]}++;
	    $fam_readct{$fields[3]} += $fields[1];
	    if (exists $annot{$fields[2]}{$fields[3]}) {
 		push @{$annot{$fields[2]}{$fields[3]}}, $fields[1];
	    }
	    else {
 		$annot{$fields[2]}{$fields[3]} = [$fields[1]];
	    }
	}
	else {
	    my $fam = $fields[5];
	    $fam =~ s/\-\d.*//;
	    $sfam_hitct{$fields[4]}++;
	    $fam_readct{$fam} += $fields[1];
	    if (exists $annot{$fields[4]}{$fam}) {
		push @{$annot{$fields[4]}{$fam}}, $fields[1];
	    }
	    else {
		$annot{$fields[4]}{$fam} = [$fields[1]];
	    }
	}
    }
    close $in;

    say $outsum join "\t", "Superfamily", "Family", "ReadCt/TotalReadCt", "PerCov";

    for my $sf (reverse sort { $sfam_hitct{$a} <=> $sfam_hitct{$b} } keys %sfam_hitct) {
	for my $f (reverse sort { $fam_readct{$a} <=> $fam_readct{$b} } keys %fam_readct) {
	    if (exists $annot{$sf}{$f}) {
		my $read_ct = sum @{$annot{$sf}{$f}};
		my $perc_cov = sprintf("%.12f",$read_ct/$total_ct);
		say $outsum join "\t", $sf, $f, $read_ct."/".$total_ct, $perc_cov;
	    }
	}
    }
    close $outsum;
}

sub blast2annot {
    my ($json, $filebase, $readct, $top_hit, $out) = @_;

    my $repeats = json2hash($json);

    for my $type (keys %$repeats) {
	if ($type eq 'pseudogene' || $type eq 'simple_repeat' || $type eq 'integrated_virus') {
	    if ($type eq 'pseudogene' && $$top_hit =~ /rrna|trna|snrna/i) {
		say $out join "\t", $filebase, $readct, $type, $$top_hit;
	    }
	    elsif ($type eq 'simple_repeat' && $$top_hit =~ /msat/i) {
		say $out join "\t", $filebase, $readct, $type, "Satellite", "MSAT", $$top_hit;
	    }
	    elsif ($type eq 'simple_repeat' && $$top_hit =~ /sat/i) {
		say $out join "\t", $filebase, $readct, $type, "Satellite", "SAT", $$top_hit;
	    }                                                                                                                                              
	    elsif ($type eq 'integrated_virus' && $$top_hit =~ /caul/i) {
		say $out join "\t", $filebase, $readct, $type, "Caulimoviridae", $$top_hit;
	    }
	    elsif ($type eq 'integrated_virus' && ($$top_hit eq 'PIVE' || $$top_hit eq 'DENSOV_HM')) {
		say $out join "\t", $filebase, $readct, $type, "DNA Virus", $$top_hit;
	    }
	    elsif ($type eq 'endogenous_retrovirus' && $$top_hit =~ /erv/i) {
		say $out join "\t", $filebase, $readct, $type, "Endogenous Retrovirus", $$top_hit;
	    }                                                                   
	    next;
	}
	for my $class (keys %{$repeats->{$type}}) {
	    while ( my ($superfam_index, $superfam) = each @{$repeats->{$type}{$class}} ) {
		for my $superfam_h (keys %$superfam) {
		    if ($superfam_h =~ /sine/i) {
			while (my ($sine_fam_index, $sine_fam_h) = each @{$superfam->{$superfam_h}}) {
			    for my $sine_fam_mem (keys %$sine_fam_h) {
				for my $sines (@{$repeats->{$type}{$class}[$superfam_index]{$superfam_h}[$sine_fam_index]{$sine_fam_mem}}) {
				    for my $sine (@$sines) {
					if ($sine =~ /$$top_hit/) {
					    say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $sine_fam_mem, $$top_hit;
					}
				    }
				}
			    }
			}
		    } # if /sine/i
		    elsif ($superfam_h =~ /gypsy/i && $$top_hit =~ /^RLG/) {
			say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $$top_hit;
		    }
		    elsif ($superfam_h =~ /copia/i && $$top_hit =~ /^RLC/) {
			say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $$top_hit;
		    }
		    else {
			for my $fam (@{$repeats->{$type}{$class}[$superfam_index]{$superfam_h}}) {
			    for my $mem (@$fam) {
				if ($mem =~ /$$top_hit/i) {
				    say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $$top_hit;
				}
				#else { say "NOMATCH $mem => $$top_hit"; }
			    }
			}
		    }
		}
	    }
	}
    }
    #close($out);
}

sub json2hash {
    my $json = shift;
   
    my $json_text;
    local $/;
    open my $in, '<', $json;
    $json_text = <$in>;
    close $in;

    my $repeats = JSON->new->utf8->space_after->decode($json_text);
    return $repeats;
}

sub parse_blast {
    my ($blast_out, $blast_file_path) = @_;
    my %blhits;

    my $top_hit;
    my $top_hit_num = 0;
    my $hit_ct = 0;

    for my $hit (@$blast_out) {
	chomp $hit;
	my ($ct, $hittype) = split /\t/, $hit;
	next unless defined $ct;
	$blhits{$hittype} = $ct;
	$hit_ct++;
    }

    #dump(%blhits);

    ## Need to create an array of colors to pass to R barplot   
    if ($hit_ct > 0) {
	open my $out, '>', $blast_file_path;
	$top_hit = (reverse sort { $blhits{$a} <=> $blhits{$b} } keys %blhits)[0];
	for my $hits (reverse sort { $blhits{$a} <=> $blhits{$b} } keys %blhits) {
	    say $out join "\t", $hits, $blhits{$hits};
	}
	close $out;
	return \$hit_ct, \$top_hit;
    }
    else {
	return undef, undef;
    }
}


