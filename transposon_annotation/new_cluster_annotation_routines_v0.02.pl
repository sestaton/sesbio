#!/usr/bin/env perl

use utf8;
use v5.12;
use strict;
use warnings;
use warnings FATAL => "utf8";
use open qw(:std :utf8);
use autodie qw(open);
use Getopt::Long;
use Capture::Tiny qw(:all);
use File::Basename;
use File::Path qw(make_path);
use File::Spec;
use File::Basename;
use POSIX qw(strftime);
use Graph::UnionFind;
use List::Util qw(sum max);
use JSON;
use Cwd;
BEGIN {
  @AnyDBM_File::ISA = qw( DB_File SQLite_File )
      unless @AnyDBM_File::ISA == 1;
}
use AnyDBM_File;                  
use vars qw( $DB_BTREE &R_DUP );  
use AnyDBM_File::Importer qw(:bdb);
use DBM::Deep;
use charnames qw(:full :short);
use Encode qw(encode decode);

## for debugging
use Data::Dump qw(dd);

# lexical vars
my $usage = "$0 -i cls_dir_with_merges -o annotations -r report -d repeat_blast_db -n seqct -j repbase1801_full.json [-e 1e-5]\n\n";
my $indir;
my $outdir;
my $report;
my $database;
my $json;
my $evalue;
my $seqct;

GetOptions(
           'i|indir=s'           => \$indir,
	   'o|outdir=s'          => \$outdir,
           'r|report=s'          => \$report,
           'd|database=s'        => \$database,
           'j|repbase_json=s'    => \$json,
           'e|evalue=f'          => \$evalue,
	   'n|seqct=i'           => \$seqct,
           );

# open the infile or die with a usage statement
if (!$indir || !$outdir || !$database || !$report || !$json || !$seqct) {
    print "\nERROR: Input not parsed correctly.\n";
    print $usage;
    exit(1);
}

## annotate_clusters
#annotate_clusters($cls_dir_path, $database, $json, $report, $outdir, $evalue); ## annotate_clusters sub as it is called from repeat_explorer_v0.14.pl
annotate_clusters($indir, $database, $json, $report, $outdir, $evalue, $seqct); 

exit; ## This is the end.

#
# subs
#
sub annotate_clusters {
    my ($cls_with_merges_dir, $database, $json, $report, $outdir, $evalue, $seqct) = @_;

    my ($rpname, $rppath, $rpsuffix) = fileparse($report, qr/\.[^.]*/);
    my $anno_rep = $rpname."_annotations.tsv";
    my $anno_summary_rep = $rpname."_annotations_summary.tsv";
    my $anno_rp_path = File::Spec->rel2abs($rppath.$anno_rep);
    my $anno_sum_rep_path = File::Spec->rel2abs($rppath.$anno_summary_rep);
    my $total_readct = 0;
    $evalue //= 10;
    my $top_hit_superfam = {};

    ## get input files
    opendir my $dir, $outdir."/".$cls_with_merges_dir || die "\nERROR: Could not open directory: $cls_with_merges_dir. Exiting.\n";
    my @clus_fas_files = grep /\.fa.*$/, readdir $dir;
    closedir $dir;
    

    if (scalar @clus_fas_files < 1) {
        say "\nERROR: Could not find any fasta files in $cls_with_merges_dir. Exiting.\n";
        exit(1);
    }

    ## set path to output dir
    my $annodir = $outdir."/".$cls_with_merges_dir."_annotations";
    my ($oname, $opath, $osuffix) = fileparse($annodir, qr/\.[^.]*/);
    my $out_path = File::Spec->rel2abs($opath.$oname);
    make_path($annodir, {verbose => 0, mode => 0711,}); # allows for recursively making paths

    my ($dname, $dpath, $dsuffix) = fileparse($database, qr/\.[^.]*/);
    my $db_path = File::Spec->rel2abs($dpath.$dname);
    my @blasts; # container for each report (hash) // need to rethink how duplicate entries will be handled
    my @superfams;

    open my $out, '>>', $anno_rp_path;
    chdir $cls_with_merges_dir;

    say $out join "\t", "Cluster", "Read_count", "Type", "Class", "Superfam", "(SINE_family; if present)","Top_hit";  
    for my $file (@clus_fas_files) {
        my $query = $outdir."/".$cls_with_merges_dir."/".$file;
        my ($fname, $fpath, $fsuffix) = fileparse($query, qr/\.[^.]*/);
        my $blast_res = $fname;
        my ($filebase, $readct) = split /\_/, $fname, 2;
        $total_readct += $readct;
        $blast_res =~ s/\.[^.]+$//;
        $blast_res .= "_blast_$evalue.tsv";
        my $blast_file_path = File::Spec->catfile($out_path, $blast_res);

        my $blastcmd = "blastn -dust no -query $query -evalue $evalue -db $db_path -outfmt 6 | ".
	               "sort -k1,1 -u | ".                       # count each read in the report only once
                       "cut -f2 | ".                             # keep only the ssids
                       "sort | ".                                # sort the list
                       "uniq -c | ".                             # reduce the list
                       "sort -bnr | ".                           # count unique items
                       "perl -lane 'print join(\"\\t\",\@F)'";   # create an easy to parse format
	#say $blastcmd;
        my @blast_out = qx($blastcmd);
	
        my ($hit_ct, $top_hit, $blhits) = parse_blast_to_top_hit(\@blast_out, $blast_file_path);
        next unless defined $top_hit && defined $hit_ct;
	#say $file, " ==> ", $blast_res, " ==> ", dd $blhits; 
	push @blasts, $blhits;
        $top_hit_superfam = blast2annot($json, $filebase, $readct, $top_hit, $out);
	push @superfams, $top_hit_superfam unless !%$top_hit_superfam;
    }
    close $out;

    clusters_annotation_to_summary($anno_rp_path, $anno_sum_rep_path, $total_readct, $seqct, \@blasts, \@superfams);
}

sub parse_blast_to_top_hit {
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

    if ($hit_ct > 0) {
        open my $out, '>', $blast_file_path;
        $top_hit = (reverse sort { $blhits{$a} <=> $blhits{$b} } keys %blhits)[0];
	keys %blhits; #reset iterator
        for my $hits (reverse sort { $blhits{$a} <=> $blhits{$b} } keys %blhits) {
            say $out join "\t", $hits, $blhits{$hits};
        }
        close $out;
        return \$hit_ct, \$top_hit, \%blhits;
    }
    else { ## if (!%blhits) {
        unlink $blast_file_path;
        return undef, undef, undef;
    }
}

sub blast2annot {
    my ($json, $filebase, $readct, $top_hit, $out) = @_;

    my %top_hit_superfam;
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
					    $top_hit_superfam{$$top_hit} = $sine_fam_mem;
                                            say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $sine_fam_mem, $$top_hit;
                                        }
                                    }
                                }
                            }
                        }
                    } 
                    elsif ($superfam_h =~ /gypsy/i && $$top_hit =~ /^RLG/) {
			$top_hit_superfam{$$top_hit} = $superfam_h;
                        say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $$top_hit;
                    }
                    elsif ($superfam_h =~ /copia/i && $$top_hit =~ /^RLC/) {
			$top_hit_superfam{$$top_hit} = $superfam_h;
                        say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $$top_hit;
                    }
                    else {
                        for my $fam (@{$repeats->{$type}{$class}[$superfam_index]{$superfam_h}}) {
                            for my $mem (@$fam) {
                                if ($mem =~ /$$top_hit/i) {
				    $top_hit_superfam{$$top_hit} = $superfam_h;
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
    return \%top_hit_superfam;
}

sub clusters_annotation_to_summary  {
    my ($anno_rp_path, $anno_sum_rep_path, $total_readct, $seqct, $blasts, $superfams) = @_;
    
    my %top_hit_superfam;
    @top_hit_superfam{keys %$_} = values %$_ for @$superfams;
   
    #dd \%top_hit_superfam;
    for my $f (keys %top_hit_superfam) {
	if ($f =~ /^((RL.\-\w+)\-\d)/) {
	    my $fam = $2;
	    $top_hit_superfam{$fam} = $top_hit_superfam{$f};
	    delete $top_hit_superfam{$f};
	}
    }
    #dd \%top_hit_superfam;
    open my $outsum, '>', $anno_sum_rep_path;
    
    my %annot;
    my %fams;
    my $total_ct = 0;
    my $hashct = scalar @$blasts;
    my $hitct;
    for my $blast (@$blasts) {
	for my $fam (keys %$blast) {
	    $total_ct += $blast->{$fam};
	    if ($fam =~ /^((RL.\-\w+)\-\d)/) {
		my $famname = $2;
		if (exists $fams{$famname}) {
		    $fams{$famname} += $blast->{$fam};
		}
		else {
		    $fams{$famname} = $blast->{$fam};
		}
	    }
	    else {
		if (exists $fams{$fam}) {
		    $fams{$fam} += $blast->{$fam};
		}
		else {
		    $fams{$fam} = $blast->{$fam};
		}
	    }
	}
    }
    #say "total_ct: $total_ct\tTotal hashct: $hashct\thitct: $hitct";
    #dd \%fams;
    my $total_gcov = 0;
                                                          ### NEED TO SIMPLIFY THE STATS REPORTED, once I know they're correct
    say $outsum join "\t", "ReadNum", "Superfamily", "Family", "ReadCt/AllReads", "ReadCt/ReadsWithHit", "HitPerc", "GPerc";
    for my $k (reverse sort { $fams{$a} <=> $fams{$b} } keys %fams) {
	$total_gcov += $fams{$k};
	if (exists $top_hit_superfam{$k}) {
	    my $gperc = sprintf("%.12f",$fams{$k}/$seqct);
	    my $with_hit_perc = sprintf("%.12f",$fams{$k}/$total_ct);
	    say $outsum join "\t", $seqct, $top_hit_superfam{$k}, $k, $fams{$k}."/".$seqct, $fams{$k}."/".$total_ct, $with_hit_perc, $gperc;
	}
    }
    #my $total_gcov_perc = sprintf("%.12f", $total_gcov/$seqct);
    #my $total_hit_perc = sprintf("%.12f", $total_gcov/$total_ct);
    #say $outsum "=" x 60;
    #say $outsum "Cumulative genomic percentage: ",$total_gcov_perc," ($total_gcov"."/"."$seqct)";
    #say $outsum "Cumulative hit percentage    : ",$total_hit_perc," ($total_gcov"."/"."$total_ct)";
    close $outsum;
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
