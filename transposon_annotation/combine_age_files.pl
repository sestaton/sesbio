#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use POSIX       qw(strftime);
use Cwd         qw(abs_path);
use Time::HiRes qw(gettimeofday);
use File::Basename;
use File::Spec;
use Sort::Naturally;
use Tephra::Annotation::Util;

my $usage = "$0 genome age_file age_file age_fas age_fas\n";
my $genome = shift or die $usage;
my $ltrage = shift or die $usage;
my $tirage = shift or die $usage;
my $ltrfas = shift or die $usage;
my $tirfas = shift or die $usage;

my @ages = ($ltrage, $tirage);
my @fastas = ($ltrfas, $tirfas);

my ($ct, $age_sum) = combine_age_files($genome, \@ages, \@fastas);
say STDERR "$ct LTR/TIR/TRIM elements logged to $age_sum";

sub combine_age_files {
    my ($genome, $ages, $fastas) = @_;

    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    my $age_sum = File::Spec->catfile( abs_path($path), $name.'_tephra_ltr-tir_age_summary.tsv' );
    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);

    my $util = Tephra::Annotation::Util->new;

    my ($famname, %families, %ages);
    for my $fasta (@$fastas) {
        open my $in, '<', $fasta or die "\nERROR: Could not open file: $fasta\n";
        
        while (my $line = <$in>) {
            chomp $line;
            if ($line =~ /^>([A-Z]{3}(?:_singleton_)?(?:_?family\d+)?)_(?:LTR|TRIM|terminal)/) {
                $famname = $1;
                $line =~ s/>//;
                push @{$families{$famname}}, $line;
            }
        }
        close $in;
    }

    my $tetype;
    for my $agefile (@$ages) { 
        open my $tab, '<', $agefile or die "\nERROR: Could not open file: $agefile\n";
        
        while (my $line = <$tab>) {
            chomp $line;
            next if $line =~ /^(?:LTR|TIR)-ID/;
            #LTR-ID Divergence Age Ts:Tv
            my ($id, $div, $age, $tstv) = split /\t/, $line;
            if ($id =~ /^([A-Z]{3}(?:_singleton_)?(?:_?family\d+)?)_(LTR|TRIM|terminal)/) {
                my $fam = $1;
		my $type = $2;
		if ($type =~ /ltr/i) {
		    $tetype = 'LTR';
		}
		elsif ($type =~ /trim/i) {
		    $tetype = 'TRIM';
		}
		elsif ($type =~ /terminal/i) {
		    $tetype = 'TIR';
		}

		if (exists $families{$fam}) {
                    my $famsize = @{$families{$fam}};
                    push @{$ages{$famsize}{$fam}}, join "||", $id, $div, $age, $tstv, $tetype;
                }
            }
        }
        close $tab;
    }

    open my $out, '>', $age_sum or die "\nERROR: Could not open file: $age_sum\n";
    say $out join "\t", 'Type', 'Superfamily', 'Family', 'Family_size', 'ElementID', 'Divergence', 'Age', 'Ts:Tv';

    my $ct = 0;
    for my $famsize (reverse sort { $a <=> $b } keys %ages) {
        for my $fam (nsort keys %{$ages{$famsize}}) { 
            for my $agestr (@{$ages{$famsize}{$fam}}) {
                $ct++;
                my ($id, $div, $age, $tstv, $type) = split /\|\|/, $agestr;
                my $sfam = $util->map_superfamily_name($fam);
                $sfam = $sfam ? $sfam : 'Unknown'; # the method above returns 0 if the superfamily is unknown
                say $out join "\t", $type, $sfam, $fam, $famsize, $id, $div, $age, $tstv;
            }
        }
    }
    close $out;
    my $age_pd = ' ' x length($ct);
    say STDERR "Command - Generating combined age files at: $st.";
    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    say STDERR "Results - Finished generating combined age files for $ct transposons at: $st.";

    return ($ct, $age_sum);
}
