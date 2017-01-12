#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use List::Util qw(sum);
use Data::Dump::Color;

my $file = shift or die $!;
open my $in, '<', $file or die $!;
my $document = do {
    local $/ = undef;
    open my $fh, "<", $file
        or die "could not open $file: $!";
    <$fh>;
};

my @struct = eval $document;
my $genome_length = 166489593;

#dd \@struct and exit;

write_masking_results(\@struct, $genome_length);

sub write_masking_results {
    my ($reports, $genome_length) = @_;
    my $outfile = 'outfile';

    my %final_rep;
    my $repeat_map = _build_repeat_map();

    ## build offsets here
    my ($classlen, $orderlen, $namelen, $masklen);
    for my $report (@$reports) {
        for my $rep_type (keys %$report) {
            my $total = sum(@{$report->{$rep_type}});
            my ($class, $order, $name) = @{$repeat_map->{$rep_type}}{qw(class order repeat_name)} ;
	    ($classlen, $orderlen, $namelen) = (length($class), length($order), length($name)); 
            $final_rep{$class}{$order}{$name} = $total;
        }
    }
    
    ($classlen,$orderlen, $namelen) = ($classlen+10, $orderlen+10, $namelen+15);
    my $totallen = $namelen+100;
    my $masked_total = 0;
    say "========================== 'Tephra maskref' masking results for: $outfile";
    printf "%-${classlen}s %-${classlen}s %-${orderlen}s %-${namelen}s\n", "Class", "Order", "Superfamily", "Percent Masked";
    #print pack '(A10)*', ("Class", "Order", "Superfamily", "Percent Masked"), "\n";
    say "-" x 80;
    for my $class (keys %final_rep) {
        for my $order (keys %{$final_rep{$class}}) {
            for my $name (keys %{$final_rep{$class}{$order}}) {
                $masked_total += $final_rep{$class}{$order}{$name};
                my $repmasked = sprintf("%.2f",($final_rep{$class}{$order}{$name}/$genome_length)*100);
                printf "%-${classlen}s %-${classlen}s %-${orderlen}s %-${namelen}s\n",
		    $class, $order, $name, "$repmasked% ($final_rep{$class}{$order}{$name}/$genome_length)";
		say length("$repmasked% ($final_rep{$class}{$order}{$name}/$genome_length)");
            }
        }
    }

    my $masked = sprintf("%.2f",($masked_total/$genome_length)*100);
    say "-" x 80;
    my $res = "$masked% ($masked_total/$genome_length)";
    my $reslen = 55 - length($res);
    say $reslen;
    printf "%-15s %-55s\n", "Total masked bases:", "$masked% ($masked_total/$genome_length)";

    return;
}

sub _build_repeat_map {
    my %repeat_map = (
        ## Class I
        # DIRS
        'RLD' => { class => 'Class I', order => 'DIRS', repeat_name => 'DIRS' },
        'RYN' => { class => 'Class I', order => 'DIRS', repeat_name => 'Ngaro' },
        'RYX' => { class => 'Class I', order => 'DIRS', repeat_name => 'Unknown DIRS' },
        'RYV' => { class => 'Class I', order => 'DIRS', repeat_name => 'VIPER' },
        # LINE 
        'RII' => { class => 'Class I', order => 'LINE', repeat_name => 'I' },
        'RIJ' => { class => 'Class I', order => 'LINE', repeat_name => 'Jockey' },
        'RIL' => { class => 'Class I', order => 'LINE', repeat_name => 'L1' },
        'RIR' => { class => 'Class I', order => 'LINE', repeat_name => 'R2' },
        'RIT' => { class => 'Class I', order => 'LINE', repeat_name => 'RTE' },
        'RIX' => { class => 'Class I', order => 'LINE', repeat_name => 'Unknown LINE' },
        'RIC' => { class => 'Class I', order => 'LINE', repeat_name => 'CR1' },
        # LTR
        'RLB' => { class => 'Class I', order => 'LTR', repeat_name => 'Bel/Pao' },
        'RLC' => { class => 'Class I', order => 'LTR', repeat_name => 'Copia' },
        'RLE' => { class => 'Class I', order => 'LTR', repeat_name => 'ERV' },
        'RLG' => { class => 'Class I', order => 'LTR', repeat_name => 'Gypsy' },
        'RLR' => { class => 'Class I', order => 'LTR', repeat_name => 'Retrovirus' },
        'RLX' => { class => 'Class I', order => 'LTR', repeat_name => 'Unknown LTR' },
        # PLE
        'RPP' => { class => 'Class I', order => 'Penelope', repeat_name => 'Penelope' },
        'RPX' => { class => 'Class I', order => 'Penelope', repeat_name => 'Unknown PLE' },
	# SINE
        'RSS' => { class => 'Class I', order => 'SINE', repeat_name => '5S' },
        'RSL' => { class => 'Class I', order => 'SINE', repeat_name => '7SL' },
        'RST' => { class => 'Class I', order => 'SINE', repeat_name => 'tRNA' },
        'RSX' => { class => 'Class I', order => 'SINE', repeat_name => 'Unknown SINE' },
        'RXX' => { class => 'Class I', order => 'SINE', repeat_name => 'Unknown retrotransposon' },
        ## Class II
        # - Subclass 1
        # Crypton
        'DYC' => { class => 'Class II', order => 'Crypton', repeat_name => 'Crypton' },
        'DYX' => { class => 'Class II', order => 'Crypton', repeat_name => 'Unknown Crypton' },
        # TIR
        'DTC' => { class => 'Class II', order => 'TIR', repeat_name => 'CACTA' },
        'DTA' => { class => 'Class II', order => 'TIR', repeat_name => 'hAT' },
        'DTE' => { class => 'Class II', order => 'TIR', repeat_name => 'Merlin' },
        'DTM' => { class => 'Class II', order => 'TIR', repeat_name => 'Mutator' },
        'DTP' => { class => 'Class II', order => 'TIR', repeat_name => 'P' },
        'DTH' => { class => 'Class II', order => 'TIR', repeat_name => 'PIF/Harbinger' },
        'DTB' => { class => 'Class II', order => 'TIR', repeat_name => 'PiggyBac' },
        'DTT' => { class => 'Class II', order => 'TIR', repeat_name => 'Tc1/Mariner' },
        'DTR' => { class => 'Class II', order => 'TIR', repeat_name => 'Transib' },
        'DTX' => { class => 'Class II', order => 'TIR', repeat_name => 'Unknown TIR' },
        'DXX' => { class => 'Class II', order => 'TIR', repeat_name => 'Unknown DNA transposon' },
        # - Subclass 2
        # Helitron
        'DHH' => { class => 'Class II', order => 'Helitron', repeat_name => 'Helitron' },
        'DHX' => { class => 'Class II', order => 'Helitron', repeat_name => 'Unknown Helitron' },
        # Maverick
        'DMM' => { class => 'Class II', order => 'Maverick', repeat_name => 'Maverick' },
        'DMX' => { class => 'Class II', order => 'Maverick', repeat_name => 'Unknown Maverick' },
        );

    return \%repeat_map;
}
