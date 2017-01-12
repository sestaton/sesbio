#!/usr/bin/env perl

use 5.020;
use warnings;
use Time::HiRes qw(gettimeofday);
use List::Util  qw(sum);
use Set::IntervalTree;
use File::Basename;
use Data::Dump::Color;

my $usage = "$0 vmat.out";
my $file = shift or die $usage;
my (%windows, %refs, %hits, %aligns, %report, %final_rep, @reports);

my $repeat_map = _build_repeat_map();

open my $in, '<', $file or die $!;
my $t0 = gettimeofday();
my $tree = Set::IntervalTree->new;
my $comm = <$in>; # discard the args to vmatch
my $line = <$in>;
# $line is the following format: 
# l(S) h(S) r(S) t l(Q) h(Q) r(Q) d e s i
# where:
# l = length
# h = sequence header
# r = relative position
# t = type (D=direct, P=palindromic)
# d = distance value (negative=hamming distance, 0=exact, positive=edit distance)
# e = E-value
# s = score value (negative=hamming score, positive=edit score)
# i = percent identity
# (S) = in Subject
# (Q) = in Query
chomp $line;
$line =~ s/^\s+//;
$line =~ s/^\s+//;
my @f = split /\s+/, $line;
my $send = $f[2] + $f[0];
$tree->insert({ id => $f[1], match => $f[5], start => $f[2], end => $send, len => $f[0] }, $f[2], $send);
$windows{$f[2]} = { start => $f[2], end => $send, len => $f[0], overlap => 0, match => $f[5] };

while (my $line = <$in>) {
    chomp $line;
    $line =~ s/^\s+//;
    my ($slen, $sid, $spos, $htype, $qlen, $hid, $qpos, $dist, $evalue, $score, $pid) = split /\s+/, $line;
    my $end = $spos + $slen;

    my $res = $tree->fetch($spos, $end);
    
    if (@$res) {
	my ($best_start, $best_end, $overl, $best_match);
	for my $overlap (@$res) {	    
	    my ($ostart, $oend, $match, $subj, $olen) = @{$overlap}{qw(start end match id len)};

	    $best_start = $ostart <= $spos ? $ostart : $spos;
	    $best_end = $oend >= $end ? $oend : $end;
	    my $oe = $best_end == $oend ? $end : $oend;
	    $overl = $best_end - $oe;

	    $tree->remove($ostart, $oend);
	}
	
	my $nlen = $best_start > 0 ? $best_end-$best_start : $best_end;
	$windows{$best_start} = { start => $best_start, end => $best_end, len => $nlen, match => $hid, overlap => $overl };
	$tree->insert({ id => $sid, match => $hid, start => $best_start, end => $best_end, len => $nlen }, $best_start, $best_end);
    }
    else {
	$tree->insert({ id => $sid, match => $hid, start => $spos, end => $end, len => $slen }, $spos, $end);
        $windows{$spos} = { start => $spos, end => $end, len => $slen, overlap => 0, match => $hid};
    }
}
#exit;
#dd \%windows and exit;

my @del;
for my $s (sort { $a <=> $b } keys %windows) {
    if (exists $windows{ $windows{$s}{end} }) {
	my $h = $windows{ $windows{$s}{end} };
	$windows{ $windows{$s}{end} } = { start   => $h->{start}+1, 
					  end     => $h->{end},
					  match   => $h->{match},
					  len     => $h->{len}-1,
					  overlap => 0 };
    }
}
#exit;
#delete $windows{$_} for @del;

#exit;
#dd \%windows and exit;
my %seen;
for my $s (sort { $a <=> $b } keys %windows) { 
    #say join q{ }, $s, $windows{$s}{end}, $windows{$s}{len};
    my ($code) = ($windows{$s}{match} =~ /^(\w{3})-?_?/);         
    if (defined $code && exists $repeat_map->{$code}) {
	push @{$report{ $code }}, $windows{$s}{len};
    }
}
#exit;


push @reports, \%report;
write_masking_results(\@reports, $t0, $repeat_map);
#
#
sub get_best_overlap {
    my ($start, $end, $len, $res) = @_;

    for my $overlap (@$res) {
	my ($ostart, $oend, $match, $subj, $olen) = @{$overlap}{qw(start end match id len)};
	my $overl = $end - $oend;
	next if $start >= $ostart && $end <= $oend;
	my $nlen = $end - $ostart + 1;
	
	#if ($end > $oend) {
	
	#}
    }
    
}

sub write_masking_results {
    my ($reports, $t0, $repeat_map) = @_;
    my $genome_length = 10e5; #119871;
    my $outfile = 'foo';

    my %final_rep;

    my ($classlen, $orderlen, $namelen, $masklen);
    for my $report (@$reports) {
        for my $rep_type (keys %$report) {
            my $total = sum(@{$report->{$rep_type}});
            my ($class, $order, $name) = @{$repeat_map->{$rep_type}}{qw(class order repeat_name)} ;
            ($classlen, $orderlen, $namelen) = (length($class), length($order), length($name)); 
            $final_rep{$class}{$order}{$name} = $total;
        }
    }
    
    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    ($classlen,$orderlen, $namelen) = ($classlen+10, $orderlen+10, $namelen+15);
    my $masked_total = 0;
    say "=================== 'Tephra maskref' finished in $final_time minutes ==================";
    printf "%-${classlen}s %-${classlen}s %-${orderlen}s %-${namelen}s\n", "Class", "Order", "Superfamily", "Percent Masked";
    #print pack '(A10)*', ("Class", "Order", "Superfamily", "Percent Masked"), "\n";
    say "-" x 80;
    for my $class (sort keys %final_rep) {
        for my $order (sort keys %{$final_rep{$class}}) {
            for my $name (sort keys %{$final_rep{$class}{$order}}) {
                $masked_total += $final_rep{$class}{$order}{$name};
                my $repmasked = sprintf("%.2f",($final_rep{$class}{$order}{$name}/$genome_length)*100);
                printf "%-${classlen}s %-${classlen}s %-${orderlen}s %-${namelen}s\n",
		$class, $order, $name, "$repmasked% ($final_rep{$class}{$order}{$name}/$genome_length)";
            }
        }
    }

    my $masked = sprintf("%.2f",($masked_total/$genome_length)*100);
    say "=" x 80;
    say "Input file:         genome";
    say "Output file:        $outfile";
    say "Total masked bases: $genome_length";
    say "Total masked bases: $masked% ($masked_total/$genome_length)";
    #my $totlen = $classlen*2 + $orderlen + $namelen;
    #my $res = "$masked% ($masked_total/$genome_length)";
    #my $reslen = length($res);
    #say $reslen;
    #printf "%-15s %51s\n", "Total masked bases:", "$masked% ($masked_total/$genome_length)";
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
