#!/usr/bin/env perl

use 5.020;
use warnings;
use Data::Printer;

local $/ = '';
my %hits;
while (my $line = <>) {
    chomp $line;
    next if $line =~ /^#/;
    my ($hits, $subj, $query) = split /\n/, $line;
    $hits =~ s/^\s+//;
    my ($slen, $sid, $spos, $htype, $qlen, $hid, $qpos, $dist, $evalue, $score, $pid) = split /\s+/, $hits;
    my ($sstring, $send) = ($subj =~ /Sbjct: (\w+)\s+(\d+)/);
    my ($qstring, $qend) = ($query =~ /Query: (\w+)\s+(\d+)/);
    push @{$hits{$hid}},
        join "||", $slen, $sid, $spos, $htype, $qlen, $hid, $qpos, $dist, $evalue, $score, $pid, $sstring, $send, $qstring, $qend;
        
}

p %hits;
