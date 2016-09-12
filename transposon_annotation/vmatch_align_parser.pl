#!/usr/bin/env perl

use 5.020;
use warnings;

while (my $line = <>) {
    chomp $line;
    next if $line =~ /^#/;
    $line =~ s/^\s+//;
    if ($line =~ /^\d+/) { 
	#simply skip alignments
        my ($slen, $sid, $spos, $htype, $qlen, $hid, $qpos, $dist, $evalue, $score, $pid) = split /\s+/, $line;
        if ($pid >= 90 && $slen >= 50 && $qlen >= 50) {
            my ($code) = ($sid =~ /^(\w{3})-?_?/);
            say $code;
        }
    }
    #my ($hits, $subj, $query) = split /\n/, $line;
    #$hits =~ s/^\s+//;
    #my ($slen, $sid, $spos, $htype, $qlen, $hid, $qpos, $dist, $evalue, $score, $pid) = split /\s+/, $hits;
    #my ($sstring, $send) = ($subj =~ /Sbjct: (\w+)\s+(\d+)/);
    #my ($qstring, $qend) = ($query =~ /Query: (\w+)\s+(\d+)/);
    #push @{$hits{$hid}},
        #join "||", $slen, $sid, $spos, $htype, $qlen, $hid, $qpos, $dist, $evalue, $score, $pid, $sstring, $send, $qstring, $qend;
        
}
