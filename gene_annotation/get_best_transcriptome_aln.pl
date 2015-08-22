#!/usr/bin/env perl

use 5.020;
use strict;
use warnings;
use autodie;
use File::Find;
use File::Basename;
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(:all);
use Try::Tiny;
use experimental 'signatures';
use Data::Dump;

my (@elite, @weedy, @wild);
my $dir           = 'scaffold_500_variant_calls_bams_assemblies';
my $blastdb       = 'epsps_lg4_lg16_complete_cds';
my $combined_orfs = 'targeted_scaffold_500_denovo_epsps_proteins_orfs.faa';
my $combined_alns = 'targeted_scaffold_500_denovo_epsps_proteins_orfs.aln';

find( sub { push @elite, $File::Find::name if -f and /^elite/ and /clean.fasta$/ }, $dir );
find( sub { push @weedy, $File::Find::name if -f and /^weedy/ and /clean.fasta$/ }, $dir );
find( sub { push @wild,  $File::Find::name if -f and /^wild/  and /clean.fasta$/ }, $dir );

#my @top_weeds = grep { ! /pasa/ } @weedy; # we just want things in the top dir
#my @top_elite = grep { ! /pasa/ } @elite;
#my @top_wild  = grep { ! /pasa/ } @wild;
#dd \@top_wild and exit;

process_assemblies(\@weedy, $blastdb, $combined_orfs);
process_assemblies(\@elite, $blastdb, $combined_orfs);
process_assemblies(\@wild,  $blastdb, $combined_orfs);
align_all($combined_orfs, $combined_alns);

#
# Methods
#
sub process_assemblies ($arr_ref, $blastdb, $combined_orfs) {
    for my $file (@$arr_ref) {
	my $orffile = get_orfs($file);
	my $best_hitid = run_blast($orffile, $blastdb);
	my $orfs = select_bestorf($orffile, $best_hitid, $combined_orfs);
	#align_all($orfs, $combined_alns);
    }
}

sub align_all ($orfs, $combined_alns) {
    my $cmd = "muscle -in $orfs -out $combined_alns";
    try {
        my ($stdout, $stderr, @res) = capture { system([0..5], $cmd); };
    }
    catch {
        say "\nERROR: $cmd failed. Here is the exception: $_\n";
    };
}

sub select_bestorf ($file, $best_hitid, $combined_orfs) {
    my ($ffile, $fdir, $fext) = fileparse($file, qr/\.[^.]*/);
    my ($classification, $acc) = ($ffile =~ /^(weedy|wild|elite)\.(\w+)_k31/);

    #say "DEBUG: $classification $acc";
    my @aux = undef;
    my ($name, $comm, $seq, $qual);
    open my $out, '>>', $combined_orfs;
    open my $in, '<', $file;
    
    while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
	if ($name eq $best_hitid) {
	    #say "Match!";
	    my $id = ">$classification"."_".$acc."_".$name;
	    say $out join "\n", $id, $seq;
	}
    }
    close $in;
    close $out;

    return $combined_orfs if -s $combined_orfs;
}

sub get_orfs ($query) {
    #my ($ffile, $fdir, $fext) = fileparse($query, qr/\.[^.]*/);
    #my $orffile = $ffile."_orfs.faa";
    my $orffile = $query;
    $orffile =~ s/\.fasta/_orfs.faa/;

    my $cmd = "hmmer2go getorf -i $query -o $orffile";
    try {
        my ($stdout, $stderr, @res) = capture { system([0..5], $cmd); };
    }
    catch {
        say "\nERROR: $cmd failed. Here is the exception: $_\n";
    };

    return $orffile;
}

sub run_blast ($query, $blastdb) {
    my $cmd = "blastall -p blastp -i $query -d $blastdb -m 8 | ";
    $cmd .= "sort -nrk 12 |sort -k1,1 -u |head -1 | cut -f1";

    my $bestid;
    try {
	my ($stdout, $stderr, @res) = capture { system([0..5], $cmd); };
	chomp $stdout;
	$bestid = $stdout;
    }
    catch {
	say "\nERROR: $cmd failed. Here is the exception: $_\n";
    };

    return $bestid;
}

sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!@$aux);
    return if ($aux->[1]);
    if (!defined($aux->[0])) {
	while (<$fh>) {
	    chomp;
	    if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
		$aux->[0] = $_;
		last;
	    }
	}
	if (!defined($aux->[0])) {
	    $aux->[1] = 1;
	    return;
	}
    }
    my ($name, $comm) = /^.(\S+)(?:\s+)(\S+)/ ? ($1, $2) :
                        /^.(\S+)/ ? ($1, '') : ('', '');
    my $seq = '';
    my $c;
    $aux->[0] = undef;
    while (<$fh>) {
	chomp;
	$c = substr($_, 0, 1);
	last if ($c eq '>' || $c eq '@' || $c eq '+');
	$seq .= $_;
    }
    $aux->[0] = $_;
    $aux->[1] = 1 if (!defined($aux->[0]));
    return ($name, $comm, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
	chomp;
	$qual .= $_;
	if (length($qual) >= length($seq)) {
	    $aux->[0] = undef;
	    return ($name, $comm, $seq, $qual);
	}
    }
    $aux->[1] = 1;
    return ($name, $seq);
}
