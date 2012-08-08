#!/usr/bin/env perl 

=head1 NAME 
                                                                       
getorf.pl - Search a muliti-fasta file and keep the longest ORFs 

=head1 SYNOPSIS    

getorf.pl -i seqs.fas -o seqs_trans.faa

=head1 DESCRIPTION
                                                                   
Translate a nucleotide multi-fasta file in all 6 frames and select 
the longest ORF for each sequence. The ORFs are reported as nucleotide
sequences by default, but translated may also be reported. The minimum 
ORF length to report can be given as an option.

=head1 DEPENDENCIES

This script uses EMBOSS, so it must be installed. 
EMBOSS v6.2+ must be installed (the latest is v6.5.7 as of this writing).

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The fasta files to be translated.

=item -o, --outfile

A file to place the translated sequences.

=back

=head1 OPTIONS

=over 2

=item -l, --orflen

The minimum length for which to report an ORF (Default: 80).
Lowering this value will not likely result in any significant hits 
from iprscan or other search programs (though there may be a reason to do so).

=item -t, --translate

Report translated ORFs instead of nucleotide sequences for each ORF.

=itme -s, --sameframe

Report all ORFs in the same (sense) frame.

=item -nm, --nomet

Do not report only those ORFs starting with Methionine (Default: Yes).

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut    

#
# Includes
#
use 5.010;
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use File::Basename;
use File::Temp;
use Bio::SeqIO;
use Pod::Usage;

#
# Vars
#
my $infile; 
my $outfile;
my $orflen;
my $sense;
my $find;
my $nomet;
my $help;
my $man;

#
# Counters
#
my $seqstot = 0;
my $orfseqstot = 0;
my $orfct = 0;
my $fcount = 0;

GetOptions(#Required
	   'i|infile=s'     => \$infile,
	   'o|outfile=s'    => \$outfile,
	   #Options
	   'l|orflen=i'     => \$orflen,
	   't|translate'    => \$find,
	   's|sameframe'    => \$sense,
	   'nm|nomet'       => \$nomet,
	   'h|help'         => \$help,
	   'm|man'          => \$man,
	  );

pod2usage( -verbose => 2 ) if $man;

#
# Check @ARGVs
#  
usage() and exit(0) if $help;

if (!$infile || !$outfile) {
    print "\nERROR: No input was given.\n";
    usage();
    exit(1);
}

if (defined $find) {
    $find = '1';
} 
else {
    $find = '3';
} 

my $getorf = find_prog("getorf");

if (-e $outfile) { 
# Because we are appending the ORFs from each sequence to the same output,
# there is the possibility to add to existing data, if the file exists. So,
# test to make sure it does not exist.
    die "\nERROR: $outfile already exists. Exiting.\n";
}
open(my $out, ">>", $outfile) or die "\nERROR: Could not open file: $outfile";

my ($fasnum, $seqhash) = seqct($infile);

$orflen //= '80';

if ($$fasnum >= 1) {
    print "\n========== Searching for ORFs with minimum length of $orflen.\n";
} else {
    die "\nERROR: No sequences were found! Check input. Exiting.\n";
}

my ($iname, $ipath, $isuffix) = fileparse($infile, qr/\.[^.]*/);

while (my ($id, $seq) = each %$seqhash) {
    $fcount++;
    my $orffile = getorf($iname,$isuffix,$fcount,$id,$seq,$find,$nomet);

    if (-s $orffile) {
	$orfseqstot++;
	my $longest_seq = largest_seq($orffile,$sense);

	while (my ($k, $v) = each %$longest_seq) {
	    if (defined $sense) {
		my ($sense_name, $sense_seq) = revcom($k,$v);
		print $out join("\n",(">".$sense_name, $sense_seq)), "\n";
	    }
	    else {
		print $out join("\n",(">".$k, $v)), "\n";
	    }
	}
    }
    unlink($orffile);
}
close($out);

print "\n========== $fcount sequences in $infile.\n";
print "\n========== $orfseqstot sequences processed with ORFs above $orflen.\n";

exit;

#
# Subs
#
sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!defined(@$aux));
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
    # SES mod 5/17/12
    my $name;
    if (/^.?(\S+\s\S+.*)/) {          # Illumina 1.8+, now more greedy 8/1 SES
        $name = $1;
    }
    elsif (/^.?(\S+)/) {              # Illumina 1.3+
        $name = $1;
    } 
    else {
        $name = '';                   # ?
    }
    #my $name = /^.(\S+)/? $1 : '';   # Heng Li's original regex
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
    return ($name, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
}

sub find_prog {
    my $prog = shift;
    my $path = qx(which $prog 2>&1 /dev/null);
    
    if ($path =~ /^which\: no getorf/) {
	say 'Couldn\'t find getorf in PATH. Will keep looking.';
	$path = "/usr/local/emboss/latest/bin/getorf";           # path at zcluster
    }

    # Instead of just testing if getorf exists and is executable 
    # we want to make sure we have permissions, so we try to 
    # invoke getorf and examine the output. 
    my $getorf_path = qx($path --help 2>&1 /dev/null);
    given ($getorf_path) {
	when (/Version\: EMBOSS/) { say 'Found it! /usr/local/emboss/latest/bin/getorf'; }
	when (/^-bash: \/usr\/local\/emboss\/bin\/getorf\: No such file or directory$/) { die "Could not find getorf. Exiting.\n"; }
	when ('') { die "Could not find getorf. Exiting.\n"; }
	default { die "Could not find getorf. Trying installing EMBOSS or adding it's location to your PATH. Exiting.\n"; }
    }	
    return($path);
}

sub seqct {
    my $f = shift;
    
    open(my $fh, "<", $f) or die "\nERROR: Could not open file: $f\n";
    my ($name, $seq, $qual);
    my @aux = undef;
    my $seqct = 0;
    my %seqhash;
    while (($name, $seq, $qual) = readfq(\*$fh, \@aux)) {
	$seqct++;
	# EMBOSS uses characters in identifiers as delimiters, which can produce some
        # unexpected renaming of sequences, so warn that it's not this script doing
        # the renaming.
	given ($name) {
	    when (/\:|\;|\||\(|\)|\.|\s/) { die "WARNING: Identifiers such as \"$name\" will produce unexpected renaming with EMBOSS."; }
	    when ('') { say 'WARNING: Sequences appear to have no identifiers. Continuing.'; }
	}
	$seqhash{$name} = $seq;
    }
    close($fh);
    return(\$seqct,\%seqhash);
}

sub getorf {
    my ($iname, $isuffix, $fcount, $id, $seq, $find, $nomet) = @_;
    my $tmpiname = $iname."_".$fcount."_XXXX";
    my $cwd = getcwd();
    my $fname = File::Temp->new( TEMPLATE => $tmpiname,
                                 DIR => $cwd,
                                 SUFFIX => $isuffix,
                                 UNLINK => 0);

    open(my $fh, ">", $fname) or die "\nERROR: Could not open file: $fname\n";

    print $fh join("\n",(">".$id, $seq)),"\n";

    close($fh);

    my $orffile = $fname."_orfs";

    my $getorfcmd;
    if (defined $nomet) {
	$getorfcmd = "$getorf -sequence $fname -outseq $orffile -minsize $orflen -find $find -nomethionine -auto";
    }
    else {
	$getorfcmd = "$getorf -sequence $fname -outseq $orffile -minsize $orflen -find $find -methionine -auto";
    }

    system($getorfcmd);
    unlink($fname);

    return($orffile);

}

sub largest_seq {
    my ($file, $sense) = shift;
 
    open(my $fh, "<", $file) or die "\nERROR: Could not open file: $file\n";
    
    my ($name, $seq, $qual);
    my @aux = undef;

    my %seqhash;
    while (($name, $seq, $qual) = readfq(\*$fh, \@aux)) {
	    $seqhash{$name} = $seq;
    }
    close($fh);

    # modified from:
    # http://stackoverflow.com/a/5958473
    my $max;
    my %hash_max;
    keys %seqhash; # reset iterator
    while(my ($key, $value) = each %seqhash) {
	if ( !defined $max || length($value) > $max ) {
	    %hash_max = ();
	    $max = length($value);
	}
	$hash_max{$key} = $value if $max == length($value);
    }
    
    return(\%hash_max);
}

sub revcom {
    my ($name, $seq) = @_;

    # If the sequence has been revcom'd 
    # we don't want the ID to say REVERSE. 
    $name =~ s/\(R.*//;   
    my $revcom = reverse $seq;
    $revcom =~ tr/ACGTacgt/TGCAtgca/;
    return($name, $revcom);
}

sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script -i infile -o outfile [-l] [-t] [-s] [-nm] [-h] [-m]

Required:
 -i|infile     :       A multifasta file. The longest ORF for each sequence will be reported.
 -o|outfile    :       A file to put the ORFs for each sequence.

Options:
 -l|orflen     :       An interger that will serve as the lower threshold
                       length for ORFs to consider prior to translating.
 -t|translate  :       If given, the longest ORF for each sequence will be translated
                       and the protein sequence will be reported.
 -s|sameframe  :       Report all ORFs in the same (sense) frame.
 -nm|nomet     :       Report all ORFs, not just those starting with Methionine (Default: Only report
                       ORFs starting with Methionine).
 -h|help       :       Print a usage statement.
 -m|man        :       Print the full documantion.

END
}
