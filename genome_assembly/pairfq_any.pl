#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
BEGIN {
  @AnyDBM_File::ISA = qw( DB_File SQLite_File )
      unless @AnyDBM_File::ISA == 1; # if loaded already, AnyDBM_File::ISA has a length of one;
}
use AnyDBM_File;
use vars qw( $DB_BTREE &R_DUP );
use lib qw(/iob_home/jmblab/statonse/apps/perlmod/AnyDBM_File-Importer-0.012/blib/lib);
use AnyDBM_File::Importer qw(:bdb);

my ($fread, $rread, $fpread, $rpread, $fsread, $rsread);

GetOptions(
           'f|forward=s'        => \$fread,
           'r|reverse=s'        => \$rread,
           'fp|forw_paired=s'   => \$fpread,
           'rp|rev_paired=s'    => \$rpread,
           'fs|forw_unpaired=s' => \$fsread,
           'rs|rev_unpaired=s'  => \$rsread,
          );

if (!$fread  || !$rread  ||
    !$fpread || !$rpread ||
    !$fsread || !$rsread) {
    print "\nERROR: Command line not parsed correctly. Check input.\n\n";
    usage();
    exit(1);
}

open(my $f, '<', $fread) or die "ERROR: Could not open file: $fread\n";
open(my $r, '<', $rread) or die "ERROR: Could not open file: $rread\n";
open(my $fp, '>', $fpread) or die "ERROR: Could not open file: $fpread\n";
open(my $rp, '>', $rpread) or die "ERROR: Could not open file: $rpread\n";
open(my $fs, '>', $fsread) or die "ERROR: Could not open file: $fsread\n";
open(my $rs, '>', $rsread) or die "ERROR: Could not open file: $rsread\n";

my %rseqhash;
$DB_BTREE->{cachesize} = 100000;
$DB_BTREE->{flags} = R_DUP;
tie( %rseqhash, 'AnyDBM_File', ':memory:', 0666, $DB_BTREE);

my @raux = undef;
my ($fname, $fseq, $fqual, $fid, $rname, $rseq, $rqual, $rid);
my ($fct, $rct, $fpct, $rpct, $pct, $fsct, $rsct, $sct) = (0, 0, 0, 0, 0, 0, 0, 0);

while (($rname, $rseq, $rqual) = readfq(\*$r, \@raux)) {
    $rct++;
    if ($rname =~ /(\/\d)$/) {
	$rname =~ s/$1//;
    } 
    elsif ($rname =~ /\s(\d\:\w\:\d\:\w+)$/) {
	$rid = $1; chomp($rid);
	$rid =~ s/^.//;
	$rname =~ s/\s.*//;
	$rname .= "|".$rid;
    } else {
	die "\nERROR: Could not determine Illumina encoding. Exiting.\n";
    }
    $rseqhash{$rname} = join("\t",$rseq, $rqual) if defined $rqual;
    $rseqhash{$rname} = $rseq if !defined $rqual;
}
close($r);

my ($forw_id, $rev_id);
my @faux = undef;
while (($fname, $fseq, $fqual) = readfq(\*$f, \@faux)) {
    $fct++;
    if ($fname =~ /(\/\d)$/) {
	$fname =~ s/\/\d//;
    }
    elsif ($fname =~ /\s(\d\:\w\:\d\:\w+)$/) {
        $fid = $1; chomp($fid);
        $fid =~ s/^.//;
	$fname =~ s/\s.*//;
	$fname .= "|".$fid;
    } else {
	die "\nERROR: Could not determine Illumina encoding. Exiting.\n";
    }
    if ($fname =~ /\|/) {
	my ($name, $id) = split(/\|/,$fname);
	$forw_id = $name." "."1".$id;
	$rev_id  = $name." "."2".$id;
    }
    if (exists $rseqhash{$fname}) {
	$fpct++; $rpct++;
	if (defined $fqual) {
	    my ($rread, $rqual) = split(/\t/,$rseqhash{$fname});
	    if ($fname =~ /\|/) {
		print $fp join("\n","@".$forw_id, $fseq, "+", $fqual),"\n";
		print $rp join("\n","@".$rev_id, $rread, "+", $rqual),"\n";
	    } else {
		print $fp join("\n","@".$fname."/1", $fseq, "+", $fqual),"\n";
                print $rp join("\n","@".$fname."/2", $rread, "+", $rqual),"\n";
	    }
	} else {
	    if ($fname =~ /\|/) {
		print $fp join("\n",">".$forw_id, $fseq),"\n";
		print $rp join("\n",">".$rev_id, $rseqhash{$fname}),"\n";
	    } else {
                print $fp join("\n",">".$fname."/1", $fseq),"\n";
                print $rp join("\n",">".$fname."/2", $rseqhash{$fname}),"\n";
            }
	}
        delete $rseqhash{$fname};
    } else {
	$fsct++;
	if (defined $fqual) {
	    if ($fname =~ /\|/) {
		print $fs join("\n","@".$forw_id, $fseq, "+", $fqual),"\n";
	    } else {
		print $fs join("\n","@".$fname."/1", $fseq, "+", $fqual),"\n";
	    }
	} else {
	    if ($fname =~ /\|/) {
		print $fs join("\n",">".$forw_id, $fseq),"\n";
	    } else {
		print $fs join("\n",">".$fname."/1", $fseq),"\n";
	    }
	}
    }
}
close($f);
close($fp);
close($rp);
close($fs);

my $rev_id_up;
while (my ($rname_up, $rseq_up) = each %rseqhash) {
    $rsct++;
    if ($rname_up =~ /\|/) {
	my ($uname, $uid) = split(/\|/,$rname_up);
	$rev_id_up .= $uname." "."2".$uid;
    }
    if ($rseq_up =~ /\t/) {
	my ($rread_up, $rqual_up) = split(/\t/,$rseq_up);
	if ($rname_up =~ /\|/) {
	    print $rs join("\n","@".$rev_id_up, $rread_up, "+", $rqual_up), "\n";
	} else {
	    print $rs join("\n","@".$rname_up."/2", $rread_up, "+", $rqual_up), "\n";
	}
    } else {
	if ($rname_up =~ /\|/) {
	    print $rs join("\n",">".$rev_id_up, $rseq_up), "\n";
	} else {
	    print $rs join("\n",">".$rname_up."/2", $rseq_up), "\n";
	}
    }
    $rev_id_up = undef;
}
close($rs);

$pct = $fpct + $rpct;
$sct = $fsct + $rsct;
untie %rseqhash;

print "\nTotal forward reads in $fread:              $fct";
print "\nTotal reverse reads in $rread:              $rct";
print "\nTotal paired reads in $fpread and $rpread:  $pct";
print "\nTotal forward paired reads in $fpread:      $fpct";
print "\nTotal reverse paired reads in $rpread:      $rpct";
print "\nTotal upaired reads in $fsread and $rsread: $sct"; 
print "\nTotal forward unpaired reads in $fsread:    $fsct";
print "\nTotal reverse unpaired reads in $rsread:    $rsct\n";

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
    ################################ SES mod 5/17/12
    my $name;
    if (/^.?(\S+\s\S+)/) {
	$name = $1;
    }
    elsif (/^.?(\S+)/) {
	$name = $1;
    } else {
	$name = '';
    }
    #my $name = /^.(\S+)/? $1 : ''; # original REGEX
    ################################
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

sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-f] [-r] [-fp] [-rp] [-fs] [-rs]

Arguments:
-f|forward        :       File of foward reads (usually with "/1" or " 1" in the header)
-r|reverse        :       File of reverse reads (usually with "/2" or " 2" in the header)
-fp|forw_paired   :       Name for the file of paired forward reads
-rp|rev_paired    :       Name for the file of paired reverse reads
-fs|forw_unpaired :       Name for the file of singleton forward reads
-rs|rev_unpaired  :       Name for the file of singleton reverse reads

EOF
}
