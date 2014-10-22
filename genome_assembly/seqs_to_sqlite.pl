#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use DBI;
use DBD::SQLite;

my $infile;
my $help;

GetOptions('i|infile=s' => \$infile, 'h|help' => \$help);

usage() and exit(0) if $help;

if (!$infile) {
    print "\nERROR: No input was given.\n";
    usage();
    exit(1);
}

my ($name, $comm, $seq, $qual, $header, $seqstr);
my @aux = undef;
my $dbfile = 'seqs.db';

my $fh = get_fh($infile);

my $dbh = DBI->connect( "DBI:SQLite:dbname=$dbfile", "", "",
			{ PrintError => 0 , RaiseError => 1 } );
$dbh->do("DROP TABLE IF EXISTS seqs");
$dbh->do("CREATE TABLE seqs(gene_name VARCHAR(50) PRIMARY KEY, sequence TEXT)");

while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
    $header = mk_key($name, $comm) if defined $comm && $comm ne '';
    $header = $name if !defined $comm || $comm eq '' ;
    $seqstr = mk_key($seq, $qual) if defined $qual;
    $seqstr = $seq if !defined $qual;
    $dbh->do("INSERT INTO genes VALUES('$header','$seqstr')")
}

my $sbh = DBI->connect( "DBI:SQLite:dbname=$dbfile", "", "" );

my $sth = $sbh->prepare("SELECT * FROM seqs");
$sth->execute();

while (my $row = $sth->fetchrow_hashref()) {
    if ($row->{gene_name} =~ /~~/ && $row->{sequence} =~ /~~/) {
	my ($id, $com) = mk_vec($row->{gene_name});
	my ($nt, $ql)  = mk_vec($row->{sequence});
	say join "\n", "@".$id.q{ }.$com, $nt, "+", $ql;
    }
    elsif ($row->{gene_name} =~ /~~/ && $row->{sequence} !~ /~~/) {
        my ($id, $com) = mk_vec($row->{gene_name});
	say join "\n", ">".$id.q{ }.$com, $row->{sequence};
    }
    elsif ($row->{gene_name} !~ /~~/ && $row->{sequence} =~ /~~/) {
	my ($nt, $ql)  = mk_vec($row->{sequence});
	say join "\n", "@".$row->{gene_name}, $nt, "+", $ql;
    }
    else {
	say join "\n", ">".$row->{gene_name}, $row->{sequence};
    }
}

$sth->finish();
$dbh->disconnect();
unlink $dbfile;

exit;
#
# methods
#
sub get_fh {
    my ($file) = @_;

    my $fh;
    if ($file =~ /\.gz$/) {
        open $fh, '-|', 'zcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /\.bz2$/) {
        open $fh, '-|', 'bzcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    else {
        open $fh, '<', $file or die "\nERROR: Could not open file: $file\n";
    }

    return $fh;
}

sub mk_key { return join "~~", @_ }

sub mk_vec { return split /\~\~/, shift }

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
    my ($name, $comm);
    defined $_ && do {
        ($name, $comm) = /^.(\S+)(?:\s+)(\S+)/ ? ($1, $2) : 
	                 /^.(\S+)/ ? ($1, '') : ('', '');
    };
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

sub usage {
    my $script = basename($0);
  print STDERR <<END
USAGE: $script -i s_1_sequence.fasta

Required:
    -i|infile   :    FastA/Q file of reads/contigs.

Options:
    -h|help     :    Print usage statement.

END
}
    
