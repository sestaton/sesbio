#!/usr/bin/perl -w

use strict;
use File::Basename;
use Getopt::Long;

my $infile;
my $outfile;

GetOptions(# Required arguments
           'i|infile=s'         => \$infile,
           'o|outfile=s'        => \$outfile,
           );

# open the infile or die with a usage statement
die usage() if !$infile or !$outfile;
open(my $in, '<', $infile) || print "\nERROR: Could not open file: $infile\n";
open(my $out, '>', $outfile) || print "\nERROR: Could not open file: $outfile\n";

#
# comments must be removed or they will be counted
#
my @contigs = map +(split "\t")[1], <$in>;
close($in);

@contigs = sort { $a <=> $b } @contigs;
my @lengths;
foreach my $len (@contigs) {
    chomp $len;
    push(@lengths,$len);
}
my %seen = ();
my @unique_lengths = grep { ! $seen{$_} ++ } @lengths;   

my $unique = @unique_lengths;
my $total = @contigs;

print "\n","There are: ", $total, " total contigs.\n";
print "\n","There are: ", $unique, " unique contig lengths.\n\n";

count_unique ( @contigs );
close($out);

exit;

#
# sub
#
sub count_unique {

    my @array = @_;
    my %count;
    map { $count{$_}++ } @array;

    map {print $out $_."\t".${count{$_}}."\n"} sort keys(%count);

}

sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $0 -i inreport -o outreport                                                                                                                                                                                                                                   
This script takes as input a list of contigs with the length of each                                                                      
contig separated by tabs such as:                                                                                                                                                                                                                                                
contig-14445    445                                                                                                                      contig-8589     589                                                                                                                      contig-13776    776                                                                                                                      
(...)                                                                                                                                                                                                                                                                             
and the count of each length is output such as (minus the headers):                                                                       
                                                                                                                                         length  count                                                                                                                             
1009      1                                                                                                                              1038      1                                                                                                                              
1134      1                                                                                                                              1152      1                                                                                                                             
1236      1                                                                                                                              1269      1                                                                                                                             
130       581                                                                                                                            131       554                                                                                                                            
132       470                                                                                                                            133       480                                                                                                                            (...)                                                                                                                                    
                                                                                                                                         
NB: the lengths will not be reported in the same order as the input file.

END
}
