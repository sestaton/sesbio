#!/usr/bin/env perl

use 5.014;
use strict;
use warnings;
use autodie qw(open);
use Data::Dump qw(dd);
use JSON;
# JSON::Parse <- has json_to_perl method

my $usage = "$0 repbase1801.json\n";
my $infile = shift or die $usage;

my $json_text;
{
    local $/;
    open(my $in, '<', $infile);
    $json_text = <$in>;
    close($in);
}

#my $json = JSON->new->utf8->space_after->decode($infile);
my $json = JSON->new->utf8->decode($json_text);
dd $json;
