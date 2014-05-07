#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::DB::EUtilities;
 
my $id = 527031;
 
my $factory = Bio::DB::EUtilities->new(-eutil => 'esummary',
                                       -email => 'mymail@foo.bar',
                                       -db    => 'taxonomy',
                                       -id    => $id );
 
my ($name) = $factory->next_DocSum->get_contents_by_name('ScientificName');
 
say $name;


