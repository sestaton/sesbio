#!/usr/bin/env perl

##NB: This script requires Graph::GEXF (https://github.com/sestaton/graph-gexf)
##TODO: add example file formats and role into tests for this library

use 5.012;
use strict;
use warnings;
use autodie;
use Data::Dump;
use Graph::GEXF;

my $usage = "$0 hitsort_file cls_file\n";
my $hitsort   = shift or die $usage;
my $clsfile   = shift or die $usage;
my $clusterid = shift;

$clusterid //= 'CL100'; 

my %nodes;
my %idmap;
my %clsstats;

open my $hs, '<', $hitsort;
while (<$hs>) {
    chomp;
    my @f = split /\t/;
    $idmap{$f[0]} = $f[1];
    $idmap{$f[1]} = $f[0];
    push @{$nodes{$f[0]}}, join "|", $f[1], $f[2];
}
close $hs;

my $edge = 0;
my $node = 0;

{
    local $/ = '>';
    
    open my $in, '<', $clsfile;   
    while (my $line = <$in>) {
	$line =~ s/>//g;
	next if !length($line);
	my ($clsid, $seqids) = split /\n/, $line;
	my ($id, $seqnum) = split /\s+/, $clsid;
	my @ids = split /\s+/, $seqids;
	my @edges;
	my $graph = Graph::GEXF->new( visualization => 0, graph_mode => 'static' );
	if ($id eq $clusterid) { # this is just to demonstrate how to make a graph for one cluster
	    for my $i (@ids) {
		my $i_j = $idmap{$i};
		if (exists $nodes{$i}) {
		    my $nodeobj= $graph->add_node($node);

		    for my $n (@{$nodes{$i}}) {
			my ($n_i, $w_i) = split /\|/, $n;

			my $edge = Graph::GEXF::Edge->new(
							  source => $i,
							  target => $n_i,
							  weight => $w_i
							  );
			$nodeobj->add_edge( $node => $edge );
			$node++;
			$edge++;
		    }
		}
		elsif (exists $nodes{$i_j}) {
		    my $nodeobj= $graph->add_node($node);

		    for my $n (@{$nodes{$i_j}}) {
			my ($n_i_j, $w_i_j) = split /\|/, $n;

			my $edge = Graph::GEXF::Edge->new(
							  source => $i,
							  target => $n_i_j,
							  weight => $w_i_j
							  );
			$nodeobj->add_edge( $node => $edge );
			$node++;
			$edge++;
		    }
		}
	    }
	    my $xml = $graph->to_xml;
	    my $graph_file = $id."_graph.gexf";
	    open my $gr, '>', $graph_file or die "\n[ERROR]: Could not open file: $graph_file\n";
	    say $gr $xml;
	    close $gr;
	    $clsstats{scalar(@ids)} = join "|", $id, $edge;
	}
	#my $xml = $graph->to_xml;
	#my $graph_file = $id."_graph.gexf";
	#open my $gr, '>', $graph_file or die "\n[ERROR]: Could not open file: $graph_file\n";
	#say $gr $xml;
	#close $gr;
	#$clsstats{scalar(@ids)} = join "|", $id, $edge;
	$edge = 0;
    }
    close $in;
}

say join "\t","Cluster","Nodes","Edges";
for my $n (reverse sort { $a <=> $b } keys %clsstats) {
    my ($k, $e) = split /\|/, $clsstats{$n};
    say join "\t", $k, $n, $e;
}
