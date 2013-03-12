#!/usr/bin/env perl

## Load JSON file to get hits for each type DONE
## push matches into array, should stay ordered

use v5.12;
use strict;
use warnings;
use autodie qw(open);
use File::Spec qw(catfile rel2abs);
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long;
use JSON;
use Data::Dump qw(dd dump);

my $usage = "$0 -i cluster_fas_files_dir -d repeat_db -j repbase.json -o annot_clusters_dir -r report.tsv\n";
my $indir; 
my $database;
my $outdir;
my $json;
my $report;

GetOptions(
	   'i|indir=s'  => \$indir,
	   'o|outdir=s' => \$outdir,
	   'd|db=s'     => \$database,
	   'j|json=s'   => \$json,
	   'r|report=s' => \$report,
	   );

die $usage if !$indir or !$outdir or !$database or !$json or !$report;

## get input files
opendir(DIR, $indir) || die "\nERROR: Could not open directory: $indir. Exiting.\n";
my @clus_fas_files = grep /\.fa.*$/, readdir(DIR);
closedir(DIR);


if (scalar(@clus_fas_files) < 1) {
    say "\nERROR: Could not find any fasta files in $indir. Exiting.\n";
    exit(1);
}

## set path to output dir
my ($oname, $opath, $osuffix) = fileparse($outdir, qr/\.[^.]*/);
my $out_path = File::Spec->rel2abs($opath.$oname);

make_path($outdir, {verbose => 0, mode => 0711,}); # allows for recursively making paths

chdir($indir);

my @colors;

open(my $out, '>', $report);                                                                                                                             
say $out join "\t", "Cluster", "Read_count", "Type", "Class", "Superfam", "SINE_family","Top_hit";  ## this is not accurate and needs to be updated

for my $file (sort { $a =~ /CL(\d+)/ <=> $b =~ /CL(\d+)/ } @clus_fas_files) {
    my ($fname, $fpath, $fsuffix) = fileparse($file, qr/\.[^.]*/);
    my $blast_res = $fname;
    my ($filebase, $readct, @blname) = split /\_/, $fname;
    $blast_res =~ s/\.[^.]+$//;
    $blast_res .= "_blast.tsv";
    my $blast_file_path = File::Spec->catfile($out_path, $blast_res);

    my $blastcmd = "blastn -dust no -query $file -evalue 1e-5 -db $database -outfmt 6 | ".
	           "cut -f2 | ".                          # keep only the ssids
		   "sort | ".                             # sort the list
		   "uniq -c | ".                          # reduce the list
		   "sort -bnr | ".                        # count unique items
		   "perl -lane 'print join(\"\\t,\\@F)'"; # create an easy to parse format
    my @blast_out = qx($blastcmd);

    my ($hit_ct, $top_hit) = parse_blast(\@blast_out, $blast_file_path);
    next unless defined $top_hit;
    my $color = blast2annot($json, $filebase, $readct, $top_hit, $out);
    #$color = defined($color) ? $color : "nomatch";
    #push @colors, $color;
    ## return to hit by frequency

    if (not defined $hit_ct) {
	#say $$hit_ct," sequences have a BLAST hit ",$$top_hit," is the top hit ($$color) in ",$file;
	say "No hits for $file";
    }
}
close($out);

sub blast2annot {
    my ($json, $filebase, $readct, $top_hit, $out) = @_;

    my $repeats = json2hash($json);
    #my $repeats_path = File::Spec->rel2abs($repeats);
    my $color; # this is for annotating the barplot

    #open(my $out, '>', $report);
    #say $out join "\t", "Cluster", "Read_count", "Type", "Class", "Superfam", "SINE_family","Top_hit";

    for my $type (keys %$repeats) {
	if ($type eq 'pseudogene' || $type eq 'simple_repeat' || $type eq 'integrated_virus') {
	    if ($type eq 'pseudogene' && $$top_hit =~ /rrna|trna|snrna/i) {
		say $out join "\t", $filebase, $readct, $type, $$top_hit;
	    }
	    elsif ($type eq 'simple_repeat' && $$top_hit =~ /msat/i) {
		say $out join "\t", $filebase, $readct, $type, "Satellite", "MSAT", $$top_hit;
	    }
	    elsif ($type eq 'simple_repeat' && $$top_hit =~ /sat/i) {
		say $out join "\t", $filebase, $readct, $type, "Satellite", "SAT", $$top_hit;
	    }                                                                                                                                              
	    elsif ($type eq 'integrated_virus' && $$top_hit =~ /caul/i) {
		say $out join "\t", $filebase, $readct, $type, "Caulimoviridae", $$top_hit;
	    }
	    elsif ($type eq 'integrated_virus' && ($$top_hit eq 'PIVE' || $$top_hit eq 'DENSOV_HM')) {
		say $out join "\t", $filebase, $readct, $type, "DNA Virus", $$top_hit;
	    }
	    elsif ($type eq 'endogenous_retrovirus' && $$top_hit =~ /erv/i) {
		say $out join "\t", $filebase, $readct, $type, "Endogenous Retrovirus", $$top_hit;
	    }                                                                   
	    next;
	}
	for my $class (keys %{$repeats->{$type}}) {
	    while ( my ($superfam_index, $superfam) = each @{$repeats->{$type}{$class}} ) {
		for my $superfam_h (keys %$superfam) {
		    if ($superfam_h =~ /sine/i) {
			while (my ($sine_fam_index, $sine_fam_h) = each @{$superfam->{$superfam_h}}) {
			    for my $sine_fam_mem (keys %$sine_fam_h) {
				for my $sines (@{$repeats->{$type}{$class}[$superfam_index]{$superfam_h}[$sine_fam_index]{$sine_fam_mem}}) {
				    for my $sine (@$sines) {
					if ($sine =~ /$$top_hit/) {
					    say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $sine_fam_mem, $$top_hit;
					    #$color = "aquamarine"; ## fix this to be distinct
					}
					#else { say "NOMATCH $sine => $$top_hit"; }
				    }
				}
			    }
			}
		    } # if /sine/i
		    elsif ($superfam_h =~ /gypsy/i && $$top_hit =~ /^RLG/) {
			say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $$top_hit;
		    }
		    elsif ($superfam_h =~ /copia/i && $$top_hit =~ /^RLC/) {
			say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $$top_hit;
		    }
		    else {
			for my $fam (@{$repeats->{$type}{$class}[$superfam_index]{$superfam_h}}) {
			    for my $mem (@$fam) {
				if ($mem =~ /$$top_hit/i) {
				    #if ($mem =~ /^RLG/ || $mem =~ /gypsy/i) { $color = "darkmarine"; }
				    #elsif ($mem =~ /^RLC/ || $mem =~ /copia/i) { $color = "aquamarine4"; }
				    #elsif ($mem =~ /rrna/ || $mem =~ /sat/i) { $color = "azure2"; }
				    #elsif ($mem =~ /helitron/i) { $color = "darkolivegree3"; }
				    #else {
				    #    for my $classii (@classII) {
				    #	$color = "chartreuse" if $classii =~ /$$top_hit/i; 
				    #    }
				    #    for my $nonltr (@non_ltr_retros) { 
				    #	$color = "aquamarine2" if $nonltr =~ /$$top_hit/i;
				    #    }
				    #}
				    say $out join "\t", $filebase, $readct, $type, $class, $superfam_h, $$top_hit;
				}
				#else { say "NOMATCH $mem => $$top_hit"; }
			    }
			}
		    }
		}
	    }
	}
    }
    #close($out);
    return \$color;
}

sub json2hash {
    my $json = shift;
   
    my $json_text;
    local $/;
    open(my $in, '<', $json);
    $json_text = <$in>;
    close($in);

    my $repeats = JSON->new->utf8->space_after->decode($json_text);
    return $repeats;
}

sub parse_blast {
    my ($blast_out, $blast_file_path) = @_;
    my %blhits;

    my $top_hit;
    my $top_hit_num = 0;
    my $hit_ct = 0;

    for my $hit (@$blast_out) {
	chomp $hit;
	my ($ct, $hittype) = split /\t/, $hit;
	next unless defined $ct;
	$blhits{$hittype} = $ct;
	$hit_ct++;
    }

    #dump(%blhits);

    ## Need to create an array of colors to pass to R barplot   
    if ($hit_ct > 0) {
	open(my $out, '>', $blast_file_path);
	$top_hit = (reverse sort { $blhits{$a} <=> $blhits{$b} } keys %blhits)[0];
	for my $hits (reverse sort { $blhits{$a} <=> $blhits{$b} } keys %blhits) {
	    say $out join "\t", $hits, $blhits{$hits};
	}
	close($out);
	return \$hit_ct, \$top_hit;
    }
    else {
	return undef, undef;
    }
}

sub write_cls_size_dist_summary {
    my ($cls_file, $seqct, $cwd) = @_;

    ## open it, write it, then unlink it

    my $cls_size_dist_plot = $cls_file;
    $cls_size_dist_plot .= ".png";
    my $cls_size_dist_rscript = $cls_file;
    $cls_size_dist_rscript .= ".rscript";
    open(my $rscript, '>', $cls_size_dist_rscript);

    ## 
    ##

    say $rscript "cls=scan(file=\"$cls_file\",what=character(),sep=\"\\n\",comment.char=\">\",quiet=T)
                  cls=strsplit(cls,split=\"[ \\t]\")

png(filename=\"$cls_size_dist_plot\",width=800,height=500)

options(scipen=999) # this turns of printing in scientific notation
#n=${seqct}
#formatC(n, format = \"d\")
#format(n, scientific = FALSE)
#setwd(${cwd})     #################### if this is necessary set it to the cwd
# plot barplot:
NinClusters=length(unlist(cls))
#NinAll=length(seqs)
NinAll=$seqct
NinSingles=NinAll-NinClusters
clsLength=sapply(cls,length)
clsLength=sort(clsLength,decreasing=T)
clsLengthAll=c(clsLength,rep(1,NinSingles))

barplot(clsLengthAll,width=clsLengthAll,space=0,ylim=c(0,max(clsLength)*1.2),ylab=\"number of reads [reads]\",xlab=\"number of reads [%]\",main=paste(NinAll, \"reads total\"))
rect(0,0,NinClusters,clsLength[[1]]*1.2,col=\"#FF000010\")
rect(NinClusters,0,NinAll,clsLength[[1]]*1.2,col=\"#00FF0010\")
text(NinClusters/2,clsLength[[1]]*1.05, labels=paste(NinClusters,\"reads in\\n\",length(cls),\"clusters\"))
text(NinClusters+NinSingles/2,clsLength[[1]]*1.05, labels=paste(NinSingles,\"singlets\"))
axis(1,at=seq(0,NinAll,length.out=11),label=seq(0,100,by=10))

#legend(60,8,c(\"Ty3\/gypsy\",\"Ty1\/copia\",\"Non-LTR\",\"Helitron\",\"Class II\",\"rRNA\"),lwd=4,col=c(\"darkgreen\",\"aquamarine4\",\"aquamarine\",\"darkolivegreen3\",\"chartreuse\",\"azure2\"),bty=\"n\")

tmp=capture.output(dev.off())";
    close($rscript);

    system("/usr/local/R/2.15.0/lib64/R/bin/R --vanilla --slave --silent < $cls_size_dist_rscript 2> /dev/null");
}
