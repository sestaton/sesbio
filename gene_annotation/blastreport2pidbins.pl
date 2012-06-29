#!/usr/bin/perl -w


use strict;
use Getopt::Long;

my $usage = "USAGE: $0 -i blasttablein -o tophitsout <--graph> <--full_report> <--annot_by_bin> <--idlist>\n\n";
my $infile;
my $outfile; 
my $graph;
my $ids;
my $full_report;
my $annotations;

GetOptions(
	   'i|infile=s'     => \$infile,
	   'o|outfile=s'    => \$outfile,
	   'g|graph'        => \$graph,
	   'f|full_report'  => \$full_report,
	   'a|annot_by_bin' => \$annotations,
	   'id|idlist'      => \$ids,
	   );

if (!$infile || !$outfile) {
    die "\nERROR: No infile/outfile was found. Exiting.\n\n",$usage;
}

open(my $in , '<', $infile) or die "\nERROR: Could not open file: $infile\n";
open(my $out , '>', $outfile) or die "\nERROR: Could not open file: $outfile\n";

my $totalhits = 0;
my $onehun2ninety = 0;
my $eightynine2eighty = 0;
my $seventynine2seventy = 0;
my $sixtynine2sixty = 0;
my $fiftynine2fifty = 0;
my $lessthanfifty = 0;

my $binned_annot = $outfile;
$binned_annot =~ s/\.[^\.]*$//;
my $onehun2ninety_annot = $binned_annot."_90-100pid_annot.txt";
my $eightynine2eighty_annot = $binned_annot."_80-89pid_annot.txt";
my $seventynine2seventy_annot = $binned_annot."_70-79pid_annot.txt";
my $sixtynine2sixty_annot = $binned_annot."_60-69pid_annot.txt";
my $fiftynine2fifty_annot = $binned_annot."_50-59pid_annot.txt";
my $lessthanfifty_annot = $binned_annot."_lt50pid_annot.txt";

my $onehun2ninety_idlist = $binned_annot."_90-100pid_idlist.txt";
my $eightynine2eighty_idlist = $binned_annot."_80-89pid_idlist.txt";
my $seventynine2seventy_idlist = $binned_annot."_70-79pid_idlist.txt";
my $sixtynine2sixty_idlist = $binned_annot."_60-69pid_idlist.txt";
my $fiftynine2fifty_idlist = $binned_annot."_50-59pid_idlist.txt";
my $lessthanfifty_idlist = $binned_annot."_lt50pid_idlist.txt";

my ($onehun2ninety_annot_hndl, $eightynine2eighty_annot_hndl, $seventynine2seventy_annot_hndl, 
    $sixtynine2sixty_annot_hndl, $fiftynine2fifty_annot_hndl, $lessthanfifty_annot_hndl);
my ($onehun2ninety_idlist_hndl, $eightynine2eighty_idlist_hndl, $seventynine2seventy_idlist_hndl,
    $sixtynine2sixty_idlist_hndl, $fiftynine2fifty_idlist_hndl, $lessthanfifty_idlist_hndl);

# add test for --annot_by_bin option; otherwise, do not create all these files
if ($annotations) {
    open($onehun2ninety_annot_hndl, '>', $onehun2ninety_annot) or die "\nERROR: Could not open file: $onehun2ninety_annot\n";
    open($eightynine2eighty_annot_hndl, '>', $eightynine2eighty_annot) or die "\nERROR: Could not open file: $eightynine2eighty_annot\n"; 
    open($seventynine2seventy_annot_hndl, '>', $seventynine2seventy_annot) or die "\nERROR: Could not open file: $seventynine2seventy_annot\n";
    open($sixtynine2sixty_annot_hndl, '>', $sixtynine2sixty_annot) or die "\nERROR: Could not open file: $sixtynine2sixty_annot\n";
    open($fiftynine2fifty_annot_hndl, '>', $fiftynine2fifty_annot)  or die "\nERROR: Could not open file: $fiftynine2fifty_annot\n";
    open($lessthanfifty_annot_hndl, '>', $lessthanfifty_annot) or die "\nERROR: Could not open file: $lessthanfifty_annot\n";
}
if ($ids) {
    open($onehun2ninety_idlist_hndl, '>', $onehun2ninety_idlist) or die "\nERROR: Could not open file: $onehun2ninety_idlist\n";
    open($eightynine2eighty_idlist_hndl, '>', $eightynine2eighty_idlist) or die "\nERROR: Could not open file: $eightynine2eighty_idlist\n";
    open($seventynine2seventy_idlist_hndl, '>', $seventynine2seventy_idlist) or die "\nERROR: Could not open file: $seventynine2seventy_idlist\n";
    open($sixtynine2sixty_idlist_hndl, '>', $sixtynine2sixty_idlist) or die "\nERROR: Could not open file: $sixtynine2sixty_idlist\n";
    open($fiftynine2fifty_idlist_hndl, '>', $fiftynine2fifty_idlist)  or die "\nERROR: Could not open file: $fiftynine2fifty_idlist\n";
    open($lessthanfifty_idlist_hndl, '>', $lessthanfifty_idlist) or die "\nERROR: Could not open file: $lessthanfifty_idlist\n";
}

while (<$in>) { 
    chomp; 
    next if /^#/ || /^Query/;
    $totalhits++;
    my @blfields = split(/\t/, $_);
    if ($blfields[2] >= 90.00 && $blfields[2] <= 100.00) {
	$onehun2ninety++;
	print $onehun2ninety_annot_hndl join("\t",(@blfields),"\n") if $annotations;
	print $onehun2ninety_idlist_hndl $blfields[0]."\n" if $ids;
    }
    elsif ($blfields[2] >= 80.00 && $blfields[2] <= 89.99) {
	$eightynine2eighty++;
	print $eightynine2eighty_annot_hndl join("\t",(@blfields),"\n") if $annotations;
	print $eightynine2eighty_idlist_hndl $blfields[0]."\n" if $ids;
    }
    elsif ($blfields[2] >= 70.00 && $blfields[2] <= 79.99) {
	$seventynine2seventy++;
	print $seventynine2seventy_annot_hndl join("\t",(@blfields),"\n") if $annotations;
	print $seventynine2seventy_idlist_hndl $blfields[0]."\n" if $ids;
    }
    elsif ($blfields[2] >= 60.00 && $blfields[2] <= 69.99) {
	$sixtynine2sixty++;
	print $sixtynine2sixty_annot_hndl join("\t",(@blfields),"\n") if $annotations;
	print $sixtynine2sixty_idlist_hndl $blfields[0]."\n" if $ids;
    }
    elsif ($blfields[2] >= 50.00 && $blfields[2] <= 59.99) {
	$fiftynine2fifty++;
	print $fiftynine2fifty_annot_hndl join("\t",(@blfields),"\n") if $annotations;
	print $fiftynine2fifty_idlist_hndl $blfields[0]."\n" if $ids;
    }
    elsif ($blfields[2] < 49.99) {
	$lessthanfifty++;
	print $lessthanfifty_annot_hndl join("\t",(@blfields),"\n") if $annotations;
	print $lessthanfifty_idlist_hndl $blfields[0]."\n" if $ids;
    } 
}

if ($annotations) {
    close($onehun2ninety_annot_hndl);
    close($eightynine2eighty_annot_hndl);
    close($seventynine2seventy_annot_hndl);
    close($sixtynine2sixty_annot_hndl);
    close($fiftynine2fifty_annot_hndl);
    close($lessthanfifty_annot_hndl);
}
if ($ids) {
    close($onehun2ninety_idlist_hndl);
    close($eightynine2eighty_idlist_hndl);
    close($seventynine2seventy_idlist_hndl);
    close($sixtynine2sixty_idlist_hndl);
    close($fiftynine2fifty_idlist_hndl);
    close($lessthanfifty_idlist_hndl);
}

print $out "#total_hits\t100-90\t89-80\t79-70\t69-60\t59-50\t<50\n";
if ($graph) {

    my $rscript = $outfile."_rscript";
    my $graph_file = $rscript.".png";
    open(my $tmprscript, '>', $rscript) or die "\nERROR: Could not open file: $rscript\n";
    print $tmprscript "arr <- c($totalhits,$onehun2ninety,$eightynine2eighty,$seventynine2seventy,$sixtynine2sixty,$fiftynine2fifty,$lessthanfifty);
for (i in 1:6) {arr.i <- arr[i+1]/arr[1]; if (i == 1) { arrbins <- arr.i } else { arrbins <- append(arrbins,arr.i)} };
png(filename=\"$graph_file\",width=1024,height=640);
barplot(arrbins,names.arg=c(\"100-90\",\"89-80\",\"79-70\",\"69-60\",\"59-50\",\"<50\"),ylab=\"Percentage of ORF alignments\",xlab=\"Amino Acid Identity\",main=\"ORF conservation\")
tmp<-dev.off();
quit();
";
    close($tmprscript);
    system("R --vanilla --slave --silent < $rscript 2>/dev/null");
    unlink($rscript);

} 

if ($full_report) {

    print $out $totalhits."\t".
               $onehun2ninety."\t".
               $eightynine2eighty."\t".
               $seventynine2seventy."\t".
               $sixtynine2sixty."\t".
               $fiftynine2fifty."\t".
               $lessthanfifty."\n";

    my $onehundred2ninety_percent = sprintf("%.2f",$onehun2ninety/$totalhits);
    my $eightynine2eighty_percent = sprintf("%.2f",$eightynine2eighty/$totalhits);
    my $seventynine2seventy_percent = sprintf("%.2f",$seventynine2seventy/$totalhits);
    my $sixtynine2sixty_percent = sprintf("%.2f",$sixtynine2sixty/$totalhits);
    my $fiftynine2fifty_percent = sprintf("%.2f",$fiftynine2fifty/$totalhits);
    my $lessthanfifty_percent = sprintf("%.2f",$lessthanfifty/$totalhits);

    print $out "#100-90_percent\t89-80_percent\t79-70_percent\t69-60_percent\t59-50_percent\t<50_percent\n";
    print $out $onehundred2ninety_percent."\t".
               $eightynine2eighty_percent."\t".
               $seventynine2seventy_percent."\t".
               $sixtynine2sixty_percent."\t".
               $fiftynine2fifty_percent."\t".
               $lessthanfifty_percent."\n";

} else {

    print $out $totalhits."\t".
	       $onehun2ninety."\t".
               $eightynine2eighty."\t".
               $seventynine2seventy."\t".
               $sixtynine2sixty."\t".
               $fiftynine2fifty."\t".
               $lessthanfifty."\n";

}

close($in);
close($out);

exit;
