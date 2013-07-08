#! /usr/bin/perl -w
# chienchi@lanl.gov
# generate some stats for Newbler or Velvet (>0.7.6) contigs.
# Assume contigs fasta file is in the assembly output folder
#   which contains other files, like Log (velvet), stats.txt(velvet),
#   and 454ContigGraph.txt(Newbler)
# 20100423
# no flag -t can generate some statistics with contig file only.
# 20100708

use strict;
use File::Basename;
use Getopt::Long;
my $plot;
my $assembly_tool="";
my $pdf;

GetOptions( "plot"     => \$plot,
            "tool=s"   => \$assembly_tool,
            "One_pdf"  => \$pdf,
            "help|?"   => sub {Usage()});
sub Usage
{
    print STDERR "perl $0 [-p -O] -t <Newbler | Velvet> <contigs.fasta>\n";
    print STDERR "     -tool       The assembly tool: Newbler or Velvet\n";
    print STDERR "     -plot       plot length histogram, GC histogram, GC vs. Depth, Len vs. Depth, Len vs. Cov
                 (need R installed)\n";
    print STDERR "                 -One_pdf    put all plots generate by -plot option in one pdf file.\n";
    print STDERR "     -help       print usage\n";
    exit;
}
if (scalar(@ARGV) != 1){&Usage;}
#if ($assembly_tool !~ /Newbler/i && $assembly_tool !~ /Velvet/i){&Usage;}

my ($file_name, $path, $suffix)=fileparse("$ARGV[0]", qr/\.[^.]*/);
my ($len,$total)=(0,0);
my @x;
my $seq_num;
my $seq;
my $reads_num;
my $kmer_cov;
my $GC_num;
my $GC_content;
my $velvet_id;
my $Newbler_id;
my $id_to_reads;
my $id_to_cov;
my $used_percent;
my $singleton;
my $exp_cov;
my ($over100k_bases,$over50k_bases,$over25k_bases,$over10k_bases,$over5k_bases,$over3k_bases,$over2k_bases,$over1k_bases)=(0,0,0,0,0,0,0,0,0);
my ($over100k_reads,$over50k_reads,$over25k_reads,$over10k_reads,$over5k_reads,$over3k_reads,$over2k_reads,$over1k_reads)=(0,0,0,0,0,0,0,0,0);

($id_to_reads,$used_percent,$singleton,$exp_cov)=&read_velvet_stats_and_log if ($assembly_tool =~ /Velvet/i);
if ($assembly_tool =~ /Newbler/i)
{
    my $assembled_reads=`grep -c "Assembled" $path/454ReadStatus.txt`;
    my $total_reads=`wc -l $path/454ReadStatus.txt | awk '{print \$1 - 1}'`;
    $singleton=`grep -c "Singleton" $path/454ReadStatus.txt`;
    chomp $singleton;
    $used_percent= $assembled_reads/$total_reads*100;
    if ($plot) {$id_to_cov=&read_Newbler_ContigGraph};
}
if ($plot){
    open (PLOT,">${file_name}_len_gc_depth_cov.txt");
    print PLOT "Len\tGC\tDepth\tCov\n";
}
while(<>){
        chomp;
	if(/^[\>\@]/){
	    $seq_num++;
		if($len>0){
                        &stats($len);
			push @x,$len;
			$GC_num = $seq =~ tr/GCgc/GCgc/;
                        $GC_content = $GC_num/$len;

                        printf PLOT ("%d\t%.4f\t%d\t%.2f\n",$len,$GC_content,$reads_num,$kmer_cov) if ($plot && $assembly_tool =~ /Velvet/i);
                        printf PLOT ("%d\t%.4f\t%d\t%.2f\n",$len,$GC_content,$reads_num,$id_to_cov->{$Newbler_id}) if ($plot && $assembly_tool =~ /Newbler/i);
		}
                if ($assembly_tool =~ /Newbler/i){
                        ($reads_num) = $_ =~ /numreads=(\d+)/;
                        ($Newbler_id) = $_ =~ /^>(\S+)/;
                }
                if ($assembly_tool =~ /Velvet/i){
                        ($velvet_id,$kmer_cov) = $_ =~ /NODE_(\d+)_length_\d+_cov_(.*)/;
                        $reads_num = $id_to_reads->{$velvet_id};
                }
		$len=0;
                $seq="";
	}
	else{
		s/\s//g;
		$len+=length($_);
                $seq.=$_;
	}
}
if ($len>0){
        &stats($len);
	push @x,$len;
	$GC_num = $seq =~ tr/GCgc/GCgc/;
        $GC_content = $GC_num/$len;
        printf PLOT ("%d\t%.4f\t%d\t%.2f\n",$len,$GC_content,$reads_num,$kmer_cov) if ($plot && $assembly_tool =~ /Velvet/i);
        printf PLOT ("%d\t%.4f\t%d\t%.2f\n",$len,$GC_content,$reads_num,$id_to_cov->{$Newbler_id}) if ($plot && $assembly_tool =~ /Newbler/i);
}
close PLOT if ($plot);
@x=sort{$b<=>$a} @x;
my $N50;
my $N90;
my ($count,$half)=(0,0);
for (my $j=0;$j<@x;$j++){
	$count+=$x[$j];
	if (($count>=$total/2)&&($half==0)){
                $N50=$x[$j];
		$half=$x[$j]
	}elsif ($count>=$total*0.9){
                $N90=$x[$j];
                last;
	}
}

my ($top10, $top20, $top40, $top100);
for (0..99)
{
   $top10+= $x[$_] if ($_<9 and $x[$_]);
   $top20+= $x[$_] if ($_<19 and $x[$_]);
   $top40+= $x[$_] if ($_<39 and $x[$_]);
   $top100+= $x[$_] if ($_<99 and $x[$_]);
}
print "Expected_coverage:\t$exp_cov\n" if ($assembly_tool =~ /Velvet/i and $exp_cov);
if ($assembly_tool =~ /velvet|newbler/i){
  printf ("Assembled_reads:\t%.2f%%\n",$used_percent);
  print "Singleton:\t$singleton\n";
}
print "Contigs_number:\t$seq_num\n";
print "N50:\t$N50\n";
print "N90:\t$N90\n";
print "Max:\t$x[0]\n";
print "Min:\t$x[-1]\n";
print "Total_bases:\t$total\n";
print "Top10_bases:\t$top10\n";
print "Top20_bases:\t$top20\n";
print "Top40_bases:\t$top40\n";
print "Top100_bases:\t$top100\n";
print ">100kb_bases:\t$over100k_bases\n";
print ">50kb_bases:\t$over50k_bases\n";
print ">25kb_bases:\t$over25k_bases\n";
print ">10kb_bases:\t$over10k_bases\n";
print ">5kb_bases:\t$over5k_bases\n";
print ">3kb_bases:\t$over3k_bases\n";
print ">2kb_bases:\t$over2k_bases\n";
print ">1kb_bases:\t$over1k_bases\n";
if ($assembly_tool =~ /velvet|newbler/i){
  print ">100kb_reads:\t$over100k_reads\n";
  print ">50kb_reads:\t$over50k_reads\n";
  print ">25kb_reads:\t$over25k_reads\n";
  print ">10kb_reads:\t$over10k_reads\n";
  print ">5kb_reads:\t$over5k_reads\n";
  print ">3kb_reads:\t$over3k_reads\n";
  print ">2kb_reads:\t$over2k_reads\n";
  print ">1kb_reads:\t$over1k_reads\n";
}

if ($plot && !$pdf){

  open (Rscript, ">Rscript$$");

  print Rscript "png(filename=\"${file_name}_len_hist.png\",width=1024,height=640);
a<-read.table(file=\"${file_name}_len_gc_depth_cov.txt\",header=TRUE);
h<-hist(a\$Len,col=\'blue\',breaks=max(a\$Len), main=\"Contig Length Distribution\",xlab=\"Contig Length\");
leg.txt<-c(\"$seq_num contigs\",\"N50 $N50\",\"N90 $N90\",\"Max $x[0]\",\"Min $x[-1]\",\"Total Base $total\");
legend(\"topright\",leg.txt);
tmp<-dev.off();

png(filename=\"${file_name}_gc_hist.png\",width=1024,height=640);
GC<-a\$GC*100;
hist(GC,breaks=c(0:100),xlab=\"GC (%)\",ylab=\"# of reads\",main=\"GC Histogram for contigs\");
GCmean<-format(mean(GC),digit=4,nsmall = 2)
GCsd<-format(sd(GC),digit=3,nsmall = 2)
legend(\"topleft\",c(paste(\"avg=\",GCmean,\"%\"),paste(\"std=\",GCsd,\"%\")),inset=0.02)
tmp<-dev.off();

png(filename=\"${file_name}_gc_vs_depth.png\",width=1024,height=640);
plot(GC,a\$Depth,xlab=\"GC (%)\",ylab=\"Coverage (# of reads in contig)\", main=\"contigs Coverages vs. GC\",pch=4);
tmp<-dev.off();

png(filename=\"${file_name}len_vs_depth.png\",width=1024,height=640);
plot(a\$Len,a\$Depth,xlab=\"Contig Length\",ylab=\"Coverage (# of reads in contig)\", main=\"contigs Coverages vs. Contig Length\",pch=4)
tmp<-dev.off();

png(filename=\"${file_name}_len_vs_cov.png\",width=1024,height=640);
plot(a\$Len,a\$Cov,xlab=\"Contig Length\",ylab=\"Mean Coverage\", main=\"Contigs Mean Coverages vs. Contig Length\",pch=4)
tmp<-dev.off();
quit();
";

  system ("R --vanilla --slave --silent < Rscript$$ 2>/dev/null");
  unlink "Rscript$$";
}

if ($plot && $pdf)
{
open (Rscript, ">Rscript$$");
print Rscript "\\documentclass[landscape]{report}
\\usepackage{fullpage}
%\\usepackage{Sweave}
\\RequirePackage[T1]{fontenc}
\\RequirePackage{graphicx,ae,fancyvrb}
\\IfFileExists{upquote.sty}{\\RequirePackage{upquote}}{}
\\setkeys{Gin}{width=7.5in,height=6in}
\\DefineVerbatimEnvironment{Sinput}{Verbatim}{fontshape=sl}
\\DefineVerbatimEnvironment{Soutput}{Verbatim}{}
\\DefineVerbatimEnvironment{Scode}{Verbatim}{fontshape=sl}
\\newenvironment{Schunk}{}{}
\\pagestyle{empty}
\\begin{document}
<<echo=FALSE>>=
a<-read.table(file=\"${file_name}_len_gc_depth_cov.txt\",header=TRUE);
@
\\begin{center}
<<fig=TRUE,echo=FALSE,>>=
h<-hist(a\$Len,col=\'blue\',breaks=max(a\$Len), main=\"Contig Length Distribution\",xlab=\"Contig Length\");
leg.txt<-c(\"$seq_num contigs\",\"N50 $N50\",\"N90 $N90\",\"Max $x[0]\",\"Min $x[-1]\",\"Total Base $total\");
legend(\"topright\",leg.txt);
@
\\end{center}

\\begin{center}
<<fig=TRUE,echo=FALSE,>>=
GC<-a\$GC*100;
hist(GC,breaks=c(0:100),xlab=\"GC (%)\",ylab=\"# of reads\",main=\"GC Histogram for contigs\");
GCmean<-format(mean(GC),digit=4,nsmall = 2)
GCsd<-format(sd(GC),digit=3,nsmall = 2)
legend(\"topleft\",c(paste(\"avg=\",GCmean,\"%\"),paste(\"std=\",GCsd,\"%\")),inset=0.02)
@
\\end{center}

\\begin{center}
<<fig=TRUE,echo=FALSE,>>=
plot(GC,a\$Depth,xlab=\"GC (%)\",ylab=\"Coverage (# of reads in contig)\", main=\"contigs Coverages vs. GC\",pch=4);
@
\\end{center}

\\begin{center}
<<fig=TRUE,echo=FALSE,>>=
plot(a\$Len,a\$Depth,xlab=\"Contig Length\",ylab=\"Coverage (# of reads in contig)\", main=\"contigs Coverages vs. Contig Length\",pch=4)
@
\\end{center}

\\begin{center}
<<fig=TRUE,echo=FALSE,>>=
plot(a\$Len,a\$Cov,xlab=\"Contig Length\",ylab=\"Mean Coverage\", main=\"Contigs Mean Coverages vs. Contig Length\",pch=4)
@
\\end{center}

\\end{document}
";
  system ("R CMD Sweave Rscript$$ 1>/dev/null");
  system ("gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=${file_name}_stats_plots.pdf -dBATCH Rscript$$*.pdf 1> /dev/null");
  #system ("mv Rplots.pdf ${file_name}_stats_plots.pdf");
  system ("rm Rscript*");
}

sub stats{
        my $len = shift;
	$total+=$len;
        $over100k_bases+=$len if ($len>100000);
        $over50k_bases+=$len if ($len>50000);
        $over25k_bases+=$len if ($len>25000);
        $over10k_bases+=$len if ($len>10000);
        $over5k_bases+=$len if ($len>5000);
        $over3k_bases+=$len if ($len>3000);
        $over2k_bases+=$len if ($len>2000);
        $over1k_bases+=$len if ($len>1000);

        if ($assembly_tool =~ /velvet|newbler/i){
         $over100k_reads+=$reads_num if ($len>100000);
         $over50k_reads+=$reads_num if ($len>50000);
         $over25k_reads+=$reads_num if ($len>25000);
         $over10k_reads+=$reads_num if ($len>10000);
         $over5k_reads+=$reads_num if ($len>5000);
         $over3k_reads+=$reads_num if ($len>3000);
         $over2k_reads+=$reads_num if ($len>2000);
         $over1k_reads+=$reads_num if ($len>1000);
        }
}

sub read_velvet_stats_and_log
{
    my @array;
    my %id_to_reads_num;
    # read stats.txt pile which contains numReads info each contigs.
    open (IN,"$path/stats.txt") || die "cannot open velvet stats.txt file";
    while (<IN>)
    {
       chomp;
       if ($_=~/^\d/)
       {
           @array=split /\t/, $_;
           $id_to_reads_num{$array[0]}=$array[10];
       }
    }
    close IN;

    # raed  Laafile which contains
    open (IN,"$path/Log") || die "cannot open velvet Log file";
    my ($used_percent,$singleton,$exp_cov);
    while(<IN>)
    {
       chomp;
        #Median coverage depth = 4.343750
        if ($_=~/Median coverage depth/)
        {
           ($exp_cov)= $_=~/(\d+\.\d+)/;
        }
        if ($_=~/^Final/){
           my ($used,$total)=$_=~/using\s+(\d+)\/(\d+)/;
           $used_percent = $used/$total*100;
           $singleton = $total - $used;
        }
    }
    close IN;
    return (\%id_to_reads_num,$used_percent,$singleton,$exp_cov);
}
sub read_Newbler_ContigGraph
{
    my @array;
    my %id_to_cov;
    open (IN,"$path/454ContigGraph.txt") ||die "cannot open 454ContigGraph.txt file";
    while (<IN>)
    {
       chomp;
       @array=split /\t/,$_;
       $id_to_cov{$array[1]}=$array[3];
    }
    return (\%id_to_cov);
}
