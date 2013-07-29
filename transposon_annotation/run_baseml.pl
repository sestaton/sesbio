#!/usr/bin/perl -w 
#_______________________________________________________________________+
#                                                                       |
# run_baseml.pl               
#_______________________________________________________________________+
#                                                                       |
# Description: 
#                                                                       |
# Author: S. Evan Staton                                                |
# Contact: statonse<at>uga.edu                                          |
# Started: 1.6.11                                                      |                                                
# Updated:                                                      
#                                                                       |
# Suggested usage:                                                      |

#_______________________________________________________________________+
# TODO: 

use strict;
use Getopt::Long;

# define vars
my $usage = "USAGE: batch_run_baseml.pl -i /path/to/dir/of/phylip_files <options>

\tRequired:
\t -i|indir      :     The input directory of interleaved phylip and *dnd tree files
\t                     from clustalw. 

\tOptions:
\t --clean       :     Remove all the output of PAML.";

my $in_dir; # = shift or die "\n$usage\n\n";
my $clean;


GetOptions(# Required arguments
	   "i|indir=s"         => \$in_dir,                    # indir (for now)
	   "clean"             => \$clean,
	  );


# NOT SURE THIS IS THE RIGHT SYNTAX FOR open() in this case
if (!$in_dir) {
    die "\n","ERROR: No input directory was given at the command line\n\n",$usage,"\n\n";
} else {
    if ($in_dir) {

	#print $in_dir,"\n";

	opendir(DIR, $in_dir) || die "ERROR: Could not open directory: $in_dir\n";
	
	my @phy_files;
	my @tree_files;
	while (my $file = readdir(DIR)) {
	    if ($file =~ m/\.phy$|\.phylip$/) {
		push(@phy_files,$file);
	    }
	    if ($file =~ m/\.dnd$/) {
		push(@tree_files,$file);
	    }
	}
	closedir(DIR);
	chdir($in_dir);
	my $div_time_files = "Divergence_time_files";
	unless ( -e $div_time_files) {
	    mkdir($div_time_files) || die "ERROR: Could not create directory: $div_time_files\n";
	}

	my $phy_count = @phy_files;
	my $tree_count = @tree_files;
	if ($phy_count == $tree_count) {
	    print "\nThere are $phy_count phylip files and $tree_count treefiles being processed...\n";
	} else {
	    if ($phy_count != $tree_count) {
		die "ERROR: Can not proceed! Any unequal number of sequence files and tree files were found.\n";	
	    }
	}
	foreach my $phy (@phy_files) {

	    foreach my $tree (@tree_files) {

		$phy =~ s/\.phy$//;
		$tree =~ s/\.dnd$//;
		if ($phy =~ $tree) {
		    
		    run_baseml($phy);
		    if ($clean) {

			# remove the PAML output but keep the summary produced
			# by this script in a separate directory.
			system("rm 2base.t rub rst* lnf rates in.basemlg *txt *out");

		    }		    
		}
	    }	    
	}
	chdir("$div_time_files");
	my $summary_file = "all_divergence_times_summary.txt";
	system("cat *txt > $summary_file");
    }
}

sub run_baseml {
    
    my $name     = shift;
    my $divfile  = $name."-divergence.txt";
    my $outfile  = $name."-paml.out";
    my $phylip   = $name.".phy";
    my $treefile = $name.".dnd";
    
    my @ctl_file = "      seqfile = $phylip 
     treefile = $treefile

      outfile = $outfile       * main result file
        noisy = 0   * 0,1,2,3: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

        model = 1   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
                    * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu
        
        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff

*        ndata = 5
        clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below
        kappa = 5  * initial or fixed kappa

    fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below
        alpha = 0.5   * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        ncatG = 5   * # of categories in the dG, AdG, or nparK models of rates
        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK 

        nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states

   Small_Diff = 7e-6
    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
*        icode = 0  * (with RateAncestor=1. try GC in data,model=4,Mgene=4)
*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 0  * Optimization method 0: simultaneous; 1: one branch a time";

    my $control_file = "baseml.ctl";
    open(TMP, ">$control_file") || die "ERROR: Could not open control file: $control_file\n";
    
    print TMP @ctl_file;
    close(TMP);
    system("baseml 2>&1 > /dev/null");
    #my $pwdivfile = "2base.t";     # using the output now to get kappa and divergence
    open(DIVIN, $outfile) || die "ERROR: Could not open divergence file: $outfile\n";
    open(DIVOUT, ">$divfile") || die "ERROR: Could not open divergence file: $divfile\n";
 
    while (<DIVIN>) {
	chomp;
	if (m/^\d\w+\s+\d\.\d+\(\s/) {
	    # 3prime_Ung         0.0269( 8.7752)
	    my ($seqid,$divergence_time,$kappa) = split(/\s+/, $_);
	    
	    $divergence_time =~ s/\($//;
	    $kappa =~ s/\)$//;
	    my $time = $divergence_time/(1e-8 * 2);   # T=k/2r, k=1.0 × 10-8

	    # alignID divergence age Ts:Tv
	    print DIVOUT join("\t",($phylip,$divergence_time,$time,$kappa)),"\n";
	}
	elsif (m/^(\d\w+)         (\d\.\d+\()(\d+\.\d+\))/) {
	    # 3prime_RL1         0.0087(999.0000)
	    
	    my $divergence_time = $2;
	    my $kappa = $3;
	    $divergence_time =~ s/\($//;
	    $kappa =~ s/\)$//;
	    my $time = $divergence_time/(1e-8 * 2);

	    # alignID divergence age Ts:Tv
	    print DIVOUT join("\t",($phylip,$divergence_time,$time,$kappa)),"\n";
	}
	elsif (m/^(\d\w+)         (\d\.\d+\()(\-\d+\.\d+\))/) {
	    # 3prime_RL1         0.0017(-0.0025)

	    my $divergence_time = $2;
	    my $kappa = $3;
	    $divergence_time =~ s/\($//;
	    $kappa =~ s/\)$//;
	    my $time = $divergence_time/(1e-8 * 2);

	    # alignID divergence age Ts:Tv
	    print DIVOUT join("\t",($phylip,$divergence_time,$time,$kappa)),"\n";
	}
    }

    close(DIVIN);
    close(DIVOUT);
    my $summary_file_path = "Divergence_time_files";
    system("cp $divfile $summary_file_path");
    unlink($control_file);           # instead of overwriting, delete so the last one is not left

}

exit;
