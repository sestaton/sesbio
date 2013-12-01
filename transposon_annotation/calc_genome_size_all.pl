#!/usr/bin/env perl

# to run:
# nohup perl calc_genome_size_all.pl 2>&1 > calc_genome_size_all.out &
#
# to calculate genome size:
# perl -pe 'while (<>) { @t = split; if ($t[1] == 60 && $t[2] == 70) { print join "\t", @t[0..2], $t[4], 2*$t[4], "\n" } }' < calc_genome_size_all.out > calc_genome_size_all_adjusted_cvals.tsv

use 5.010;
use strict;
use warnings;
use IPC::System::Simple qw(system);

my $target = q{combined_aster_est_assembly_over300bp.fasta};

my $phe_fafile = q{Phoebanthus_5m_interl.fasta};
my $han_fafile = q{Ann1238_GGCTAC_L004_5m_interl.fasta};
my $har_fafile = q{Harg_5m_interl.fasta};
my $hte_fafile = q{Hteph_5m_interl.fasta};
my $hpo_fafile = q{Hport_5m_interl.fasta};
my $hve_fafile = q{Hvert_5m_interl.fasta};
my $sen_fafile = q{Senecio_5m_interl.fasta};
my $saf_fafile = q{Saff_5m_interl.fasta};
my $ger_fafile = q{Gerb_5m_interl.fasta};
my $tks_fafile = q{TKS_5m_interl.fasta};
my $cp_fafile  = q{CP_5m_interl.fasta};
my $age_fafile = q{Ageratina_5m_interl.fasta};
my $cal_fafile = q{Calyc_5m_interl.fasta};
my $das_fafile = q{Dasyphyllum_5m_interl.fasta};
my $gna_fafile = q{Gnaph_5m_interl.fasta};

my $phe_blfile = q{combined_est_assemb_over300bp_Phoebanthus_5m.bln};
my $han_blfile = q{combined_est_assemb_over300bp_Ann1238_5m.bln};
my $har_blfile = q{combined_est_assemb_over300bp_Harg_5m.bln};
my $hte_blfile = q{combined_est_assemb_over300bp_Hteph_5m.bln};
my $hpo_blfile = q{combined_est_assemb_over300bp_Hport_5m.bln};
my $hve_blfile = q{combined_est_assemb_over300bp_Hvert_5m.bln};
my $sen_blfile = q{combined_est_assemb_over300bp_Senecio_5m.bln};
my $saf_blfile = q{combined_est_assemb_over300bp_Saff_5m.bln};
my $ger_blfile = q{combined_est_assemb_over300bp_Gerb_5m.bln};
my $tks_blfile = q{combined_est_assemb_over300bp_TKS_5m.bln};
my $cp_blfile  = q{combined_est_assemb_over300bp_CP_5m.bln};
my $age_blfile = q{combined_est_assemb_over300bp_Ageratina_5m.bln};
my $cal_blfile = q{combined_est_assemb_over300bp_Calyc_5m.bln};
my $das_blfile = q{combined_est_assemb_over300bp_Dasyphyllum_5m.bln};
my $gna_blfile = q{combined_est_assemb_over300bp_Gnaph_5m.bln};

say STDERR join "\t", "blast_file", "filter_len", "PID", "cval_by_mean", "cval_by_trimmed_mean_15p", 
    "cval_by_trimmed_mean_10p", "cval_by_trimmed_mean_5p", "cval_by_geom_mean", "cval_by_median";

for my $fl (qw(50 60 70 80 90 100)) {
    for my $pid (qw(70 80 90 95)) {
	system([0..5], "perl calc_genome_size_v0.03.pl -i $phe_blfile -q $target -t $phe_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $han_blfile -q $target -t $han_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $har_blfile -q $target -t $har_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $hte_blfile -q $target -t $hte_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $hpo_blfile -q $target -t $hpo_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $hve_blfile -q $target -t $hve_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $sen_blfile -q $target -t $sen_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $ger_blfile -q $target -t $ger_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $saf_blfile -q $target -t $saf_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $tks_blfile -q $target -t $tks_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $cp_blfile -q $target -t $cp_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $age_blfile -q $target -t $age_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $cal_blfile -q $target -t $cal_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $das_blfile -q $target -t $das_fafile -fl $fl -pid $pid");
	system([0..5], "perl calc_genome_size_v0.03.pl -i $gna_blfile -q $target -t $gna_fafile -fl $fl -pid $pid");
    }
}

