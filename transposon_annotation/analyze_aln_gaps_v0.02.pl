#!/usr/bin/perl -w

# TODO: return col matrix for each seq (for multifasta) DONE

#       extract the sequence around the gap to look at direct repeats flanking indels DONE
#       (this is not a general procedures so maybe leave in a dedicated script)

#       May not be calculating stats correctly for a MSE DONE 

#       Write flanking gap seqs to separate file for each input seq in alignment

#       Take search space as option as well as PID for flanking gap
#        sequence matches

#
#

use strict;
use Cwd;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SearchIO;
use lib qw(/iob_home/jmblab/statonse/apps/perlmod/bioperl-run/lib); # Using my own copy of the Run package  
use Bio::Tools::Run::StandAloneBlast;
use Data::Dumper;
use Statistics::Descriptive;
use Getopt::Long;
use File::Basename;
use File::Temp;
use List::MoreUtils qw(uniq);

my $aln_file;
my $outfile;
my $statsfile;
my $dr_pid;
my $gap_stats;

my $pos = 0;
my $gap = 0;
my $del = 0;
my @indels;
my @flanking_seqs;

GetOptions(
	   'i|infile=s'     => \$aln_file,
	   'o|outfile=s'    => \$outfile,
	   's|statsfile=s'  => \$statsfile,
	   'p|repeat_pid=i' => \$dr_pid,
	   );

if (!$aln_file) {
    &usage();
    exit(0);
}

$dr_pid = defined($dr_pid) ? $dr_pid : '10';
my $cwd = getcwd();
open(my $out, '>>', $outfile) or die "\nERROR: Could not open file: $outfile\n";
open(my $stats_out, '>>', $statsfile) or die "\nERROR: Could not open file: $statsfile\n";

my ($seqs_in_aln,$count) = split_aln($aln_file);

foreach my $fas (@$seqs_in_aln) {

    my $aln_in = Bio::AlignIO->new(-fh => \*$fas,       # we are passing an open filehandle, not a filename
				   -format => 'fasta');

    my ($fname, $fpath, $fsuffix) = fileparse($fas, qr/\.[^.]*/);
    my $seq_out = $fname;
    $seq_out .= "_gap_flanking_sequences.fasta";
    open(my $each_out, '>>', $seq_out) or die "\nERROR: Could not open file: $seq_out\n";
    print $each_out join("\t",("Query_ID","Hit_ID","HSP_len","Hit_start","Hit_stop","Query_start","Query_stop","HSP_PID")),"\n";

    while ( my $aln = $aln_in->next_aln() ) {
	my $aln_len = $aln->length;
	foreach my $hash (@{$aln->gap_col_matrix}) {
	    foreach my $key (keys %$hash) {
		$pos++; 
		if ($hash->{$key} eq '1') { # gap columns are coded as '1'
		    $gap++;
		    push(@indels,$pos); 
		}   
	    }
	}
	
	my ($indel_lengths, $indel_ranges) = get_indel_range(@indels);   
	@indels = ();

	$gap_stats = get_stats($fas,$pos,$gap,$indel_lengths,$fname) if $statsfile;
	
    	foreach my $line (@$indel_ranges) {
	    $del++;
	    my ($indel_len, $indel_spos, $indel_epos) = split(/\t/,$line); 

	    next unless $indel_len >= 10; # Ma et al. 2004, Genome Research
	    
	    # What we want to do is search 20 bp around the gap. 
	    # this is based on the 1-15 direct repeats found in Arabidopsis... (Devos et al. 2002)
	    # sunflower may be different
	    my $upstream_spos = $indel_spos - 20;
	    my $upstream_epos = $indel_spos - 1;
	    my $downstream_spos = $indel_epos + 1;
	    my $downstream_epos = $indel_epos + 20;
        
	    if ($upstream_spos < 1 || $upstream_epos > $aln_len) {                                                    
		print "Deletion $del has a flanking repeat out of bounds.\n";
		next;
	    }
	    if ($downstream_epos > $aln_len) {
		print "Deletion $del has the downstream flanking repeat out of bounds.\n";
		last;
	    }
	    
	    foreach my $seq ($aln->each_seq) { # seq "is a" a Bio::LocatableSeq object (which is a part of Bio::PrimarySeq)

		my $upstr_seq = $seq->subseq($upstream_spos,$upstream_epos);
		my $downstr_seq = $seq->subseq($downstream_spos,$downstream_epos);
		
		my $upstr_id = $seq->id."_upstr-del-".$del."_".$upstream_spos."-".$upstream_epos;
		my $downstr_id = $seq->id."_downstr-del-".$del."_".$downstream_spos."-".$downstream_epos;
		
		my $upstream_seqobj = Bio::Seq->new(-seq => $upstr_seq,
						    -id => $upstr_id,
						    -alphabet => 'dna');
		
		my $downstream_seqobj = Bio::Seq->new(-seq => $downstr_seq,
						      -id => $downstr_id,
						      -alphabet =>'dna');
		
		# Here is where we do the comparison
		blast_compare($indel_len,$upstream_seqobj,$downstream_seqobj,$cwd,$each_out);

	    }       
	    
	}
	
    }
    $pos = 0;
    $del = 0;
    close($each_out);
    collate($seq_out,$out);
    collate($gap_stats,$stats_out);
    unlink($fas);
    unlink($seq_out);
    unlink($gap_stats);
}

close($out);
close($stats_out);

exit;
#
# Subroutines
#
sub split_aln {

    my ($input) = @_;

    my ($iname, $ipath, $isuffix) = fileparse($input, qr/\.[^.]*/);

    my $seq_in  = Bio::SeqIO->new(-file  => $input,
                                  -format => 'fasta');
    my $count = 0;
    my $fcount = 1;
    my @split_files;
    $iname =~ s/\.fa.*//;     # clean up file name like seqs.fasta.1
    
    my $tmpiname = $iname."_".$fcount."_XXXX";
    my $fname = File::Temp->new( TEMPLATE => $tmpiname,
                                 DIR => $cwd,
                                 UNLINK => 0,
				 SUFFIX => ".fasta");
    
    my $seq_out = Bio::SeqIO->new(-file   => ">$fname", 
                                  -format => 'fasta');

    push(@split_files,$fname);
    while (my $seq = $seq_in->next_seq) {
        if ($count % 1 == 0 && $count > 0) {
            $fcount++;
            $tmpiname = $iname."_".$fcount."_XXXX";
            $fname = File::Temp->new( TEMPLATE => $tmpiname,
                                      DIR => $cwd,
                                      UNLINK => 0,
				      SUFFIX => ".fasta");

            $seq_out = Bio::SeqIO->new(-file   => ">$fname", 
                                       -format => 'fasta');

            push(@split_files,$fname);
        }
        $seq_out->write_seq($seq);
        $count++;
    }

    my @unique_split_files = uniq(@split_files);
    return (\@unique_split_files,$count);
}

sub get_indel_range {

    my @indels = @_;
    my @indel_ranges;
    my @indel_lengths;

    # The algorithm below is based on 
    # something I found on stackoverflow
    # for converting an array with sequential 
    # numbers into an array with ranges
    # http://stackoverflow.com/questions/117691
    my $gap_range = [ $indels[0] ];

    for my $indel (@indels[1..$#indels]) {
        if ($indel == $gap_range->[-1] + 1) {
            push(@$gap_range,$indel);
        }
        else {
            my $gap_length = ($gap_range->[-1] - $gap_range->[0]) + 1;
            push(@indel_lengths,$gap_length);
            push(@indel_ranges,$gap_length."\t".$gap_range->[-1]."\t".$gap_range->[0]."\n");
            $gap_range = [ $indel ];
        }
    }
    return(\@indel_lengths,\@indel_ranges);

}

sub blast_compare {

    my ($indel_len,$upstream_seqobj,$downstream_seqobj,$cwd,$each_out) = @_;
	
    my $bl2_out = File::Temp->new(TEMPLATE => "bl2seq_out_XXXX",
				  DIR => $cwd);
    my $bl2_out_fname = $bl2_out->filename;

    my @params = (-F => 'F', -p => 'blastn', -W => 4, -o => $bl2_out_fname);
    my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
								    		
    $factory->bl2seq($upstream_seqobj,$downstream_seqobj);
    
    my $bl2seq_report = Bio::SearchIO->new(-file => $bl2_out_fname,
					   -format => 'blast');
    
    while( my $result = $bl2seq_report->next_result ) {
	
	my $query      = $result->query_name();
	my $qlen       = $result->query_length();
	
	while( my $hit = $result->next_hit ) {
	    
	    my $hitid    = $hit->name();
	    
	    while( my $hsp = $hit->next_hsp ) {
		
		my $hsplen    = $hsp->length('total');
		my $hstart    = $hsp->start('hit');
		my $hstop     = $hsp->end('hit'); chomp($hstart);
		my $qstart    = $hsp->start('query'); chomp($qstart);
		my $qstop     = $hsp->end('query');
		my $qstring   = $hsp->query_string;
		my $hstring   = $hsp->hit_string;
		my $hpid      = $hsp->percent_identity;
		
		if( $hsplen > 2 && $hpid >= $dr_pid ) {
		    my $hitseq = $downstream_seqobj->seq;
		    next unless $hstring =~ /$hitseq/gi;
		    print $each_out join("\t",($query,$hitid,$hsplen,$hstart,$hstop,$qstart,$qstop,$hpid)),"\n\n";     
		    my $hit_gap_char =()= $hstring =~ /\-/gi;
		    my $hit_pad = $hit_gap_char + $hstart;
		    #my $five_pr_upst_len = length(length($upstream_seqobj->seq) - $qstart);
		    #my $leading_upstr_string = substr($upstream_seqobj->seq,1,$five_pr_upst_len);
		    #my $upstr_gap = 0; my $downstr_gap;
		    #while ($upstream_seqobj->seq =~ /\-/g) { $upstr_gap++; }
		    #while ($downstream_seqobj->seq =~ /\-/g) { $downstr_gap++; }
		    #$hstart += $upstr_gap; $qstart += $downstr_gap;
		    #my $qstart_filled = sprintf("%${qstart}s", uc($qstring));
		    #my $hstart_filled = sprintf("%${hstart}s", uc($hstring));
		    $qstring =~ s/^(.*)/' ' x ($qstart) . $1/mge; #http://stackoverflow.com/questions/670693
		    $hstring =~ s/^(.*)/' ' x ($hit_pad) . $1/mge;
		    
		    print $each_out "Gap length        : $indel_len\n";
		    print $each_out "Query string      : ",$upstream_seqobj->seq,"\n";
		    print $each_out "Query match string:",uc($qstring),"\n"; 
		    print $each_out "Hit match string  :",uc($hstring),"\n";
		    print $each_out "Hit string        : ",$downstream_seqobj->seq,"\n\n";
		    
		}
		
	    }
	    
	}
	
    }

    
}

sub get_stats {

    my ($fas,$pos,$gap,$indel_lengths,$fname) = @_;

    my $gap_stats = $fname."_gap_stats.txt";
    open(my $statsout, '>', $gap_stats) or die "\nERROR: Could not open file: $gap_stats\n";
    my $stat = Statistics::Descriptive::Full->new;

    $stat->add_data(@$indel_lengths);

    my ($fasname, $faspath, $fassuffix) = fileparse($fas, qr/\.[^.]*/);
    my $count = $stat->count;
    my $gap_percent = sprintf("%.2f",$gap/$pos);
    my $mean = sprintf("%.2f",$stat->mean);
    my $min = $stat->min;
    my $max = $stat->max;

    print $statsout "$fasname\t$gap\t$count\t$gap_percent\t$mean\t$min\t$max\n";
    
    close($statsout);
    return($gap_stats);

}

sub collate {

    my ($file_in, $fh_out) = @_;
    open(my $fh_in, '<', $file_in) or die "\nERROR: Could not open file: $file_in\n";
    while(<$fh_in>) {
	print $fh_out $_;
    }
    close($fh_in);

}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -i seqs.fas -o blast_result 

Required:
    -i|infile        :    An alignment file in fasta format.
    -o|outfile       :    File name to write the extracted sequences to.

Options:
    -s|statsfile     :    A file to write alignment stats to.
    -p|repeat_pid    :    The percent identity threshold for retaining repeats that flank gaps. (Default: 10, which is a 2 bp match).

END
}
