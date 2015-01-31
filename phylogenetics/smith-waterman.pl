#!/usr/bin/env perl

#use strict; use warnings;

die "usage: $0 <seq 1> <seq 2>\\n" unless @ARGV == 2;

my ($seq1, $seq2) = @ARGV;

my $MATCH    =  1;
my $MISMATCH = -1;
my $GAP      =  1;

my @matrix;

$matrix[0][0]{score}   = 0;
$matrix[0][0]{pointer} = "none";

for (my $j = 1; $j <= length($seq1); $j++) {
    $matrix[0][0]{score}    = 0;
    $matrix[0][$j]{pointer} = "none";
}

for (my $i = 1; $i <= length($seq2); $i++) {
    $matrix[$i][0]{score}   = 0;
    $matrix[$i][0]{pointer} = "none";
}

# fill
my $max_i     = 0;
my $max_j     = 0;
my $max_score = 0;

for (my $i = 1; $i <= length($seq2); $i++) {
    for (my $j = 1; $j <= length($seq1); $j++) {
        my ($diagonal_score, $left_score, $up_score);
        
        # calculate match score
        my $letter1 = substr($seq1, $j-1, 1);
        my $letter2 = substr($seq2, $i-1, 1);       
        if ($letter1 eq $letter2) {
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH; #TODO: this needs to be initialized prior to incrementing
        }
        else {
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
        }
        
        # calculate gap scores
        $up_score   = $matrix[$i-1][$j]{score} + $GAP;              #TODO: same as above, need to initialize this value
        $left_score = $matrix[$i][$j-1]{score} + $GAP;
        
        if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
            $matrix[$i][$j]{score}   = 0;
            $matrix[$i][$j]{pointer} = "none";
            next; # terminate this iteration of the loop
        }
        
        # choose best score
        if ($diagonal_score >= $up_score) {
            if ($diagonal_score >= $left_score) {
                $matrix[$i][$j]{score}   = $diagonal_score;
                $matrix[$i][$j]{pointer} = "diagonal";
            }
            else {
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
            }
        } 
	else {
            if ($up_score >= $left_score) {
                $matrix[$i][$j]{score}   = $up_score;
                $matrix[$i][$j]{pointer} = "up";
            }
            else {
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
            }
        }
        
        # set maximum score
        if ($matrix[$i][$j]{score} > $max_score) {
            $max_i     = $i;
            $max_j     = $j;
            $max_score = $matrix[$i][$j]{score};
        }
    }
}

# trace-back
my $align1 = "";
my $align2 = "";

my $j = $max_j;
my $i = $max_i;

while (1) {
    last if $matrix[$i][$j]{pointer} eq "none";
    
    if ($matrix[$i][$j]{pointer} eq "diagonal") {
        $align1 .= substr($seq1, $j-1, 1);
        $align2 .= substr($seq2, $i-1, 1);
        $i--; $j--;
    }
    elsif ($matrix[$i][$j]{pointer} eq "left") {
        $align1 .= substr($seq1, $j-1, 1);
        $align2 .= "-";
        $j--;
    }
    elsif ($matrix[$i][$j]{pointer} eq "up") {
        $align1 .= "-";
        $align2 .= substr($seq2, $i-1, 1);
        $i--;
    }   
}

$align1 = reverse $align1;
$align2 = reverse $align2;

print "$align1\n";
print "$align2\n";
