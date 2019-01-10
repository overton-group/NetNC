#! /usr/bin/perl -w

# this script estimates the proportion of neutrally bound genes in the dataset 
# based on the local FDR, estimated from global FDR values
# local FDR is summated in order to estimate of the proportion of neutral binding across
# the entire dataset

# Consider a toy example: if the top-ranked (best scoring) gene global fdr = 0.1 implies 10% chance of thes gene being 'noise' 
# if the threhold value for the second gene global fdr = 0.2 implies 20% chance of noise globally but 
# 0.3 chance of the second gene being noise ((0.3 + 0.1) / 2 = 0.2), therefore the second gene 
# has local FDR of 0.3. The threshold for a third gene could have global FDR = 1, and the third gene would 
# have a local FDR of 1. The normalised summated lcFDR for these three genes would therefore be 1.4 / 3 = 0.466667
# however, the global FDR is 1 for the same genes.

# @author Ian Overton (i.overton@qub.ac.uk)

#    estimateLocalFDR_conservativeSmoothing.pl is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    estimateLocalFDR_conservativeSmoothing.pl is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with estimateLocalFDR_conservativeSmoothing.pl in the file LICENSE.txt
#   If not, see <http://www.gnu.org/licenses/>.
#
# v1.1

use strict;

my $fd = $ARGV[0]; # node FDR file generated by NodeFDR_forNetNC.pl
die "\nusage: estimateLocalFDR_conservativeSmoothing.pl [NodeFDR_file] > sum_lcFDR_output_file.txt\n\n## this script estimates the proportion of neutrally bound genes in the dataset using local FDR,\n## estimated from node-centric pFDR values produced by NodeFDR_forNetNC.pl\n## The output is the sum of local FDR across the input dataset\n" unless ($fd);

open (A, $fd) or die "cant open $fd $!\n";

my %s;
my $n = 0;
my $totSignal = 0;

while (<A>) {
	chomp();
	my @a = split (/\t/, $_);
	my $fdr = $a[4];
	++$s{$fdr};
}

# loop through the fdr values, analysing per datapoint

# expanding out the formula localFDR = ((n x FDR2) - ((n - X) x FDR1)) / X
# gives: ((n x FDR1 + n x FDR_diff) - (n x FDR1 - (X * FDR1) ) / X
# = ((n x FDR_diff) - (- (X * FDR1))  / X
# = ((n x FDR_diff + (X * FDR1)) / X
# = ((n x FDR_diff) / X) + FDR1

my $prev_fdr = 0;
my $prev_lc = 0;
foreach my $f (sort { $a <=> $b } keys %s) {
	$n += $s{$f};
	my $diff_fdr = $f - $prev_fdr;
	
	my $nFDR_diff = $n * $diff_fdr;
	my $nFDR_diff_divX = $nFDR_diff / $s{$f};
	
	my $lcFDR = $nFDR_diff_divX + $prev_fdr;
	$lcFDR = $prev_lc if ($prev_lc > $lcFDR);
	$lcFDR = 1 if ($lcFDR > 1);
	$totSignal += ($lcFDR * $s{$f});
	$prev_fdr = $f;
	$prev_lc = $lcFDR;
	
}

my $proportion = $totSignal / $n;
print "$fd: $proportion\n";