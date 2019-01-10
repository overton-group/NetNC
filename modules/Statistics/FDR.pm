#    Copyright (C) 2018 Queen's University Belfast
#
#
#    FDR.pm is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    FDR.pm is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with FDR.pm in the file LICENSE.txt
#    If not, see <http://www.gnu.org/licenses/>.
#

package Statistics::FDR;

#
# Functions for FDR calculation and thresholding in the NetNC software
#

#@author Ian Overton - i.overton@qub.ac.uk
# v1.1

use strict;

######
#
# new
#
######
#
# constructor method
#

sub new {

my $inv = shift;
my $class = ref($inv) || $inv;
my $self = { @_ };
bless ($self, $class);
return $self;

}

#####################
#
# calculate FDR
#
#####################

#
# Calculates global FDR by looping through the edge scores in the FG dataset (from high to low)
# and counts the number of BG edges above the FG netNC score (by looping through the BG edges, inside the FG loop, 
# keeping running totals to save compute)


sub calculate_FDR {

	my ($var, $fg, $bg) = @_;
	my @fg_keys = @{$fg->{summary}};
	my @bg_keys = @{$bg->{summary}};	
	my $fg_all_count = scalar(@fg_keys);
	my $bg_all_count = scalar(@bg_keys);
	my $all_count = $bg_all_count + $fg_all_count;
	my $pi_zero = $bg_all_count / $all_count;
	my $fdr = Statistics::FDR->new; # object to hold q-values
	my %fg_runTot; # running total for FG
	my %bg_runTot; # running total for BG
	my $prev_fgk; # key to access running total value at next step for FG
	my $prev_bgk; # key to access running total value at next step for BG
	my $corr = 1 / $var->{resample_iterations}; # factor to correct for multiple resampling in the negative dataset
	# make a sorted BG array
	my @sorted_bgks = sort { $b <=> $a } @bg_keys; # descending numeric sort (high to low)
	
	foreach my $fgk (sort { $b <=> $a } @fg_keys) { # descending numeric sort order (high to low)
		if ($prev_fgk) {
			my $rt_fg = $fg_runTot{$prev_fgk};
			$rt_fg += $fg->{$fgk};
			$prev_fgk = $fgk;
			$fg_runTot{$prev_fgk} = $rt_fg;
			
			# code to loop through the BG array until the score is greater than the FG score
			# and remove the array values 
			
			my $rt_bg = $bg_runTot{$prev_bgk};
			my $count_elements = 0;
			
			BGCH:foreach my $bgk (@sorted_bgks) {
				if ($bgk >= $fgk) {
					$rt_bg += $bg->{$bgk};
					++$count_elements;
				
				} else {
					$bg_runTot{$bgk} = $rt_bg; # this key was not included in summing the rt, so the value is for the keys prior to the key considered				
					$prev_bgk = $bgk;
					last BGCH;
				
				}
			}
			
			# Assign FDR

			# calculate FDR for the FG score, not storing a summary, because this is included in the $fg object 
			if ($rt_bg == 0) {
				$fdr->{$fgk} = 0;
				
			} else {
				$fdr->{$fgk} = ($rt_bg * $corr) / $rt_fg; # use the BG to estimate the proportion of negative examples (ie false positives) in the FG data
			
			}
			
			# remove the values already counted from the bg array
			splice(@sorted_bgks, 0, $count_elements);
			
		} else { # no previous fgk so no previous bgk
			$prev_fgk = $fgk;
			my $rt_fg = $fg->{$fgk};
			my $rt_bg = 0;
			
			$fg_runTot{$prev_fgk} = $rt_fg;
			
			# code to loop through the BG array until the score is greater than the FG score
			# and remove the array values 
			my $count_elements = 0;
			BGCH:foreach my $bgk (@sorted_bgks) {
				if ($bgk >= $fgk) {
					$rt_bg += $bg->{$bgk};
					++$count_elements;
					
				} else {
					$bg_runTot{$bgk} = $rt_bg;				
					$prev_bgk = $bgk;
					last BGCH;
				
				}
			}
			
			# calculate FDR for the FG score, not storing a summary, because this is included in the $fg object 
			if ($rt_bg == 0) {
				$fdr->{$fgk} = 0;
				
			} else {
				$fdr->{$fgk} = ($rt_bg * $corr) / $rt_fg; # use the BG to estimate the proportion of  negative examples (ie false positives) in the FG data
			
			}
			
			# remove the values already counted from the bg array
			splice(@sorted_bgks, 0, $count_elements);

		}
		
		
	
	}
	
	# loop through the FDR values to ensure monotonicity (from low to high ensuring fdr_n >= fdr_n+1)
	my @sorted_fgks = (sort { $a <=> $b } @fg_keys); # ascending numeric sort - identify minimum pFDR associated with a given score (ie q-value) Storey & tibrishani 2003 PNAS, Storey 2002
	my $prev_fdr;
	foreach my $fgk (@sorted_fgks) {
		my $fdr_val = $fdr->{$fgk};
		
		if (defined($prev_fdr)) {
			if ($fdr_val > $prev_fdr) {
				$fdr->{$fgk} = $prev_fdr;
			
			} else {
				$prev_fdr = $fdr_val;
			
			}
			
		} else {
			$prev_fdr = $fdr_val;
		}
		
		$fdr->{$fgk} = 1 if ($prev_fdr > 1); # assess agains previous fdr
	}
	
	return $fdr;
	
}


#################
#
# FDR_threshold
#
#################
#
# Prints FDR and nLP values, thresholds edges at FDR value
# 

sub FDR_threshold {

	my ($var, $fdr, $fg) = @_;
	my $outfile = $var->{output_location_root};
	my $fdr_threshold = $var->{fdr_threshold};
	my (%fdr, $minSignif_nlP, $actual_fdr);
	my @fg_keys = @{$fg->{summary}};
	open (FDR, ">$outfile.FDR.txt")  or die "cant open FDR $outfile.FDR.txt $!";
	
	# assign FDR threshold to NLP value and print out the FDR assignemnt for NLPs
	foreach my $fgk (sort { $b <=> $a } @fg_keys) {
		if ($fdr_threshold >= $fdr->{$fgk}) {
			$minSignif_nlP = $fgk;
			$actual_fdr = $fdr->{$fgk};
			print FDR sprintf ("%.5e %s %.10f\n", "$fgk" , "\t", $fdr->{$fgk});
		
		} else {
			print FDR sprintf ("%.5e %s %.10f\n", "$fgk" , "\t", $fdr->{$fgk});
		
		}
	}
	
	if (defined ($minSignif_nlP)) {
		## Write out the significant interactors
		open (P, ">$outfile.FDRthresholded_pairs.txt") or die "cant open $outfile.PEPthresholded_pairs.txt $!";
		open (I, "$outfile.FG.id1_id2_nlp.txt") or die "cant open I $outfile.FG.id1_id2_nlp.txt $!";
		open (M, ">$outfile.minSignif_nlP.log");
		print M "The -log p-value threshold value was set to $minSignif_nlP - which was the closest match to FDR <= $fdr_threshold \nThe actual FDR at this threshold was estimated to be $actual_fdr\n";
		print M "Resample iterations = " . $var->{resample_iterations} . "\n";
		close M;
		
        # write out the significant pairs in tab-delimited (weighted) format
		while (<I>) {
			chomp();
			my @data = split (/\s+/, $_);	
			if ($data[2] >= $minSignif_nlP) {
				print P $data[0] . "\t" . $data[1] . "\t" . $data[2] . "\n";
				
            } 
        }
	} else {
		open (M, ">$outfile.minSignif_nlP.log");
		print M "No -log p-value threshold value met the criterion FDR <= $fdr_threshold \nTherefore no threshold value was assigned\n";
		print M "Resample iterations = " . $var->{resample_iterations} . "\n";
		close M;

	}
	
	close I;
	close P;

}


1;
