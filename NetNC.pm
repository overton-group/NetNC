#    Copyright (C) 2018 Queen's University Belfast
#
#
#    NetNC.pm is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    NetNC.pm is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with NetNC.pm in the file LICENSE.txt
#    If not, see <http://www.gnu.org/licenses/>.
#

package NetNC;

use strict;
use Statistics::FishersExactTest;
use Math::BigFloat;

#@author Ian Overton - i.overton@qub.ac.uk
# Version 1.5

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

###################
#
# make_netnc_vars
#
###################
#
# establishes vars object

sub make_netnc_vars { 
	my $optref = shift;
	my %opt = %$optref;
	my $var = NetNC->new; 
	$var->{threshold} = 0;
	$var->{track_ids} = 1;
	
	# command line options
	$var->{output_location_root} = ($opt{"o"}) ? $opt{"o"} : 0;
	$var->{network_location} = ($opt{"n"}) ? $opt{"n"} : 0;
	$var->{test_nodes} = ($opt{"i"}) ? $opt{"i"} : 0;
	$var->{pvalue_file} = ($opt{"p"}) ? $opt{"p"} : 0;
	$var->{fdr_threshold} = ($opt{"t"}) ? $opt{"t"} : 0.1;
	$var->{resample_file} = ($opt{"s"}) ? $opt{"s"} : 0;
	$var->{no_bg_pval_write} = ($opt{"x"}) ? $opt{"x"} : 0;
	$var->{background_genelist} = ($opt{"l"}) ? $opt{"l"} : 0;
	$var->{mixmodel} = ($opt{"M"}) ? $opt{"M"} : 0;
	$var->{mincut} = ($opt{"c"}) ? $opt{"c"} : 0;
	$var->{mincut_density_threshold} = ($opt{"d"}) ? $opt{"d"} : 0.306;
	$var->{functional_target_identify} = ($opt{"F"}) ? $opt{"F"} : 0;
	$var->{pathway_identify} = ($opt{"E"}) ? $opt{"E"} : 0;
	$var->{remove_doublets} = 0 ? $opt{"B"} : 1;
	
	if ($var->{pvalue_file}) { 
		warn "\nWARNING: Assuming resample_iterations=100 (default value) for FDR calculation based on the background p-value file provided by the -p option.\nIf the default settings were not used, please restart NetNC including the -z option with a value that matches the number of resample iterations used to generate the supplied background p-value file; otherwise, FDR estimation will be compromised. For example, if a background p-value file is generated with -z 500 and then given to a subsequent run of NetNC without a value in the -z option the FDR calculation will be out by a factor of 5.\nThe resample_iterations value for NetNC runs is recorded in the file with extension .minSignif_nlP.log\n" unless ($opt{"z"});
		$var->{resample_iterations} = ($opt{"z"}) ? $opt{"z"} : 100; # used to determine the number of resamples for pFDR calculation
		
	} else {
		$var->{resample_iterations} = ($opt{"z"}) ? $opt{"z"} : 100; # used to determine the number of resamples for pFDR calculation
	}
	
	if ($var->{functional_target_identify}) {
		die "ERROR: options -E (pathway identification) and -F (functional target identification) are mutually exclusive\n" if ($var->{pathway_identify});
		## set paramters for FTI - these are the 'noise-robust' parameters
		$var->{mincut} = 1;
		$var->{mincut_density_threshold} = 0.306;
		$var->{fdr_threshold} = 0.12025;
		print STDOUT "\n NetNC running in Functional Target Identification Mode\n";
		
	} elsif ($var->{pathway_identify}) {
		die "ERROR: options -E (pathway identification) and -F (functional target identification) are mutually exclusive\n" if ($var->{functional_target_identify});
		## set paramters for PID - these are the 'noise-robust' parameters
		$var->{mincut} = 1;
		$var->{mincut_density_threshold} = 0.5044;
		$var->{fdr_threshold} = 0.106;
		print STDOUT "\nNetNC running in Pathway Identification Mode\n";
	}
	
	# print help message if insufficient command-line options are provided
	die "USAGE: NetNC_v2pt2.pl -n [network] -i [test_NodeList] -o [output file prefix]\n\nOPTIONAL ARGUMENTS:\n-E [Flag that specifies 'Pathway Identification' mode, including iterative minimum cut]\n-F [Flag that specifies 'Functional Target Identification' mode, including iterative minimum cut]\n-M [Flag to enable node-centric mixture modelling analysis]\n\n-t [NetNC pFDR threshold (default 0.1)]\n-p [precomputed list of resampled negative log pvalues (serving as H0 for FDR calculation)]\n-z [number of resamples for H0 generation, default 100]\n-s [previous resampling (file of node lists)]\n-x [flag to avoid writing resampled negative log p-values to disk]\n-l [background list (e.g. detected genes from microarray experiment)]\n-B [flag to prevent removal of doublet components (i.e. a graph of two nodes connected by a single edge) during FDR thresholding]\n\n-c [flag to enable iterative minimum cut]\n-d [minimum cut density threshold, default 0.306]\n\n-h [print help message]\n" unless ((($var->{network_location}) and ($var->{test_nodes})) and ($var->{output_location_root})); 
	return $var;
}

#################
#
# read_list
#
#################
#
# reads a list of genes (or any data list)
# identifiers given at beginning of the line and no whitespace is included

sub read_list {

	my $var = shift;
	my $list = $var->{test_nodes};
	die "no list file given to NetNC::read_list $!\n" unless ($list);
	open (IN, $list) or  die "cant open $list in NetNC::read_list $!\n";
	
	my $obj = NetNC->new;
	$obj->{structure} = "obj->{summary} = tabdelimited summary of all identifiers\nobj->{identifier} = 1\n";
	
	my %seen;
	while (<IN>) {
		chomp();
		if (/^(\S+)/) {
			my $a = $1;
			$obj->{summary} .= "$a\t" unless ($seen{$a});		
			$obj->{$a} = 1;
			$seen{$a} = 1;	
		}	
	}

	$var->{total_test_nodes} = scalar(keys %seen);
	return $obj;
	close IN;


}

#################
#
# read_network_lite 
#
#################
#
# reads an undirected network in the format: ID1 ID2 edgeWeight
# only reads half the possible matrix
#
# Input is network file
# returns an object with the information
# object->{structure} gives information on what is in the object

sub read_network_lite {

	my $var = shift;
	my $net = $var->{network_location};
	die "no network file given to NetNC::read_network $!\n" unless ($net);
	open (IN, $net) or  die "cant open $net in NetNC::read_network $!\n";
	
	my $obj = NetNC->new;
	$obj->{structure} = "obj->{node_summary} = tabdelimited summary of all nodes\nobj->{edge_values}{node1_id,node2_id} = values for edges\n";
	
	my %seen;
	while (<IN>) {
		chomp();
		my @a = split(/\t/, $_);
		$obj->{node_summary} .= $a[0] . "\t" unless ($seen{$a[0]});
		$obj->{node_summary} .= $a[1] . "\t" unless ($seen{$a[1]});
		$obj->{edge_values}{$a[0]  . "," . $a[1]} = $a[2] unless (($seen{$a[1] . "," . $a[0]}) or ($seen{$a[0]
		. "," . $a[1]}));

		$seen{$a[1] . "," . $a[0]} = 1;
		$seen{$a[1]} = 1;
		$seen{$a[0]} = 1;		
		
	}

	close IN;
	return $obj;


}

########################
#
# extract_subnetwork
#
########################
#
# Pulls out the connections between a subset of nodes
# with an optional threshold value
#
# Inputs: object from read_list 
#         object from read_network


sub extract_subnetwork {
	my ($net, $list, $var) = @_;

	# set the threshold
	my $threshold;
	$var->{threshold} = 0 unless ($var->{threshold});
	$threshold = $var->{threshold};
	
	my $obj = NetNC->new;
	$obj->{structure} = "\$obj->{summary} = summary of all edges\nobj->{edge_values}{\$genepair} = edge values for the gene pairs\n";
	my @nodes = split (/\t/, $list->{summary});
	
	my %seen;

	$obj->{structure} .= "\$obj->{connections}{\$node} = all nodes that connect to that node"; 
	foreach my $n (@nodes) {
		foreach my $nn (@nodes) {
			next if ($n eq $nn);
			next if (($seen{$n . "," . $nn}) or ($seen{$nn . "," . $n}));
		
			if ($net->{edge_values}{$n  . "," . $nn}) {
				next unless ($net->{edge_values}{$n  . "," . $nn} >= $threshold);
				$obj->{summary} .= "$n;;;$nn\t";
				$obj->{edge_values}{$n  . ";;;" . $nn} = $net->{edge_values}{$n  . "," . $nn};
				$obj->{connections}{$n} .= "$nn\t";
				$obj->{connections}{$nn} .= "$n\t";
				
			#	print "$n, $nn, " . $obj->{connections}{$n}  . "\n";
			#	exit;
				
			} elsif ($net->{edge_values}{$nn  . "," . $n}) {
				next unless ($net->{edge_values}{$nn  . "," . $n} >= $threshold);
				$obj->{summary} .= "$nn;;;$n\t";
				$obj->{edge_values}{$nn  . ";;;" . $n} = $net->{edge_values}{$nn  . "," . $n};
				$obj->{connections}{$n} .= "$nn\t";
				$obj->{connections}{$nn} .= "$n\t";
			}
			
			$seen{$n . "," . $nn} = 1;
		}
	}
	
	return $obj;
}

###########################
#
# subset_network_contingency
#
###########################
#
# generates contingency information for the node in a subset network
# only making tables for genes that have direct interactions
#
# Input is the network and the list object used to generate the network
# Modified to include total number of nodes to be the total that have interactions in the dataset
# Modified to return A B C D values for new FET code


sub subset_network_contingency {
	# set up variables, make objects
	my ($net, $var) = @_;
	my $contingency = NetNC->new;
	my $cointeractor_list = NetNC->new;
	my $pair_int = NetNC->new; # object to receive co-interaction counts
	
	$contingency->{keepIDs} = $var->{track_ids};
	
	my %reset; # for flushing hashes from memory
		
	# identify nodes that have interactions with other subset nodes
	my @edges = split(/\t/, $net->{summary});
	my %subset_interactor_nodes;
	foreach my $ed (@edges) {
		my ($n1, $n2) = split (/;;;/, $ed);
		$subset_interactor_nodes{$n1} = 1;
		$subset_interactor_nodes{$n2} = 1;
	
	}
	
	my $total_subset_interactor_nodes = scalar(keys %subset_interactor_nodes);
	my %seen_pair;
	my %edge_count;
	
	# identify cointeractors & determine edge count values	
	foreach my $id (keys %subset_interactor_nodes) {
		my @connected = split(/\t/, $net->{connections}{$id});
		$edge_count{$id} = scalar(@connected);
		
		foreach my $c (@connected) {	
			# only analyse half of the matrix
			next if ($seen_pair{$id . ",$c"} or $seen_pair{$c . ",$id"});
			$seen_pair{$id . ",$c"} = 1;
			
			my @con_to_connected = split(/\t/, $net->{connections}{$c});
			my %seen;
			foreach my $gene (@con_to_connected) {
				next if ($gene eq $id); # same ids aren't a co-interaction
				next if ($seen{$gene});
				$seen{$gene} = 1;
				if (($net->{edge_values}{"$id;;;$gene"}) or ($net->{edge_values}{"$gene;;;$id"})) {
				
					++$pair_int->{$id . "\t" . $c};
					$cointeractor_list->{$id . "\t" . $c} .= "$gene\t";
				
				}
			
			}
		}
	}
	
	# make contingency tables
	
	my $total_nodes_lesspair = $total_subset_interactor_nodes - 2;
	%seen_pair = %reset;
	
	if ($contingency->{keepIDs}) { # retain IDS (for FG dataset)
	
		foreach my $id (keys %subset_interactor_nodes) {		
			next unless ($net->{connections}{$id}); # must connect to at least one of the genes of interest
			my @connected = split(/\t/, $net->{connections}{$id});
			foreach my $c (@connected) {
				next if ($seen_pair{$id . "\t" . $c});
				my $edge_count_c = $edge_count{$c};
				my $edge_count_id = $edge_count{$id};
				# to account for the interaction between the node pair
				if (($net->{edge_values}{"$c;;;$id"}) or ($net->{edge_values}{"$id;;;$c"})) { # there is an interaction between the node pair
					--$edge_count_id;
					--$edge_count_c; 
				} 
							
				if ($pair_int->{$id . "\t" . $c}) {
					my $cointeract = $pair_int->{$id . "\t" . $c};
					my $cont_b = $edge_count_id - $cointeract;
					my $cont_c = $edge_count_c - $cointeract;
					my $cont_d = $total_nodes_lesspair - ($cont_c + $cont_b + $cointeract);
					$contingency->{perl_fisher} .= "$cointeract\t$cont_b\t$cont_c\t$cont_d\n";
					$contingency->{ids}{"$cointeract\t$cont_b\t$cont_c\t$cont_d"} .= "$id\t$c;;;";
					$seen_pair{$id . "\t" . $c} = 1;
					$seen_pair{$c . "\t" . $id} = 1;
					
				} elsif ($pair_int->{$c . "\t" . $id}) {
				
					my $cointeract = $pair_int->{$c . "\t" . $id};
					my $cont_b = $edge_count_id - $cointeract;
					my $cont_c = $edge_count_c - $cointeract;
					my $cont_d = $total_nodes_lesspair - ($cont_c + $cont_b + $cointeract);
					$contingency->{perl_fisher} .= "$cointeract\t$cont_b\t$cont_c\t$cont_d\n";
					$contingency->{ids}{"$cointeract\t$cont_b\t$cont_c\t$cont_d"} .= "$id\t$c;;;";
					$seen_pair{$id . "\t" . $c} = 1;
					$seen_pair{$c . "\t" . $id} = 1;
				
				} else { # no cointeractors
					my $cointeract = 0;
					my $cont_b = $edge_count_id - $cointeract;
					my $cont_c = $edge_count_c - $cointeract;
					my $cont_d = $total_nodes_lesspair - ($cont_c + $cont_b + $cointeract);
					$contingency->{perl_fisher} .= "0\t$cont_b\t$cont_c\t$cont_d\n";
					$contingency->{ids}{"0\t$cont_b\t$cont_c\t$cont_d"} .= "$id\t$c;;;";
					$seen_pair{$id . "\t" . $c} = 1;
					$seen_pair{$c . "\t" . $id} = 1;
				}
			}
		}
		
	} else {
		
		foreach my $id (keys %subset_interactor_nodes) {
			
			next unless ($net->{connections}{$id}); # must connect to at least one of the genes of interest
			my @connected = split(/\t/, $net->{connections}{$id});
			foreach my $c (@connected) {
				next if ($seen_pair{$id . "\t" . $c});
				my $edge_count_c = $edge_count{$c};
				my $edge_count_id = $edge_count{$id};
				# to account for the interaction between the node pair
				if (($net->{edge_values}{"$c;;;$id"}) or ($net->{edge_values}{"$id;;;$c"})) { # there is an interaction between the node pair
					--$edge_count_id;
					--$edge_count_c; 
				} 
							
				if ($pair_int->{$id . "\t" . $c}) {
				
					my $cointeract = $pair_int->{$id . "\t" . $c};
					my $cont_b = $edge_count_id - $cointeract;
					my $cont_c = $edge_count_c - $cointeract;
					my $cont_d = $total_nodes_lesspair - ($cont_c + $cont_b + $cointeract);
					$contingency->{perl_fisher} .= "$cointeract\t$cont_b\t$cont_c\t$cont_d\n";
					$seen_pair{$id . "\t" . $c} = 1;
					$seen_pair{$c . "\t" . $id} = 1;
					
				} elsif ($pair_int->{$c . "\t" . $id}) {
				
					my $cointeract = $pair_int->{$c . "\t" . $id};
					my $cont_b = $edge_count_id - $cointeract;
					my $cont_c = $edge_count_c - $cointeract;
					my $cont_d = $total_nodes_lesspair - ($cont_c + $cont_b + $cointeract);
					$contingency->{perl_fisher} .= "$cointeract\t$cont_b\t$cont_c\t$cont_d\n";
					$seen_pair{$id . "\t" . $c} = 1;
					$seen_pair{$c . "\t" . $id} = 1;
				
				} else { # no cointeractors
					my $cointeract = 0;
					my $cont_b = $edge_count_id - $cointeract;
					my $cont_c = $edge_count_c - $cointeract;
					my $cont_d = $total_nodes_lesspair - ($cont_c + $cont_b + $cointeract);
					$contingency->{perl_fisher} .= "0\t$cont_b\t$cont_c\t$cont_d\n";
					$seen_pair{$id . "\t" . $c} = 1;
					$seen_pair{$c . "\t" . $id} = 1;
				}
			}
		}
		
	
	}
	return ($contingency, $cointeractor_list);
}


#######################
#
# unique_contingency
#
########################
#
# Input is a subset_contingency object
# adds the attributes uniq_contingency
# and unseen_contingency

sub unique_contingency {
	my $contingency = shift;
	my @contingency_data = split (/\n/, $contingency->{perl_fisher});
	my $tmp_cont = shift;
	my $iter = shift;

	foreach my $key (@contingency_data) {
		if ($contingency->{uniq_contingency}{$key}) {
			$contingency->{uniq_contingency}{$key}++;
			
		} else {
			$contingency->{unseen_contingency} .= "$key\n"; 
			$contingency->{uniq_contingency}{$key} = 1;
			
		}
	}
}

#######################
#
# fisher_contingency_uniq
#
########################
#
# Input subset_contingency object and vars object
# calculates fisher test for the contingency object
# using Statistics::FishersExactTest
# this achieves arbitrary precision using Math::Pari
# prints the data to file

sub fisher_contingency_uniq {
	my $contingency = shift;
	my $var = shift;
	my @data = split (/\n/, $contingency->{unseen_contingency});
	my $p_vals = NetNC->new; # object to keep the p-values for FDR calculation 
	my @p_vals; # array to store pvalues
	$p_vals->{summary} = \@p_vals; # store array reference, used by FDR.pm
	my %seen; # hash to prevent storing of repeat p-values
	if ($contingency->{keepIDs}) {
		my $out = $var->{output_location};
		open (OO, ">$out.id1_id2_nlp.txt") or die "cant open output file $out.id1_id2_nlp.txt $!\n";
		foreach my $co (@data) {
			my ($n11, $n1p, $np1, $npp) = split (/\t/, $co);
			my $p_value = Statistics::FishersExactTest::fishers_exact($n11, $n1p, $np1, $npp);
			my $nlp = -log($p_value);
			$nlp = 0 if ($nlp =~ /0\.E/);
			$nlp = 0 if ($p_value >= 1);
			unless ($seen{$nlp}) {
				push(@p_vals, $nlp);
				$seen{$nlp} = 1;
			}
						
			my @ids = split(/;;;/, $contingency->{ids}{$co});
			foreach my $id (@ids) {
				print OO sprintf ("%s %.5e", "$id\t", "$nlp") . "\n";
				++$p_vals->{$nlp};
			}

		}
		
		close OO;
	
	} else { 
		my $out = $var->{output_location};
		if ($var->{no_bg_pval_write}) { # option to avoid writing resampled p-values to file
			foreach my $co (@data) {
				my ($n11, $n1p, $np1, $npp) = split (/\t/, $co);
				my $p_value = Statistics::FishersExactTest::fishers_exact($n11, $n1p, $np1, $npp);
				my $nlp = -log($p_value);
				$nlp = 0 if ($nlp =~ /0\.E/);
				$nlp = 0 if ($p_value >= 1);
				unless ($seen{$nlp}) {
					push(@p_vals, $nlp);
					$seen{$nlp} = 1;
				}
				
				for (my $n = 0; $n < $contingency->{uniq_contingency}{$co}; ++$n) {
					++$p_vals->{$nlp};
				}
			}
		} else {
			open (O, ">$out.nlPonly.txt");	
			foreach my $co (@data) {
				my ($n11, $n1p, $np1, $npp) = split (/\t/, $co);
				my $p_value = Statistics::FishersExactTest::fishers_exact($n11, $n1p, $np1, $npp);
				my $nlp = -log($p_value);
				$nlp = 0 if ($nlp =~ /0\.E/);
				$nlp = 0 if ($p_value >= 1);
				unless ($seen{$nlp}) {
					push(@p_vals, $nlp);
					$seen{$nlp} = 1;
				}
				
				for (my $n = 0; $n < $contingency->{uniq_contingency}{$co}; ++$n) {
					print O "$nlp\n";
					++$p_vals->{$nlp};
				}
			}
			
			close O;
		}
	}
	return $p_vals;
}

####################
#
# resample_network 
#
###################
#
# Generates resampled lists from a network
# takes vars object and network object as inputs
#

sub resample_network {

	my ($network, $var) = @_;
	my $resample = NetNC->new;
	RN: {	
		# no need to resample if there is already a p-value file provided
		last RN if ($var->{pvalue_file});
		
		# if a precalculated sampling file is provided then read it in
		if ($var->{resample_file}) {
			my $sample_file = $var->{resample_file};
			open (A, $sample_file) or die "cant open $sample_file $!\n";
			my $iteration;
			while (<A>) {
				chomp();
				if (/###\s+([0-9]+)/) {
					$iteration = $1;
					$resample->{max_iteration} = $1; ## added this
					
				} elsif (/^(\S+)/) {
					my $node = $1;			
					$resample->{$iteration} .= "$node\t";
				}
			}
			$var->{resample_iterations} = $resample->{max_iteration} + 1;
			return $resample;
	
		}
		
		my @gene = split(/\t/, $network->{node_summary});
		if ($var->{background_genelist}) {
			my (%bg_genes, @filtered_genes);
			my $bg_list = $var->{background_genelist};
			open (BG, "$bg_list") or die "cant open file $bg_list $!\n";
			while (<BG>){ 
				chomp();
				my $bg_gene = $1 if (/(\S+)/);
				$bg_genes{$bg_gene} = $bg_gene; # hash containing IDs for the background
			}
			
			close BG;
			foreach my $node (@gene) {
				push (@filtered_genes, $node) if ($bg_genes{$node});
	
			}
	
			@gene = @filtered_genes;
			
		}
		
		my $iteration = 0; 
		my $max_iterate = $var->{resample_iterations} ? $var->{resample_iterations} : 0; # if false, has been set to 0 by user (default is 100)
		my $totgene = scalar(@gene);
		my $output = $var->{output_location_root} . "_$max_iterate" . "_resamples.txt";
		my $tot_test_nodes = $var->{total_test_nodes};
		open (O, ">$output") or die "cant open output file $output in NetNC::resample_network $!\n";
		die "The node list for resampling is smaller than the input node list, resampling will not complete\n" if ($tot_test_nodes > $totgene);	
		until ($iteration >= $max_iterate) { 
			my $rand = 0; 
			my %seenrand; 
			my (%randlist, %seenout); 
			print O "### $iteration\n"; 
			
			until ($rand >= $tot_test_nodes) { 
				my $r = int(rand $totgene);
				next if ($seenrand{$r});
				print O "$gene[$r]\n";
				$resample->{$iteration} .= "$gene[$r]\t";
				++$rand; 
				$seenrand{$r} = 1; 
			}  
			
			++$iteration; 
		}
		
		$resample->{max_iteration} = $max_iterate - 1; # counting from 0
		close O;
	}
	return $resample;

}

#####################
#
# resample_contingency 
#
####################
#
# calculates the contingency tables for the resampled data
# 
#
sub resample_contingency {

	my ($resample, $var, $network)  = @_;
	# set some variables
	$var->{output_location} =  $var->{output_location_root} . ".BG";
	$var->{track_ids} = 0; # don't keep track of IDs for background data to save memory
	my $max_iterat = 0;
	$max_iterat = $resample->{max_iteration} unless ($var->{pvalue_file});
	
	# Make the contingency tables for each of the background samples,
	my $contingency_bg = NetNC->new;
	BAK:for (my $iter = 0; $iter <= $max_iterat; ++$iter) {
		last BAK if ($var->{pvalue_file});
		print STDERR "# processing resample $iter\n";
		my $bg = NetNC->new;
		$bg->{summary} = $resample->{$iter};
		my $bg_subset_net = NetNC::extract_subnetwork($network, $bg, $var);
		my ($contingency_tmp, $coint_bg) = NetNC::subset_network_contingency($bg_subset_net, $var); ## track IDs if doing combined P-values
		# minimise memory footprint
		$contingency_bg->{perl_fisher} = $contingency_tmp->{perl_fisher};
		$contingency_bg->unique_contingency();
		
	}
	
	# calculate -log p-values for resampled list
	my ($bg_p, %seen);
	if ($var->{pvalue_file}) {
		my @p_vals;
		$bg_p->{summary} = \@p_vals;
		open (IN, $var->{pvalue_file}) or die "cant open file " . $var->{pvalue_file} . "$! \n";
		while (<IN>) {
			chomp();
			if (/^(\S+)/) {
				my $val = $1;
				my $p = Math::BigFloat->new($val); # high precision is required for reproducibility
				++$bg_p->{$p};
			
				unless ($seen{$p}) {
					push(@p_vals, $p);
					$seen{$p} = 1;
				}
			}
		
		}
		close IN;
	
	} else {
		$bg_p = $contingency_bg->fisher_contingency_uniq($var);
	}

	
	return $bg_p;
}


#############
#
# mixmod_mincut
#
#############
#
# Runs node-centric
# mixture modelling analysis
# and/or edge-centric Iterative 
# minimum cut. Also (optionally)
# removes doublet components if 
# mincut is called

sub mixmod_mincut {

	my ($var, $NetNC_home, $path_to_R) = @_;
	
	# mixture modelling
	if ($var->{mixmodel}) { 
		print STDOUT "\nRunning node-centric analysis\n";
		warn "ERROR in mixture modelling: failed to run FCS.pl $!\n" if (system("$NetNC_home/FCS.pl " . $var->{test_nodes} . " " . $var->{output_location_root} . ".FG.id1_id2_nlp.txt " . $var->{output_location_root} . ".FCS.txt") == 1);
		warn "ERROR in mixture modelling: failed to run NetNCmixmodel.R $!\n" if (system("$path_to_R --slave --no-save --args " . $var->{output_location_root} . ".FCS.txt " . $var->{output_location_root} . "_NFCS-mixturemodel $NetNC_home < $NetNC_home/mixtureModelling/NetNCmixmodel.R") == 1);

	# write network with coherent nodes where edges have bonferroni-corrected p<=0.05
		my $out_root = $var->{output_location_root};
		my %coherent;
		open (CO, $out_root . "_NFCS-mixturemodel.coherentNodes") or warn "cant read coherent nodes $!\n";
		while (<CO>) {
			chomp();
			my $node = $1 if (/^(\S+)/);
			$coherent{$node} = 1;
		
		}
		close CO;
		
		my %candidate_network;
		my $hypothesis_count;
		open (FG, $out_root . ".FG.id1_id2_nlp.txt") or warn "cant read $out_root" . ".FG.id1_id2_nlp.txt $!\n";
		while (<FG>) {
			chomp();
			my @data = split(/\t/, $_);
			$candidate_network{$data[0] . ";;;" . $data[1]} = $data[2] if (($coherent{$data[0]}) and ($coherent{$data[1]}));
			++$hypothesis_count if (defined($data[2]));
		}
		close FG;
		# Bonferroni-corrected edge significance threshold
		my $signif_thresh = -log(0.05 / $hypothesis_count);
		
		open (MNET, ">$out_root" . "_NFCS-mixturemodel_coherentNet.txt") or warn "cant open output file $out_root" . "_NFCS-mixturemodel_coherentNet.txt $!\n";
		foreach my $key (keys %candidate_network) {
			if ($candidate_network{$key} >= $signif_thresh) { ## edge passes statistical significance threshold
			
			my $ids = $key;
			$ids =~ s/;;;/\t/g;
			
			print MNET "$ids\t$candidate_network{$key}\n";
			
			}
				
		}
		close MNET;
	}

	# Iterative mimimum cut
	if ($var->{mincut}) {
		my $mc_dens = $var->{mincut_density_threshold};
		my $mc = $mc_dens; 
		$mc =~ s/\./pt/g;
		print STDOUT "\nRunning iterative minimum cut\n";
		warn "ERROR in iterative minimum cut: failed to run itercut.py $!\n" if (system("python $NetNC_home/mincut/itercut.py -i " . $var->{output_location_root} . ".FDRthresholded_pairs.txt -o " . $var->{output_location_root} . ".FDRthresholded_mincutDensity$mc.txt -t $mc_dens")  == 1);
	
        $var->remove_doublets(); ## remove doublet components
	
	}


}

##########
#
# remove_doublets
#
##########
#
# eliminates components formed from
# only two nodes and writes out the NetNC output
#

sub remove_doublets {
	my $var = shift;
	if ($var->{remove_doublets}) {
		my (%edge_count, %result, %total_connector);
		
		if (($var->{functional_target_identify}) or ($var->{pathway_identify})) {
			my $mc_dens = $var->{mincut_density_threshold};
			my $mc = $mc_dens; 
			$mc =~ s/\./pt/g;
			open (I, $var->{output_location_root} . ".FDRthresholded_mincutDensity$mc.txt") or die "cant open I " . $var->{output_location_root} . ".FDRthresholded_mincutDensity$mc.txt $!";
			open (OUT, ">" . $var->{output_location_root} . ".FDRthresholded_mincutDensity$mc" . "_noDoublets.txt") or die "cant open OUT " . $var->{output_location_root} . ".FDRthresholded_mincutDensity$mc" . "_noDoublets.txt $!";
		
			while (<I>) { ## read through the results and count the interactors
				chomp();
				my @data = split (/\s+/, $_);
				$result{$data[0] . ";;;" . $data[1]} = $_;
				$edge_count{$data[0]} .= "$data[1]\t";
				$edge_count{$data[1]} .= "$data[0]\t";
			}
		
			foreach my $k (keys %edge_count) { ## count the number of unique edges for each node with an edge passing the fdr threshold
				my @connectors = split(/\t/, $edge_count{$k});
				my (@uniq, %seen_connector);
			
				foreach my $connector (@connectors) {
					push (@uniq, $connector) unless $seen_connector{$connector}++;
				}
				$total_connector{$k} = scalar(@uniq);
			
			}
	
			foreach my $k (keys %result) {
				my @data = split (/\s+/, $result{$k});
				next if (($total_connector{$data[0]} == 1) and ($total_connector{$data[1]} == 1));
				print OUT $data[0] . "\t" . $data[1] . "\t" . $data[2] . "\n"; 
			}
	
			close I;
			close OUT;
			print "NetNC Functional Target Identification (FTI) output at:" . $var->{output_location_root} . ".FDRthresholded_mincutDensity$mc" . "_noDoublets.txt" if ($var->{functional_target_identify});
        	print "NetNC Pathway Identification (PID) output at:" . $var->{output_location_root} . ".FDRthresholded_mincutDensity$mc" . "_noDoublets.txt" if ($var->{pathway_identify});
		
    	}
    
    }

}




#########
#
# destroy
#

sub destroy {
	my $o = shift;
	my $d;
	$o = $d;

}




1;

