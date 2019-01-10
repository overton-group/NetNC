#!/usr/bin/perl -w
#
# netNC defines focus networks based on network neigbourhood clustering
# of induced subgraph (eg defined by ChIP-seq peaks or differential expression). 
# P-value distributions are calculated for the induced subgraph and 
# resampled sets, providing for edge-centric false discovery rate estimation
# node-centric degree-normalised values are calculated and Gaussian mixture
# modelling provides for separation pf noise (low-scoring) from signal
# (high-scoring) nodes 

# @author Ian Overton (i.overton@qub.ac.uk)
# version 2.2

#    Copyright (C) 2018 Queen's University Belfast
#
#
#    netNC_v2pt2.pl is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    netNC_v2pt2.pl is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with netNC_v2pt2.pl in the file LICENSE.txt
#    If not, see <http://www.gnu.org/licenses/>.
#

use lib '/path/to/directory/netNC/modules/'; # EDIT THIS LINE TO INDICATE THE PATH TO NetNC MODULES DIRECTORY
use lib "/path/to/the/perl/modules/directory/";  # EDIT THIS LINE TO INDICATE THE PATH TO THE PERL MODULES DIRECTORY CONTAINING Math::Pari (e.g. /usr/lib/perl5/site_perl/5.18.1/x86_64-linux-thread-multi/)
my $path_to_R = "/path/to/R"; ## EDIT THIS LINE TO INDICATE THE PATH TO YOUR R INSTALLATION (e.g. /usr/bin/R)
my $NetNC_home = "/path/to/directory/netNC/"; ## EDIT THIS LINE TO PROVIDE THE PATH TO THE DIRECTORY WHERE YOU HAVE THIS NETNC SCRIPT AND FCS.pl

use strict;
use NetNC;
use Statistics::FDR;
use Getopt::Std;

# parse command-line options
my %opt;
getopts ('n:i:t:o:p:z:s:l:d:EFMBxhc', \%opt);

my $help_text = "## NetNC defines focus networks based on network neigbourhood clustering
## of an induced subgraph (e.g. defined by ChIP-seq peaks or differential expression). 
## P-value distributions are calculated for the induced subgraph and 
## resampled sets, providing for edge-centric false discovery rate estimation.
## Optionally, node-centric degree-normalised values are calculated and Gaussian mixture
## modelling provides for unsupervised separation of noise nodes (low-scoring) from
## coherent nodes (high-scoring). Iterative minimum cut may also be run to improve
## identification of functionally coherent subgraphs ('Functional Target Identification' mode) 
## or to provide densely connected individual clusters ('Pathway Identification') mode.\n\n";

my $example_usage = "EXAMPLE USAGE:

# Functional Target Identification mode
NetNC_v2pt2.pl -n MyNetwork.txt -i NodeList.txt -o /path/to/my/outdir/fileprefix -F

# Pathway Identification mode and Node-Centric Analysis, also specifying a background list for resampling (e.g. genes detected from a microarray experiment)
NetNC_v2pt2.pl -n MyNetwork.txt -i NodeList.txt -o /path/to/my/outdir/fileprefix -E -M -l /path/to/backgroundGenelist.txt

# Node centric analysis using background negative log p-values generated from a previous run
NetNC_v2pt2.pl -n MyNetwork.txt -i NodeList.txt -o /path/to/my/outdir/fileprefix -M -p /path/to/my/outdir/fileprefix.BG.nlPonly.txt\n\nPlease note that the \'node list\' does not have to be a strict subset of nodes in the network provided using the -n option\n";

# print help message if needed
die "$help_text" . "USAGE: NetNC_v2pt2.pl -n [network] -i [test_NodeList] -o [output file prefix]\n\nOPTIONAL ARGUMENTS:\n-E [Flag that specifies 'Pathway Identification' mode, including iterative minimum cut]\n-F [Flag that specifies 'Functional Target Identification' mode, including iterative minimum cut]\n-M [Flag to enable node-centric mixture modelling analysis]\n\n-t [NetNC pFDR threshold (default 0.1)]\n-p [precomputed list of resampled negative log pvalues (serving as H0 for FDR calculation)]\n-z [number of resamples for H0 generation, default 100]\n-s [previous resampling (file of node lists)]\n-x [flag to avoid writing resampled negative log p-values to disk]\n-l [background list (e.g. detected genes from microarray experiment)]\n-B [flag to prevent removal of doublet components (i.e. a graph of two nodes connected by a single edge) during FDR thresholding]\n\n-c [flag to enable iterative minimum cut]\n-d [minimum cut density threshold, default 0.306]\n\n-h [print this message]\n\nNetwork is required in tab-delimited format:\nnode_id1\\tnode_id2\\tweight\\n\n(Note: \\t indicates a tab, \\n a newline)\n\nBackground and test node lists are required in newline-delimited format:\nNode_id\\n\n\n$example_usage" if ($opt{"h"});

# establish vars object, including command-line options
my $var = NetNC::make_netnc_vars(\%opt);

# read in the test nodes (e.g. from ChIP-seq)
my $fg = $var->read_list();

# read the network
my $network = $var->read_network_lite();

# Extract the network corresponding to genes of interest
my $subset_network = NetNC::extract_subnetwork($network, $fg, $var);

# Make contingency tables for the subset network
my ($contingency_fg, $cointeractors) = $subset_network->subset_network_contingency($var);
$contingency_fg->unique_contingency();

# calculate hypergeometric p-values (FET) for network neighbourhood connectivity
$var->{output_location} = $var->{output_location_root} . ".FG"; # indicating foreground (test) data
my $fg_p = $contingency_fg->fisher_contingency_uniq($var);
$contingency_fg->destroy;

# Generate resampled data for null distribution
my $resample = $network->resample_network($var); # make background draws from network

# calculate hypergeometric p-values (FET) for network neighbourhood connectivity on resampled data
my $bg_p = $resample->resample_contingency($var, $network);

# calculate FDR values
my $fdr = $var->Statistics::FDR::calculate_FDR($fg_p, $bg_p);

## Apply the appropriate FDR threshold and write out thresholded datafiles including sif format
$var->Statistics::FDR::FDR_threshold($fdr, $fg_p);

# Optionally run mixture modelling and iterative minimum cut
$var->mixmod_mincut($NetNC_home, $path_to_R);

exit; 



