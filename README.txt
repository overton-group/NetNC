##
## README for Network Neighbourhood Clustering (NetNC) software
##

## Licensed under the GNU General Public License (GPL) version 3
## you should have recieved a copy of the GNU General Public License
## with the NetNC software in the file LICENSE.txt
## if not, see <http://www.gnu.org/licenses/>.


##
## INTRODUCTION
##

In this README, the usage and testing sections provide examples of 
some ways to run NetNC and its components. The testing section also 
includes description of output files. 
The final two sections summarise contributors and how to cite NetNC.

The NetNC distribution includes three main components:
1) Network Neighbourhood Clustering - NetNC_v2pt2.pl

This program runs Network Neighbourhood Clustering and optionally
calls the components below according to the command-line options
supplied. Requires the perl module Math::Pari

After installation (see below), run 'perl NetNC_v2pt2.pl -h' for 
further details.

Files in the modules/ directory provide functions for NetNC_v2pt2.pl

2) Iterative minimum cut - itercut.py (see directory mincut/ )
The script mincut/itercut.py iteratively cuts edges using weighted min-cut max-flow 
algorithm with a graph density stopping criterion.

Further details are given in mincut/MinCut_README.txt

3) Node centric analysis with Gaussian mixture modelling
The script FCS.pl generates a Node Functional Coherence Scores (NFCS)
which are the sum of node edge weight(s) normalised by degree.

The script mixtureModelling/NetNCmixmodel.R runs Gaussian mixture modelling on
NFCS values, then processes the results to infer functionally coherent nodes
and estimate the proportion of noise; mixtureModelling/netNCmix.R
provides functions for the NetNCmixmodel.R script. A little more detail is given
in the usage section, below.



##
## INSTALLATION
##

1. Network Neighbourhood clustering - NetNC_v2pt2.pl

a) The arbitary precision calculations require the perl module Math::Pari.
Install this from cpan (www.cpan.org). For example from the command line:
cpan
install Math::Pari

b) In the NetNC_v2pt2.pl edit the 'use lib' statements to include:
i) the path to Math::Pari (typically this will be where cpan installs
perl modules) and ii) the full path to the NetNC modules/ (this is the 'modules' 
directory created when you unzipped the NetNC distribution).

Also edit the $path_to_R and $NetNC_home variables to specify 
the path to your R installation and the NetNC home directory (this is the location of 
NetNC_v2pt2.pl and FCS.pl scripts)

The lines to edit in NetNC_v2pt2.pl are indicated with comments near to the start of the 
script. NetNC has been tested with perl v5.16, v5.18 and v5.20

2. Iterative Minimum cut - itercut.py

Install the dependencies: networkx 1.8 and numpy 

For example: easy_install networkx==1.8
	     easy_install numpy

(itercut.py is verified to work with networkx 1.8. Earlier versions are 
not compatible)

More information is given in mincut/MinCut_README.txt

3. The Gaussian mixture modelling requires R to be installed and has been 
tested with R versions 3.0.2, 3.2.2, 3.3.0, 3.3.1 on Linux, but is expected 
to work on Mac OSX, Windows and with other recent versions of R.



##
## USAGE
##

1. NetNC

NetNC can be run with many possible combinations of options, here are some
illustrative examples:
-- Functional Target Identification mode (calls itercut.py):
NetNC_v2pt2.pl -n MyNetwork.txt -i NodeList.txt -o /path/to/my/outdir/fileprefix -F

-- Pathway Identification mode (calls itercut.py) and Node-Centric Analysis
(including Gaussian Mixture Modelling) also specifying a background list for 
resampling (e.g. detected genes from a microarray experiment):
NetNC_v2pt2.pl -n MyNetwork.txt -i NodeList.txt -o /path/to/my/outdir/fileprefix -E -M -l /path/to/backgroundGenelist.txt

-- Node centric analysis (including Gaussian Mixture Modelling) using background 
negative log p-values generated from a previous run of NetNC:
NetNC_v2pt2.pl -n MyNetwork.txt -i NodeList.txt -o /path/to/my/outdir/fileprefix -M -p /path/to/my/outdir/fileprefix.BG.nlPonly.txt -z [resample_number]

Where -z [resample_number] indicates the number of resamples used in generating 
the list of -log p-values provided using the -p option (with suffix '.BG.nlPonly.txt')

A help message is produced by running:
NetNC_v2pt2.pl -h


2. Iterative minimum cut

The iterative minimum cut is (optionally) called by NetNC_v2pt2.pl
If you wish to run the iterative minimum cut as a stand-alone program, the basic
usage is:
python minCut/itercut.py -i graphfile.txt -o cutgraph.txt -t 0.3

Note that -t gives the density stopping criterion, a value other than 0.3 may
be used. The itercut.py help message is produced by running:
python minCut/itercut.py 


3. Node-centric analysis with Gaussian mixture modelling

Node-centric analysis is intended to be (optionally) called by NetNC_v2pt2.pl
This subsection describes components of the analysis that may be called separately. 

The script FCS.pl takes as input the genelist and the file written by NetNC with 
the suffix ".FG.id1_id2_nlp.txt". For example:
FCS.pl input_genelist NetNC_results.FG.id1_id2_nlp.txt FCS_outputFile.txt

FCS.pl output has the format: Node edgePvalueSum degree degree-normalised_edgePvalueSum

Gaussian mixture modelling can be run with the NetNCmixmodel.R script as follows:

R --slave --no-save --args FCS_output.txt output_file NetNC_home_directory < NetNC_home_directory/mixtureModelling/NetNCmixmodel.R

The 'NetNC_home_directory' is the directory where NetNC_v2pt2.pl and FCS.pl are located. 

NetNCmixmodel.R runs mixture modelling on the NFCS scores written by FCS.pl
and seeks to identify a NFCS threshold to define 'noise' nodes. If a unimodal 
model is returned, a tiny Gaussian noise component is added to NFCS=0 values 
and mixture modelling is rerun (the .log file details which analyses were run). 
The output file with suffix '.coherentNodes' includes nodes scoring above the 
NFCS threshold value determined by this procedure (please the manuscript methods
section for further details). Of course functions in netNCmix.R may also be 
called directly.



##
## TESTING THE INSTALLATION AND EXAMPLE OUTPUT
##

The directory 'test/' includes data to enable you to test that NetNC
is working correctly. 

1. The tests may be run as follows (from the netNC/ directory):

a) Network Neighbourhood Clustering (NNC)

NetNC_v2pt2.pl -i test/test_genelist.txt -n test/network/test_net.txt -o test/output/NNConly/NNCz10 -z 10

This usually runs in less than one minute.

Please note: the above command runs 10 resamples (-z 10); however, at 
least 100 resamples are recommended for general use (more samples enable
greater precision for estimation of pFDR).

b) Node-centric Analysis mode:

NetNC_v2pt2.pl -n test/network/test_net.txt -i test/test_genelist.txt -o test/output/NodeCentric/NodeCent_z10 -M -p test/output/NNConly/NNCz10.BG.nlPonly.txt -z 10

Typically takes less than three minutes to complete.

This command is also an example of using the background p-value distribution generated
in test a) above, with the -p option; however it is not neccessary to do so for 
Node-centric analysis. If using the -p option it is also important to specify the 
number of resamples taken (-z  option) so that pFDR estimation is correct; if -z is not 
specified, the default number of resamples (100) is assumed and a warning message is 
issued.

c) Functional Target Identification mode:
NetNC_v2pt2.pl -n test/network/test_net.txt -i test/test_genelist.txt -o test/output/FTI/FTIz10 -z 10 -F

Runtime is around 15 minutes. This mode calculates NNC and runs the iterative minimum
cut, applying pFDR and network density thresholds optimised for the Functional Target 
Identification (FTI) task.

d) Pathway Identification mode and node-centric analysis:
NetNC_v2pt2.pl -n test/network/test_net.txt -i test/test_genelist.txt -o test/output/PID/PID_NodeCent_z10 -z10 -E -M -l test/test_background_genelist.txt

Usually takes about 15 minutes to compute. This mode calculates NNC and runs the iterative
minimum cut, applying pFDR and network density thresholds optimised for the Pathway 
Identification (PID) task. This command also includes an example of using the -l option
to specify a background node list for the resampling procedure (e.g. detected genes from a 
microarray experiment). 
Please note that to ensure correct estimation of pFDR the node list given to NetNC with 
the -i option should be a subset of the background list provided in the -l option; 
At present, no warning is given by NetNC if the -i (input) node list is not a subset of 
the -l (background) node list.

e) Iterative minimum cut in standalone mode:

python mincut/itercut.py -i test/output/NNConly/NNCz10.FDRthresholded_pairs.txt -o test/output/mincutOnly/NNCz10_FDR0pt1_mincutThresh0pt1.txt -t 0.1

Runtime is usually 10 minutes. This command depends on test a) above having
successfully completed.

2. Example output for the tests a) through e), above, are given in subdirectories
of test/exampleOutput. Please note that there might be differences between the
example output and the output of the tests below - due to variation introduced by
the resampling procedure. Subsections below summarise files provided in the example
output directories:

a) test/exampleOutput/NNConly

NNCz10.FG.id1_id2_nlp.txt - The Network Neighbourhood Clustering p-values for the
			    subnetwork induced by test_genelist.txt in test_net.txt
			    Format: nodeID1 nodeID2 -log(p-value)

NNCz10_10_resamples.txt  -  Resampled nodes taken for empirical estimation of pFDR.
			    Each resample begins with the following header:
			    '### [resample number]'
			    Where [resample number] is an integer.
			    The -z 10 option will produce 10 headers (numbered from
			    0 to 9). The list of resampled nodes is given
			    below each header.
			     
NNCz10.BG.nlPonly.txt   -   List of -log p-values for the networks induced by
			    all resampling from test_net.txt, which provides
			    the empirical null distribution for pFDR estimation.
			    (No headers are given in this file).

NNCz10.minSignif_nlP.log -  Reports the NNC -log(p-value) estimated to match
			    the pFDR threshold. Also records the number of resamples 
			    taken for pFDR estimation.

NNCz10.FDRthresholded_pairs.txt - Network of edges that meet the pFDR threshold (default 
				  0.1). Format: nodeID1 nodeID2 -log(p-value)

NNCz10.FDR.txt		 -  The mapping between NNC score and pFDR.
			    Format: NNC_score pFDR
		

b) test/exampleOutput/NodeCentric

NodeCent_z10.FG.id1_id2_nlp.txt - The Network Neighbourhood Clustering p-values for the
			    	  subnetwork induced by test_genelist.txt in test_net.txt
			    	  Format: nodeID1 nodeID2 -log(p-value)

NodeCent_z10.minSignif_nlP.log - Reports the NNC -log(p-value) estimated to 
				 match the pFDR threshold. Also records the number of 
			    	 resamples taken for pFDR estimation.

NodeCent_z10.FDRthresholded_pairs.txt - Network of edges that meet the pFDR threshold 
					(default 0.1). 
					Format: nodeID1 nodeID2 -log(p-value)

NodeCent_z10.FDR.txt     	-  The mapping between NNC score and pFDR.
			    	   Format: NNC_score pFDR
				   
NodeCent_z10.FCS.txt		- Node Functional Coherence Score (NFCS) values
				  calculated by FCS.pl.
				  Format: Node edgePvalueSum degree NFCS
				  (NFCS is the sum of -log(p-values) normalised by degree)
				  
NodeCent_z10_NFCS-mixturemodel.txt - Describes the Gaussian mixture model fitted to the 
				     NFCS distribution using Expectation-Maximisation
				     and model selection with Bayesian Information 
				     Criterion regularisation.
				    
NodeCent_z10_NFCS-mixturemodel.log - Indicates the NFCS threshold determined from analysis
				     of the Gaussian mixture modelling, and summarises the
				     steps taken in determining the threshold value.
				     
NodeCent_z10_NFCS-mixturemodel.coherentNodes - List of nodes that pass the NFCS threshold
					       and so classed as functionally coherent.
					       
NodeCent_z10_NFCS-mixturemodel_coherentNet.txt - Network where nodes are coherent in 
						the mixture modelling analysis and edges 
						have bonferroni-corrected -log(p)<=0.05					       
					       

c) test/exampleOutput/FTI/
					       
FTIz10.FG.id1_id2_nlp.txt - The Network Neighbourhood Clustering p-values for the
			    subnetwork induced by test_genelist.txt in test_net.txt
			    Format: nodeID1 nodeID2 -log(p-value)

FTIz10_10_resamples.txt  -  Resampled nodes taken for estimation of pFDR.
			    Each resample begins with the header:
			    '### [resample number]'
			    Where [resample number] is an integer.
			    The -z 10 option will produce 10 headers (numbered from
			    0 to 9). The list of resampled nodes is given
			    below each header.
			     
FTIz10.BG.nlPonly.txt   -   List of -log p-values for the networks induced by
			    all resampling from test_net.txt, which provides
			    the null distribution for empirical estimation of pFDR. 
			    There are no headers in this file.

FTIz10.minSignif_nlP.log -  Reports the NNC -log(p-value) estimated to match
			    the pFDR threshold. Also records the number of resamples 
			    taken for pFDR estimation.

FTIz10.FDRthresholded_pairs.txt - Network of edges that meet the pFDR threshold 
				  (default 0.1). Format: nodeID1 nodeID2 -log(p-value)

FTIz10.FDR.txt		 -  The mapping between NNC score and pFDR.
			    Format: NNC_score pFDR
			    
FTIz10.FDRthresholded_mincutDensity0pt306.txt - Edges passing the pFDR and minimum cut
						density thresholds for Functional Target 
						Identification.	
						Format: nodeID1 nodeID2 -log(p-value)		    

FTIz10.FDRthresholded_mincutDensity0pt306_noDoublets.txt - Final NetNC output corresponding
							   to the data used in benchmarking.
							   Edges pass the pFDR and minimum
							   cut density thresholds for FTI
							   and any two-node components are
							   deleted.

d) test/exampleOutput/PID/

PID_NodeCent_z10.FG.id1_id2_nlp.txt - Network Neighbourhood Clustering p-values for the
			   	      subnetwork induced by test_genelist.txt in 
				      test_net.txt.
			    	      Format: nodeID1 nodeID2 -log(p-value)
				  
PID_NodeCent_z10_10_resamples.txt  -  Resampled nodes taken for estimation of pFDR.
			    	      Each resample begins with the header:
			    	      '### [resample number]'
			    	      Where [resample number] is an integer.
			    	      The -z 10 option will produce 10 headers (numbered 
				      from 0 to 9). The list of resampled nodes is given
	  			      below each header.
			     
PID_NodeCent_z10.BG.nlPonly.txt   -   List of -log p-values for the networks induced by
			    	      all resampling from test_net.txt, which provides
			    	      the null distribution for empirical estimation of
				      pFDR. There are no headers in this file.
			    
PID_NodeCent_z10.minSignif_nlP.log - Reports the -log(p-value) estimated to match the 
				     NNC pFDR threshold. Also records the number 
				     of resamples taken for pFDR estimation.

PID_NodeCent_z10.FDRthresholded_pairs.txt - Network of edges that meet the pFDR threshold 
					    (default 0.1). 
					    Format: nodeID1 nodeID2 -log(p-value)

PID_NodeCent_z10.FDR.txt     	-  The mapping between NNC score and pFDR.
			    	   Format: NNC_score pFDR
				   
PID_NodeCent_z10.FCS.txt	- Node Functional Coherence Score (NFCS) values
				  calculated by FCS.pl.
				  Format: Node edgePvalueSum degree NFCS
				  (NFCS is the sum of -log(p-values) normalised by degree)
				  
PID_NodeCent_z10_NFCS-mixturemodel.txt - Describes the Gaussian mixture model fitted to the 
				         NFCS distribution using Expectation-Maximisation
				         and model selection with Bayesian Information 
				         Criterion regularisation.
				    
PID_NodeCent_z10_NFCS-mixturemodel.log - Indicates the NFCS threshold determined from 
					 analysis of the Gaussian mixture modelling, and 
					 summarises the steps taken in determining the 
					 threshold value.
				     
PID_NodeCent_z10_NFCS-mixturemodel.coherentNodes - List of nodes that pass the 
						   NFCS threshold and so classed as 
						   functionally coherent.
						   
PID_NodeCent_z10_NFCS-mixturemodel_coherentNet.txt - Network where nodes are coherent 
						     in the mixture modelling analysis 
						     and edges have bonferroni-corrected 
						     -log(p)<=0.05					       
							   
					       
PID_NodeCent_z10.FDRthresholded_mincutDensity0pt5044.txt - Edges passing the pFDR and 
							   minimum cut density thresholds 
							   for Pathway Identification.	
							   File format: 
							   nodeID1 nodeID2 -log(p-value)		    

PID_NodeCent_z10.FDRthresholded_mincutDensity0pt5044_noDoublets.txt - Edges passing the pFDR and
	                                                              minimum cut density thresholds
            	                                                      for Pathway Identification;
								      two node components are delete.
                                                          	      File format: 
	                                                              nodeID1 nodeID2 -log(p-value)

e) test/exampleOutput/mincutOnly/ 					       		

NNCz10_FDR0pt1_mincutThresh0pt1.txt - The network of edges passing the minimum cut density
				      stopping criterion (0.1 in this example). These edges
				      also pass the NNC pFDR threshold of 0.1 which was
				      applied by NetNC prior to running the minimum cut 
				      algorithm in stand-alone mode.
				      Format: nodeID1 nodeID2 -log(p-value)
				      (-log(p-value) previously calculated by NetNC).


##
## Contributors
##

The NetNC software distribution was developed by Ian Overton, Jeremy Owen (iterative 
minimum cut) and Alex Lubbock (Gaussian Mixture Modelling).


##
## Citing NetNC
##

If you use any of the code in this NetNC software distribution please cite:
Overton IM, Sims A, Owen JA, Heale B, Ford M, Lubbock ALR, Pairo-Castineira E, Essafi E (2018).
"Functionally Coherent Transcription Factor Target Networks Illuminate Control of Epithelial 
Remodelling and Oncogenic Notch". BioRxiv 2018 https://doi.org/10.1101/455709.
