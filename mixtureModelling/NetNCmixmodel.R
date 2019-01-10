#    Copyright (C) 2016 The University of Edinburgh
#
#
#    NetNCmixmodel.R is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    NetNCmixmodel.R is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with NetNCmixmodel.R.  If not, see <http://www.gnu.org/licenses/>.
#

#
# This script performs Mixture modelling on Node Functional Coherence Scores (NFCS)
# in order to infer functionally coherent nodes and estimate the proportion of
# 'noise' in the input nodes
#

# Authors: Ian Overton (i.overton@qub.ac.uk) and Alex Lubbock (code@lubbock.com)
# v1.1

# read command-line arguments
co <- commandArgs()
nfcs_file = co[[5]]
out = co[[6]]
path = co[[7]]
netnc_mix =  paste(path, "/mixtureModelling/netNCmix.R", sep="")

# load mixture modelling functions
source(netnc_mix)

# Read in the data with a list that reflects the structure of the FCS.pl output
l = list(id="", nlp=0, deg=0, normP=0)
nfcs_data = scan(file=nfcs_file, what=l)

# run mixture modelling on the nfcs scores
sel_nfcs = select_model(nfcs_data$normP)

mixmodel_file = paste(out, ".txt", sep="") 
# write the model to file at the results
header = paste("Gaussian mixture model on NFCS scores from ", nfcs_file, ":", sep="")
write(header, mixmodel_file)
# print the model
x = sel_nfcs
mixData = paste(sprintf("Fit:\n  Component(s):\t%d\t\tData points: %d\n  Log-likelihood: %.3f\t(AIC,BIC):   (%.2f,%.2f)\n",x$k,x$n,x$loglik,x$aic,x$bic), sprintf("  Initialisation: %s\tModel type:  %s\n",x$initialisation,x$model.type), "Components:\n", sep="")
write(mixData, mixmodel_file, append="TRUE")
for(i in 1:x$k) write(sprintf("  %d. Proportion: %.2f\tMean: %.2f\tStd Dev: %.2f",i,x$mix.props[i],x$means[i],x$std.dev[i]), mixmodel_file, append="TRUE")

# identify functionally coherent nodes
nfcs_coherent = coherent_nodes(nfcs_data, sel_nfcs, out)

# check to see if a unimodal model was returned and if so add Gaussian noise to the 0 values, and then rerun the mixture modelling

if (nfcs_coherent[1] == "NA") { # unimodal model returned

# make an object with gaussian noise added to 0 values
	nfcs_gn = nfcs_data
	nfcs_gn$normP[nfcs_gn$normP==0] = add_gauss_0(nfcs_gn$normP)

# mixture modelling
	sel_nfcs_GN = select_model(nfcs_gn$normP)
# output the results
	header = paste ("Gaussian mixture model on NFCS scores with Gaussian noise at NFCS==0 for ", nfcs_file, ":", sep="")
	write(header, mixmodel_file)
	x = sel_nfcs_GN
	mixData = paste(sprintf("Fit:\n  Component(s):\t%d\t\tData points: %d\n  Log-likelihood: %.3f\t(AIC,BIC):   (%.2f,%.2f)\n",x$k,x$n,x$loglik,x$aic,x$bic), sprintf("  Initialisation: %s\tModel type:  %s\n",x$initialisation,x$model.type), "Components:\n", sep="")
write(mixData, mixmodel_file, append="TRUE")
for(i in 1:x$k) write(sprintf("  %d. Proportion: %.2f\tMean: %.2f\tStd Dev: %.2f",i,x$mix.props[i],x$means[i],x$std.dev[i]), mixmodel_file, append="TRUE")
	
# identify functionally coherent nodes
	gn_data_out = paste(out, "_GN", sep="") 
	nfcs_gn_coherent = coherent_nodes(nfcs_gn, sel_nfcs_GN, gn_data_out)
	
# check to see if unimodal model returned
	if (nfcs_gn_coherent[1] == "NA") { 
		print("WARNING: Only unimodal models were returned from mixture modelling on NFCS scores")
		
	}
}

# save an image of the session (commented out)
# R_session_out = paste(out, ".RData", sep="") 
# save.image(file = R_session_out)
