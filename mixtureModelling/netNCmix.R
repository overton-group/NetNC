#    Copyright (C) 2012 The University of Edinburgh
#
#
#    netNCmix.R is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    netNCmix.R is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with netNCmix.R.  If not, see <http://www.gnu.org/licenses/>.
#
# Authors: Alex Lubbock (code@lubbock.com) and Ian Overton (i.overton@qub.ac.uk)


# Mixture modelling and clustering
##
## Univariate mixture of Gaussians; k-means++; quantile discretisation

##' Implementation of k-means++ - k-means with improved seeding
##'
##' @title K-means++ clustering algorithm
##' @param x A numeric vector of data to cluster
##' @param k The number of clusters (must be >=1)
##' @param iter.max Maximum number of iteration for k-means algorithm
##' @param algorithm Choice of k-means algorithm
##' @return A list of centres, and data points belonging to those centres
##' @seealso \R's \code{kmeans} for k-means without the improved seeding
##' @references Arthur, D. & Vassilvitskii, S. k-means++: the advantages of careful seeding. 1027-1035 (2007).
##' @keywords clustering
##' @export
##' @examples
##' mix <- c(rnorm(100,mean=0,sd=2),rnorm(50,mean=5,sd=1))
##' cl <- kmeansplusplus(mix,2)
##' plot(mix,col=tma.lcols[cl$cluster],pch=cl$cluster,main="K-means++ clustering",ylab="Value")
kmeansplusplus <- function(x,k,iter.max = 10, algorithm=c("Hartigan-Wong","Lloyd","Forgy","MacQueen")) {
  if(k<1) stop("k must be >=1")
  if(is.null(dim(x))) dim(x) <- c(length(x),1)
  if(k>nrow(x)) stop("k must be >n")
  centres <- matrix(0,k,ncol(x))
  centres[1,] <- x[sample(1:nrow(x),1),]
  if(k>1) {
    for(cn in 2:k) {
      dist.nearest.centroid <- apply(x,1,function(xeach){min(apply(t(t(centres)-xeach)^2,1,sum))})
      # Avoid divide by zero
      dist.nearest.centroid[dist.nearest.centroid < .Machine$double.eps] <- .Machine$double.eps
      centres[cn,] <- x[sample(1:nrow(x),1,prob=dist.nearest.centroid/sum(dist.nearest.centroid)),]
    }
  } else {
    centres <- 1
  }
  
  return(kmeans(x,centres,iter.max=iter.max,algorithm=algorithm))
}

##' Discretises a \code{tma} or vector by either quantiles or using a mixture of
##' Gaussians model.
##'
##' Missing data values are removed automatically and silently.
##' 
##' Discretisation means to convert a range of continous values into a fixed
##' number of values. An example is using tertiles (that is, three quantiles):
##' The bottom third of the continuous values will be converted to 1, the next
##' third will be converted to 2, and the top third will be converted to 3.
##'
##' If a whole \code{tma} is supplied rather than a numeric vector, discretisation
##' is done on a per-protein basis.
##'
##' @title Discretise a \code{tma} object or numeric vector
##' @param dataset A \code{tma} object, or a vector of numeric values
##' @param method Either 'quantile' for quantiles, or 'mog' for mixture of Gaussians
##' @param quantiles The number of quantiles to use; only used for 'quantile' method
##' @param model A mixture model to use with 'mog', or NULL to learn a mixture model
##' from the data
##' @param \dots Arguments to \code{\link{EM.findk}}; only used for 'mog' method
##' @return The discretised \code{tma} or vector
##' @keywords clustering mixturemodels
##' @export
##' @examples
##' mix <- c(rnorm(100,mean=0,sd=2),rnorm(50,mean=5,sd=1))
##' cl <- discretise(mix,quantiles=2)
##' plot(mix,col=tma.lcols[cl],pch=cl,main="Quantiles",ylab="Value")
##'
discretise <- function(dataset,method="quantile",quantiles=4,model=NULL,...) {
  if(class(dataset)=="tma") dataset <- dataset$score.matrix
  if(method!="mog" && !is.null(model)) stop("Cannot use a model unless method is 'mog'")
  if(!is.null(model) && class(model)!="mixmodel") stop("Model is not a valid 'mixmodel'")
  drop.dims <- FALSE
  if(is.null(dim(dataset))) {
    nm <- names(dataset)
    drop.dims <- TRUE
    dim(dataset) <- c(1,length(dataset))
    colnames(dataset) <- nm
  }
  if(any(is.na(dataset))) dataset <- omit.na.patients(dataset)
  res <- matrix(NA,nrow=nrow(dataset),ncol=ncol(dataset),dimnames=dimnames(dataset))
  for(i in 1:nrow(dataset)) {
    data.vec <- as.numeric(dataset[i,])
    switch(method,
           quantile={
             res[i,] <- quantcut(data.vec,q=seq(0,1,1/quantiles))
             levels(res[i,]) <- 1:quantiles
           },
           mog={
             if(is.null(model)) model <- EM.findk(data.vec,...)
             if(model$k==1) {
               res[i,] <- rep(1,ncol(dataset)) ##unimodal (set all to 1)
             } else {
               res[i,] <- mog.density(data.vec,model)$membership
             }
           },
           stop("recognised values for method are 'quantile' or 'mog' (mixture of Gaussians)")
           )
  }
  
  ##'Drop unimodal markers
  if(method=="mog") {
    unimodal <- apply(res==1,1,all)
    if(any(unimodal)) cat("Dropping unimodal markers:",paste(names(unimodal)[unimodal],collapse=","),"\n")
    res <- res[!unimodal,]
    attr(res,"unimodal") <- ifelse(any(unimodal),which(unimodal),NA)
  }
  attr(res,"na.action") <- attr(dataset,"na.action")

  if(drop.dims && !is.null(nrow(res)) && nrow(res)==1) {
    res <- c(res)
    names(res) <- colnames(dataset)
  }
  if(!is.null(attr(dataset,"protein")))
    attr(res,"protein") <- attr(dataset,"protein")
  attr(res,"discretisation") <- ifelse(method=="mog","Mixture of Gaussians","Quantiles")
  
  return(res)
}

##' Log likelihood of data given a mixture of Gaussian model
##'
##' @title Get the log likelihood of data given a mixture model
##' @param data.vec A numeric vector
##' @param mixmodel An object of class mixmodel
##' @return The log likelihood of the data given the model
##' @keywords internal
log.likelihood <- function(data.vec,mixmodel) {
  data.vec <- as.numeric(data.vec)
  probs <- repmat(t(mixmodel$mix.props),length(data.vec),1) * sapply(1:mixmodel$k,function(i) { dnorm(data.vec,mixmodel$means[i],mixmodel$std.dev[i]) })
  sum(log(rowSums(probs)))
}

##' Initialise the expectation maximisation (EM) algorithm
##'
##' @title Initialise EM
##' @param data.vec A numeric vector
##' @param k The number of mixture components to use in the model
##' @param model.type 'E' contrains variances to be equal across components, 'V' estimates variances independently
##' @param init.method The method to initialise the parameters for EM - quantiles or K-means++
##' @return An object of class "mixmodel" ready for updating
##' @keywords internal
init.em <- function(data.vec,k,model.type="V",init.method="kmeans++") {
  if(init.method=="quantiles")
    q <- discretise(data.vec,method="quantile",quantiles=k)
  if(init.method=="kmeans++")
    q <- kmeansplusplus(data.vec,k)$cluster
  mix.props <- as.numeric(repmat(1/k,1,k))
  means <- sapply(unique(q),function(i) mean(data.vec[i==q]))
  std.dev <- if(model.type=="V") sapply(unique(q),function(i) sd(data.vec[i==q]))
             else if(model.type=="E") rep(sd(data.vec),k)
             else stop("Unknown model type: ",model.type)
  mixmodel <- list(initialisation=init.method,k=k,means=means,std.dev=std.dev,mix.props=mix.props,model.type=model.type,n=length(data.vec))
  class(mixmodel) <- "mixmodel"
  return(mixmodel)
}

##' Compute membership probabilities for each data point for each cluster
##'
##' @title 'E' step of EM algorithm
##' @param data.vec A numeric vector
##' @param mixmodel A mixture model to update
##' @return The mixture model with updated probabilities
##' @keywords internal
e.step <- function(data.vec,mixmodel) {
  if(class(mixmodel)!="mixmodel") stop("Mixture parameters should be of mixmodel class")
  data.vec <- as.numeric(data.vec)
  probs <- repmat(t(mixmodel$mix.props),length(data.vec),1) * sapply(1:mixmodel$k,function(i) { dnorm(data.vec,mixmodel$means[i],mixmodel$std.dev[i]) })
  probs/rowSums(probs)
}

##' Update parameters using maximum likelihood
##'
##' @title 'M' step of EM algorithm
##' @param data.vec A numeric vector
##' @param mem.probs The mixture membership probabilities of each of the values in data.vec
##' @param mixmodel A mixture model to update
##' @return The mixture model with updated parameters
##' @keywords internal
m.step <- function(data.vec,mem.probs,mixmodel) {
  data.vec <- as.numeric(data.vec)
  mixmodel$mix.props <- colSums(mem.probs)/length(data.vec)
  mixmodel$means <- colSums(data.vec*mem.probs)/colSums(mem.probs)
  data.vec.mean.centred <- data.vec - repmat(t(mixmodel$means),length(data.vec),1)
  mixmodel$std.dev <- if(mixmodel$model.type=="V") sqrt(colSums(mem.probs*(data.vec.mean.centred^2))/colSums(mem.probs)) else if(mixmodel$model.type=="E") rep(sqrt(sum(mem.probs*(data.vec.mean.centred^2))/sum(mem.probs)),mixmodel$k) else stop("Unknown model type: ",mixmodel$model.type)
  mixmodel
}

##' The expectation maximisation (EM) algorithm creates a mixture of Gaussians model
##' from a supplied vector of data, for example a set of protein scores.
##'
##' A mixture of Gaussian model consists of a summation of multiple Gaussian (normal)
##' distributions. It can often be a good approximation for multi-modal distributions
##' (loosely, where a protein's density looks like multiple Gaussians have been
##' superimposed).
##'
##' This function requires the number of Gaussians used, \code{k}, to specified in
##' advance. To find the best \code{k} from the data, use \code{\link{EM.findk}}.
##'
##' @title Expectation maximisation (EM) algorithm for mixture of Gaussians models
##' @param data.vec A numeric vector of data to use
##' @param k The number of Gaussians to use in the mixture
##' @param model.type Constrain the variances of the Gaussians to be the same ('E')
##' or not the same ('V')
##' @param it.max Maximum number of iterations to attempt convergence before aborting
##' @param lhood.tol Maximum change in likelihood between iterations to consider as convergence
##' @param init.method Initialisation method. Currently supported methods are 'kmeans++'
##' or 'quantiles'.
##' @param max.restarts Maximum number of restarts to perform if model doesn't fit successfully
##' @param always.restart If TRUE, keep refitting models until max.restarts is reached,
##' selecting the best model by likelihood (slow, but may avoid local optima). If FALSE,
##' the process terminates when the first model is fitted with the required parameter
##' constraints (min.min.diff and min.sd) and likelihood tolerance (lhood.tol)
##' @param min.mean.diff Try to fit a model with the component means at least as far apart
##' as this parameter specifies, restarting the fitting process if necessary
##' @param min.sd Minimum standard deviation of each component. Useful to avoid very narrow
##' Gaussians, which are often noise spikes. Tries to fit a model meeting this criterion,
##' restarting the fitting process if necessary
##' @param suppress.warnings Suppresses warnings about the model not meeting the required
##' tolerance or parameter constraints after max.restarts is reached
##' @return An object of class \code{mixmodel} giving the mixture model's parameter set
##' @seealso \code{\link{EM.findk}} to build a mixture model using AIC or BIC, where
##' the number of Gaussians \code{k} is not known in advance\cr
##' \code{\link{scores}} to obtain a numeric vector of protein scores
##' @keywords mixturemodels
##' @export
##' @examples
##' mix <- c(rnorm(100,mean=0,sd=2),rnorm(50,mean=5,sd=1))
##' mog <- EM(mix,2)
##' print(mog)
##' cl <- mog.density(mix,mog)
##' plot(mix,col=tma.lcols[cl$membership],pch=cl$membership,
##'      main="Mixture of Gaussians using EM",ylab="Value")
##' plot(mog,mix,min.x=min(mix))
EM <- function(data.vec,k,model.type="V",it.max=1e3,lhood.tol=1e-4,init.method="kmeans++",max.restarts=5,always.restart=FALSE,min.mean.diff=0.01,min.sd=0.01,suppress.warnings=FALSE) {
  if(length(unique(data.vec))<k) stop("Need at least k=",k,"distinct data points")
  ##' Remove missing values
  if(any(is.na(data.vec))) {
    marker.name <- attr(data.vec,"protein")
    data.vec <- na.omit(data.vec)
    attr(data.vec,"protein") <- marker.name
  }

  mixmodel <- NULL
  n.restart <- 0
  n.unconverged <- 0
  while(n.restart==0 || ((always.restart || any(diff(sort(mixmodel.this$means))<min.mean.diff,na.rm=T) || any(mixmodel.this$std.dev<min.sd,na.rm=T)) && n.restart<max.restarts)) {
    ##' Initialise
    mixmodel.this <- init.em(data.vec,k,model.type=model.type,init.method=init.method)
    
    lhood.last <- log.likelihood(data.vec,mixmodel.this)
  
    i <- 1
    repeat {
      ##Read off the new membership probabilities
      mem.probs <- e.step(data.vec,mixmodel.this)
      ##Maxmisation step - adapt parameters
      mixmodel.this <- m.step(data.vec,mem.probs,mixmodel.this)
      ##Compute the likelihood
      lhood.this <- log.likelihood(data.vec,mixmodel.this)
      ##Stop if unavailable to compute likelihood
      if(is.na(lhood.this)) break
      ##Stop on convergence
      if(abs((lhood.last-lhood.this)/lhood.last)<lhood.tol) break
      ## Stop if maximum number of iterations reached
      if(i>=it.max) {
        n.unconverged <- n.unconverged+1
        break
      }
      lhood.last <- lhood.this
      i <- i+1
    }
    if(is.null(mixmodel) || is.null(mixmodel$loglik) || is.na(mixmodel$loglik) || mixmodel$loglik>mixmodel.this$loglik) mixmodel <- mixmodel.this
    n.restart <- n.restart+1
  }

  ##'Sort the MoG modes in ascending order
  mode.order <- order(mixmodel$means)
  mixmodel$means <- mixmodel$means[mode.order]
  mixmodel$std.dev <- mixmodel$std.dev[mode.order]
  mixmodel$mix.props <- mixmodel$mix.props[mode.order]

  if(!suppress.warnings) {
    if(any(is.na(mixmodel$mix.props))) warning("k=",k,": ",sum(is.na(mixmodel$mix.props))," empty model component(s) (decrease k or increase max.restarts)")
    if(n.unconverged==max.restarts) warning("k=",k,": Model did not converge (increase it.max, lhood.tol or max.restarts)")
    else if(any(diff(mixmodel$means)<min.mean.diff,na.rm=T)) warning("k=",k,": Some component means closer than requested (decrease k or increase max.restarts)")
    else if(any(mixmodel$std.dev<min.sd,na.rm=T)) warning("k=",k,": Some component standard deviations lower than requested (decrease k or increase max.restarts)")
  }
  
  mixmodel$loglik <- lhood.this
  mixmodel$bic <- tmap.bic(lhood.this,mixmodel,length(data.vec))
  mixmodel$aic <- tmap.aic(lhood.this,mixmodel)
  
  attr(mixmodel,"protein") <- attr(data.vec,"protein")
  class(mixmodel) <- "mixmodel"
  return(mixmodel)
}

tmap.bic <- function(loglik,mixmodel,n) {
  k <- ifelse(mixmodel$model.type=="V",2*mixmodel$k,mixmodel$k+1)
  ##'The mixture components count as k-1 free parameters
  if(mixmodel$k>1) k <- k+mixmodel$k-1
  -(2*loglik)+(k*log(n))
}

tmap.aic <- function(loglik,mixmodel) {
  k <- ifelse(mixmodel$model.type=="V",2*mixmodel$k,mixmodel$k+1)
  ##'The mixture components count as k-1 free parameters
  if(mixmodel$k>1) k <- k+mixmodel$k-1
  -(2*loglik)+(2*k)
}

##' Returns the probably density, or maximum likelihood Gaussian,
##' at each of the points in data.vec using the supplied
##' Mixture of Gaussian model
##'
##' Returns the pointwise Mixture of Gaussian density across the
##' supplied data points, and the pointwise maximum likelihood
##' cluster (i.e. which Gaussian does each of the supplied data
##' points most likely come from under the model).
##' 
##' @title Probability density of a Mixture of Gaussian model
##' @param data.vec Numeric vector of data points
##' @param mixmodel Mixture of Gaussian model of class \code{mixmodel}
##' @param sort.data Sort data.vec in ascending order before estimating
##' if \code{TRUE}
##' @return A list containing the following entries:
##' \itemize{
##'   \item x - contents of data.vec (optionally sorted)
##'   \item y - density at each point in x
##'   \item prob - probability of each data point belonging to each Gaussian
##'   \item membership - which Gaussian each data point is most likely
##'   to have come from under the model
##' }
##' @seealso \code{\link{EM}} to fit a mixture model with a specific number of Gaussians \code{k}\cr
##' \code{\link{EM.findk}} to fit a mixture model and find \code{k}
##' using AIC or BIC.
##' @keywords mixturemodels
##' @export
##' @examples
##' mix <- c(rnorm(100,mean=0,sd=2),rnorm(50,mean=5,sd=1))
##' mog <- EM(mix,2)
##' cl <- mog.density(mix,mog)
##' # See EM function for plotting examples
mog.density <- function(data.vec,mixmodel,sort.data=FALSE) {
  if(any(is.na(data.vec))) stop("Please remove NA values from data.vec (use na.omit)")
  if(sort.data) data.vec <- sort(data.vec)
  y <- matrix(NA,mixmodel$k,length(data.vec))
  for(i in 1:mixmodel$k) y[i,] <- mixmodel$mix.props[i]*dnorm(as.numeric(data.vec),mean=mixmodel$means[i],sd=mixmodel$std.dev[i])
  prob <- apply(y,2,function(a){a/sum(a)})
  if(is.null(dim(prob))) dim(prob) <- c(1,length(prob))
  list(x=data.vec,y=apply(y,2,sum),prob=t(prob),membership=apply(y,2,which.max))
}

##' Creates a mixture of Gaussians model with automatic selection of the number of
##' Gaussian distributions.
##' 
##' Mixture models are fitted by expectation maximisation (EM) and the number of
##' Gaussians is selected by Akaike Information Criterion (AIC) or Bayesian Information
##' Criterion (BIC). BIC penalises more heavily than AIC for extra parameters,
##' so BIC models will contain fewer Gaussians on average.
##'
##' Missing data values are removed automatically and silently.
##'
##' @title Expectation maximisation (EM) for mixture of Gaussians models with automatic
##' selection of the number of Gaussians
##' @param data.vec Numeric vector of data points
##' @param num.gaussians Range of number of component Gaussians to select from.
##' Default is 1..9.
##' @param model.types Model types to search. 'E' - constrain Gaussians to have
##' equal variance to reduce number of parameters, 'V' - estimate Gaussian
##' variances separately. Default is to try both.
##' @param select.by Model likelihood/number of parameter tradeoff to apply.
##' Options are 'aic' or 'bic' - default is 'bic'.
##' @param suppress.warnings  Suppresses warnings from \code{\link{EM}} about the final
##' model not meeting the required constraints (likehood tolerance, parameter constraints)
##' @param ... Other arguments to \code{\link{EM}}
##' @return An object of class mixmodel giving the mixture model's parameter set
##' @seealso \code{\link{EM}} to build a mixture model where the number of Gaussians
##' \code{k} is known in advance. Plotting examples for mixture models are under that
##' function\cr
##' \code{\link{scores}} can be used to  obtain a numeric vector of protein scores.
##' @keywords mixturemodels
##' @export
##' @examples
##' mix <- c(rnorm(100,mean=0,sd=2),rnorm(50,mean=5,sd=1))
##' mog <- EM.findk(mix,num.gaussians=1:4,model.types="V")
##' print(mog)
##' # See the EM function for plotting examples
EM.findk <- function(data.vec,num.gaussians=1:9,model.types=c("E","V"),select.by="bic",suppress.warnings=TRUE,...) {
  stopifnot(select.by %in% c("aic","bic"))
  model.bic <- list()
  mixmodels <- list()

  best.k <- NA
  best.type <- NA
  best.bic <- Inf
  
  for(k in num.gaussians) {
    for(model.type in model.types) {
      mixid <- paste(k,model.type,sep="")
      mixmodels[[mixid]] <- EM(data.vec,k=k,model.type=model.type,suppress.warnings=suppress.warnings,...)
      model.bic[[mixid]] <- mixmodels[[mixid]][[select.by]]
      if(!is.na(model.bic[[mixid]]) && model.bic[[mixid]]<best.bic) {
        best.bic <- model.bic[[mixid]]
        best.k <- k
        best.type <- model.type
      }
    }
  }
  mixmodels[[paste(best.k,best.type,sep="")]]
}

##' Group scores based on a supplied discretisation scheme
##'
##' Returns a list of scores in groups defined by the membership
##' vector.
##'
##' @title Group scores based on a supplied discretisation scheme.
##' @param marker.scores A numeric vector of scores
##' @param marker.membership A vector of membership
##' @return A list of grouping marker.scores by marker.membership
##' @keywords internal
group.scores <- function(marker.scores,marker.membership) {
  if(length(marker.scores)!=length(marker.membership)) stop("Membership and scores must be the same length")
  marker.grouped <- list()
  for(s in sort(unique(marker.membership))) {
    marker.grouped[[as.character(s)]] <- marker.scores[marker.membership==s]
  }
  marker.grouped
}

##' Function to plot mixture models, optionally overlaid with source data.
##'
##' @param x An object of class \code{mixmodel}
##' @param \dots Optional arguments:
##' \itemize{
##'  \item \code{data.vec}: A numeric vector to overlay on the mixture model plot,
##' e.g. the data the model was fitted with.
##'  \item \code{show.modes}: Plot vertical lines at the mode of each mixture component
##' if \code{TRUE} (default)
##'  \item \code{breaks}: The method for calculating breaks on the histogram plot (only
##' used when \code{data.vec} is supplied). Defaults to 'Sturges'. See \code{?histogram}
##' for options.
##'  \item \code{xmin}: The minimum x value on the plot (typically 0) (only used when
##' \code{data.vec} is supplied)
##' }
##' @return (Invisibly) The mixture model, which was supplied as an argument or created
##' @seealso \code{\link{EM.findk}} to fit a mixture model to data
##' @keywords mixturemodels
##' @method plot mixmodel
##' @export
##' @examples
##' mix <- c(rnorm(100,mean=0,sd=2),rnorm(50,mean=5,sd=1))
##' mog <- EM(mix,2)
##' plot(mog,data.vec=mix)
##'
##' # Or we could just plot a mixture model without overlaying
##' # the data
##' plot(mog)
plot.mixmodel <- function(x,...) {
  params <- list(...)
  if(!"data.vec" %in% names(params)) {
    xmin <- min(x$means)-3*x$std.dev[which.min(x$means)]
    xmax <- max(x$means)+3*x$std.dev[which.max(x$means)]
    mixdens <- mog.density(seq(xmin,xmax,length.out=1000),x)
    plot(mixdens, col=tma.lcols[10], type="l", lty=tma.ltys[2], lwd=3, xlab="Score", ylab="Density",main=attr(x,"protein"))
  } else {
    data.vec <- sort(params[["data.vec"]])
    breaks <- ifelse("breaks" %in% names(params),params$breaks,"Sturges")
    if(is.null(x)) x <- EM.findk(data.vec,model.types="V") ##only non-fixed variance models
    h.max <- max(hist(data.vec, breaks=breaks, plot=FALSE)$density)

    d <- variable.bw.kde(data.vec,output.domain=seq(min(data.vec,na.rm=T),max(data.vec,na.rm=T),length.out=512))
    mixdens <- mog.density(d$x,x)
    xmin <- ifelse("xmin" %in% names(params), params[["xmin"]], 0)
    if(is.null(xmin)) xmin <- min(d$x)
    hist(data.vec,freq=FALSE, breaks=breaks, xlim=c(xmin,max(d$x)), ylim=c(0,max(c(mixdens$y,h.max,d$y))), xlab="Score", ylab="Density",main=attr(x,"protein"))
    lines(d,lwd=3,col=tma.lcols[9],lty=tma.ltys[11])
    lines(mixdens, col=tma.lcols[10], type="l", lty=tma.ltys[2], lwd=3)
  }
  
  minor.tick(ny=1,nx=2,tick.ratio=1)
  minor.tick(ny=1,nx=10,tick.ratio=0.5)
  minor.tick(ny=5,nx=1,tick.ratio=0.5)

  ##plot modes
  if(!"show.modes" %in% names(params) || params[["show.modes"]]==TRUE) for(mode.num in 1:x$k) abline(v=x$means[mode.num],lty=tma.ltys[3],col=tma.lcols[mode.num],lwd=2)

  if(!"data.vec" %in% names(params))
    legend(x="topright",legend=c("Mixture Model"),lty=c(tma.ltys[2]),lwd=c(3),col=tma.lcols[c(10,1:x$k)])
  else
     legend(x="topright",legend=c("Scores (KDE)","Mixture Model"),lty=c(tma.ltys[11],tma.ltys[2]),lwd=c(2,3),col=tma.lcols[c(9,10,1:x$k)])     

  invisible(x)
}

##' Print information about a Gaussian mixture model and its parameters.
##'
##' @title Print a Gaussian mixture model
##' @param x Gaussian mixture model of class \code{mixmodel}
##' @param \dots Additional arguments (not used, kept for compatibility)
##' @seealso \code{\link{EM}} to fit a mixture model with a specific number of
##' Gaussians \code{k}\cr
##' \code{\link{EM.findk}} to fit a mixture model and find \code{k} using AIC or
##' BIC to avoid overfitting.
##' @keywords mixturemodels
##' @method print mixmodel
##' @export
##' @examples
##' mix <- c(rnorm(100,mean=0,sd=2),rnorm(50,mean=5,sd=1))
##' mog <- EM(mix,2)
##' print(mog)
print.mixmodel <- function(x,...) {
  if(!is.null(attr(x,"protein"))) {
    cat(sprintf("Gaussian mixture model for %s\n",attr(x,"protein")))
  } else {
    cat("Gaussian mixture model\n")
  }
  cat(sprintf("Fit:\n  Component(s):\t%d\t\tData points: %d\n  Log-likelihood: %.3f\t(AIC,BIC):   (%.2f,%.2f)\n",x$k,x$n,x$loglik,x$aic,x$bic))
  cat(sprintf("  Initialisation: %s\tModel type:  %s\n",x$initialisation,x$model.type))
  cat("Components:\n")
  for(i in 1:x$k) cat(sprintf("  %d. Proportion: %.2f\tMean: %.2f\tStd Dev: %.2f\n",i,x$mix.props[i],x$means[i],x$std.dev[i]))
}

##' Summarise a Gaussian mixture model, including its parameter set and formula.
##'
##' @title Summarise a Gaussian mixture model
##' @param object Gaussian mixture model of class \code{mixmodel}
##' @param \dots Additional arguments to \code{summary} (not used, kept for compatibility)
##' @return Object of class \code{summary.mixmodel} with summary attributes
##' @seealso \code{\link{EM}} to fit a mixture model with a specific number
##' of Gaussians \code{k}\cr
##' \code{link{EM.findk}} to fit a mixture model and find \code{k} using AIC
##' or BIC
##' @keywords mixturemodels
##' @method summary mixmodel
##' @export
##' @examples
##' mix <- c(rnorm(100,mean=0,sd=2),rnorm(50,mean=5,sd=1))
##' mog <- EM(mix,2)
##' summary(mog)
summary.mixmodel <- function(object,...) {
  smry <- object
  smry$formula <- sprintf("f(x) =")
  if(smry$k==1)
    smry$formula <- paste(smry$formula,sprintf("G(%.2f,%.2f)",smry$means,smry$std.dev))
  else {
    for(i in 1:smry$k) smry$formula <- paste(smry$formula,sprintf("%.2f*G(%.2f,%.2f) +",
                                                                      smry$mix.props[i],
                                                                      smry$means[i],
                                                                      smry$std.dev[i]))
    smry$formula <- substr(smry$formula,1,nchar(smry$formula)-2)
  }
  class(smry) <- "summary.mixmodel"
  smry
}


#
# Function to avoid local optima in model selection
#

select_model = function(data, bic=0) {
	for (i in 1:20) {
		if (bic) {
			new_mod  = EM.findk(data, model.type="V")
			if (bic > new_mod$bic) {
				model = new_mod
				bic = new_mod$bic
			}
	
		} else {
			model  = EM.findk(data, model.type="V")

		}
	}
	return(model)

}


#
# Add Gaussian noise to 0 values only
#

add_gauss_0 = function(data) {
	len = length(data[data==0])
	rnorm(len[1], mean=0, sd=0.1)
}

#
# Select and write functionally coherent nodes
#

coherent_nodes = function(nfcs, mixmodel, output) {
## extract the mixture model density
	density = mog.density(nfcs$normP, mixmodel, sort.data=T)

# how many modes in the model?
	modes = length(mixmodel$means)

## Do the thresholding
	if (modes == 2) { ## bimodal model always remove the first mode only
		thresh = max(density$x[density$membership <2])	
		data_out = paste(output, ".coherentNodes", sep="")
		coherent = nfcs$id[nfcs$normP > thresh]
		write(coherent, data_out, sep = "\t", ncolumns=1)
		log_file = paste(output, ".log", sep="")
		write(c("Biomodal model found, removed first mode at NFCS threshold: ", thresh), file=log_file, append="TRUE")
	
	} else {
		if (modes == 1) { # unimodal model so no thresholding
			thresh = "NA"
			log_file = paste(output, ".log", sep="")
			write("Unimodal model found, no NFCS thresholding done", file= log_file, append="TRUE")
			coherent = "NA"
		
		} else { # model with >2 modes
		# check if mode 1 is formed from isolated nodes (mean < 0.05)
			if (mixmodel$means[1] < 0.05) { # threshold first two modes
				thresh = max(density$x[density$membership <3])	
				data_out = paste(output, ".coherentNodes", sep="")
				coherent = nfcs$id[nfcs$normP > thresh]
				write(coherent, data_out, sep = "\t", ncolumns=1)
				log_file = paste(output, ".log", sep="")
				write(c("More than 2 modes found. First mode NFCS <0.05, therefore removed first two modes at NFCS threshold: ", thresh), file=log_file, append="TRUE")
		
		
			} else { # threshold first mode
				thresh = max(density$x[density$membership <2])	
				data_out = paste(output, ".coherentNodes", sep="")
				coherent = nfcs$id[nfcs$normP > thresh]
				write(coherent, data_out, sep = "\t", ncolumns=1)
				log_file = paste(output, ".log", sep="")
				write(c("More than 2 modes found. Removed first mode at NFCS threshold: ", thresh), file=log_file, append="TRUE")
			}
			
		}	
	}
	return(coherent)
}




#####
##### Generic functions and variables
#####


##
## ##### PLOT SETTINGS #####
##

##' Line types, for plot()
##' @keywords settings
##' @export
tma.ltys <- c("22", "44", "13", "1343", "73", "2262", "12223242", "F282", "F4448444", "224282F2", "F1")

##' Colours, for plot()
##' @keywords settings
##' @export
tma.lcols <- c(rgb(27, 158, 119, max=255), rgb(217, 95, 2, max=255), rgb(117, 112, 179, max=255), rgb(231, 41, 138, max=255), rgb(102, 166, 30, max=255),colors()[30], colors()[376], colors()[567], colors()[477], colors()[619], colors()[556])

##
## ##### MATRIX MANIPULATION #####
##

##' Implements Matlab's repmat function to replicate matrices. Matrix a is
##' replicated (tiled) n times vertically and m times horizontally.
##'
##' @title Matlab-style repmat
##' @param a A matrix or vector to replicate
##' @param n Number of times to replicate vertically
##' @param m Number of times to replicate horizontally
##' @return A tiled matrix of n by m copies of a
##' @keywords internal
repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}

##' Converts NA values to 0 in a matrix or vector
##'
##' @title Covert NA values in a vector to 0
##' @param x A numeric vector
##' @return The input vector with NAs converted to 0
##' @keywords internal
na.to.zero <- function(x) {x[is.na(x)] <- 0;x}

##' Finds the closest value in the second argument to each
##' of the values in the first argument.
##'
##' In the case of ties, the first matched value in y is
##' chosen.
##'
##' @title Find the closest value in a vector
##' @param x The source values
##' @param y The vector of values to match against
##' @param index.only if true, only return the index rather than
##' the actual value
##' @return The closest element to x in y
##' @keywords internal
closest.value <- function(x,y,index.only=FALSE) {
  index <- sapply(x,function(xeach){which(abs(y-xeach)==min(abs(y-xeach)))[1]})
  ifelse(index.only,return(index),return(y[index]))
}

##' Takes an input matrix and returns a permutation of the values as
##' a new matrix of the same size.
##'
##' @title Randomly permute a matrix
##' @param mat Input matrix
##' @param randomise.within.columns Only randomise the matrix values
##' within columns if \code{TRUE}, else shuffle the whole matrix
##' @return The matrix with the ordering of its values randomised
##' @export
##' @examples
##' m <- matrix(1:16,4,4) # create 4x4 matrix with numbers 1..16
##' randomise.matrix(m)
randomise.matrix <- function(mat,randomise.within.columns=FALSE) {
  if(randomise.within.columns) {
    for(i in 1:ncol(mat)) mat[,i] <- sample(x=mat[,i])
  } else {
    mat[,] <- sample(as.vector(as.matrix(mat)))
  }
  mat
}

##' Extracts values from a matrix (column by column) and convert
##' to a vector. Missing values (\code{NA}) are removed if
##' \code{na.omit} is \code{TRUE}.
##'
##' @title Extract values from a matrix
##' @param mat A matrix of numeric values
##' @param na.omit Omit missing values if \code{TRUE}
##' @return A vector containing the values from the matrix
##' @export
##' @examples
##' m <- matrix(1:16,4,4) # create 4x4 matrix with numbers 1..16
##' extract.values(m)
extract.values <- function(mat,na.omit=TRUE) {
  v <- c(as.matrix(mat))
  if(na.omit) v <- na.omit(v)
  v
}

##
## ##### KERNEL DENSITY ESTIMATION #####
##

##' Calculates the geometric mean of a vector. Used for variable bandwidth kernel density estimation.
##'
##' @title Geometric mean
##' @param x A numeric vector
##' @return The geometric mean of x
##' @keywords internal
geometric.mean <- function(x) exp(mean(log(x)))

##' Calculates KDE for a set of points exactly, rather than an approximation
##' as per \R's \code{density} function.
##'
##' Only tractable for around 10,000 data points or less - otherwise consider
##' using \R's \code{density} function for a close approximation.
##'
##' \R's \code{density} approximation is normally a very good approximation,
##' but some small values close to zero may become zero rather than just
##' very small. This makes it less suitable for mutual information
##' estimation.
##'
##' @title Exact kernel density esimation
##' @param x A numeric vector of values
##' @param bw The bandwidth to use - either a single value,
##' or a vector of values the same length as \code{x} if using adaptive
##' bandwidth estimation (with each value giving the bandwidth at the
##' corresponding data point).
##' @param output.domain The domain of values over which to estimate the
##' density. Defaults to \code{x}. To use the same domain of \code{x}
##' values as \R's \code{density}, set to \code{NULL}.
##' @param na.rm Remove missing values if \code{TRUE}
##' @return The exact kernel density estimate as a \code{density} object,
##' compatible with \R's \code{density} function.
##' @seealso \code{\link{variable.bw.kde}} for (exact) variable bandwidth KDE\cr
##' \code{density} for \R's implementation of approximate KDE
##' @keywords kde
##' @export
##' @examples
##' x <- c(rnorm(100,sd=2),runif(50),4*rnorm(75,mean=2,sd=4))
##' plot(exact.kde(x,bw.nrd0(x),output.domain=NULL),xlab="x",
##'      col="black",lty=2,main="Exact vs approximate KDE")
##' lines(density(x),col="red",lty=1)
##' legend("topright",c("Exact KDE","Approx KDE"),lty=2:1,col=c("black","red"))
exact.kde <- function(x,bw,output.domain=x,na.rm=FALSE) {
  dens <- list(call=match.call(),data.name=deparse(substitute(x)),
               has.na=any(is.na(x)),bw=bw)
  if(na.rm) {
    x <- na.omit(x)
    output.domain <- na.omit(output.domain)
  }
  if(is.null(output.domain)) output.domain <- seq(min(x)-(mean(bw)*3),max(x)+(mean(bw)*3),length.out=512)
  if(length(bw)==1) bw <- rep(bw,length(output.domain))
  dens$x <- output.domain[order(output.domain)]
  dens$y <- sapply(output.domain,function(point) sum(dnorm(point, mean=x, sd=bw)))/length(output.domain)
  dens$y <- dens$y[order(output.domain)]
  dens$n <- length(x)
  class(dens) <- "density"
  dens
}

##' Calculates variable bandwidth KDE using Abramson's two stage estimator.
##'
##' Bandwidth is first calculated using Silverman's estimator, then refined
##' in a second stage to allow local bandwidth variations in the data based
##' on the initial estimate.
##'
##' @title Variable bandwidth Kernel Density Estimation
##' @param x A numeric vector of values for estimating density
##' @param output.domain The domain of values over which to estimate the
##' density. Defaults to \code{x}. To use the same domain of \code{x}
##' values as \R's \code{density}, set to \code{NULL}.
##' @param na.rm Remove missing values if TRUE
##' @return The kernel density estimate as a \code{density} object,
##' compatible with \R's \code{density} function.
##' @references Abramson, I. S. On Bandwidth Variation in Kernel Estimates-A
##' Square Root Law. Ann. Statist. 10, 1217-1223 (1982).
##' @seealso \code{\link{exact.kde}} for exact KDE (used by this function)
##' @keywords kde
##' @export
##' @examples
##' x <- rnorm(100)
##' plot(exact.kde(x,bw=bw.nrd0(x),output.domain=NULL),col="black",lty=2,
##'      xlab="x",main="Fixed vs variable bandwidth KDE")
##' lines(variable.bw.kde(x),col="red",lty=1)
##' legend("topright",c("Fixed BW KDE","Variable BW KDE"),lty=2:1,col=c("black","red"))
variable.bw.kde <- function(x,output.domain=x,na.rm=FALSE) {
  if(na.rm) {
    x <- na.omit(x)
    output.domain <- na.omit(output.domain)
  }
  base.bw <- bw.nrd0(x)
  ##same x axis as R's density() uses by default
  if(is.null(output.domain)) output.domain <- seq(min(x)-(base.bw*3),max(x)+(base.bw*3),length.out=512)

  d.pilot <- exact.kde(x,base.bw,output.domain=output.domain)
  #prevent divide by zero
  if(any(d.pilot$y < .Machine$double.eps))
    d.pilot$y[d.pilot$y < .Machine$double.eps] <- .Machine$double.eps
  bw.adjust.factor <- (geometric.mean(d.pilot$y)/d.pilot$y)^0.5
  d <- exact.kde(x,base.bw*bw.adjust.factor,output.domain=d.pilot$x)
  d
}

##' Adapted from minor.tick function from Hmisc.
##' Included to remove dependency on Hmisc based on a single function.
##' 
##' @author Frank E. Harrell Jr.
##' @keywords internal
minor.tick <- function (nx = 2, ny = 2, tick.ratio = 0.5) 
{
    ax <- function(w, n, tick.ratio) {
        range <- par("usr")[if (w == "x") 
            1:2
        else 3:4]
        tick.pos <- if (w == "x") 
            par("xaxp")
        else par("yaxp")
        distance.between.minor <- (tick.pos[2] - tick.pos[1])/tick.pos[3]/n
        possible.minors <- tick.pos[1] - (0:100) * distance.between.minor
        low.minor <- min(possible.minors[possible.minors >= range[1]])
        if (is.na(low.minor)) 
            low.minor <- tick.pos[1]
        possible.minors <- tick.pos[2] + (0:100) * distance.between.minor
        hi.minor <- max(possible.minors[possible.minors <= range[2]])
        if (is.na(hi.minor)) 
            hi.minor <- tick.pos[2]
        axis(if (w == "x") 
             1
        else 2, seq(low.minor, hi.minor, by = distance.between.minor), 
             labels = FALSE, tcl = par("tcl") * tick.ratio)
    }
    if (nx > 1) 
        ax("x", nx, tick.ratio = tick.ratio)
    if (ny > 1) 
        ax("y", ny, tick.ratio = tick.ratio)
    invisible()
}







