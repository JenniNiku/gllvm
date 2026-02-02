#' @title Predict species richness from a gllvm object
#' @description Obtains predictions from a fitted generalized linear latent variable model object.
#'
#' @param object an object of class 'gllvm'.
#' @param SR integer, defaults to NULL. If omitted, returns a prediction for all integers up to the number of observed species.
#' @param se.fit integer or logical, defaults to 10.000. Number of simulations for confidence interval. No confidence interval is returned when set to \code{FALSE}.
#' @param type character, defaults to \code{direct}. For 'direct', simulates confidence intervals by simulating species' probabilities, and evaluating the Poisson-Binomial PMF. Alternatively, type 'empirical' simulates from the Bernoulli distribution of species, and calculates empirical probabilities instead from sums of Bernoulli simulations.
#' @param alpha numeric between 0 and 1, defaults to 0.95. Confidence level of the prediction.
#' @param seed numeric, defaults to 42. Seed for the simulation of the confidence interval.
#' @param ... arguments passed to \code{"\link{predict.gllvm}"}.
#'
#' @return A matrix of size sites by length(SR).
#' 
#' @seealso \code{\link{predict.gllvm}}
#' 
#' @author Bert van der Veen
#'
#' @examples
#' \donttest{
#'# Load a dataset from the mvabund package
#'data(spider, package = "mvabund")
#'y <- as.matrix(spider$abund)
#'X <- scale(spider$x)
#'# Fit gllvm model
#'fit <- gllvm(y = y, X, formula = ~soil.dry, family = poisson())
#'# Fitted values
#'newX = data.frame(soil.dry = seq(min(X[,"soil.dry"]),max(X[,"soil.dry"]),length.out=100))
#'predSR <- predictSR(fit, newX = newX, level = 0)
#'
#'# Visualize the results
#' par(mfrow = c(4,4))
#' for(i in 0:ncol(fit$y)){
#' plot(x = newX$soil.dry, y = predSR[,i+1], xlab = "Soil dry matter content", ylab = bquote(p(SR == .(i))), type = "l", ylim = c(0,1))
#' }
#'
#'}
#'
#'@author Bert van der Veen
#'
#'@aliases predictSR predictSR.gllvm
#'@method predictSR gllvm
#'@export
#'@export predictSR.gllvm

predictSR.gllvm <- function(object, SR = NULL, type = "direct", se.fit = 10000, alpha = 0.95, seed = 42,...){
  
  set.seed(seed)
  preds <- predict(object, type = "response", alpha = NULL, se.fit = se.fit, ...)
  fit <- preds
  
  if(is.list(preds)){
    R = dim(preds$ci.sim)[1]
    n = dim(preds$ci.sim)[2]
    
    fit <- preds$fit
    if(any(object$family == "ordinal")){
      n = nrow(fit[1,,])
    }
    probs <- gllvm.presence.prob(fit, object)
    
    if(is.null(SR))SR <-0:ncol(probs)

    if(type == "direct"){
    # it's slow, but we can do our best
    TMBcores <- TMB::openmp(DLL="gllvm")[1]
    PoisBinCores <- get_omp_threads()
    
    if(PoisBinCores == 1 && TMBcores > 1)set_omp_threads(TMBcores)
      
    predSR  <-  poisbinom(probs)
    colnames(predSR) <- paste0("SR_", SR)
      
    SR.ci <- array(dim=c(2,n,length(SR)))
    predSR.sim = matrix(0, nrow = nrow(preds$ci.sim), ncol = length(SR))
    dists = numeric(R)
    center = numeric(length(SR))

    for(i in 1:n){
      if(!any(object$family == "ordinal")){
        probs.sim <- gllvm.presence.prob(preds$ci.sim[,i,],object)
      }else{
        probs.sim <- gllvm.presence.prob(preds$ci.sim[,,i,],object)
      }
    predSR.sim <- poisbinom(probs.sim)
    
    center = predSR[i,]         # point estimate for site i
    dists = hilbert_to_provided_center(pmin(pmax(predSR.sim,.Machine$double.eps),1-.Machine$double.eps), pmin(pmax(center, .Machine$double.eps), 1-.Machine$double.eps))
    
    threshold = quantile(dists, alpha)
    
    # Keep points inside the simultaneous confidence region
    SR.ci[,i,] <- apply(predSR.sim[dists <= threshold, ],2,range)
    }
    out <- list(fit = predSR, lower = SR.ci[1,,], upper = SR.ci[2,,])
    
  }else if(type == "empirical"){
    p <- ncol(object$y)
    R <- se.fit
    n <- nrow(fit)
    if(any(object$family == "ordinal"))n <- nrow(fit[1,,])
  
    if(!is.numical(se.fit))se.fit <- 1e3
    
    predSR  <-  array(rbinom(R*nrow(fit)*ncol(fit), size = 1, prob = array(rep(probs, each = R), dim = c(R, n, p))), dim = c(R, n, p))
    
    predSR <- rowSums(predSR, dims=2)
    predSR <- pmin(pmax(t(table(factor(predSR,levels=0:p), col(predSR))/R), 1e-15),1-1e-15)
    colnames(predSR) <- paste0("SR_", SR)
    
    SR.ci <- array(dim=c(2,n,length(SR)))
    predSR.sim = array(0, dim=c(R,n,p))
    predSR.sim.mat <- matrix(0, nrow = n, ncol = p)
    dists = numeric(R)
    center = numeric(length(SR))
    
    pred.sim.mat <- matrix(0, ncol = p+1)
    for(i in 1:n){
      probs.sim <- gllvm.presence.prob(preds$ci.sim[,i,],object)
      # bad idea to fit this into memory
      # predSR.sim  <-  array(rbinom(R*R*ncol(fit), size = 1, prob = array(rep(probs.sim, each = R), dim = c(R, R, p))), dim = c(R, R, p))
      # predSR.sim.mat <- rowSums(predSR.sim, dims = 2)
      predSR.sim.mat <-  t(replicate(R, table(factor(replicate(R, sum(rbinom(p, 1, probs.sim[r, ]))), levels = 0:p))/R))
      predSR.sim.mat <- pmin(pmax(predSR.sim.mat, 1e-15),1-1e-15)

      
      center = predSR[i,] # point estimate for site i
      dists = hilbert_to_provided_center(predSR.sim.mat, center)
      
      threshold = quantile(dists, alpha)
      
      # Keep points inside the simultaneous confidence region
      SR.ci[,i,] <- apply(predSR.sim.mat[dists <= threshold, ],2,range)
    }
    out <- list(fit = predSR, lower = SR.ci[1,,], upper = SR.ci[2,,])
    
    
    }
    return(out)
  }else if(!is.list(preds)){
    probs <- gllvm.presence.prob(fit, object)  
    
    if(is.null(SR))SR <-0:ncol(probs)
    predSR  <-  poisbinom(probs)
    colnames(predSR) <- paste0("SR_", SR)
    out <- predSR
    
    # if(type == "class"){
    #   out <- apply(predSR, 1, which.max)-1
    # }
    
    return(out)
  }
}

#'@export
predictSR <- function(object, ...) {
  UseMethod("predictSR")
}

hilbert_to_provided_center <- function(mat, center) {
  # Ensure inputs are valid
  if (any(mat <= 0) || any(center <= 0)) {
    stop("Hilbert metric requires strictly positive values (interior of the simplex).")
  }
  
  if (ncol(mat) != length(center)) {
    stop("Dimension mismatch: 'center' must have the same length as the number of columns in 'mat'.")
  }
  
  # Calculate the ratio of each sample component to the corresponding center component
  # ratios[i, j] = mat[i, j] / center[j]
  ratios <- sweep(mat, 2, center, "/")
  
  # Hilbert distance d_H(x, c) = log( max(x_i/c_i) / min(x_i/c_i) )
  row_max <- apply(ratios, 1, max)
  row_min <- apply(ratios, 1, min)
  
  return(log(row_max / row_min))
}

#' @title Predict co-occurrence from a gllvm object
#' @description Obtains the probability of co-occurrnece for species from a fitted generalized linear latent variable model object.
#'
#' @param object an object of class 'gllvm'.
#' @param spp a vector of length 2, defaults to NULL. If NULL, returns the probability of co-occurrence for all species in the data.
#' @param se.fit logical, defaults to \code{FALSE}. If set to \code{TRUE}, returns a simulated confidence interval of the prediction.
#' @param alpha numeric between 0 and 1, defaults to 0.95. Confidence level of the prediction.
#' @param seed numeric, defaults to 42. Seed for the simulation of the confidence interval.
#' @param ... arguments passed to \code{"\link{predict.gllvm}"}.
#'
#' @details
#' The probability of co-occurrence is simply given by the joint probability of getting two presences for any pair of species. As such,
#' it is simply calculated as the product of the probability of occurrence for any pair of species, so that this function mostly serves
#' as convenience wrapper around \code{\link{predict.gllvm}}.
#' 
#' @return A matrix of size np(p-1)/2 by 2.
#' 
#' @seealso \code{\link{predict.gllvm}}, \code{\link{predictSR.gllvm}}
#' 
#' @author Bert van der Veen
#'
#' @examples
#' 
#' \donttest{
#'# Load a dataset from the mvabund package
#'data(spider, package = "mvabund")
#'y <- as.matrix(spider$abund)
#'X <- scale(spider$x)
#'# Fit gllvm model
#'fit <- gllvm(y = y, X, formula = ~soil.dry, family = poisson())
#'# fitted values
#'newX = data.frame(soil.dry = seq(min(X[,"soil.dry"]),max(X[,"soil.dry"]),length.out=100))
#'predPair <- predictPairwise(fit, spp = c(1,2), newX = newX, level = 0)
#'plot(x = newX$soil.dry, y = predPair$prob, type = "l", xlab = "Soil dry moisture content", ylab = "Probability of co-occurrence", ylim = c(0,1))
#'
#'}
#'
#'@author Bert van der Veen
#'
#'@aliases predictPairwise predictPairwise.gllvm
#'@method predictPairwise gllvm
#'@export
#'@export predictPairwise.gllvm
predictPairwise.gllvm <- function(object, spp = NULL, se.fit = FALSE, alpha = 0.95, seed = 42, ...){

  preds <- predict(object, type = "response", se.fit = se.fit, alpha = NULL, ...)
  fit <- preds
  
  if(is.null(spp)){spp <- 1:p}else{spp <- sort(spp)}
  spp.pairs <- t(combn(spp, 2))
  
  if(is.list(preds)){
   fit <- preds$fit
   
   probs.fit <- gllvm.presence.prob(fit, object, spp)
   n <- nrow(probs.fit)
   p <- ncol(probs.fit)

   colnames(probs.fit) <- spp
   idx <- matrix(match(spp.pairs, colnames(probs.fit)), ncol = 2)
   
   probs.fit <- cbind(
     probs.fit[, idx[,1 ], drop = FALSE],
     probs.fit[, idx[,2 ], drop = FALSE]
   )
   probs.fit <- matrix(probs.fit, ncol = 2)
   
   predPair <- matrix(exp(rowSums(log(probs.fit))), ncol = 1) # rowProd via exp-sum-log
   
   predPair.sim <- matrix(nrow=dim(preds$ci.sim)[1], ncol =nrow(predPair))
   for(r in 1:dim(preds$ci.sim)[1]){
     probs.fit.sim <- gllvm.presence.prob(preds$ci.sim[r,,], object, spp)
     colnames(probs.fit.sim) <- spp
     probs.fit.sim <- cbind(
       probs.fit.sim[, idx[,1 ], drop = FALSE],
       probs.fit.sim[, idx[,2 ], drop = FALSE]
     )
     probs.fit.sim <- matrix(probs.fit.sim, ncol = 2)
     predPair.sim[r,] <-  exp(rowSums(log(probs.fit.sim))) # rowProd via exp-sum-log

   }
   predPair.ci <- apply(predPair.sim, 2, quantile, prob = c((1-alpha)/2, 1-(1-alpha)/2))
   out <- data.frame(prob = predPair,  site = rep(seq_len(n), times = nrow(spp.pairs)), spp1 = rep(spp.pairs[,1], each = n), spp2 = rep(spp.pairs[,2], each =n), lower = predPair.ci[1,], upper = predPair.ci[2,])
  }else{
    probs.fit <- gllvm.presence.prob(fit, object, spp)
    n <- nrow(probs.fit)
    
    colnames(probs.fit) <- spp
    idx <- matrix(match(spp.pairs, colnames(probs.fit)), ncol = 2)
    
    probs.fit <- cbind(
      probs.fit[, idx[,1 ], drop = FALSE],
      probs.fit[, idx[,2 ], drop = FALSE]
    )
    probs.fit <- matrix(probs.fit, ncol = 2)
    
    predPair <- matrix(exp(rowSums(log(probs.fit))), ncol = 1)
  
    out <- data.frame(prob = predPair,  site = rep(seq_len(n), times = nrow(spp.pairs)), spp1 = rep(spp.pairs[,1], each = n),spp2 = rep(spp.pairs[,2], each =n))
  }
  return(out)
}

#'@export
predictPairwise <- function(object, ...) {
  UseMethod("predictPairwise")
}

# prediction of binomial-based distributions can be used as-is (including orderedBeta)
# not implemented because not based (indirectly) on binary data: gaussian, beta, betaH, gamma, exponential
# families implemented: poisson, NB, NB1, ZIP, ZINB, Tweedie, ordinal
gllvm.presence.prob <- function(fit, object, spp = NULL) {

  family <- object$family
  n <- nrow(fit)
  if(any(object$family == "ordinal")) n <- nrow(fit[1,,])
  if(is.null(spp))spp <- 1:ncol(object$y)
  probs <- matrix(0, nrow = n, ncol = length(spp))
    
  for(j in spp){
  if(family[j] %in% c("binomial", "ZIB", "ZNIB"))probs[,j] <- fit[, j]
  if(family[j] == "poisson") {
    probs[,j] <- ppois(matrix(0, nrow = n, ncol = 1), lambda = fit[,j], lower.tail = FALSE)
  } else if(family[j] == "negative.binomial") {
    size <- 1 / rep(object$params$phi[j], n)
    probs[,j] <- pnbinom(rep(0,n), mu = fit[,j], size = size, lower.tail = FALSE)
  } else if(family[j] == "negative.binomial1") {
    size <- fit[,j] * object$params$phi[j]
    probs[,j] <- pnbinom(rep(0, n), mu = fit, size = size, lower.tail = FALSE)
  } else if(family[j] == "betaH") {
    probs[,j] <- 1 - fit[,-seq_len(ncol(object$y)), drop = FALSE][,j]
  } else if(family[j] == "orderedBeta"){
    probs[,j] <- 1 - binomial(link = object$link)$linkinv(binomial(link = object$link)$linkfun(fit[,j]) - object$params$zeta[j,1])
  }else if(family[j] == "ZIP") {
    probs[,j] <- 1 - pzip(rep(0, n), mu = fit[,j], sigma = rep(object$params$phi[j], each = n))
  } else if(family[j] == "ZINB") {
    probs[,j] <- 1 - pzinb(rep(0, n), mu = fit[,j],
                              p = rep(object$params$phi[j], each = n),
                              sigma = rep(object$params$ZINB.phi[j], each = n))
  } else if(family[j] == "tweedie") {
    probs[,j] <- 1 - fishMod::pTweedie(rep(0, n), mu = fit[,j],
                                   phi = rep(object$params$phi[j], each = n),
                                   p = object$Power)
  } else if(family[j] == "ordinal") {
    probs[,j] <- 1 - fit[1,,j]
  }
  }
  
  return(probs)
}

get_omp_threads <- function() {
  .Call("C_get_omp_threads", PACKAGE = "gllvm")
}

set_omp_threads <- function(t) {
  t <- as.integer(t)
  .Call("C_set_omp_threads", t, PACKAGE = "gllvm")
}

poisbinom <- function(prob) {
  .Call("C_poisbinom", prob, PACKAGE = "gllvm")
}