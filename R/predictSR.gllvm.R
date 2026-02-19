#' @title Predict species richness from a gllvm object
#' @description Obtains predictions from a fitted generalized linear latent variable model object.
#'
#' @param object an object of class 'gllvm'.
#' @param SR integer, defaults to NULL. If omitted, returns a prediction for all integers up to the number of observed species.
#' @param se.fit integer or logical, defaults to 10.000. Number of simulations for confidence interval. No confidence interval is returned when set to \code{FALSE}.
#' @param alpha numeric between 0 and 1, defaults to 0.95. Confidence level of the prediction.
#' @param seed numeric, defaults to 42. Seed for the simulation of the confidence interval.
#' @param return.pred logical, defaults to \code{FALSE}. Returns the results from \code{"\link{predict.gllvm}"}.
#' @param ... arguments passed to \code{"\link{predict.gllvm}"}.
#'
#' @details This function returns probabilities of richness for species richness from gllvm objects, by first calling \code{"\link{predict.gllvm}"}, and following by predicting from a Poisson-Binomial distribution for richness. The distribution for species richness follows from a sum of binary responses (i.e., occurrence), and is naturally extended to non-normal data types as the probability to get a non-zero observation. Especially with many species and rows, the calculation may take a while.
#' @note The point estimate ("fit") can sometimes be outside of the simulated intervals when the model fit is poor, the standard errors of parameter estimates are on very different scales, or the number of simulations is too low. Try either 1) scaling and centering covariates in the model, 2) repeatedly refitting the model (see the n.init argument in \code{"\link{gllvm}"}) to find a better set of starting values and improve the fit, or  3) increasing the number of simulations via 'n.init'.
#'
#' @return list with entries "predicted" and "expected". "predicted" includes the prediction from the Poisson-Binomial distribution, returning a matrix of size sites by length(SR), with statistical uncertianty if se.fit = TRUE. "expected" includes the expected species' richness, i.e., a vector of size n, with statistical uncertianty if se.fit = TRUE.
#' 
#' @seealso \code{\link{predict.gllvm}}, \code{\link{residuals.predictSR.gllvm}}
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
#' plot(x = newX$soil.dry, y = predSR$predicted$fit[,i+1], xlab = "Soil dry matter content", ylab = bquote(p(SR == .(i))), type = "l", ylim = c(0,1))
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

predictSR.gllvm <- function(object, SR = NULL, se.fit = 10000, alpha = 0.95, seed = 42, return.pred = FALSE, ...){
  
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
    
    eSR <- rowSums(probs)
    eSR.CI <- matrix(nrow = n, ncol = 2)
    
    if(is.null(SR))SR <-0:ncol(probs)

    if(TRUE){#type == "direct"){
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
        probs.sim <- gllvm.presence.prob(aperm(preds$ci.sim[,,i,],c(2,1,3)),object)
      }
    eSR.CI[i,]  <- quantile(rowSums(probs.sim), prob = c((1-alpha)/2, 1-(1-alpha)/2))
      
    predSR.sim <- poisbinom(probs.sim)
    
    center = predSR[i,]         # point estimate for site i
    dists = hilbert_to_provided_center(pmin(pmax(predSR.sim,.Machine$double.eps),1-.Machine$double.eps), pmin(pmax(center, .Machine$double.eps), 1-.Machine$double.eps))
    
    threshold = quantile(dists, alpha)
    
    # Keep points inside the simultaneous confidence region
    SR.ci[,i,] <- apply(predSR.sim[dists <= threshold, ],2,range)
    }
    out <- list(predicted = list(fit = predSR, lower = SR.ci[1,,], upper = SR.ci[2,,]), expected = list(fit = eSR, lower = eSR.CI[,1], upper = eSR.CI[,2]))
    
  }
    # else if(type == "empirical"){
    # p <- ncol(object$y)
    # R <- se.fit
    # n <- nrow(fit)
    # if(any(object$family == "ordinal"))n <- nrow(fit[1,,])
    # 
    # if(!is.numical(se.fit))se.fit <- 1e3
    # 
    # predSR  <-  array(rbinom(R*nrow(fit)*ncol(fit), size = 1, prob = array(rep(probs, each = R), dim = c(R, n, p))), dim = c(R, n, p))
    # 
    # predSR <- rowSums(predSR, dims=2)
    # predSR <- pmin(pmax(t(table(factor(predSR,levels=0:p), col(predSR))/R), 1e-15),1-1e-15)
    # colnames(predSR) <- paste0("SR_", SR)
    # 
    # SR.ci <- array(dim=c(2,n,length(SR)))
    # predSR.sim = array(0, dim=c(R,n,p))
    # predSR.sim.mat <- matrix(0, nrow = n, ncol = p)
    # dists = numeric(R)
    # center = numeric(length(SR))
    # 
    # pred.sim.mat <- matrix(0, ncol = p+1)
    # for(i in 1:n){
    #   probs.sim <- gllvm.presence.prob(preds$ci.sim[,i,],object)
    #   # bad idea to fit this into memory
    #   # predSR.sim  <-  array(rbinom(R*R*ncol(fit), size = 1, prob = array(rep(probs.sim, each = R), dim = c(R, R, p))), dim = c(R, R, p))
    #   # predSR.sim.mat <- rowSums(predSR.sim, dims = 2)
    #   predSR.sim.mat <-  t(replicate(R, table(factor(replicate(R, sum(rbinom(p, 1, probs.sim[r, ]))), levels = 0:p))/R))
    #   predSR.sim.mat <- pmin(pmax(predSR.sim.mat, 1e-15),1-1e-15)
    # 
    #   
    #   center = predSR[i,] # point estimate for site i
    #   dists = hilbert_to_provided_center(predSR.sim.mat, center)
    #   
    #   threshold = quantile(dists, alpha)
    #   
    #   # Keep points inside the simultaneous confidence region
    #   SR.ci[,i,] <- apply(predSR.sim.mat[dists <= threshold, ],2,range)
    # }
    # out <- list(predicted = list(fit = predSR, lower = SR.ci[1,,], upper = SR.ci[2,,]), expected = list(eSR, lower = eSR.CI[,1], upper = eSR.CI[,2]))
    # 
    # }
  }else if(!is.list(preds)){
    n <- nrow(fit)
    
    probs <- gllvm.presence.prob(fit, object)  
    
    if(is.null(SR))SR <-0:ncol(probs)
    predSR  <-  poisbinom(probs)
    colnames(predSR) <- paste0("SR_", SR)
    out <- list(predicted = list(fit = predSR), expected  = list(fit = rowSums(probs)))
    
    # if(type == "class"){
    #   out <- apply(predSR, 1, which.max)-1
    # }
  }
  
  if(return.pred){
    # code from predict.gllvm
    if(is.list(preds)){
    if(!any(object$family == "ordinal")){
      ci <- apply(preds$ci.sim, 2:3, quantile, prob = c((1-alpha)/2, 1-(1-alpha)/2))
      out$predict.gllvm <- list(fit = preds$fit, lower = ci[1,,], upper = ci[2,,])
  }else if(any(object$family == "ordinal")){
    ci <- array(dim=c(2,nrow(preds$fit), nrow(object$y),ncol(object$y)))
    for(j in 1:ncol(object$y)){
      if(object$family[j] != "ordinal"){
        ci[,1,,j] <- apply(preds$ci.sim[,,j], 2, prob = c((1-alpha)/2, 1-(1-alpha)/2)) 
      }else{
        # here we need to be a bit more careful because of the simplex constraint
        # we take the same approach as in predictSR.gllvm.R
        # we measure the distance to the center of the simplex (which we will consider the prediction on the point estimates)
        # and threshold over the distances
        for(i in 1:nrow(object$y)){
          center = pmin(pmax(preds$fit[,i,j], .Machine$double.eps), 1-.Machine$double.eps)
          dists = suppressWarnings(hilbert_to_provided_center(pmin(pmax(preds$ci.sim[,,i,j], .Machine$double.eps),1-.Machine$double.eps), center))
          threshold = quantile(dists, alpha, na.rm = TRUE)
          ci[,,i,j] <- suppressWarnings(apply(preds$ci.sim[,,i,j][dists <= threshold, ],2,range, na.rm = TRUE))
        }
        if(any(!is.finite(ci)))ci[!is.finite(ci)] <- NA # guaranteed to happen on zeta.struc = "species"
      }
    }
    out$predict.gllvm <- list(fit = preds$fit, lower = ci[1,,,], upper = ci[2,,,])
  }
    }else{
      out$predict.gllvm <- list(fit = fit)
    }
  }
  
  class(out) <- "predictSR.gllvm"
  return(out)
  
}

#'@export
predictSR <- function(object, ...) {
  UseMethod("predictSR")
}

hilbert_to_provided_center <- function(mat, center) {
  # Ensure inputs are valid
  if (any(mat[!is.na(mat)] <= 0) || any(center[!is.na(center)] <= 0)) {
    stop("Hilbert metric requires strictly positive values (interior of the simplex).")
  }
  
  if (ncol(mat) != length(center)) {
    stop("Dimension mismatch: 'center' must have the same length as the number of columns in 'mat'.")
  }
  
  ratios <- sweep(mat, 2, center, "/")
  
  row_max <- apply(ratios, 1, max, na.rm = TRUE)
  row_min <- apply(ratios, 1, min, na.rm = TRUE)
  
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
#' as convenience wrapper around \code{"\link{predict.gllvm}"}.
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


#' @title Dunn-Smyth residuals for species richness from a gllvm object
#' @description Calculates Dunn-Smyth residuals for species richness.
#'
#' @param object an object of class 'gllvm'.
#' @param predSR object returned by \code{"\link{predictSR.gllvm}"}
#' @param ... not used.
#'
#' @details
#' See \code{"\link{residuals.gllvm}"} for details.
#' 
#' @return A list of length 2 with 1) the fitted and 2) the dunn-smyth residuals.
#' 
#' @seealso \code{\link{predict.gllvm}}, \code{\link{predictSR.gllvm}}, \code{\link{residuals.gllvm}}
#' 
#' @author Bert van der Veen
#'
#'@author Bert van der Veen
#'
#'@aliases residuals.predictSR.gllvm
#'@method residuals predictSR.gllvm
#'@export
#'@export residuals.predictSR.gllvm
residuals.predictSR.gllvm <- function(predSR, object, ...){
  
  # ensure ordinal starts at 0
  if(any(object$family == "ordinal"))object$y[,object$family=="ordinal"] <- object$y[,object$family=="ordinal"] - apply(object$y[,object$family=="ordinal"],2,min)
  
  y = rowSums(ifelse(object$y==0,0,1))
  n = length(y)
  pmf = predSR$predicted$fit
  
  K <- ncol(pmf) - 1
  
  cdf <- t(apply(pmf, 1, cumsum))
  
  b <- cdf[cbind(seq_len(n), pmin(y, K) + 1)]
  
  a <- numeric(n)
  nz <- y > 0
  a[nz] <- cdf[cbind(which(nz), y[nz])]
  
  u <- a + (b - a) * runif(n)
  
  if (any(u == 1, na.rm = TRUE))
    u[u == 1] <- 1 - 1e-16
  
  if (any(u == 0, na.rm = TRUE))
    u[u == 0] <- 1e-16
  
    list(fitted = predSR$expected$fit, residuals = qnorm(u))
}

#'@export
plot.predictSR.gllvm <- function(predSR, object, ...){
  res <- residuals.predictSR.gllvm(predSR, object)
  
  plot(xlab = "Expected species richness", ylab = "Dunn-Smyth residuals", x = res$fitted, y = res$residuals)
}

