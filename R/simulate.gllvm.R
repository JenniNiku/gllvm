#' @title Simulate data from gllvm fit
#' @description Generate new data using the fitted values of the parameters
#'
#' @param object an object of class 'gllvm'.
#' @param nsim an optional positive integer specifying the number of simulated datasets. Defaults to 1.
#' @param seed an optional integer to set seed number, passed to set.seed. Defaults to a random seed number.
#' @param conditional if \code{conditional = FALSE} simulates marginally over the latent variables.
#' @param ... not used.
#'
#' @details
#' simulate function for gllvm objects. 
#' 
#' @return A matrix containing generated data.
#' @author David Warton, Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @examples
#'  \donttest{
#'# Load a dataset from the mvabund package
#'data(antTraits, package = "mvabund")
#'y <- as.matrix(antTraits$abund)
#'X <- scale(antTraits$env[, 1:3])
#'# Fit gllvm model
#'fit <- gllvm(y = y, X, family = poisson())
#'# Simulate data
#'newdata <- simulate(fit)
#'}
#'@aliases simulate simulate.gllvm
#'@method simulate gllvm
#'@export
#'@export simulate.gllvm

simulate.gllvm = function (object, nsim = 1, seed = NULL, conditional = FALSE, ...) 
{
  # code chunk from simulate.lm to sort out the seed thing:
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  if (is.null(seed)) 
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  nRows = dim(object$y)[1]
  nCols = dim(object$y)[2]
  if((object$num.lv+object$num.lv.c)==0){
    lvsNew <- NULL
  }else if(conditional == FALSE){
    # generate new latent variables
    lvsNew = matrix(rnorm(nsim*nRows*(object$num.lv+object$num.lv.c)),ncol=(object$num.lv+object$num.lv.c))
  }else{
    lvsNew = object$lvs[rep(1:nRows,nsim),,drop=FALSE]
  }
  if(is.null(object$X))   {
    prs = predict.gllvm(object,newLV = lvsNew,type="response")
  }  else if(is.null(object$TR)){ 
    Xnew <- object$X[rep(1:nRows,nsim),,drop=FALSE]; colnames(Xnew) <- colnames(object$X)
    prs = predict.gllvm(object,newX=Xnew, newLV = lvsNew,type="response")
  } else {
    Xnew <- object$X[rep(1:nRows,nsim),,drop=FALSE]; colnames(Xnew) <- colnames(object$X)
    prs = predict.gllvm(object,newX=Xnew, newLV = lvsNew,type="response")
  }
  # generate new data
  nTot = nsim*nRows*nCols # total number of values to generate
  if(object$family=="negative.binomial")
    invPhis = matrix(rep(object$params$inv.phi,each=nsim*nRows), ncol=nCols)
  if(object$family=="ZINB")
    invPhis = matrix(rep(object$params$ZINB.inv.phi,each=nsim*nRows), ncol=nCols)
  if(object$family=="tweedie")
    phis = matrix(rep(object$params$phi, each = nsim*nRows), ncol = nCols)
  if(object$family %in% c("gaussian", "gamma", "beta","ZIP","ZINB","ZIB"))
    phis = matrix(rep(object$params$phi, each = nsim*nRows), ncol = nCols)
  if(object$family == "ordinal"){
    if(object$zeta.struc=="species"){
      sims = matrix(0, nrow = nsim * nRows, ncol = nCols)
      for(j in 1:nCols){
        k <- sort(unique(object$y[,j]))
        for(i in 1:(nsim * nRows)){
          sims[i,j] <- sample(k,1,prob=prs[,i,j][!is.na(prs[,i,j])])
        }
      }
      dimnames(prs)[[3]] <- colnames(object$y)
      dimnames(prs)[[2]] <- 1:(nsim * nRows)
      prs <- prs[1,,]
    }else{
      sims = matrix(0, nrow = nsim * nRows, ncol = nCols)
      k <- sort(unique(c(object$y)))
      for(j in 1:nCols){
        for(i in 1:(nsim * nRows)){
          sims[i,j] <- sample(k,1,prob=prs[,i,j][!is.na(prs[,i,j])])
        }
      }
      dimnames(prs)[[3]] <- colnames(object$y)
      dimnames(prs)[[2]] <- 1:(nsim * nRows)
      prs <- prs[1,,]
    }
    
    
  }
  
  newDat = switch(object$family, "binomial"=rbinom(nTot, size = rep(object$Ntrials,each=nsim*nRows), prob = prs),
                  "poisson" = rpois(nTot, prs),
                  "negative.binomial" = rnbinom(nTot, size = invPhis, mu = prs),
                  "gaussian" = rnorm(nTot, mean = prs, sd = phis),
                  "gamma" = rgamma(nTot, shape = phis, scale = prs/phis),
                  "exponential" = rexp(nTot, rate = 1/prs),
                  "tweedie" = fishMod::rTweedie(nTot, mu = c(prs), phi = c(phis), p = object$Power),
                  "ordinal" = sims,
                  "beta" = rbeta(nTot, shape1 = phis*prs, shape2 = phis*(1-prs)),
                  "ZIP" = ifelse(rbinom(nTot, size = 1, prob = phis) > 0, 0, rpois(nTot, lambda = prs)),
                  "ZINB" = ifelse(rbinom(nTot, size = 1, prob = phis) > 0, 0, rnbinom(nTot, size = invPhis, mu = prs)),
                  "ZIB" = ifelse(rbinom(nTot, size = 1, prob = phis) > 0, 0, rbinom(nTot, size = rep(object$Ntrials,each=nsim*nRows), prob = prs)),
                  stop(gettextf("family '%s' not implemented ", object$family), domain = NA))
  # reformat as data frame with the appropriate labels
  newDat = as.data.frame(matrix(newDat,ncol=nCols))
  try(colnames(newDat) <- colnames(prs), silent = TRUE)
  try(rownames(newDat) <- rownames(prs), silent = TRUE)
  try(dimnames(newDat) <- dimnames(prs), silent = TRUE)
  return(newDat)
}


#' @export simulate
simulate <- function(object, ...)
{
  UseMethod(generic = "simulate")
}