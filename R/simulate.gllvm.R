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
simulate.gllvm = function (object, nsim = 1, seed = NULL, conditional = FALSE, newX = NULL, ...) 
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
  if(!is.null(newX)){
    if(conditional)stop("Cannot simulate with 'conditional = FALSE' and new covariates.")
    nRows = nrow(newX)
  }
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
    if(is.null(newX)){
    Xnew <- object$X[rep(1:nRows,nsim),,drop=FALSE]; colnames(Xnew) <- colnames(object$X)
    }else{
      Xnew <- newX[rep(1:nRows,nsim),,drop=FALSE]
    }
    prs = predict.gllvm(object,newX=Xnew, newLV = lvsNew,type="response")
  } else {
    if(is.null(newX)){
    Xnew <- object$X[rep(1:nRows,nsim),,drop=FALSE]; colnames(Xnew) <- colnames(object$X)
    }else{
      Xnew <- newX[rep(1:nRows,nsim),,drop=FALSE]
    }
    prs = predict.gllvm(object,newX=Xnew, newLV = lvsNew,type="response")
  }
  # generate new data
  nTot = nsim*nRows*nCols # total number of values to generate
  if(any(object$family %in% c("negative.binomial", "negative.binomial1", "ZINB"))){
    invPhis <- matrix(NA,nsim*nRows,nCols)
    if(any(object$family %in% c("negative.binomial", "negative.binomial1")))
      invPhis[,object$family %in% c("negative.binomial", "negative.binomial1")] = matrix(rep(object$params$inv.phi[object$family %in% c("negative.binomial", "negative.binomial1")],each=nsim*nRows), ncol=sum(object$family %in% c("negative.binomial", "negative.binomial1")))
    if(any(object$family=="ZINB"))
      invPhis[,object$family %in% c("ZINB")] = matrix(rep(object$params$ZINB.inv.phi[object$family %in% c("ZINB")],each=nsim*nRows), ncol=sum(object$family %in% c("ZINB")))
  }
  if(any(object$family %in% c("gaussian", "gamma", "beta","ZIP","ZINB","ZIB", "ZNIB", "tweedie")))
    phis = matrix(rep(object$params$phi, each = nsim*nRows), ncol = nCols)

  if(any(object$family == "ordinal")){
    ordi_ind <- c(1:nCols)[object$family == "ordinal"]
    sims = matrix(0, nrow = nsim * nRows, ncol = length(ordi_ind))
    if(object$zeta.struc=="species"){
      for(j in ordi_ind){
        k <- sort(unique(object$y[,j]))
        for(i in 1:(nsim * nRows)){
          sims[i,ordi_ind==j] <- sample(k,1,prob=prs[!is.na(prs[,i,j]),i,j, drop=FALSE])
        }
      }
      dimnames(prs)[[3]] <- colnames(object$y)
      dimnames(prs)[[2]] <- 1:(nsim * nRows)
      prs <- prs[1,,]
    }else{
      k <- sort(unique(c(object$y[,object$family == "ordinal"])))
      for(j in ordi_ind){
        for(i in 1:(nsim * nRows)){
          sims[i,ordi_ind==j] <- sample(k,1,prob=prs[!is.na(prs[,i,j]),i,j, drop=FALSE])
        }
      }
      dimnames(prs)[[3]] <- colnames(object$y)
      dimnames(prs)[[2]] <- 1:(nsim * nRows)
      prs <- prs[1,,]
    }
    
    
  }
  if(any(object$family=="ZNIB")){
    phis2 = matrix(rep(object$params$ZINB.phi, each = nsim*nRows), ncol = nCols)
    
    phis1 = phis/(1+phis + phis2);
    phis2 = phis2/(1+phis + phis2);
    phis3 = phis1+phis2
  }
  
  newDat <- matrix(NA,nsim*nRows, nCols)
  # Need to be fixed, now works only if same family
  for(fam in unique(object$family)){
    prsfam = prs[,object$family == fam, drop=FALSE]
    nTotfam = nsim*nRows*sum(object$family == fam) # total number of values to generate for given family fam
    newDat[,object$family == fam] = matrix(
                switch(fam, #object$family[1], 
                  "binomial"=rbinom(nTotfam, size = rep(object$Ntrials,each=nsim*nRows), prob = prsfam),
                  "poisson" = rpois(nTotfam, prsfam),
                  "negative.binomial" = rnbinom(nTotfam, size = invPhis[,object$family == fam], mu = prsfam),
                  "negative.binomial1" = rnbinom(nTotfam, size = prsfam/(invPhis[,object$family == fam]), mu = prsfam),
                  "gaussian" = rnorm(nTotfam, mean = prsfam, sd = phis[,object$family == fam]),
                  "gamma" = rgamma(nTotfam, shape = phis[,object$family == fam], scale = prsfam/(phis[,object$family == fam])),
                  "exponential" = rexp(nTotfam, rate = 1/prsfam),
                  "tweedie" = fishMod::rTweedie(nTotfam, mu = c(prsfam), phi = c(phis[,object$family == fam]), p = object$Power),
                  "ordinal" = sims,
                  "beta" = rbeta(nTotfam, shape1 = phis[,object$family == fam]*prsfam, shape2 = phis[,object$family == fam]*(1-prsfam)),
                  "ZIP" = ifelse(rbinom(nTotfam, size = 1, prob = phis[,object$family == fam]) > 0, 0, rpois(nTotfam, lambda = prsfam)),
                  "ZINB" = ifelse(rbinom(nTotfam, size = 1, prob = phis[,object$family == fam]) > 0, 0, rnbinom(nTotfam, size = invPhis[,object$family == fam], mu = prsfam)),
                  "ZIB" = ifelse(rbinom(nTotfam, size = 1, prob = phis[,object$family == fam]) > 0, 0, rbinom(nTotfam, size = rep(object$Ntrials[object$family == fam],each=nsim*nRows), prob = prsfam)),
                  "ZNIB" = {
                    z <- mapply(function(prob) sample(3, 1, prob = prob, replace = TRUE), prob = as.list(data.frame(rbind(c(phis1[,object$family == fam]),c(phis2[,object$family == fam]),c(1-phis3[,object$family == fam])))))
                    zeros <- rep(0, nTotfam)
                    Ntrials <- rep(object$Ntrials[object$family == fam], each = nsim * nRows)
                    sim1 <- rbinom(nTotfam, size = Ntrials, prob = prsfam)
                    (z==1)*zeros + (z==2)*Ntrials + (z==3)*sim1
                  },
                  stop(gettextf("family '%s' not implemented ", object$family), domain = NA))
                ,ncol=sum(object$family == fam))
  # newDat = as.data.frame(matrix(newDat,ncol=nCols))
  }
  
  # reformat as data frame with the appropriate labels
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