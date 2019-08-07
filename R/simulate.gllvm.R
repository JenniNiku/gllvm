#' @title Simulate data from gllvm fit
#' @description Generate new data Using the fitted values of the parameters
#'
#' @param object an object of class 'gllvm'.
#' @param nsim an optional positive integer specifying the number of simulated datasets. Defaults to 1.
#' @param seed an optional integer to set seed number, passed to set.seed. Defaults to a random seed number.
#' @param ... not used.
#'
#' @details
#' simulate function for gllvm objects.  Note this simulates marginally over LVs.
#' an option is to add a conditional argument, which would fix the LVs.
#' David Warton
#' 
#' @return A matrix containing generated data.
#' @author David Warton, Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @examples
#'  \donttest{
#'# Load a dataset from the mvabund package
#'data(antTraits)
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

simulate.gllvm = function (object, nsim = 1, seed = NULL, ...) 
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
  
  nRows = dim(object$lvs)[1]
  nCols = dim(object$params$theta)[1]
  # generate new latent variables
  lvsNew = matrix(rnorm(nsim*nRows*object$num.lv),ncol=object$num.lv)
  if(is.null(object$X)) 
  {
    prs = predict.gllvm(object,newLV = lvsNew,type="response")
  }
  else
    prs = predict.gllvm(object,newX=object$X[rep(1:nRows,nsim),], newLV = lvsNew,type="response")
    
  # generate new data
  nTot = nsim*nRows*nCols # total number of values to generate
  if(object$family=="negative.binomial")
    invPhis = matrix(rep(object$params$inv.phi,each=nsim*nRows), ncol=nCols)
  if(object$family=="tweedie")
    phis = matrix(rep(object$params$phi, each = nsim*nRows), ncol = nCols)
  if(object$family=="gaussian")
    phis = matrix(rep(object$params$phi, each = nsim*nRows), ncol = nCols)
  newDat = switch(object$family, "binomial"=rbinom(nTot, size = 1, prob = prs),
                  "poisson" = rpois(nTot, prs),
                  "negative.binomial" = rnbinom(nTot, size = invPhis, mu = prs),
                  "gaussian" = rnorm(nTot, mean = prs, sd = phis),
                  "tweedie" = tweedie::rtweedie(nTot, mu = prs, phi = phis, power = object$Power),
                  stop(gettextf("family '%s' not implemented ", object$family), domain = NA))
  # reformat as data frame with the appropriate labels
  newDat = as.data.frame(matrix(newDat,ncol=nCols))
  dimnames(newDat)=dimnames(prs)
  return(newDat)
}


#' @export simulate
simulate <- function(object, ...)
{
  UseMethod(generic = "simulate")
}
