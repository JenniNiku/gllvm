#' @title Log-likelihood of gllvm
#' @description  Extracts Log-likelihood from 'gllvm' objects.
#'
#' @param object   an object of class 'gllvm'.
#' @param ...	not used.
#'
#' @author David I. Warton, Jenni Niku
#'
#' @examples
#' \dontrun{
#'## Load a dataset from the mvabund package
#'data(antTraits, package = "mvabund")
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# log-Likelihood:
#'logLik(fit)
#'}
#'@aliases logLik logLik.gllvm
#'@method logLik gllvm
#'@importFrom stats logLik
#'
#'@export
#'@export logLik.gllvm
#'
logLik.gllvm <- function(object, ...)
{
  if(is.null(object$randomB)) object$randomB = FALSE
  if(is.null(object$num.lvcor)) object$num.lvcor = 0
  logL <- object$logL
  if(is.finite(logL)){
  if (!is.null(object$params$inv.phi)) {
    object$params$inv.phi <- NULL
  }
    if(object$family=="ordinal"){
      if(object$zeta.struc=="species")object$params$zeta <-object$params$zeta[,-1]
      if(object$zeta.struc=="common")object$params$zeta <-object$params$zeta[-1]
    }
  
    # backward compatibility
    
    if(is.null(object$params$row.params.random) && !inherits(object$row.eff, "formula") && object$row.eff == "random")object$params$row.params.random <- object$params$row.params
    if(is.null(object$beta0com))object$beta0com<-FALSE 
    if(is.null(object$col.eff$col.eff))object$col.eff$col.eff <- FALSE
    
    # end backward compatibility
    
  if (!is.null(object$params$row.params.random))
    object$params$row.params.random <- NULL
  if(object$beta0com) object$params$beta0 <- 1
  if (!is.null(object$randomX) || object$col.eff$col.eff == "random"){
    object$params$Br <- NULL
    object$params$sigmaB <- object$params$sigmaB[lower.tri(object$params$sigmaB, diag = TRUE)]
    object$params$sigmaB <- unique(object$params$sigmaB[object$params$sigmaB!=0])
  }
  if(object$randomB!=FALSE){
    object$params$LvXcoef <- NULL
  }else if(object$randomB==F&(object$num.RR+object$num.lv.c)>0){
    #correct nr. df. given orthogonality constraints
    object$params$LvXcoef[upper.tri(object$params$LvXcoef)] <- NA
  }
  if(object$quadratic=="LV"){
    object$params$theta[-1,-c(1:(object$num.lv+object$num.RR+object$num.lv.c))]<-NA
  }
  if(object$num.lv>0|object$num.RR>0|object$num.lv.c>0){
    object$params$theta[object$params$theta==0|object$params$theta==1]<-NA  
  }
  attributes(logL)$df <- length(unlist(object$params)[!is.na(unlist(object$params))])
  if(length(unique(object$disp.group))<ncol(object$y)){
    attributes(logL)$df <- attributes(logL)$df + length(unique(object$disp.group)) - ncol(object$y)
  }
  attributes(logL)$nobs <- prod(dim(object$y))
  }
  class(logL) <- "logLik"
  return(logL)
}
