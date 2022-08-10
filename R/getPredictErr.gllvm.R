#' @title Extract prediction errors for latent variables from gllvm object
#' @description  Calculates the prediction errors for latent variables and random effects for gllvm model.
#'
#' @param object   an object of class 'gllvm'.
#' @param CMSEP logical, if \code{TRUE} conditional mean squared errors for predictions are calculated. If \code{FALSE}, prediction errors are based on covariances of the variational distributions for \code{method ="VA"} and \code{method ="EVA"}.
#' @param cov if \code{TRUE}, return as covariances/variances of predictions. Otherwise \code{FALSE} (default) return as standard errors of predictions.
#' @param ...	 not used
#'
#' @details 
#' Calculates conditional mean squared errors for predictions.
#' If variational approximation is used, prediction errors can be based on covariances 
#' of the variational distributions, and therefore they do not take into account 
#' the uncertainty in the estimation of (fixed) parameters. 
#'
#' @return Function returns following components:
#'  \item{lvs }{prediction errors for latent variables}
#'  \item{row.effects }{prediction errors for random row effects if included}
#'
#' @author Francis K.C. Hui, Jenni Niku, David I. Warton
#'
#' @examples
#'\dontrun{
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# prediction errors for latent variables:
#'getPredictErr(fit)
#'}
#'
#'@aliases getPredictErr getPredictErr.gllvm
#'@method getPredictErr gllvm
#'@export
#'@export getPredictErr.gllvm
getPredictErr.gllvm = function(object, CMSEP = TRUE, cov = FALSE, ...)
{
  n <- nrow(object$y)
  num.lv <- object$num.lv
  num.lv.c <- object$num.lv.c
  num.RR <- object$num.RR
  out <- list()
  if(object$method == "LA"){
    if((num.lv+num.lv.c+num.RR)>0) out$lvs <- sqrt(apply(object$prediction.errors$lvs,1,diag))
    if(object$row.eff == "random") out$row.effects <- sqrt(abs(object$prediction.errors$row.params))
  }
  
  if((object$method %in% c("VA", "EVA"))){
    if(object$row.eff == "random"){
      if(is.null(object$Ar)){
        if(dim(object$A)[3]>ncol(object$lvs)){
          object$Ar<-object$A[,1,1]
          object$A<-object$A[,-1,-1]
        }
      }
    }
    if(CMSEP) {
      sdb <- CMSEPf(object)
      # sdb<-sdA(object)

      if(num.RR>0){
        #variational covariances but add 0s for RRR
        A <- array(0,dim=c(n,num.lv.c+num.RR+num.lv,num.lv.c+num.RR+num.lv))
        A[,-c((num.lv.c+1):(num.lv.c+num.RR)),-c((num.lv.c+1):(num.lv.c+num.RR))] <- object$A
      } else if((object$num.lvcor > 1) && (object$Lambda.struc %in% c("diag_unstructured","NN_unstructured","unstructured_unstructured"))) {
          A<-array(diag(object$A[,,1]), dim = c(nrow(object$A[,,1]), object$num.lvcor,object$num.lvcor))
          for (i in 1:dim(A)[1]) {
            A[i,,]<-A[i,,]*object$AQ
          }
      } else if((object$num.lvcor > 0) & (object$corP$cstruc !="diag")) {
        A<-array(0, dim = c(nrow(object$A[,,1]), object$num.lvcor,object$num.lvcor))
        for (i in 1:object$num.lvcor) {
          A[,i,i]<- diag(object$A[,,i])
        }
      } else{A<-object$A}
      if(object$row.eff == "random"){
        object$Ar<-sdb$Ar+object$Ar
      }
      if((num.lv+num.lv.c)>0){ object$A<-sdb$A+A} else{object$A <- sdb$A}
      # if(!is.null(object$randomX)) object$Ab<-sdb$Ab+object$Ab
    }
      r=0
      if(cov){
        if(object$row.eff=="random"){
          # r=1
          out$row.effects <- (object$Ar)
        }
        if(length(dim(object$A))==2){
          out$lvs <- (object$A[,1:(num.lv+num.lv.c+num.RR)+r])
        } else {
          if((num.lv+num.lv.c+num.RR) ==1) {
            out$lvs <- (as.matrix(object$A[,1:(num.lv+num.lv.c+num.RR)+r,1:(num.lv+num.lv.c+num.RR)+r]))
          } else {
            out$lvs <- (object$A[,1:(num.lv+num.lv.c+num.RR)+r,1:(num.lv+num.lv.c+num.RR)+r])
          }
        }
      } else {
        if(object$row.eff=="random"){
          # r=1
          out$row.effects <- sqrt(object$Ar)
          }
        if(length(dim(object$A))==2){
          out$lvs <- sqrt(object$A[,1:(num.lv+num.lv.c+num.RR)+r])
        } else {
          if((num.lv+num.lv.c+num.RR) ==1) {
            out$lvs <- sqrt(as.matrix(object$A[,1:(num.lv+num.lv.c+num.RR)+r,1:(num.lv+num.lv.c+num.RR)+r]))
          } else {
            out$lvs <- sqrt(apply((object$A[,1:(num.lv+num.lv.c+num.RR)+r,1:(num.lv+num.lv.c+num.RR)+r]),1,diag))
          }
        }
      }
  }
  if(((num.lv+num.lv.c) > 1) & is.matrix(out$lvs)) out$lvs <- t(out$lvs)
  return(out)
}

#'@export getPredictErr
getPredictErr <- function(object, ...)
{
  UseMethod(generic = "getPredictErr")
}
