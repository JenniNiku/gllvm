#' @title Extract latent variables
#' @description  Extract latent variables from gllvm object.
#' 
#' @param object an object of class 'gllvm'.
#' @param type type of latent variable scores to retrieve from a gllvm object. For models with unconstrained latent variables, defaults to "residual". For models with constrained latent variables, defaults to conditional. Alternatively, "marginal" returns linear combination scores without residual error.
#' @param ... not used
#' 
#' @details
#' Function retrieves the site scores for a GLLVM. Each type corresponds to a separate term of the model. For a GLLVM with unconstrained latent variables
#' the default is "residual".
#' 
#' For GLLVMs with constrained latent variables, "conditional" returns the complete site scores, due to both fixed- and latent effects,
#' where the latent effect is always scaled by the diagonal of the species loadings so that it can be small relative to the fixed-effects.
#' 
#' Type "marginal"  returns linear combination scores, i.e. the site scores only due to fixed-effects.
#' 
#' If both unconstrained and constrained latent variables are included in the model, type "marginal" returns LC scores for constrained latent variables
#' but "residual" scores for unconstrained latent variables.
#'@aliases getLV getLV.gllvm
#'@method getLV gllvm
#'@export
#'@export getLV.gllvm

getLV.gllvm <- function(object, type = NULL, ...)
{
  if(!is.null(type)){
    if(!type%in%c("residual","conditional","marginal")){
    stop("Type should be one of: residual, conditional, or marginal")
    }
    if(type=="conditional"&object$num.lv.c==0&object$num.RR==0){
      stop("Cannot retrieve constrained latent variables for an unconstrained latent variable model.")
    }else if(type=="residual"&object$num.lv.c==0&object$num.lv==0){
      stop("Cannot retrieve residual latent variables for a model without latent effect.")
    }
  }
  
  if(is.null(type)){
  if(object$num.lv.c==0&object$num.RR==0){
    type = "residual"
  }else{
    if((object$num.lv.c+object$num.lv)>0){
      type = "conditional"
    }else{
      type = "marginal"
    }
    
  }
  } 

  
  n <- nrow(object$y)
  if(type == "residual"){
    lvs <- object$lvs
  }else if(type=="conditional"){
      lvs <- object$lvs
      lvs[,1:object$num.lv.c] <- t(t(lvs[,1:object$num.lv.c,drop=F])*object$params$sigma.lv[1:object$num.lv.c])
      lvs[,1:object$num.lv.c] <- lvs[,1:object$num.lv.c]+object$lv.X%*%object$params$LvXcoef[,1:object$num.lv.c,drop=F]
  }else if(type=="marginal"&object$num.lv.c>0){
    lvs <-object$lv.X%*%object$params$LvXcoef[,1:object$num.lv.c,drop=F]
  }
  
  if(object$num.RR>0&object$num.lv>0&object$num.lv.c>0){
    lvs <- cbind(lvs[,1:object$num.lv.c,drop=F],object$lv.X%*%object$params$LvXcoef[,-c(1:object$num.lv.c),drop=F], lvs[,-c(1:object$num.lv.c),drop=F])
  }else if(object$num.RR>0&object$num.lv>0&object$num.lv.c==0){
    lvs <- cbind(object$lv.X%*%object$params$LvXcoef[,-c(1:object$num.lv.c),drop=F], lvs)
  }else if(object$num.RR>0&object$num.lv==0&object$num.lv.c>0){
    lvs <- cbind(lvs[,1:object$num.lv.c,drop=F],object$lv.X%*%object$params$LvXcoef[,-c(1:object$num.lv.c),drop=F])
  }else if(object$num.lv.c>0&object$num.lv>0&type=="marginal"){
    lvs <- cbind(lvs, object$lvs[,-c(1:object$num.lv.c)])
  }
  
  
  if((object$num.lv.c+object$num.RR)>0){
    if((object$num.lv.c+object$num.RR)>0&object$num.lv>0)colnames(lvs)<-c(paste("CLV",1:(object$num.lv.c+object$num.RR),sep=""),paste("LV",1:object$num.lv,sep=""))
    if((object$num.lv.c+object$num.RR)>0&object$num.lv==0)colnames(lvs)<-paste("CLV",1:(object$num.lv.c+object$num.RR),sep="")
  }
  
  return(lvs)
}

#'@export getLV
getLV <- function(object, ...)
{
  UseMethod(generic = "getLV")
}