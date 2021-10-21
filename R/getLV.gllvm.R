#' @title Extract latent variables
#' @description  Extract latent variables from gllvm object.
#' 
#' @param object an object of class 'gllvm'.
#' @param type type of latent variable scores to retrieve from a gllvm object. For models with unconstrained latent variables, defaults to "unconstrained". For models with constrained latent variables, defaults to "constrained". A third option is "scaled", which returns latent variables multiplied with scale parameter of the loadings. Alternatively, "LC" returns linear combination scores without residual error.
#' @param ... not used
#' 
#' @details
#' Function retrieves the site scores for a GLLVM. Each type corresponds to a separate term of the model. For a GLLVM with unconstrained latent variables
#' the default is "unconstrained". The type "scaled" returns site scores scaled with the diagonal of the species loadings.
#' 
#' For GLLVMs with constrained latent variables, "constrained" returns the complete site scores, due to both fixed- and latent effects,
#' where the latent effect is always scaled by the diagonal of the species loadings so that it can be small relative to the fixed-effects.
#' 
#' Here, "scaled" returns the same as for a GLLVM with unconstrained latent variables, namely the latent effect scaled with the diagonal
#' of the species loadings. Type "LC" (linear combination scores) returns the site scores only due to fixed-effects.
#' 
#' If both unconstrained and constrained latent variables are included in the model, type "LC" returns LC scores for constrained latent variables
#' but "scaled" scores for unconstrained latent variables.
#'@seealso  \code{\link{getLV.gllvm}}.

#'@aliases getLV getLV.gllvm
#'@method getLV gllvm
#'@export
#'@export getLV.gllvm

getLV.gllvm <- function(object, type = NULL, ...)
{
  if(!is.null(type)){
    if(!type%in%c("unconstrained","constrained","scaled","LC")){
    stop("Type should be one of: unconstrained, constrained, scaled or LC")
    }
    if(type=="constrained"&object$num.lv.c==0&object$num.RR==0){
      stop("Cannot retrieve constrained latent variables for an unconstrained latent variable model.")
    }else if(type=="scaled"&object$num.lv.c==0&object$num.lv==0){
      stop("Cannot retrieve scaled latent variables for a model without latent effect.")
    }else if(type=="unconstrained"&object$num.lv.c==0&object$num.lv==0){
      stop("Cannot retrieve unconstrained latent variables for a model without latent effect.")
    }
  }
  
  if(is.null(type)){
  if(object$num.lv.c==0&object$num.RR==0){
    type = "unconstrained"
  }else{
    if((object$num.lv.c+object$num.lv)>0){
      type = "constrained"
    }else{
      type = "LC"
    }
    
  }
  } 

  
  n <- nrow(object$y)
  if(type == "unconstrained"){
    lvs <- object$lvs
  }else if(type=="scaled"){
    lvs <- t(t(object$lvs)*object$params$sigma.lv)
  }else if(type=="constrained"){
      lvs <- t(t(object$lvs)*object$params$sigma.lv)
      lvs[,1:object$num.lv.c] <- lvs[,1:object$num.lv.c]+object$lv.X%*%object$params$LvXcoef[,1:object$num.lv.c,drop=F]
  }else if(type=="LC"&object$num.lv.c>0){
    lvs <-object$lv.X%*%object$params$LvXcoef[,1:object$num.lv.c,drop=F]
  }
  
  if(object$num.RR>0&object$num.lv>0&object$num.lv.c>0){
    lvs <- cbind(lvs[,1:object$num.lv.c,drop=F],object$lv.X%*%object$params$LvXcoef[,-c(1:object$num.lv.c),drop=F], lvs[,-c(1:object$num.lv.c),drop=F])
  }else if(object$num.RR>0&object$num.lv>0&object$num.lv.c==0){
    lvs <- cbind(object$lv.X%*%object$params$LvXcoef[,-c(1:object$num.lv.c),drop=F], lvs)
  }else if(object$num.RR>0&object$num.lv==0&object$num.lv.c>0){
    lvs <- cbind(lvs[,1:object$num.lv.c,drop=F],object$lv.X%*%object$params$LvXcoef[,-c(1:object$num.lv.c),drop=F])
  }else if(object$num.lv.c>0&object$num.lv>0&type=="LC"){
    lvs <- cbind(lvs, t(t(object$lvs[,-c(1:object$num.lv.c)])*object$params$sigma.lv))
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