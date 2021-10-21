#' @title Extract latent variables
#' @description  Extract latent variables from gllvm object.
#' 
#' @param object an object of class 'gllvm'.
#' @param type type of latent variable scores to retrieve from a gllvm object. For models with unconstrained latent variables, defaults to "unconstrained". For models with constrained latent variables, defaults to "constrained". A third option is "scaled", which returns latent variables multiplied with scale parameter of the loadings. Alternatively, "LC" returns linear combination scores without residual error.
#' @param ... not used
#' 
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
      stop("Cannot retrieve scaled latent variables for a model without scale parameter.")
    }else if(type=="unconstrained"&object$num.lv.c==0&object$num.lv==0){
      stop("Cannot retrieve unconstrained latent variables, for a model without latent effect.")
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
  }else if(type=="LC"){
    lvs <- object$lv.X%*%object$params$LvXcoef
  }
  
  if(type!="LC"){
  if(object$num.RR>0&object$num.lv>0&object$num.lv.c>0){
    lvs <- cbind(lvs[,1:object$num.lv.c,drop=F],object$lv.X%*%object$params$LvXcoef[,-c(1:object$num.lv.c),drop=F], lvs[,-c(1:object$num.lv.c),drop=F])
  }else if(object$num.RR>0&object$num.lv>0&object$num.lv.c==0){
    lvs <- cbind(object$lv.X%*%object$params$LvXcoef[,-c(1:object$num.lv.c),drop=F], lvs)
  }else if(object$num.RR>0&object$num.lv==0&object$num.lv.c>0){
    lvs <- cbind(lvs[,1:object$num.lv.c,drop=F],object$lv.X%*%object$params$LvXcoef[,-c(1:object$num.lv.c),drop=F])
  }
  }
  
  if((object$num.lv.c+object$num.RR)>0){
    if((object$num.lv.c+object$num.RR)>0&object$num.lv>0&type!="LC")colnames(lvs)<-c(paste("CLV",1:(object$num.lv.c+object$num.RR),sep=""),paste("LV",1:object$num.lv,sep=""))
    if((object$num.lv.c+object$num.RR)>0&object$num.lv==0|type=="LC")colnames(lvs)<-paste("CLV",1:(object$num.lv.c+object$num.RR),sep="")
  }
  
  return(lvs)
}

#'@export getLV
getLV <- function(object, ...)
{
  UseMethod(generic = "getLV")
}