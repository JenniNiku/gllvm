#' @title Functions to extract ecological quantities from GLLVM where species are a quadratic function of the latnt variables.
#' @description Extracts species optima, tolerances or gradient lengths, potentially with standard errors (derived with the delta method).
#'
#' @param object   an object of class 'gllvm'.
#' @param sd.errors logical. If \code{TRUE}, also returns standard errors.
#'
#' 
#' @author Bert van der Veen
#'
#'@aliases optima optima.gllvm.quadratic tolerances tolerances.gllvm.quadratic gradient.length gradient.length.gllvm.quadratic
#'@export
#'@export optima tolerances gradient.length
optima.gllvm <- function(object,sd.errors = T) {
  if(!inherits(object,"gllvm.quadratic")){
    stop("Optima can only be extracted for a GLLVM where species are a quadratic function of the latent variables.")
  }
    num.lv <- object$num.lv
    p <- ncol(object$y)
    opt<-object$params$theta[,1:num.lv]/(2*object$params$theta[,-c(1:num.lv)])
    if(is.null(object$sd)){
      cat("Standard errors not present in model, calculating...")
      object$Hess<-se.gllvm(object)$Hess
    }
    V <- object$Hess$cov.mat.mod
    idx <- names(object$TMBfn$par[object$Hess$incl])%in%c("lambda","lambda2")
    V <- V[idx,idx,drop=F]
    colnames(V) <- row.names(V) <- names(object$TMBfn$par[object$Hess$incl])[idx]
    
    
      for (q in 1:num.lv) {
        for (j in 1:p) {
          if (q > 1 & j < q) {
            # add zeros where necessary
            V <- cbind(cbind(V[, 1:c(p * q - 1)], 0), V[, (p * q ):ncol(V)])
            V <- rbind(rbind(V[1:c(p * q - 1), ], 0), V[(p * q):nrow(V), ])
          }
        }
      }
    
    colnames(V)[colnames(V)==""]<-"lambda"
      
      opt.sd<-matrix(0,nrow=p,ncol=num.lv)
      
      for (j in 1:p) {
          idx <- colnames(V) == "lambda" | colnames(V) == "lambda2"
          V.theta <- V[idx, idx]
          if (object$quadratic == "LV") {
            idx <- c((c(1:num.lv) - 1) * p + j, p * num.lv + 1:num.lv)
          } else {
            idx <- c((c(1:num.lv) - 1) * p + j, ((1 + p * num.lv + (num.lv * (j - 1))):(p * num.lv + (num.lv * (j - 1)) + num.lv)))
          }
      
      
      V.theta2 <- V.theta[idx, idx]
      for (i in 1:num.lv) {
        du <- c((2 * object$params$theta[j, num.lv + i, drop = F])^-1, 2 * (object$params$theta[j, i, drop = F] / (2 * object$params$theta[j, num.lv + i, drop = F])^2))
        opt.sd[j, i] <- sqrt(abs(t(du) %*% V.theta2[c(i, num.lv + i), c(i, num.lv + i)] %*% du))
      }
      }
    
      if(num.lv>1)colnames(opt) <- colnames(opt.sd) <- paste("LV",1:num.lv,sep="")
      if(num.lv>1){
        row.names(opt) <- row.names(opt.sd) <- colnames(object$y)
      }else{
        names(opt) <- names(opt.sd) <- colnames(object$y)
      }
      
 return(list(optima=opt,sd=opt.sd)) 
}

tolerances.gllvm <- function(object,sd.errors = T) {
  if(!inherits(object,"gllvm.quadratic")){
    stop("Tolerances can only be extracted for a GLLVM where species are a quadratic function of the latent variables.")
  }
  num.lv <- object$num.lv
  p <- ncol(object$y)
  tol<-1/sqrt(-2*object$params$theta[,-c(1:num.lv)])
  if(is.null(object$sd)){
    cat("Standard errors not present in model, calculating...")
    object$Hess<-se.gllvm(object)$Hess
  }
  V <- object$Hess$cov.mat.mod
  idx <- names(object$TMBfn$par[object$Hess$incl])%in%c("lambda","lambda2")
  V <- V[idx,idx,drop=F]
  colnames(V) <- row.names(V) <- names(object$TMBfn$par[object$Hess$incl])[idx]
  
  
  for (q in 1:num.lv) {
    for (j in 1:p) {
      if (q > 1 & j < q) {
        # add zeros where necessary
        V <- cbind(cbind(V[, 1:c(p * q - 1)], 0), V[, (p * q ):ncol(V)])
        V <- rbind(rbind(V[1:c(p * q - 1), ], 0), V[(p * q):nrow(V), ])
      }
    }
  }
  
  colnames(V)[colnames(V)==""]<-"lambda"
  
  tol.sd<-matrix(0,nrow=p,ncol=num.lv)
  
  for (j in 1:p) {
    idx <- colnames(V) == "lambda" | colnames(V) == "lambda2"
    V.theta <- V[idx, idx]
    if (object$quadratic == "LV") {
      idx <- c((c(1:num.lv) - 1) * p + j, p * num.lv + 1:num.lv)
    } else {
      idx <- c((c(1:num.lv) - 1) * p + j, ((1 + p * num.lv + (num.lv * (j - 1))):(p * num.lv + (num.lv * (j - 1)) + num.lv)))
    }
    
    
    V.theta2 <- V.theta[idx, idx]
    for (i in 1:num.lv) {
      dt <- 1 / (2 * object$params$theta[, -c(1:num.lv),drop=F][j, i] * (sqrt(-2 * object$params$theta[, -c(1:num.lv),drop=F][j, i]))) # need to be calculated with covariance of gamma3 if gamma2>0..that also requires subtracting theta3 from theta2
      tol.sd[j, i] <- sqrt(abs(V.theta2[-c(1:num.lv), -c(1:num.lv),drop=F][i, i] * dt^2))
    }
  }
  if(num.lv>1)colnames(tol) <- colnames(tol.sd) <- paste("LV",1:num.lv,sep="")
  if(num.lv>1){
    row.names(tol) <- row.names(tol.sd) <- colnames(object$y)
  }else{
    names(tol) <- names(tol.sd) <- colnames(object$y)
  }
  return(list(tolerances=tol,sd=tol.sd)) 
}

gradient.length.gllvm <- function(object,sd.errors = T) {
  if(!inherits(object,"gllvm.quadratic")){
    stop("Gradient length can only be extracted for a GLLVM where species are a quadratic function of the latent variables.")
  }
  num.lv <- object$num.lv
  p <- ncol(object$y)
  quadratic <- object$quadratic
  tol<-1/sqrt(-2*object$params$theta[,-c(1:num.lv),drop=F])
  
if(is.null(object$sd)){
    cat("Standard errors not present in model, calculating...")
    object$Hess<-se.gllvm(object)$Hess
  }
  V <- object$Hess$cov.mat.mod
  idx <- names(object$TMBfn$par[object$Hess$incl])%in%c("lambda2")
  V <- V[idx,idx,drop=F]
  colnames(V) <- row.names(V) <- names(object$TMBfn$par[object$Hess$incl])[idx]
  
  gradSD <- NULL
  for(q in 1:num.lv){
    #covariance matrix for all quadratic coefficients
    if(quadratic==TRUE){
      covmat <- V[colnames(V)=="lambda2",colnames(V)=="lambda2"][(p*(q-1)+1):(p*q),(p*(q-1)+1):(p*q)] 
    }else{
      covmat <- V[colnames(V)=="lambda2",colnames(V)=="lambda2"][q,q]
    }
    #selec the parameters that represent the median, middle 2 if even otherwise middle 1
    if(p%%2==0){
      #even
      idx<-order(1/sqrt(2*-object$params$theta[,num.lv+q]),decreasing=T)[(p/2):((p/2)+1)]
    }else{
      #odd
      idx<-order(1/sqrt(2*-object$params$theta[,num.lv+q]),decreasing=T)[p-round(p/2)]
    }
    if(quadratic=="LV"|length(idx)==1){
      if(quadratic=="LV"){
        #only a function of one parameter, so this is derivative of gradient length (since there is no median)
        grad<-2*abs(object$params$theta[1,num.lv+q])^-0.5
        gradSD <- c(gradSD,grad^2*covmat)
      }else{
        grad<-2*abs(object$params$theta[idx,num.lv+q])^-0.5#Only needs the right idx, no transformations
        gradSD <- c(gradSD,grad^2*covmat[idx,idx])
      }
      
    }else{
      #this one is a little more, due to the presence of the median on the tolerance scale.
      grad<-4*sum(abs(object$params$theta[idx,num.lv+q])^-0.5)^-2*abs(object$params$theta[idx,num.lv+q])^-1.5#4*sqrt(0.5)abs(out$params$theta[idx,num.lv+i])^3 
      gradSD <- c(gradSD,grad%*%covmat[idx,idx]%*%grad)
    }
    
  }
  grad.length.sd <- sqrt(gradSD)
  
  if(quadratic==TRUE){
    grad.length <- 4*sqrt(0.5)*1/apply(tol,2,median)
  }else if(quadratic=="LV"){
    grad.length <- 4*sqrt(0.5)*1/tol[1,]
  }
  
  names(grad.length) <- names(grad.length.sd) <- paste("LV",1:num.lv,sep="")
  
  return(list(gradient.length=grad.length,sd=grad.length.sd)) 
}


#'@export optima
optima <- function(object, ...)
{
  UseMethod(generic = "optima")
}

#'@export tolerances
tolerances <- function(object, ...)
{
  UseMethod(generic = "tolerances")
}
#'@export gradient.length
gradient.length <- function(object, ...)
{
  UseMethod(generic = "gradient.length")
}



