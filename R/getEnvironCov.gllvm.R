#' @title Extract species covariances due to environmental random effects from gllvm object
#' @description  Calculates the species environment covariance matrix for a gllvm object.
#'
#' @param object an object of class 'gllvm'.
#' @param x (optional) vector of covariate values to calculate the covariance for. Defaults to a vector of 1s. If both 'randomX' and random species effects are present in the model this should be a list of length two.
#' @param relative defaults to TRUE. For \code{\link{getEnvironCor}}, if there are residual latent effects in the model (such as due to num.lv or num.lv.c), this scales the covariance matrix by the sum of the environmental and residual variances instead, so that the diagonal of the correlation matrix is smaller than one. This matrix then represents the contribution of the environment to the overall species associations in the model.
#' @param ... arguments passed on to \code{\link{getResidualCov.gllvm}}
#' 
#' @return Function returns the following components:
#'  \item{cov}{species covariances due to covariates}
#'  \item{trace.randomB}{ trace of the covariance matrix due to random canonical coefficients}
#'  \item{trace.randomB.quad}{ trace of the covariance matrix components due to squared model terms}
#'  \item{trace.col.eff}{ trace of the covariance matrix due to random column (species) effects}
#' @details
#' Species covariance matrix due to the environment is calculated.
#' 
#' Covariances due to the covariates can only be calculated when random effects are included in the model, and are thus limited to reduced rank models (those including constrained and concurrent ordinations) fitted with random slopes, models fitted with random effects via the formula interface, or the fourth corner model fitted with random slopes.
#' For full rank models with random slopes, i.e., with the formula interface or the fourth corner model, the covariances of species are formulated as:
#' 
#'  \deqn{\Sigma_e = kronecker(C\rho + (1-\rho)I_p, R),}
#'  
#' where \eqn{C} is a correlation matrix for the columns in the response (e.g., a Phylogenetic matrix), \eqn{\rho} the signal parameter, and R the covariance matrix for the random effects. Here, \deqn{I = kronecker(I_p, x)}, with x a vector of covariate values for each of the random effects, which defaults to a vector of 1s.
#' when there are covariate-specific phylogenetic signal parameters in the model, this is instead:
#'
#' \deqn{\Sigma_e = kronecker(x_i', I_m)*bdiag(L_k)*kronecker(\Sigma_r, I_m)*bdiag(L_k')*kronecker(x_i, I_m),}
#' 
#' where \eqn{bdiag(L_k)} is a block-diagonal lower triangular matrix, and each \eqn{L_k} the lower triangular matrix of the covariance matrix for each covariate.
#'
#' For reduced rank models, the covariance is separately defined for the different variance structures of the canonical coefficients in the package. With LV-specific variances and without correlations, we have:
#' 
#'  \deqn{\Sigma_e = \Theta*S*\Theta',}
#' 
#' where \eqn{\Theta} is the matrix of loadings, and S the (diagonal) covariance matrix for the canonical coefficients. With predictor-specific variances, we instead have:
#' 
#' \deqn{\Sigma_e = \sum^K_{k=1} \Theta(I_d*\sigma_k^2)\Theta',}
#'
#' with I_d an identity matrix for the number of constrained and informed latent variables, and \eqn{\sigma_k^2} the variance per predictor for the canonical coefficients. When correlations are included, we have:
#' 
#'  \deqn{\Sigma_e = kronecker(x_i', I_m)kronecker(\Sigma,\Theta\Theta')kronecker(x_i, I_m),}
#'
#' and similarly for LV-specific variances with correlations. Expressions for the quadratic models in the package are determined similarly but not documented here for brevity.
#' 
#' @author Bert van der Veen
#'
#'@seealso  \code{\link{getEnvironCor}} ,\code{\link{getResidualCov.gllvm}}, \code{\link{getResidualCor.gllvm}},.
#' @examples
#' \dontrun{
#'# Example with the spider dataset
#'data(eSpider)
#'y = eSpider$abund[eSpider$nonNA,]
#'X = eSpider$X[eSpider$nonNA,]
#'fit <- gllvm(y, X = scale(X), num.RR = 2, 
#'             randomB = "P", family = "negative.binomial")
#'envcov <- getEnvironCov(fit)
#'envcov$trace.randomB
#'# As proportion of variance in the model
#'envcov$trace.randomB/sum(envcov$trace.randomB)
#'}
#'@aliases getEnvironCov getEnvironCor getEnvironCor.gllvm
#'@method getEnvironCov gllvm
#'@export
#'@export getEnvironCov.gllvm

getEnvironCov.gllvm <- function(object, x = NULL, ...){
  
  if(isFALSE(object$randomB) & isFALSE(object$col.eff$col.eff)& is.null(object$randomX)){
    stop("Canot calculate correlations without random effects for covariates in the model.")
  }
x.list <- NULL

if(is.list(x) && length(x) == 2){
  x.list <- x
  x <- x.list[[1]]
}else if(object$col.eff$col.eff=="random" && !isFALSE(object$randomB) && is.null(x)){
  x.list <- list(rep(1,ncol(object$col.eff$spdr)), rep(1, ncol(object$lv.X.design)))
  x <- x.list[[1]]
}else if(object$col.eff$col.eff=="random" && !isFALSE(object$randomB)){
 stop("x should either be NULL, or a list of length 2 for this model.")
}
if(object$col.eff$col.eff=="random" || !is.null(object$randomX)){
  if(is.null(x)){
    x <- rep(1,nrow(object$params$Br))
  }else if(length(x)!=nrow(object$params$Br)){
    stop("Supplied 'x' of  incorrect length.")
  }
  C <- kronecker(t(x), as(diag(ncol(object$y)),"TsparseMatrix"))
  if(is.null(object$params$rho.sp)){
    cov.environ.col <- diag(c(x%*%object$params$sigmaB%*%x), ncol(object$y))
  }else if(length(object$params$rho.sp)==1){
    object$col.eff$colMat <- as.matrix(object$col.eff$colMat*object$params$rho.sp+(1-object$params$rho.sp)*diag(ncol(object$y)))
    cov.environ.col <- as.matrix(C%*%kronecker(object$col.eff$colMat,object$params$sigmaB)%*%Matrix::t(C))
  }else if(length(object$params$rho.sp)>1){
    Ls<-sapply(object$params$rho.sp,function(rho)t(chol(object$col.eff$colMat*rho+diag(1-rho,ncol(object$col.eff$colMat)))),simplify=FALSE)
    cov.environ.col = as.matrix(C%*%Matrix::bdiag(Ls)%*%kronecker(object$params$sigmaB,diag(ncol(object$col.eff$colMat)))%*%Matrix::t(Matrix::bdiag(Ls))%*%Matrix::t(C))
  }
  trace.environ.col <- sum(diag(cov.environ.col))
}
  if(!is.null(object$lv.X) && is.null(object$lv.X.design))object$lv.X.design <- object$lv.X #for backward compatibility
  if(!isFALSE(object$randomB)){
    if(!is.null(x.list)){
      x <- x.list[[2]]
    }
    if(is.null(x)){
      x <- rep(1,nrow(object$params$LvXcoef))
    }else if(length(x)!=ncol(object$lv.X.design)){
      stop("Supplied 'x' of incorrect length.")
    }
    
    if(object$randomB %in% c("P","LV")){
      # generally speaking, Beta_rr ~ MN(0, U, V)
      # with V given by species loadings and U by the (albeit present or not) covariance matrix among predictors
      # with randomB = "LV" we have LV-specific variances in V, and U is a correlation matrix
      # with randomB = "P" we have a covariance matrix U, and LV-specific variance
      U <- as(diag(ncol(object$lv.X.design)), "TsparseMatrix")
      if(!is.null(object$params$corsLvXcoef)){
        U <- object$params$corsLvXcoef
      }
      if(object$randomB == "P"){
        sigma1 <- c(1, tail(object$params$sigmaLvXcoef, object$num.lv.c+object$num.RR-1))#first LV serves as "reference"
        sigma2 <- head(object$params$sigmaLvXcoef, ncol(object$lv.X.design))
        U = diag(sigma2)%*%U%*%diag(sigma2)
      }
      if(object$randomB=="LV"){
        sigma1 <- object$params$sigmaLvXcoef
      }
      
      # U <- diag(sigma1)%*%U%*%diag(sigma1)
      # V <- object$params$theta[,1:(object$num.lv.c+object$num.RR)]%*%diag(sigma2^2)%*%t(object$params$theta[,1:(object$num.lv.c+object$num.RR)])
      # sum(diag(C%*%kronecker(V,U)%*%t(C)))
      
      C = kronecker(t(x), as(diag(ncol(object$y)),"TsparseMatrix")) # sum over the predictors
      cov.environ.randomB <- sapply(1:(object$num.lv.c+object$num.RR),function(q)as.matrix(C%*%kronecker(U, object$params$theta[,q,drop=FALSE]%*%t(object$params$theta[,q,drop=FALSE])*sigma1[q]^2)%*%t(C)),simplify=FALSE)

      trace.environ.randomB <- lapply(cov.environ.randomB, function(x)sum(diag(x)))

      if(!isFALSE(object$quadratic)){
        # add tr(D_jSigma_zD_kSigma_z) for z = B^t x
        I<-matrix(0,ncol=object$num.lv.c+object$num.RR,nrow=(object$num.lv.c+object$num.RR)*ncol(object$lv.X.design))
        for(q in 1:(object$num.lv.c+object$num.RR)){
          I[1:ncol(object$lv.X.design)+(q-1)*ncol(object$lv.X.design),q] <- x
        }
        Sigmaz <- t(I)%*%kronecker(diag(sigma1^2),U)%*%(I)
        theta <- object$params$theta[,-c(1:(object$num.lv.c+object$num.RR+object$num.lv))][,1:(object$num.lv.c+object$num.RR)]
        cov.environ.randomB.quad <- sapply(1:(object$num.lv.c+object$num.RR), function(q)2*abs(theta)[,q,drop=F]%*%t(abs(theta[,q,drop=F]))*Sigmaz[q,q],simplify=F)
        trace.environ.randomB.quad <- lapply(cov.environ.randomB.quad, function(x)sum(diag(x)))
      }
    }else if(object$randomB=="single" | object$randomB=="iid"){
      cov.environ.randomB <- sapply(1:(object$num.lv.c+object$num.RR), function(q)object$params$theta[,q,drop=F]%*%t(object$params$theta[,q,drop=F])*object$params$sigmaLvXcoef^2,simplify=F)
      cov.environ.randomB.trace <- lapply(cov.environ.randomB, function(x)sum(diag(x)))  
      if(!isFALSE(object$quadratic)){
        Sigmaz <- Reduce("+", sapply(x,function(xk)xk^2*object$params$sigmaLvXcoef,simplify=F))
        # add tr(D_jSigma_zD_kSigma_z) for z = B^t x
        theta <- object$params$theta[,-c(1:(object$num.lv.c+object$num.RR+object$num.lv))][,1:(object$num.lv.c+object$num.RR)]
        cov.environ.randomB.quad <- sapply(1:(object$num.lv.c+object$num.RR), function(q)2*abs(theta)[,q,drop=F]%*%t(abs(theta[,q,drop=F]))*Sigmaz^2,simplify=F)
        trace.environ.randomB.quad <- lapply(cov.environ.randomB.quad, function(x)sum(diag(x)))
      }
    }
  }
  covMat <- 0
  out <- list()
  if(object$col.eff$col.eff=="random" || !is.null(object$randomX)){
  covMat <- cov.environ.col
  out$trace.col.eff <- trace.environ.col
  }
  
  if(!isFALSE(object$randomB)){
  if(is.list(cov.environ.randomB)){
    covMat <- covMat + Reduce("+",cov.environ.randomB)
    out$trace.randomB <- unlist(trace.environ.randomB)
  }else{
    covMat <- cov.environ.randomB
    out$trace.randomB <- trace.environ.randomB
  }
  
    names(out$trace.randomB) <- colnames(object$params$theta[,1:(object$num.RR+object$num.lv.c),drop=F])

  
    if(!isFALSE(object$quadratic)){
      covMat <- covMat + Reduce("+", cov.environ.randomB.quad)
      out$trace.randomB.quad <- unlist(trace.environ.randomB.quad)
      names(out$trace.randomB.quad) <- colnames(object$params$theta[,1:(object$num.RR+object$num.lv.c),drop=F])
    }
  }
  out$cov <- covMat
  
  colnames(out$cov) <- row.names(out$cov) <- colnames(object$y)
  
return(out)
}

#'@export getEnvironCor
getEnvironCor <- function(object, ...)
{
  UseMethod(generic = "getEnvironCor")
}

#'@export
#'@export getEnvironCor.gllvm
getEnvironCor.gllvm <- function(object, x = NULL, relative = TRUE, ...)
{
  ResCov <- matrix(0, ncol = ncol(object$y), nrow = nrow(object$y))
  if(relative)ResCov <- getResidualCov(object, ...)$cov
  EnvCov <- getEnvironCov.gllvm(object, x)$cov
  sds <- sqrt(diag(ResCov)+diag(EnvCov))
  EnvCors <- sweep(EnvCov, 1, sds, "/")
  EnvCors <- sweep(EnvCors, 2, sds, "/")
  
  colnames(EnvCors) <- row.names(EnvCors) <- colnames(object$y)
  
  return(EnvCors)
}

#'@export getEnvironCov
getEnvironCov <- function(object, ...)
{
  UseMethod(generic = "getEnvironCov")
}
