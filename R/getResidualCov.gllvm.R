#' @title Extract residual covariance matrix from gllvm object
#' @description  Calculates the residual covariance matrix for gllvm model.
#'
#' @param object an object of class 'gllvm'.
#' @param adjust The type of adjustment used for  negative binomial, binomial and normal distribution when computing residual correlation matrix. Options are 0 (no adjustment), 1 (the default adjustment) and 2 (alternative adjustment for NB distribution), see details.
#' @param x (optional) vector of covariate values to calculate the covariance for, when applicable.
#' @param ...  not used.
#'
#' @return Function returns following components:
#'  \item{cov }{residual covariance matrix}
#'  \item{trace }{trace of the residual covariance matrix, the total variance explained}
#'  \item{var.q }{trace of the residual covariance matrix per latent variable, variance explained per latent variable}
#'  \item{var.q2 }{trace of the squared term of the residual covariance matrix per latent variable, for quadratic responses. Variance explained per latent variable by the quadratic term}

#' @details 
#' Residual covariance matrix, storing information on species co-occurrence that is not explained by the environmental variables (if included), is calculated using the matrix of latent variables loadings, that is, \eqn{\Theta\Theta'}, and the dispersion parameter related to the distribution of choice, is applicable (e.g. in the case of negative-binomial distributed responses).
#' 
#' When the responses are modelled using the negative binomial distribution, the residual variances for each species must be adjusted for overdispersion. The two possible adjustment terms are \eqn{log(\phi_j + 1)} (\code{adjust = 1}) and \eqn{\psi^{(1)}(1/\phi_j)} (\code{adjust = 2}), where \eqn{\psi^{(1)}} is the trigamma function.
#' 
#' The negative binomial model can be written using different parameterizations. 
#' The residual covariance with \code{adjust = 1} can be obtained using the lognormal-Poisson parametrization, that is,
#' \deqn{Y_{ij} \sim Poisson(\mu_{ij} \lambda_j),}
#' where \eqn{\lambda_j \sim lognormal(-\sigma^2/2, \sigma^2)} and \eqn{\sigma^2 = log(\phi_j + 1)} and \eqn{log(\mu_{ij}) = \eta_{ij}}. Now \eqn{E[Y_{ij}] = \mu_{ij}} and variance \eqn{V(\mu_{ij}) = \mu_{ij} + \mu_{ij}^2 (exp(\sigma^2) - 1) = \mu_{ij} + \mu_{ij}^2 \phi_j}, which are the same as for the NB distribution.
#' Therefore, on linear predictor scale, we have the variance 
#' \deqn{V(log(\mu_{ij} \lambda_j)) = V(log\mu_{ij}) + V(log\lambda_j) = V(u_i'\theta_j) + \sigma^2 = \theta_j'\theta_j + log(\phi_j + 1).}
#' which leads to the residual covariance matrix \eqn{\Theta \Theta' + \Psi}, where \eqn{\Psi} is the diagonal matrix with \eqn{log(\phi_j + 1)} as diagonal elements (\code{adjust = 1}).
#' 
#' Or, for a GLLVM where species are a quadratic function of the latent variables, we instead have
#' \deqn{V(log(\mu_{ij} \lambda_j)) = V(log\mu_{ij}) + V(log\lambda_j) = V(u_i'\theta_j-u_i' D_j u_i) + \sigma^2 }
#' \deqn{ = \theta_j'\theta_j + 2diag(D_j)'diag(D_j)log(\phi_j + 1).}
#' which leads to the residual covariance matrix \eqn{\Theta \Theta' + 2 \Gamma_j \Gamma_j' + diag(\Phi)}, where \eqn{\Gamma_j} holds the quadratic coefficients.
#' Since the quadratic coefficients are constrained to be positive, the residual covariance in the latter case is, given the same coefficients on the linear term, equal or more positive than in the linear case.
#' 
#' The residual covariance matrix with \code{adjust = 2} can be obtained by using Poisson-Gamma parametri-zation
#' \deqn{Y_{ij} \sim Poisson(\mu_{ij} \lambda_j),}
#' where \eqn{\lambda_j \sim Gamma(1/\phi_j, 1/\phi_j)} and \eqn{\mu_{ij}} is as above. The mean and the variance are of similar form as above and we have that
#' \deqn{V(log(\mu_{ij} \lambda_j)) = V(log\mu_{ij}) + V(log\lambda_j) = \theta_j'\theta_j + \psi^{(1)}(1/\phi_j),}
#' where \eqn{\psi^{(1)}} is the trigamma function.
#' 
#' In the case of binomial distribution, the adjustment terms (\code{adjust = 1}) are 1 for probit link and \eqn{\pi^2/3} for logit link.
#' These are obtained by treating binomial model as latent variable model. Assume
#' \deqn{Y^*_{ij} = \eta_{ij} + e_{ij},}
#' where \eqn{e_{ij} \sim N(0, 1)} for probit model, and \eqn{e_{ij} \sim logistic(0, 1)} for logit model.
#' Then binary response is defined as \eqn{Y_{ij} = 1}, if \eqn{Y^*_{ij} > 0} and 0 otherwise.
#' Now we have that \eqn{\mu_{ij} = P(Y_{ij} = 1) = P(Y^*_{ij} > 0) = P(\eta_{ij} > -e_{ij}) = P(e_{ij} <= \eta_{ij})} which leads to probit and logit models.
#' On linear predictor scale we then have that
#' \deqn{V(\eta_{ij} + e_{ij}) = V(\eta_{ij}) + V(e_{ij}).}
#' For the probit model, the residual covariance matrix is then \eqn{\Theta\Theta' + I_m}, and for the logit model \eqn{\Theta\Theta' + \pi^2/3 I_m}.
#' Similarly as above, for a GLLVM where species are a quadratic function of the latent variables, the term \eqn{2\Gamma_j\Gamma_j'} is added to the residual covariance matrix.
#' 
#' For normal distribution, we can write
#' \deqn{Y_{ij} = \eta_{ij} + e_{ij},}
#' where \eqn{e_{ij} \sim N(0, \phi_j^2)} and thus we have that
#' \deqn{V(\eta_{ij} + e_{ij}) = V(\eta_{ij}) + V(e_{ij}).}
#' For the gaussian model, the residual covariance matrix is then \eqn{\Theta\Theta' + diag(\Phi^2)}.
#' 
#' 
#' @author Francis K.C. Hui, Jenni Niku, David I. Warton, Bert van der Veen
#'
#' @examples
#' \dontrun{
#'# Load a dataset from the mvabund package
#'data(antTraits, package = "mvabund")
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# residual covariance:
#'rescov <- getResidualCov(fit)
#'rescov$cov
#'# Trace of the covariance matrix
#'rescov$trace
#'# Variance explained per latent variable
#'rescov$var.q
#'}
#'@aliases getResidualCov getResidualCov.gllvm
#'@method getResidualCov gllvm
#'@export
#'@export getResidualCov.gllvm
getResidualCov.gllvm = function(object, adjust = 1, x = NULL, ...)
{
  #backward compatibility
  opts <- list(...)
  if("site.index"%in%names(opts)){
    site.index <- opts$site.index
    x <- object$lv.X.design[site.index[1],]
  }
  
  if(object$quadratic!=FALSE&&is.null(x)&object$num.lv.c>0){
    x <- rep(1,nrow(object$params$LvXcoef))
  }else if(object$quadratic!=FALSE&&!is.null(x)){
    if(length(x)!=ncol(object$lv.X.design)){
    stop("Supplied 'x' of incorrect length.")
    }
  }
    
  if(!is.null(object$lv.X) && is.null(object$lv.X.design))object$lv.X.design <- object$lv.X #for backward compatibility
  
  if((object$num.lv+object$num.lv.c)==0){
    stop("No latent variables present in model.")
  }
  
  # Remove Reduced Rank parameters if present without residual term
  if(object$num.RR>0){
    if(object$quadratic==FALSE){
    if(object$num.lv.c>0){
      object$params$theta <- object$params$theta[,-c((object$num.lv.c+1):(object$num.lv.c+object$num.RR)),drop=F]
      object$params$LvXcoef <- object$params$LvXcoef[,1:object$num.lv.c,drop=F]
    }else{
      object$params$theta <- object$params$theta[,-c(1:object$num.RR),drop=F]
    }
    }else{
        theta <- object$params$theta[,1:(object$num.lv+object$num.lv.c+object$num.RR),drop=F]
        theta2 <- object$params$theta[,-c(1:(object$num.lv+object$num.lv.c+object$num.RR)),drop=F]
        theta <- theta[,-c((object$num.lv.c+1):(object$num.lv.c+object$num.RR)),drop=F]
        theta2 <- theta2[,-c((object$num.lv.c+1):(object$num.lv.c+object$num.RR)),drop=F]
        object$params$theta <- cbind(theta,theta2)
        if(object$num.lv.c>0){
          object$params$LvXcoef <- object$params$LvXcoef[,1:object$num.lv.c,drop=F]
        }
    }
  }
  
  if((object$num.lv+object$num.lv.c)>1){
    Sigma <- diag(object$params$sigma.lv)  
  }else{
    Sigma <- object$params$sigma.lv
  } 

  ResCov <- matrix(0,ncol=ncol(object$y),nrow=ncol(object$y))
  if(inherits(object,"gllvm.quadratic")){
    ResCov <- ResCov + object$params$theta[, 1:(object$num.lv+object$num.lv.c), drop = F]%*% Sigma %*% t(object$params$theta[, 1:(object$num.lv+object$num.lv.c), drop = F] %*% Sigma) + 2 * object$params$theta[, -c(1:(object$num.lv+object$num.lv.c)), drop = F] %*% Sigma^2 %*% t(object$params$theta[, -c(1:(object$num.lv+object$num.lv.c)), drop = F] %*% Sigma^2)
    ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) object$params$sigma.lv[q]^2*object$params$theta[, q] %*% t(object$params$theta[, q]), simplify = F)
    ResCov.q2 <- sapply(1:(object$num.lv+object$num.lv.c), function(q) 2*object$params$sigma.lv[q]^4*object$params$theta[, q+(object$num.lv+object$num.lv.c)] %*% t(object$params$theta[, q+(object$num.lv+object$num.lv.c)]), simplify = F)
    if(object$num.lv.c>0 && isFALSE(object$randomB)){
      ResCov <- ResCov + Reduce("+",sapply(1:object$num.lv.c,function(q)4*c(x%*%object$params$LvXcoef[,q,drop=F]*x%*%object$params$LvXcoef[,q,drop=F])*(object$params$theta[,-c(1:(object$num.lv+object$num.lv.c)),drop=F][,q,drop=F] *object$params$sigma.lv[q]^2)%*%t(object$params$theta[,-c(1:(object$num.lv+object$num.lv.c)),drop=F][,q,drop=F])-2*c(x%*%object$params$LvXcoef[,q,drop=F]*object$params$sigma.lv[q]^2+x%*%object$params$LvXcoef[,q,drop=F]*object$params$sigma.lv[q]^2)*object$params$theta[,-c(1:(object$num.lv+object$num.lv.c)),drop=F][,q,drop=F]%*%t(object$params$theta[,1:(object$num.lv+object$num.lv.c),drop=F][,q,drop=F]),simplify=F))
      ResCov.q2[1:object$num.lv.c] <- sapply(1:object$num.lv.c, function(q)ResCov.q2[[q]]+4*c(x%*%object$params$LvXcoef[,q,drop=F]*x%*%object$params$LvXcoef[,q,drop=F])*(object$params$theta[,-c(1:(object$num.lv+object$num.lv.c)),drop=F][,q,drop=F] *object$params$sigma.lv[q]^2)%*%t(object$params$theta[,-c(1:(object$num.lv+object$num.lv.c)),drop=F][,q,drop=F]),simplify=F)#what follows is a covariance term..-2*c(x%*%object$params$LvXcoef[,q,drop=F]*object$params$sigma.lv[q]^2+x%*%object$params$LvXcoef[,q,drop=F]*object$params$sigma.lv[q]^2)*object$params$theta[,-c(1:(object$num.lv+object$num.lv.c)),drop=F][,q,drop=F]%*%t(object$params$theta[,1:(object$num.lv+object$num.lv.c),drop=F][,q,drop=F]),simplify=F)
    }
    #if(object$num.lv.c>0)ResCov <- ResCov - Reduce("+",sapply(1:object$num.lv.c,function(q)2*(c(object$lv.X.design[site.index[1],,drop=F]%*%object$params$LvXcoef[,q,drop=F])*(abs(object$params$theta[,-c(1:(object$num.lv+object$num.lv.c)),drop=F][,q,drop=F])%*%t(object$params$theta[,q,drop=F]))+c(object$lv.X.design[site.index[2],,drop=F]%*%object$params$LvXcoef[,q,drop=F])*(abs(object$params$theta[,-c(1:(object$num.lv+object$num.lv.c)),drop=F][,q,drop=F])%*%t(object$params$theta[,q,drop=F]))),simplify=F))
  }else{
    ResCov <- object$params$theta[, 1:(object$num.lv+object$num.lv.c), drop = F]%*% Sigma %*% t(object$params$theta[, 1:(object$num.lv+object$num.lv.c), drop = F] %*% Sigma)
    ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) (object$params$theta[, q, drop = F]* object$params$sigma.lv[q]) %*% t(object$params$theta[, q, drop = F] * object$params$sigma.lv[q]), simplify = F)
  }
  
  
  if(adjust > 0 && object$family %in% c("negative.binomial", "binomial", "gaussian")){
  if(object$family == "negative.binomial"){ 
    if(adjust == 1) {
        ResCov <- ResCov + diag(log(object$params$phi + 1), ncol=ncol(ResCov))
      }else if(adjust == 2){
        ResCov <- ResCov + diag(trigamma(1/object$params$phi), ncol=ncol(ResCov))
     }
    
  }
    if(object$family == "binomial"){ 
      if(object$link == "probit"){
        ResCov <- ResCov + diag(ncol(object$y))
      } 
      if(object$link == "logit"){
        ResCov <- ResCov + diag(ncol(object$y))*pi^2/3
      } 
    }
    if(object$family == "gaussian"){
          ResCov <- ResCov + diag((object$params$phi^2))
    }
  }
  ResCov.q <- sapply(1:(object$num.lv+object$num.lv.c), function(q) sum(diag(ResCov.q[[q]])))
  if(object$num.lv.c==0)names(ResCov.q) <- paste("LV", 1:object$num.lv, sep = "")
  if(object$num.lv==0)names(ResCov.q) <- paste("CLV", 1:object$num.lv.c, sep = "")
  if(object$num.lv.c>0&object$num.lv>0)names(ResCov.q) <-c( paste("CLV", 1:object$num.lv.c, sep = ""), paste("LV", 1:object$num.lv, sep = ""))
  if(inherits(object,"gllvm.quadratic")){
    ResCov.q2 <- sapply(1:(object$num.lv+object$num.lv.c), function(q) sum(diag(ResCov.q2[[q]])))
    if(object$num.lv.c==0)names(ResCov.q2) <- paste("LV", 1:object$num.lv, "^2",sep = "")
    if(object$num.lv==0)names(ResCov.q2) <- paste("CLV", 1:object$num.lv.c, "^2",sep = "")
    if(object$num.lv>0&object$num.lv.c>0)names(ResCov.q2) <- c(paste("CLV", 1:object$num.lv.c, "^2",sep = ""),paste("LV", 1:object$num.lv, "^2",sep = ""))
  }
  
  if(object$quadratic==F){
    if(object$num.lv.c>0){
      if(any(ResCov.q[1:object$num.lv.c]<0.1)){
        warning("The residual variance of ",paste(colnames(object$lvs)[which(ResCov.q[1:object$num.lv.c]<0.01)],collapse=", and ")," is very small. This might indicate that the latent variable is nearly perfectly represented by covariates alone. \n")
      }
    }  
  }else{
    if(object$num.lv.c>0){
      if(any((ResCov.q+ResCov.q2)[1:object$num.lv.c]<0.01)){
        warning("The residual variance related to ",paste(colnames(object$lvs)[which((ResCov.q+ResCov.q2)[1:object$num.lv.c]<0.01)],collapse=", and ")," is very small. This might indicate that the latent variable is nearly perfectly represented by covariates alone. \n")
      }
    }  
  }
  
  
  colnames(ResCov) <- colnames(object$y)
  rownames(ResCov) <- colnames(object$y)
  if(inherits(object,"gllvm.quadratic")){
    out <- list(cov = ResCov, trace = sum(diag(ResCov)), var.q = ResCov.q, var.q2 = ResCov.q2)
  }else{
    out <- list(cov = ResCov, trace = sum(diag(ResCov)), var.q = ResCov.q)  
  }
  return(out)
}

#'@export getResidualCov
getResidualCov <- function(object, ...)
{
  UseMethod(generic = "getResidualCov")
}
