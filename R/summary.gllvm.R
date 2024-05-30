#' @title Summarizing gllvm model fits
#' @description A summary of the fitted 'gllvm' object, including function call, distribution family and model parameters.
#' 
#' @details Various options are available to include extra parameter estimates in the summary, which have been excluded by default, for readability.
#'
#' @param object an object of class 'gllvm'
#' @param x a summary object
#' @param by By = "all" (default) will return a Wald statistics per predictor and LV if the ordination includes predictors, by = "terms" will return a multivariate Wald statistic per predictor (displayed at first LV), and by = "LV" will do the same but per dimension (displayed at first predictors).
#' @param digits the number of significant digits to use when printing
#' @param signif.stars If \code{TRUE}, significance stars are printed for each coefficient, defaults to \code{TRUE}
#' @param dispersion option to return dispersion parameters, defaults to \code{FALSE}
#' @param spp.intercepts option to return species intercepts, defaults to \code{FALSE}
#' @param row.intercepts option to return row intercepts, defaults to \code{FALSE} 
#' @param Lvcoefs option to return species scores in the ordination, defaults to \code{FALSE}. Returns species optima for quadratic model.
#' @param rotate defaults to \code{TRUE}. If \code{TRUE} rotates the output of the latent variables to principal direction, so that it coincides with the ordiplot results. If both unconstrained and constrained latent variables are included, predictor slopes are not rotated.
#' @param type to match "type" in \code{\link{ordiplot.gllvm}}
#' @param ...	 not used.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>, Bert van der Veen
#'
#' @examples
#' \dontrun{
#'## Load a dataset from the mvabund package
#' data(antTraits)
#' y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#' fit <- gllvm(y = y, family = poisson())
#' summary(fit)
#'}
#'@export
#'@export print.summary.gllvm 

summary.gllvm <- function(object, by = "all", digits = max(3L, getOption("digits") - 3L),
                          signif.stars = getOption("show.signif.stars"), dispersion = FALSE, spp.intercepts = FALSE, row.intercepts = FALSE, Lvcoefs = FALSE, rotate = TRUE, type = NULL,
                          ...) {
  
  n <- NROW(object$y)
  p <- NCOL(object$y)
  nX <- dim(object$X)[2]
  nTR <- dim(object$TR)[2]
  num.lv <- object$num.lv
  num.lv.c <- object$num.lv.c
  num.RR <- object$num.RR
  quadratic <- object$quadratic
  family <- object$family
  
  if(!by%in%c("all","terms","LV"))stop("'by' must be one of 'all', 'terms', or 'LV'.")
  #calculate rotation matrix
  if(rotate && (num.lv+num.lv.c+num.RR)>1){
    lv <- getLV(object, type = type)
    
    do_svd <- svd(lv)
    svd_rotmat_sites <- do_svd$v
    if(num.lv.c>0|(num.RR+num.lv)>0){
      do_svd$v <- svd(getLV(object))$v
    }
  }else{
    svd_rotmat_sites <- diag(num.lv.c+num.RR+num.lv)
  }
  if((num.RR+num.lv+num.lv.c)>0 && Lvcoefs){
    if(quadratic==FALSE){
      M <- cbind(object$params$beta0, object$params$theta[,1:(num.lv.c+num.RR+num.lv)]%*%svd_rotmat_sites)  
    }else{
      M <- cbind(object$params$beta0, optima(object,sd.errors=F)%*%svd_rotmat_sites)  
    }  
  }else{
    M <- cbind(object$params$beta0)
  }
  
  
  sumry <- list()
  sumry$digits <- digits
  sumry$signif.stars <- signif.stars
  sumry$dispersion <- dispersion
  sumry$row.intercepts <- row.intercepts
  sumry$Lvcoefs <- Lvcoefs
  sumry$num.lv <- num.lv
  sumry$num.lv.c <- num.lv.c
  sumry$num.RR <- num.RR
  sumry$quadratic <- quadratic
  sumry$formula <- object$formula
  if(is.null(object$col.eff$col.eff))object$col.eff$col.eff <- FALSE # backward compatibility
  if(object$col.eff$col.eff=="random")sumry$formula <- object$call$formula
  sumry$lv.formula <- object$lv.formula
  sumry$'log-likelihood' <- object$logL
  crit <- inf.criteria(object)
  sumry$df <- crit$k
  sumry$AIC <- crit$AIC
  sumry$AICc <- crit$AICc
  sumry$BIC <- crit$BIC
  
  crit <-
    newnams <- c("Intercept")
  
  if((num.lv+num.lv.c+num.RR)>0 && Lvcoefs){
    newnams <- c(newnams, dimnames(object$params$theta)[[2]][1:(num.lv+num.lv.c+num.RR)])
  }
  if(object$col.eff$col.eff=="random" || !is.null(object$randomX)){
    REcovs <- data.frame(Name = colnames(object$params$sigmaB), Variance = format(round(diag(object$params$sigmaB), digits), nsmall = digits), Std.Dev = format(round(sqrt(diag(object$params$sigmaB)), digits), nsmall = digits))
    if(!is.null(object$params$rho.sp)){
      REcovs <- do.call(cbind, list(data.frame(Name = REcovs$Name), data.frame(Signal = format(round(object$params$rho.sp, digits), nsmall = digits), row.names = NULL), REcovs[,-1]))
    }
    if(!all(object$params$sigmaB[row(object$params$sigmaB)!=col(object$params$sigmaB)]==0)){
      cors <- format(round(cov2cor(object$params$sigmaB), digits), nsmall = digits)
      cors[upper.tri(cors, diag = TRUE)] <- ""
      REcovs <- cbind(REcovs, cors, deparse.level = 0L)
      colnames(REcovs)[tail(1:ncol(REcovs), ncol(cors))] <- c("Corr", rep("", ncol(cors) - 1))
    }
    sumry$REcovs <- REcovs
  }
  
  if (!is.logical(object$sd)&!is.null(object$X)&is.null(object$TR)) {
    pars <- c(object$params$Xcoef)
    se <- c(object$sd$Xcoef)
    zval <- pars/se
    pvalue <- 2 * pnorm(-abs(zval))
    coef.table <- cbind(pars, se, zval, pvalue)
    dimnames(coef.table) <- list(paste(rep(colnames(object$X.design),each=ncol(object$y)),colnames(object$y),sep=":"), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    if(object$col.eff$col.eff == "random"){
      coef.table<- coef.table[!duplicated(coef.table),,drop=FALSE]
      row.names(coef.table)[grepl("RE_mean_", row.names(coef.table))] <- sub("RE_mean_", "RE mean:",grep("RE_mean_",colnames(object$X.design), value = TRUE))
    }
  }else if(!is.logical(object$sd)&!is.null(object$X)){
    pars <- c(object$params$B)
    se <- c(object$sd$B)
    nX <- dim(object$X)[2]
    cnx <- rep(colnames(object$X.design), each = ifelse(!is.null(object$TR),1,p))
    
    if(!is.null(object$params$Xcoef)){
      rnc <- rep(rownames(object$params$Xcoef), nX)
      newnam <- paste(cnx, rnc, sep = ":")
    }else if(!is.null(object$X.design)){
      newnam <- colnames(object$X.design)
    }
    zval <- pars/se
    pvalue <- 2 * pnorm(-abs(zval))
    coef.table <- cbind(pars, se, zval, pvalue)
    dimnames(coef.table) <- list(newnam, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  }else if(!is.logical(object$sd)&object$col.eff$col.eff=="random"){
  pars <- c(object$params$B)
  se <- c(object$sd$B)
  nX <- dim(object$col.eff$Xt)[2]
  newnam <- names(object$params$B)
  
  zval <- pars/se
  pvalue <- 2 * pnorm(-abs(zval))
  coef.table <- cbind(pars, se, zval, pvalue)
  dimnames(coef.table) <- list(newnam, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  }else{
    coef.table <- NULL
  }
  
  if(spp.intercepts&!is.logical(object$sd)){
      pars <- unique(object$params$beta0)
      se <- unique(object$sd$beta0)
      zval <- pars/se
      pvalue <- 2 * pnorm(-abs(zval))
      coef.table.int <- cbind(pars, se, zval, pvalue)
      if(!object$beta0com){
        newnam <- colnames(object$y)
      }else if(object$beta0com){
        newnam <- "Community intercept"
      }
      dimnames(coef.table.int) <- list(newnam, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

      coef.table <- rbind(coef.table.int, coef.table)
  }
  
  if(!is.null(object$lv.X) && is.null(object$lv.X.design))object$lv.X.design <- object$lv.X #for backward compatibility
  if (!is.logical(object$sd)&!is.null(object$lv.X.design)&object$randomB==FALSE) {
    if(!rotate|num.lv>0&(num.lv.c+num.RR)>0){
      pars <- c(object$params$LvXcoef)
      se <- c(object$sd$LvXcoef)
      zval <- pvalue <- rep(NA,length(pars))
      if(by=="all"){
      zval <- pars/se
      pvalue <- 2 * pnorm(-abs(zval))
      coef.table.constrained <- cbind(pars, se, zval, pvalue)
      dimnames(coef.table.constrained) <- list(paste(colnames(object$lv.X.design),"(CLV",rep(1:(object$num.lv.c+object$num.RR),each=ncol(object$lv.X.design)),")",sep=""), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
      }else if(by=="terms"){
        covB <- object$Hess$cov.mat.mod
        colnames(covB) <- row.names(covB) <- names(object$TMBfn$par)[object$Hess$incl]
        covB <- covB[row.names(covB)=="b_lv",colnames(covB)=="b_lv"]
        for(i in 1:ncol(object$lv.X.design)){
          idx <- seq(i,ncol(object$lv.X.design)*(object$num.lv.c+object$num.RR),ncol(object$lv.X.design))
          b <- object$params$LvXcoef[i,]
          zval[i] <- b%*%MASS::ginv(covB[idx,idx])%*%b
          pvalue[i] <- 1-pchisq(zval[i],object$num.lv.c+object$num.RR)
        }
        coef.table.constrained <- cbind(pars, se, zval, pvalue)
        dimnames(coef.table.constrained) <- list(paste(colnames(object$lv.X.design),"(CLV",rep(1:(object$num.lv.c+object$num.RR),each=ncol(object$lv.X.design)),")",sep=""), c("Estimate", "Std. Error", "X2 value", "Pr(>X2)"))
      }else if(by=="LV"){
        covB <- object$Hess$cov.mat.mod
        colnames(covB) <- row.names(covB) <- names(object$TMBfn$par)[object$Hess$incl]
        covB <- covB[row.names(covB)=="b_lv",colnames(covB)=="b_lv"]
        for(i in 1:(object$num.RR+object$num.lv.c)){
          b <- object$params$LvXcoef[,i]
          zval[1+ncol(object$lv.X.design)*(i-1)] <- b%*%MASS::ginv(covB[(1:ncol(object$lv.X.design))+ncol(object$lv.X.design)*(i-1),(1:ncol(object$lv.X.design))+ncol(object$lv.X.design)*(i-1)])%*%b
          pvalue[1+ncol(object$lv.X.design)*(i-1)] <- 1-pchisq(zval[1+ncol(object$lv.X.design)*(i-1)],ncol(object$lv.X.design))
        }
        coef.table.constrained <- cbind(pars, se, zval, pvalue)
        dimnames(coef.table.constrained) <- list(paste(colnames(object$lv.X.design),"(CLV",rep(1:(object$num.lv.c+object$num.RR),each=ncol(object$lv.X.design)),")",sep=""), c("Estimate", "Std. Error", "X2 value", "Pr(>X2)"))
      }
    }else{
      LVcoef <- (object$params$LvXcoef%*%svd_rotmat_sites)
      pars <- c(LVcoef)
      covB <- object$Hess$cov.mat.mod
      colnames(covB) <- row.names(covB) <- names(object$TMBfn$par)[object$Hess$incl]
      covB <- covB[row.names(covB)=="b_lv",colnames(covB)=="b_lv", drop=FALSE]
      zval <- pvalue <- rep(NA,length(pars))
      rotSD <- matrix(0,ncol=num.RR+num.lv.c,nrow=ncol(object$lv.X.design)) 
      for(i in 1:ncol(object$lv.X.design)){
        rotSD[i,] <- sqrt(abs(diag(t(svd_rotmat_sites[1:(num.lv.c+num.RR),1:(num.lv.c+num.RR)])%*%covB[seq(i,(num.RR+num.lv.c)*ncol(object$lv.X.design),by=ncol(object$lv.X.design)),seq(i,(num.RR+num.lv.c)*ncol(object$lv.X.design),by=ncol(object$lv.X.design))]%*%svd_rotmat_sites[1:(num.lv.c+num.RR),1:(num.lv.c+num.RR)])))
      }
      se <- c(rotSD)
      if(by=="all"){
        pars <- c(LVcoef)
        zval <- pars/se
        pvalue <- 2 * pnorm(-abs(zval))
        coef.table.constrained <- cbind(pars, se, zval, pvalue)
        dimnames(coef.table.constrained) <- list(paste(colnames(object$lv.X.design),"(CLV",rep(1:(object$num.lv.c+object$num.RR),each=ncol(object$lv.X.design)),")",sep=""), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
      }else if(by=="terms"){
        # This test is invariant to rotation, but the rotation matrix is included for posterity
        for(i in 1:ncol(object$lv.X.design)){
          idx <- seq(i,ncol(object$lv.X.design)*(object$num.lv.c+object$num.RR),ncol(object$lv.X.design))
          b <- LVcoef[i,]
          zval[i] <- b%*%MASS::ginv(t(svd_rotmat_sites[1:(num.lv.c+num.RR),1:(num.lv.c+num.RR)])%*%covB[idx,idx]%*%svd_rotmat_sites[1:(num.lv.c+num.RR),1:(num.lv.c+num.RR)])%*%b
          pvalue[i] <- 1-pchisq(zval[i],object$num.lv.c+object$num.RR)
        }
        coef.table.constrained <- cbind(pars, se, zval, pvalue)
        dimnames(coef.table.constrained) <- list(paste(colnames(object$lv.X.design),"(CLV",rep(1:(object$num.lv.c+object$num.RR),each=ncol(object$lv.X.design)),")",sep=""), c("Estimate", "Std. Error", "X2 value", "Pr(>X2)"))
      }else if(by=="LV"){
        covB <- object$Hess$cov.mat.mod
        colnames(covB) <- row.names(covB) <- names(object$TMBfn$par)[object$Hess$incl]
        covB <- covB[row.names(covB)=="b_lv",colnames(covB)=="b_lv"]

        for(q in 1:(object$num.RR+object$num.lv.c)){
         # Rotated cov matrix for the qth LV
          covBnew <- matrix(0,ncol=ncol(object$lv.X.design),nrow=ncol(object$lv.X.design))        
        for(k in 1:ncol(object$lv.X.design)){
          for(k2 in 1:ncol(object$lv.X.design)){
            for(q2 in 1:(object$num.RR+object$num.lv.c)){
              for(q3 in 1:(object$num.RR+object$num.lv.c)){
                covBnew[k,k2] <- covBnew[k,k2]+svd_rotmat_sites[1:(object$num.RR+object$num.lv.c),1:(object$num.RR+object$num.lv.c)][q2,q]*svd_rotmat_sites[1:(object$num.RR+object$num.lv.c),1:(object$num.RR+object$num.lv.c)][q3,q]*covB[(object$num.RR+object$num.lv.c)*(k-1)+q2,(object$num.RR+object$num.lv.c)*(k2-1)+q3]
              }
            }
          }
        }
          b <- LVcoef[,q]
          zval[1+ncol(object$lv.X.design)*(q-1)] <- b%*%solve(covBnew)%*%b
          pvalue[1+ncol(object$lv.X.design)*(q-1)] <- 1-pchisq(zval[1+ncol(object$lv.X.design)*(q-1)],ncol(object$lv.X.design))
        }
        coef.table.constrained <- cbind(pars, se, zval, pvalue)
        dimnames(coef.table.constrained) <- list(paste(colnames(object$lv.X.design),"(CLV",rep(1:(object$num.lv.c+object$num.RR),each=ncol(object$lv.X.design)),")",sep=""), c("Estimate", "Std. Error", "X2 value", "Pr(>X2)"))
        }
      }
  }else{
    coef.table.constrained <- NULL
  }
  if(!is.null(type)){
    if(type=="residual"){
      coef.table.constrained <- NULL
    }
  }
  
  colnames(M) <- newnams
  rownames(M) <- colnames(object$y)
  sumry$Call <- object$call
  sumry$family <- object$family
  sumry$Coefficients <- M
  
  if (!is.null(object$TR)) {
    # if (!is.null(object$X)) {
    #   sumry$'Covariate coefficients' <- object$params$B
    # }
  } else {
    # if (!is.null(object$X)) {
    #   sumry$'Environmental coefficients' <- object$params$Xcoef
    # }
  }
  if (!is.null(object$params$row.params)) {
    sumry$'Row intercepts' <- object$params$row.params
  }
  sumry$rstruc <- object$rstruc
  if (object$row.eff == "random") {
    object$params$sigma2 = object$params$sigma[1] ^ 2
    names(object$params$sigma2) = "sigma^2"
    sumry$'Variance of random row intercepts' <- object$params$sigma2
  }
  
  if (object$family == "negative.binomial") {
    sumry$'Dispersion parameters' <- object$params$phi
  }
  if (object$family == "gamma") {
    sumry$'Shape parameters' <- object$params$phi
  }
  if (object$family == "tweedie") {
    sumry$'Dispersion parameters' <- object$params$phi
  }
  if (object$family == "ZIP") {
    sumry$'Zero inflation p' <- object$params$phi
  }
  if(object$family == "gaussian"){
    sumry$'Standard deviations' <- object$params$phi
  }
  if(!is.null(object$X) | spp.intercepts | object$col.eff$col.eff == "random"){
    sumry$'Coef.tableX' <- coef.table
  }
  if((num.lv+num.lv.c)>0){
    
    #add zeros for num.rr
    object$params$sigma.lv <- c(object$params$sigma.lv[1:num.lv.c],rep(0,num.RR), if(num.lv>0)object$params$sigma.lv[-c(1:(num.lv.c+num.RR))])
    
    if(rotate){
      object$params$sigma.lv <- sqrt(diag(t(svd_rotmat_sites)%*%diag(object$params$sigma.lv^2, ncol = num.lv+num.lv.c, nrow = num.lv+num.lv.c)%*%svd_rotmat_sites))
    }
    sumry$sigma.lv <- object$params$sigma.lv
  }
  
  if((object$num.lv.c+object$num.RR)>0&!is.null(coef.table.constrained)){
    sumry$'Coef.tableLV' <- coef.table.constrained
  }
  class(sumry) <- "summary.gllvm"
  return(sumry)
}

#'@export
#'@rdname summary.gllvm 
print.summary.gllvm <- function (x, ...) 
{
  cat("\nCall:\n", paste(deparse(x$Call), sep = "\n", 
                         collapse = "\n"), "\n\n", sep = "")
  cat("Family: ", x$family, "\n\n")
  
  AIC <- round(x$AIC,x$digits)
  BIC <- round(x$BIC,x$digits)
  AICc <- round(x$AICc,x$digits)
  
  cat("AIC: ", AIC, "AICc: ", AICc, "BIC: ", BIC, "LL: ", zapsmall(x$`log-likelihood`, x$digits), "df: ", x$df, "\n\n")
  
  cat("Informed LVs: ", x$num.lv.c, "\n")
  cat("Constrained LVs: ", x$num.RR,"\n")
  cat("Unconstrained LVs: ", x$num.lv, "\n")
  
  #this scenario we don't want the SD from num.lv as it is meaningless
  if(x$num.lv>0&(x$num.RR+x$num.lv.c)>0 & isFALSE(x$quadratic))x$sigma.lv <- x$sigma.lv[1:(x$num.lv.c+x$num.RR)]
  
  #only print SD from LV if model is quadratic or if (hybrid) concurrent
  if((x$num.lv.c)>0|!isFALSE(x$quadratic)){cat("Residual standard deviation of LVs: ", zapsmall(x$sigma.lv,x$digits),"\n\n")}else{cat("\n")}
  
  cat("Formula: ", paste(x$formula, collapse = ""), "\n")
  cat("LV formula: ", ifelse(is.null(x$lv.formula),"~ 0", paste(x$lv.formula,collapse="")), "\n")
  
  if(!is.null(x$REcovs)){
    cat("\nRandom effects:\n")
    print(x$REcovs, row.names = FALSE, right = FALSE)
  }
  
  df <- x[["df"]]
  if(!is.null(x$Coef.tableX)){
    cat("\nCoefficients predictors:\n")
    coefs <- x$Coef.tableX
    
    printCoefmat(coefs, digits = x$digits, signif.stars = x$signif.stars, 
                 na.print = "NA")
  }
  if(x$Lvcoefs){
    if((x$num.lv+x$num.lv.c+x$num.RR)>0){
      if(x$quadratic==F){
        cat("\nCoefficients LVs: \n")  
      }else{
        cat("\nOptima LVs: \n")
      }
      
      print(x$Coefficients[,-1,drop=F])
    }
  }
  
  if(!is.null(x$Coef.tableLV)){
    cat("\nCoefficients LV predictors:\n")
    coefs <- x$Coef.tableLV
    
    printCoefmat(coefs, digits = x$digits, signif.stars = x$signif.stars, 
                 na.print = "")
  }
  if(x$row.intercepts){
    if(!is.null(x$`Row intercepts`)){
      cat("\n Row intercepts with variance", zapsmall(x$'Variance of random row intercepts',x$digits), ":\n")
      print(zapsmall(x$`Row intercepts`,x$digits))
    }
  }
  if(x$dispersion){
    
    if (x$family == "negative.binomial") {
      phi <- x$'Dispersion parameters'
    }
    if (x$family == "gamma") {
      phi <- x$'Shape parameters'
    }
    if (x$family == "tweedie") {
      phi <- x$'Dispersion parameters'
    }
    if (x$family == "ZIP") {
      phi <- x$'Zero inflation p'
    }
    if(x$family == "gaussian"){
      phi <- x$'Standard deviations'
    }
    if(x$family%in%c("negative.binomial","gamma","tweedie","ZIP","gaussian")){
      # names(phi) <- row.names(x$Coefficients)
      cat("\n(Dispersion estimates for ", x$family, ":\n")
      print(phi)
    }
  }
  
  invisible(x)
}

#'@export
#'@rdname plot.summary.gllvm 
plot.summary.gllvm <- function (x, component = NULL, ...) 
{
  args <- list(...)
  
  if(!is.null(component)){
    component <- match.arg(component, c("main", "LV"))
  }else if(!is.null(x$Coef.tableX)){
    component <- "main"
  }else if(!is.null(x$Coef.tableLV)){
    component <- "LV"
  }

  if(component == "main" && !is.null(x$Coef.tableX)){
    coefs <- x$Coef.tableX
  }else if(component == "LV" && !is.null(x$Coef.tableLV)){
    coefs <- x$Coef.tableLV
  }else{
    stop("Either nothing to plot, or forgot to select 'component' of either 1) 'main' or 2) 'LV'.")
  }
  
  if(!"mar"%in%names(args)){
    par(mar = c(4,7,2,1))
  }
  
  upper = coefs[,1]+ qnorm(0.95)*coefs[,2]
  lower = coefs[,1]+ qnorm(1-0.95)*coefs[,2]
  
  plot(x = coefs[,1], y = 1:nrow(coefs), yaxt = "n", ylab = "", xlab = "Estimate", pch = "x", xlim = c(min(lower), max(upper)), ...)
  segments(x0 = lower, y0 = 1:nrow(coefs), x1 = upper, y1 = 1:nrow(coefs))
  axis( 2, at = 1:nrow(coefs), labels = row.names(coefs), las = 1, ...)
  abline(v = 0, lty = 1)
}