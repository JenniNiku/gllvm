#' @title Summarizing gllvm model fits
#' @description A summary of the fitted 'gllvm' object, including function call, distribution family and model parameters.
#' 
#' Various options are available to include extra parameter estimates in the summary, which have been excluded by default, for readability.
#'
#' @param object an object of class 'gllvm'
#' @param x a summary object
#' @param digits the number of significant digits to use when printing
#' @param signif.stars If \code{TRUE}, significance stars are printed for each coefficient, defaults to \code{TRUE}
#' @param dispersion option to return dispersion parameters, defaults to \code{FALSE}
#' @param spp.intercepts option to return species intercepts, defaults to \code{FALSE}
#' @param row.intercepts option to return row intercepts, defaults to \code{FALSE} 
#' @param theta option to return species scores in the ordination, defaults to \code{FALSE}
#' @param ...	not used.
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

summary.gllvm <- function(object, digits = max(3L, getOption("digits") - 3L),
                          signif.stars = getOption("show.signif.stars"), dispersion = FALSE, spp.intercepts = FALSE, row.intercepts = FALSE, theta = FALSE,
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

  M <- cbind(object$params$beta0, object$params$theta)
  sumry <- list()
  sumry$digits <- digits
  sumry$signif.stars <- signif.stars
  sumry$dispersion <- dispersion
  sumry$spp.intercepts <- spp.intercepts
  sumry$row.intercepts <- row.intercepts
  sumry$theta <- theta
  sumry$num.lv <- num.lv
  sumry$num.lv.c <- num.lv.c
  sumry$num.RR <- num.RR
  sumry$formula <- object$formula
  sumry$lv.formula <- object$lv.formula
  sumry$'log-likelihood' <- object$logL
  crit <- inf.criteria(object)
  sumry$df <- crit$k
  sumry$AIC <- crit$AIC
  sumry$AICc <- crit$AICc
  sumry$BIC <- crit$BIC

  crit <-
    newnams <- c("Intercept")

  if (num.lv > 0){
    if(quadratic != FALSE){
      if((num.lv.c+num.RR)==0)newnams <- c(newnams, paste("theta.LV", 1:num.lv, sep = ""), paste("theta.LV", 1:num.lv, "^2" ,sep = ""))
      if(num.lv.c+num.RR>0)newnams <- c(newnams, paste("theta.CLV", 1:(num.lv.c+num.RR), sep = ""), paste("theta.LV", 1:num.lv, sep = ""), paste("theta.CLV", 1:(num.lv.c+num.RR), "^2" ,sep = ""),paste("theta.LV", 1:num.lv, "^2" ,sep = ""))
    }
    if(quadratic == FALSE){
      if((num.lv.c+num.RR)==0)newnams <- c(newnams, paste("theta.LV", 1:num.lv, sep = ""))
    if((num.lv.c+num.RR)>0)newnams <- c(newnams, paste("theta.CLV", 1:(num.lv.c+num.RR), sep = ""), paste("theta.LV", 1:num.lv, sep = ""))
    }
  }else if(num.lv==0&(num.lv.c+num.RR)>0){
    if(quadratic != FALSE){
      if((num.lv.c+num.RR)>0)newnams <- c(newnams, paste("theta.CLV", 1:(num.lv.c+num.RR), sep = ""), paste("theta.CLV", 1:(num.lv.c+num.RR), "^2" ,sep = ""))
    }
    if(quadratic == FALSE){
      if((num.lv.c+num.RR)>0)newnams <- c(newnams, paste("theta.CLV", 1:(num.lv.c+num.RR), sep = ""))
    }
  }
  
  if (!is.logical(object$sd)&!is.null(object$X)&is.null(object$TR)) {
    pars <- c(object$params$Xcoef)
    se <- c(object$sd$Xcoef)
    
    zval <- pars/se
    pvalue <- 2 * pnorm(-abs(zval))
    coef.table <- cbind(pars, se, zval, pvalue)
    dimnames(coef.table) <- list(paste(rep(colnames(object$X.design),each=ncol(object$y)),colnames(object$y),sep=":"), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  }else if(!is.logical(object$sd)&!is.null(object$X)){
    pars <- c(object$params$B)
    se <- c(object$sd$B)
    nX <- dim(object$X)[2]
    cnx <- rep(colnames(object$X.design), each = p)
    rnc <- rep(rownames(object$params$Xcoef), nX)
    newnam <- paste(cnx, rnc, sep = ":")
    zval <- pars/se
    pvalue <- 2 * pnorm(-abs(zval))
    coef.table <- cbind(pars, se, zval, pvalue)
    dimnames(coef.table) <- list(newnam, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  }else{
    coef.table <- NULL
  }
  
  if (!is.logical(object$sd)&!is.null(object$lv.X)) {
    pars <- c(object$params$LvXcoef)
    se <- c(object$sd$LvXcoef)
    
    zval <- pars/se
    pvalue <- 2 * pnorm(-abs(zval))
    coef.table.constrained <- cbind(pars, se, zval, pvalue)
    dimnames(coef.table.constrained) <- list(paste(rep(colnames(object$lv.X),2),"(LV",rep(1:(object$num.lv.c+object$num.RR),each=ncol(object$lv.X)),")",sep=""), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  }else{
    coef.table.constrained <- NULL
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
  if(!is.null(object$X)){
  sumry$'Coef.tableX' <- coef.table
  }
  if((num.lv+num.lv.c)>0){
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
  
  cat("Constrained LVs: ", x$num.lv.c, "\n")
  cat("Reduced Ranks: ", x$num.RR,"\n")
  cat("Unconstrained LVs: ", x$num.lv, "\n")
  if((x$num.lv+x$num.lv.c)>0){cat("Standard deviation of LVs: ", zapsmall(x$sigma.lv,x$digits),"\n\n")}else{cat("\n")}
  
  cat("Formula: ", paste(x$formula,collapse=""), "\n")
  cat("LV formula: ", ifelse(is.null(x$lv.formula),"~0", paste(x$lv.formula,collapse="")), "\n")
  
  df <- x[["df"]]
  if(!is.null(x$Coef.tableX)){
    cat("\nCoefficients predictors:\n")
    coefs <- x$Coef.tableX
    
    printCoefmat(coefs, digits = x$digits, signif.stars = ifelse(!is.null(x$Coef.tableLV),F,x$signif.stars), 
                 na.print = "NA")
  }
  if(x$theta){
  if((x$num.lv+x$num.lv.c)>0){
    cat("\nCoefficients LVs: \n")
    
    print(x$Coefficients[,-1,drop=F])
  }
  }
  
  if(!is.null(x$Coef.tableLV)){
    cat("\nCoefficients LV predictors:\n")
    coefs <- x$Coef.tableLV
    
    printCoefmat(coefs, digits = x$digits, signif.stars = x$signif.stars, 
                 na.print = "NA")
  }
  if(x$spp.intercepts){
  cat("\n Species Intercepts: \n")
  print(zapsmall(x$Coefficients[,1],x$digits))
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
    names(phi) <- row.names(x$Coefficients)
    cat("\n(Dispersion estimates for ", x$family, ":\n")
    print(phi)
  }
  }
  
  cat("\n")
  invisible(x)
}