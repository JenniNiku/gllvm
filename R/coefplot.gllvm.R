#' @title Plot covariate coefficients and confidence intervals
#' @description Plots covariate coefficients and their confidence intervals.
#'
#' @param fit an object of class 'gllvm'.
#' @param y.label logical, if \code{TRUE} (default) colnames of y with respect to coefficients are added to plot.
#' @param which.Xcoef vector indicating which X-coefficients will be plotted. Can be vector of covariate names or numbers. Default is \code{NULL} when all covariate coefficients are plotted.
#' @param cex.ylab the magnification to be used for axis annotation relative to the current setting of cex.
#' @param mfrow same as \code{mfrow} in \code{par}. If \code{NULL} (default) it is determined automatically.
#' @param mar vector of length 4, which defines the margin sizes: \code{c(bottom, left, top, right)}. Defaults to \code{c(4,5,2,1)}.
#' @param ...	additional graphical arguments.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>, Francis K.C. Hui, Sara Taskinen
#'
#' @examples
#'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'X <- as.matrix(antTraits$env)
#'TR <- antTraits$traits
#'# Fit model with environmental covariates
#'fit <- gllvm(y, X, formula = ~ Bare.ground + Shrub.cover,
#'             family = "poisson")
#'coefplot.gllvm(fit)
#'
#'\donttest{
#'# Fit model with all environmental covariates
#'fitx <- gllvm(y, X, family = "negative.binomial")
#'coefplot.gllvm(fitx, mfrow = c(3,2))
#'coefplot.gllvm(fitx, which.Xcoef = 1:2)
#'
#'# Fit gllvm model with environmental and trait covariates
#'TR <- antTraits$traits
#'fitT <- gllvm(y = y, X = X, TR = TR, family = "negative.binomial")
#'coefplot.gllvm(fitT)
#'}
#'@export


coefplot.gllvm <- function(fit, y.label = TRUE, which.Xcoef=NULL,cex.ylab=0.5, mfrow=NULL, mar=c(4,6,2,1),...)
{
  if(is.null(fit$X)) stop("No X covariates in the model.");
  if(is.null(fit$TR)){
    if(is.null(which.Xcoef)) which.Xcoef <- c(1:NCOL(fit$params$Xcoef))
    Xcoef <-as.matrix(fit$params$Xcoef[,which.Xcoef])
  cnames <- colnames(fit$params$Xcoef)[which.Xcoef]
  k <- length(cnames)
  labely <- rownames(Xcoef)
  m <- length(labely)
  Xc <- Xcoef
  if(is.null(mfrow) && k>1) mfrow=c(1,k);
  if(!is.null(mfrow)) par(mfrow = mfrow, mar = mar)
  if(is.null(mfrow)) par(mar = mar)
  for (i in 1:k) {
    Xc <- Xcoef[,i]
    sdXcoef<-as.matrix(fit$sd$Xcoef[,which.Xcoef])
    lower <- Xc - 1.96 * sdXcoef[,i]
    upper <- Xc + 1.96 * sdXcoef[,i]
    Xc <- sort(Xc)
    lower <- lower[names(Xc)]
    upper <- upper[names(Xc)]

    col.seq <- rep("black", m)
    col.seq[lower < 0 & upper > 0] <- "grey"

    At.y <- seq(1,m)
    plot(x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = c(min(lower), max(upper)), pch = "x",cex.lab=1.3,...)

    segments(x0 = lower, y0 = At.y, x1 = upper, y1 = At.y, col = col.seq)
    abline(v=0, lty=1)
    if (y.label) axis(2, at=At.y, labels=names(Xc), las=1,cex.axis=cex.ylab) }
  } else{
    Xcoef <-fit$params$B

    xnames <- names(fit$params$B)
    m <- length(xnames); k=1
    #if(fit$get.trait) k=2;
    if(is.null(mfrow)) mfrow=c(1,k);
    par(mfrow = mfrow, mar = mar)
    Xc <- Xcoef
    sdXcoef<-fit$sd$B
    lower <- Xc - 1.96 * sdXcoef
    upper <- Xc + 1.96 * sdXcoef
    Xc <- sort(Xc)
    lower <- lower[names(Xc)]
    upper <- upper[names(Xc)]

    col.seq <- rep("black", m)
    col.seq[lower < 0 & upper > 0] <- "grey"

    At.y <- seq(1,m)
    plot(x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = "Coefficients", xlim = c(min(lower), max(upper)), pch = "x",cex.lab=1.3,...)

    segments(x0 = lower, y0 = At.y, x1 = upper, y1 = At.y, col = col.seq)
    abline(v=0, lty=1)
    if (y.label) axis(2, at=At.y, labels=names(Xc), las=1,cex.axis=cex.ylab)

  }

}
