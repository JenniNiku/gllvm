#' @title Plot covariate coefficients and confidence intervals
#' @description Plot covariate coefficients and their confidence intervals
#'
#' @param fit An object of class 'gllvm'
#' @param y.label Logical, if \code{TRUE} (default) colnames of y with respect to coefficients are added to plot.
#' @param cex.ylab The magnification to be used for axis annotation relative to the current setting of cex.
#' @param mfrow Same as \code{mfrow} in \code{par}. If \code{NULL} (default) it is determined automatically.
#' @param mar vector of length 4, which defines the margin sizes: \code{c(bottom, left, top, right)}. Defaults to \code{c(4,5,2,1)}.
#' @param which.Xcoef vector indicating which X-coefficients will be plotted. Can be vector of covariate names or numbers. Default is \code{NULL} when all covariate coefficients are plotted.
#' @param ...	Additional graphical arguments.
#' @details
#' Plot covariate coefficients and their confidence intervals
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>, Francis K.C. Hui, Sara Taskinen
#'
#' @examples
#' \dontrun{
#' library(mvabund) ## Load a dataset from the mvabund package
#' data("antTraits")
#' y <- antTraits$abund
#' X <- antTraits$env
#' TR <- antTraits$traits
#' X <- scale(X)
#'
#' fitx <- gllvm(y = y, X = X, family = "negative.binomial")
#' coefplot.gllvm(fitx,mfrow=c(3,2), mar = c(4,8,2,1),cex.ylab=0.8)
#' coefplot.gllvm(fitx, which.Xcoef = c("Bare.ground","Shrub.cover"),
#'           mar = c(4,8,2,1), cex.ylab = 0.8)
#'
#' fitT <- gllvm(y = y, X = X, TR = TR, family = "negative.binomial")
#' coefplot.gllvm(fitT, mar = c(4,8,2,1))
#'}
#'@export


coefplot.gllvm <- function(fit, y.label = TRUE, which.Xcoef=NULL,cex.ylab=1, mfrow=NULL, mar=c(4,5,2,1),...)
{
  if(is.null(fit$X)) stop("No X covariates in the model.");
  if(is.null(fit$TR)){
    if(is.null(which.Xcoef)) which.Xcoef <- c(1:NCOL(fit$params$Xcoef))
    Xcoef <-fit$params$Xcoef[,which.Xcoef]
  cnames <- colnames(Xcoef)
  k <- length(cnames)
  labely <- rownames(Xcoef)
  m <- length(labely)
  Xc <- Xcoef
  if(is.null(mfrow)) mfrow=c(1,k);
  par(mfrow = mfrow, mar = mar)
  for (i in 1:k) {
    Xc <- Xcoef[,i]
    sdXcoef<-fit$sd$Xcoef[,which.Xcoef]
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
    if(is.null(which.Xcoef)) which.Xcoef <- c(1:length(fit$params$Xcoef))
    Xcoef <-fit$params$Xcoef[which.Xcoef]

    xnames <- names(Xcoef)
    m <- length(xnames); k=1
    if(fit$get.trait) k=2;
    if(is.null(mfrow)) mfrow=c(1,k);
    par(mfrow = mfrow, mar = mar)
    Xc <- Xcoef
    sdXcoef<-fit$sd$Xcoef[which.Xcoef]
    lower <- Xc - 1.96 * sdXcoef
    upper <- Xc + 1.96 * sdXcoef
    Xc <- sort(Xc)
    lower <- lower[names(Xc)]
    upper <- upper[names(Xc)]

    col.seq <- rep("black", m)
    col.seq[lower < 0 & upper > 0] <- "grey"

    At.y <- seq(1,m)
    plot(x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = "X-coefficients", xlim = c(min(lower), max(upper)), pch = "x",cex.lab=1.3,...)

    segments(x0 = lower, y0 = At.y, x1 = upper, y1 = At.y, col = col.seq)
    abline(v=0, lty=1)
    if (y.label) axis(2, at=At.y, labels=names(Xc), las=1,cex.axis=cex.ylab)
    if(fit$get.trait) {
    tnames <- names(fit$params$Tcoef);
    m <- length(tnames);
    Tc <- fit$params$Tcoef
    lower <- Tc - 1.96 * fit$sd$Tcoef
    upper <- Tc + 1.96 * fit$sd$Tcoef
    Tc <- sort(Tc)
    lower <- lower[names(Tc)]
    upper <- upper[names(Tc)]

    col.seq <- rep("black", m)
    col.seq[lower < 0 & upper > 0] <- "grey"

    At.y <- seq(1,m)
    plot(x = Tc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = "Trait coefficients", xlim = c(min(lower), max(upper)), pch = "x",cex.lab=1.3,...)

    segments(x0 = lower, y0 = At.y, x1 = upper, y1 = At.y, col = col.seq)
    abline(v=0, lty=1)
    axis(2, at=At.y, labels=names(Tc), las=1,cex.axis=cex.ylab)
    }

  }

}
