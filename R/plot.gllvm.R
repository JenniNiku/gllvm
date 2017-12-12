#' @title Plot Diagnostics for an gllvm Object
#' @description Four plots (selectable by which) are currently available: a plot of residuals against
#' linear predictors of fitted values, a Normal Q-Q plot of residuals, residuals against row index and residuals against column index.
#'
#' @param object an object of class 'gllvm'.
#' @param which if a subset of the plots is required, specify a subset of the numbers 1:4, see caption below.
#' @param caption captions to appear above the plots.
#' @param var.colors colors for responses, vector with length of number of response variables or 1. Defaults to NULL, when different responses have different colors.
#' @param ...	additional graphical arguments.
#'
#' @details
#' plot.gllvm is used for model diagnostics. Dunn-Smyth residuals or randomized quantile residuals (Dunn and Smyth, 1996) are used in plots. Colors indicate different species.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @references
#'
#' Dunn, P. K., and Smyth, G. K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics, 5, 236-244.
#'
#' Hui, F. K. C., Taskinen, S., Pledger, S., Foster, S. D., and Warton, D. I. (2015).  Model-based approaches to unconstrained ordination. Methods in Ecology and Evolution, 6:399-411.
#'
#'@seealso \code{\link{gllvm}}, \code{\link{residuals.gllvm}}
#' @examples
#' \dontrun{
#'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = "negative.binomial")
#'# Plot residuals
#'plot(fit, mfrow = c(2,2))
#'plot(fit, which = 1:2, mfrow = c(1,2))
#'}
#'@export


plot.gllvm <- function(object, which=1:4, caption=c("Residuals vs linear predictors", "Normal Q-Q","Residuals vs row index", "Residuals vs column index"),var.colors=NULL, ...) {
  n <- NROW(object$y)
  p <- NCOL(object$y)

  mains=rep("",4)
  mains[which]=caption[which]

  res=residuals.gllvm(object)
  ds.res <- res$residuals
  eta.mat <- res$linpred
  csum=order(colSums(as.matrix(object$y)))
  if(!is.null(var.colors)){
    col=rep(1,p)
    col[1:p]=(var.colors)
  } else
  if(p<8)
    col=(1:p)[csum]
  else
    col = (rainbow(p+1)[2:(p+1)])[csum]

par(...)

    if(1 %in% which) {plot(eta.mat,ds.res,xlab="linear predictors",ylab="Dunn-Smyth-residuals", col=rep(col,each=n), main = mains[1],...); abline(0,0,col=2)}
    if(2 %in% which) {qqnorm(c(ds.res),main = mains[2],ylab="Dunn-Smyth residuals",col=rep(col,each=n)); qqline(c(res$residuals), col = 2)}
    if(3 %in% which) {
      plot(rep(1:n,p),ds.res,xlab="site.index",ylab="Dunn-Smyth-residuals", col=rep(col,each=n),main = mains[3],...); abline(0,0,col=2)}
    if(4 %in% which) {
      plot(rep(1:p,each=n),ds.res,xlab="spp.index",ylab="Dunn-Smyth-residuals", col=rep(col[csum],each=n), main = mains[4],...); abline(0,0,col=2)}
}
