#' @title Plot Diagnostics for an gllvm Object
#' @description Five plots (selectable by which) are currently available: a plot of residuals against
#' linear predictors of fitted values, a Normal Q-Q plot of residuals with a simulated point-wise 95\% confidence interval envelope, residuals against row index and column index and scale location plot.
#'
#' @param x an object of class 'gllvm'.
#' @param which if a subset of the plots is required, specify a subset of the numbers 1:5, see caption below.
#' @param caption captions to appear above the plots.
#' @param var.colors colors for responses, vector with length of number of response variables or 1. Defaults to NULL, when different responses have different colors.
#' @param add.smooth	logical indicating if a smoother should be added.
#' @param envelopes logical, indicating if simulated point-wise confidence interval envelope will be added to Q-Q plot, defaults to \code{TRUE}
#' @param reps number of replications when simulating confidence envelopes for normal Q-Q plot
#' @param envelope.col colors for envelopes, vector with length of two
#' @param n.plot number of species (response variables) to be plotted. Defaults to \code{NULL} when all response variables are plotted. Might be useful when data is very high dimensional.
#' @param ...	additional graphical arguments.
#'
#' @details
#' plot.gllvm is used for model diagnostics. Dunn-Smyth residuals (randomized quantile residuals) (Dunn and Smyth, 1996) are used in plots. Colors indicate different species.
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
#'data(antTraits, package = "mvabund")
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model with Poisson family
#'fit <- gllvm(y, family = poisson())
#'# Plot residuals
#'plot(fit, mfrow = c(3,2))
#'
#' \donttest{
#'# Fit gllvm model with negative binomial family
#'fitnb <- gllvm(y = y, family = "negative.binomial")
#'# Plot residuals
#'plot(fitnb, mfrow = c(3,2))
#'# Plot only two first plots
#'plot(fitnb, which = 1:2, mfrow = c(1,2))
#'}
#'}
#'@export


plot.gllvm <- function(x, which = 1:5, caption = c("Residuals vs linear predictors", "Normal Q-Q",
                                                   "Residuals vs row", "Residuals vs column", "Scale-Location"), 
                       var.colors = NULL, add.smooth = TRUE, envelopes = TRUE, reps = 150, 
                       envelope.col = c("blue","lightblue"), n.plot = NULL, ...) {
  n <- NROW(x$y)
  p <- NCOL(x$y)
  
  sppind=1:p
  if(!is.null(n.plot)) {
    sppind <- sort(sample(1:p, n.plot))
    p <- n.plot
  }
  
  mains <- rep("", 4)
  mains[which] <- caption[which]
  
  res <- residuals(x)
  ds.res <- res$residuals[,sppind]

  eta.mat <- res$linpred[,sppind]
  xxx <- boxplot(c(eta.mat), outline = FALSE,plot = FALSE)$stats
  
  csum = order(colSums(as.matrix(x$y))[sppind])
  if (!is.null(var.colors)) {
    col <- rep(1, p)
    col[1:p] <- (var.colors)
  } else
    if (p < 8)
      col <- (1:p)[csum]
  else
    col <- (grDevices::rainbow(p + 1)[2:(p + 1)])[csum]
  
  gr.pars <- list(...)
  par(...)
  
  if(1 %in% which) {
    if(is.null(gr.pars$xlim)) {
      plot(eta.mat[!is.na(ds.res)], ds.res[!is.na(ds.res)], xlab = "linear predictors", ylab = "Dunn-Smyth-residuals",
           type = "n", col = rep(col[csum], each = n)[!is.na(ds.res)], main = mains[1], xlim = c(min(xxx), max(xxx))); abline(0, 0, col = "grey", lty = 3)
    } else {
      plot(eta.mat[!is.na(ds.res)], ds.res[!is.na(ds.res)], xlab = "linear predictors", ylab = "Dunn-Smyth-residuals", type =
             "n", col = rep(col[csum], each = n), main = mains[1], ...); abline(0, 0, col = "grey", lty = 3)
    }
    
    if(add.smooth) gamEnvelope(eta.mat[!is.na(ds.res)], ds.res[!is.na(ds.res)], col = rep(col[csum], each = n)[!is.na(ds.res)], envelopes = envelopes, envelope.col = envelope.col, ...)
    #      panel(eta.mat, ds.res, col = rep(col, each = n), cex = 1, cex.lab = 1, cex.axis = 1, lwd = 1)
  }
  if(2 %in% which) {
    qq.x<-qqnorm(c(ds.res)[!is.na(ds.res)], main = mains[2], ylab = "Dunn-Smyth residuals", col = rep(col[csum], each = n)[!is.na(ds.res)], cex = 0.5, xlab = "theoretical quantiles");
    qqline(c(res$residuals)[!is.na(ds.res)], col = envelope.col[1])
    if(envelopes){
      K <- reps
      yy <- quantile(ds.res[!is.na(ds.res)], c(0.25, 0.75), names = FALSE, type = 7, na.rm = TRUE)
      xx <- qnorm(c(0.25, 0.75))
      slope <- diff(yy) / diff(xx)
      int <- yy[1L] - slope * xx[1L]
      Xm <- Ym <- NULL
      for (i in 1:K) {
        ri <- (rnorm(n * p-sum(is.na(ds.res)), int, sd = slope))
        Ym <- cbind(Ym, sort(ri))
      }
      Xm <- sort(qq.x$x)
      cis <- apply(Ym, 1, quantile, probs = c(0.025, 0.975))
      
      n.obs <- n*p-sum(is.na(ds.res))
      polygon(Xm[c(1:n.obs,n.obs:1)], c(cis[1, ],cis[2, n.obs:1]), col = envelope.col[2], border = NA)
      points(qq.x, col = rep(col[csum], each = n)[!is.na(ds.res)], cex = 0.5)
      qqline(c(res$residuals), col = envelope.col[1])
    }
  }
  if(3 %in% which) {
    plot(rep(1:n, p)[!is.na(ds.res)], ds.res[!is.na(ds.res)], xlab = "site index", ylab = "Dunn-Smyth-residuals", col =
           rep(col[csum], each = n)[!is.na(ds.res)], main = mains[3], ...);
    abline(0, 0, col = "grey", lty = 3)
    if(add.smooth) panel.smooth(rep(1:n, p)[!is.na(ds.res)], ds.res[!is.na(ds.res)], col = rep(col[csum], each = n)[!is.na(ds.res)], col.smooth = envelope.col[1], cex = NA, ...)
    #panel(rep(1:n, p), ds.res, col = rep(col, each = n), cex = 1, cex.lab = 1, cex.axis = 1, lwd = 1)
  }
  if(4 %in% which) {
    plot(rep(1:p, each = n)[!is.na(ds.res)], ds.res[!is.na(ds.res)], xlab = "species index", ylab = "Dunn-Smyth-residuals", col =
           rep(col[csum], each = n)[!is.na(ds.res)], main = mains[4], ...);  abline(0, 0, col = "grey", lty = 3)
    if(add.smooth) panel.smooth(rep(1:p, each = n)[!is.na(ds.res)], ds.res[!is.na(ds.res)], col = rep(col[csum], each = n)[!is.na(ds.res)], col.smooth = envelope.col[1], cex = NA, ...)
    #panel(rep(1:p, each = n), ds.res, col = rep(col[csum], each = n), cex = 1, cex.lab = 1, cex.axis = 1, lwd = 1)
  }
  if(5 %in% which) {
    sqres <- sqrt(abs(ds.res[!is.na(ds.res)]))
    yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name("Dunn-Smyth-residuals"))))
    if(is.null(gr.pars$xlim)) {
      plot(eta.mat[!is.na(ds.res)], sqres, xlab = "linear predictors", ylab = yl, col = rep(col[csum], each = n)[!is.na(ds.res)],
           main = mains[5], xlim = c(min(xxx), max(xxx)), ...);
    } else {
      plot(eta.mat[!is.na(ds.res)], sqres, xlab = "linear predictors", ylab = yl, col = rep(col[csum], each = n)[!is.na(ds.res)], main = mains[5], ...);
    }
    if(add.smooth) panel.smooth(eta.mat[!is.na(ds.res)], sqres, col = rep(col[csum], each = n)[!is.na(ds.res)], col.smooth = envelope.col[1], cex = NA, ...)
    #panel(eta.mat, sqres, col = rep(col, each = n), cex = 1, cex.lab = 1, cex.axis = 1, lwd = 1)
  }
  
}
