#' @title Extract residual correlations from gllvm object
#' @description  Calculates the residual covariance and correlation matrices for gllvm model.
#'
#' @param object   an object of class 'gllvm'.
#'
#' @return Function returns following components:
#'  \item{cor }{residual correlation matrix}
#'  \item{cov }{residual covariance matrix}
#'  \item{trace }{trace of residual covariance matrix}
#'
#' @author Francis K.C. Hui, Jenni Niku, David I. Warton
#'
#' @examples
#'# Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#'fit <- gllvm(y = y, family = poisson())
#'# residual correlations:
#'cr <- getResidualCor(fit)
#'
#'\dontrun{
#'# Plot residual correlations:
#'install.packages("corrplot", "gclus")
#'library(corrplot)
#'library(gclus)
#'rbPal <- colorRampPalette(c('darkblue', 'white', 'darkred'))
#'breaks <- seq(min(cr$cor), max(cr$cor), length.out = 40)
#'Colors <- rbPal(100)[as.numeric(cut(cr$cor, breaks = breaks))]
#'corrplot(cr$cor[order.single(cr$cor), order.single(cr$cor)], diag = F,
#'   type = "lower", method = "square", tl.cex = 0.8, tl.srt = 45, tl.col = "red")
#'   }
#'
#'@aliases getResidualCor getResidualCor.gllvm
#'@method getResidualCor gllvm
#'@export
#'@export getResidualCor.gllvm
getResidualCor.gllvm = function(object)
{
  ResCov <- object$params$theta %*% t(object$params$theta)
  Res.sd <- 1 / sqrt(diag(ResCov))
  Res.Cor <- diag(Res.sd) %*% ResCov %*% diag(Res.sd) * 0.99999
  colnames(Res.Cor) <- colnames(object$y)
  rownames(Res.Cor) <- colnames(object$y)
  out <- list(cor = Res.Cor, cov = ResCov, trace = sum(diag(ResCov)))
  return(out)
}

#'@export getResidualCor
getResidualCor <- function(object)
{
  UseMethod(generic = "getResidualCor")
}
