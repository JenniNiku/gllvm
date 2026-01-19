#'@export

print.gllvm <- function(x, ...) {
  cat("Call: \n")
  print(x$call)
  cat("family: \n")
  print(unique(x$family))
  cat("method: \n")
  print(x$method)
  cat("\n")
  cat("log-likelihood: ", x$logL, "\n")
  if(!is.null(x$params$inv.phi)){ x$params$inv.phi <- NULL; }
  crit <- inf.criteria(x)
  df <- crit$k
  cat("Residual degrees of freedom: ", length(x$y) - df, "\n")
  cat("AIC: ", crit$AIC, "\n")
  cat("AICc: ", crit$AICc, "\n")
  cat("BIC: ", crit$BIC, "\n")
  invisible(crit)
}
