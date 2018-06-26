#'@export

print.gllvm <- function(x, ...) {
  cat("Call: \n")
  print(x$call)
  cat("family: \n")
  print(x$family)
  cat("method: \n")
  print(x$method)
  cat("\n")
  cat("log-likelihood: ",x$logL,"\n")
  if(!is.null(x$params$inv.phi)){ x$params$inv.phi<-NULL; }
  crit=inf.criteria(x)
  df=crit$k
  cat("Degrees of freedom: ",df,"\n")
  cat("BIC: ",crit$BIC,"\n")
  cat("AIC: ",crit$AIC,"\n")
  cat("AICc: ",crit$AICc,"\n")
}
