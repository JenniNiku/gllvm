#' @title Analysis Of Deviance for gllvm
#' @description  Computes an analysis of deviance table for two generalized linear latent variable model fits.
#'
#' @param object   an object of class 'gllvm'.
#' @param which either "multi" or "uni". If "uni", performs anova for each species separately.
#' @param method method used to adjust p-values for multiple testing when \code{which="uni"}. One of "holm" (default), "hochberg", "hommel", "bonferonni", "BH", BY", "fdr", or "none". See \code{\link{p.adjust}} for more information.
#' @param ...   one or more objects of class 'gllvm'
#'
#' @details
#' Computes likelihood-ratio test for two or more gllvm models.
#' Test results makes sense only for nested models.
#' Notice also that this test is not designed for testing models which have
#' degrees of freedom difference larger than 20. For such models
#' the P-value should be treated as very approximate.
#'
#' @author Jenni Niku, Bert van der Veen
#'
#' @examples
#' \donttest{
#'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- antTraits$abund
#'X <- antTraits$env
#'TR <- antTraits$traits
#'# Fit gllvm model
#'fit1 <- gllvm(y, X, TR, formula = ~ Bare.ground + Shrub.cover, family = poisson())
#'fit2 <- gllvm(y, X, TR, formula = ~ Bare.ground + Shrub.cover +
#'              (Bare.ground + Shrub.cover) : Webers.length, family = poisson())
#'# Test if the model with fourth corner interaction terms is significantly
#'# better using likelihood-ratio test:
#'anova(fit1, fit2)
#'}
#'@export

anova.gllvm <- function(object, ... ,which="multi",method="holm") {
  objects <- list(object, ...)
  if (length(objects) < 2)
    stop("At least two objects are needed for tests.")
  if (any(!(sapply(objects, function(x)inherits(x,"gllvm")))))
    stop("The function 'anova.gllvm' can only be used for a gllvm object.")

  tt <- sapply(objects, function(x)
    x$method)
  if (!(all(tt == "VA") == !all(tt == "LA")))
    stop("The objects are not comparable when they are fitted using different methods.")

  y <- object$y
  n <- NROW(y)
  p <- NCOL(y)
  diff <- sapply(objects, function(x) sum(x$y - y))
  if (any(!(diff == 0)))
    stop("The objects can not be compared")

  if(which=="multi"){
    df.list <- sapply(objects, function(x)
      attr(logLik.gllvm(x), "df"))
    objects_order <- objects[order(df.list)]
    formulas <- sapply(objects_order, function(x)
      formula(x$terms))
    
    #if(test=="LR"){
    df.list <- sapply(objects_order, function(x) attr(logLik.gllvm(x), "df"))
    ll.list <- sapply(objects_order, logLik.gllvm)
    
    D <- 2 * (ll.list[-1] - ll.list[1:(length(df.list) - 1)])
    df.chisq <- (df.list[-1] - df.list[1:(length(df.list) - 1)])
    Pval <- 1 - pchisq(D, df.chisq)
    paste("Model", 1:length(objects_order))
    result <- data.frame( Resid.Df = n * p - df.list, D = c(0, D), Df.diff = c(0, df.chisq), P.value = c("", signif(Pval)) )
    #}
    if (any(result$Df > 20))
      warning( "This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.\n")
    for (i in 1:length(objects_order)) {
      formchar <- as.character(formulas[[i]])
      if (length(formchar) == 3)
        formchar <- formchar[c(2, 1, 3)]
      cat("Model ", i, ": ", formchar, "\n")
    }
    return(result)
  }else if(which=="uni"){
    if(object$method=="LA"){
      stop("Species-specific anova has not yet been implemented for the Laplace approximation.")
    }
    df.list <- sapply(objects, function(x)
      attr(logLik.gllvm(x), "df"))
    objects_order <- objects[order(df.list)]
    formulas <- sapply(objects_order, function(x)
      formula(x$terms))
    
    #if(test=="LR"){
    df.list <- sapply(objects_order, function(x) {
      df <- attr(logLik.gllvm(x), "df")
      
      df <- df - p*x$num.lv +  x$num.lv * (x$num.lv - 1) / 2
      
      if(length(x$params$beta0)==1){
        df <- df-1
      }
 
    if(x$row.eff=="fixed"){
      df <- df-n-1
    }else if(x$row.eff=="random"){
      df<-df-length(x$params$sigma)
      if(!is.null(x$params$rho)) df<-df-length(x$params$rho)
    }
      if(x$family=="ordinal"){
    if(x$zeta.struc=="common"){
      df<-(df-length(unlist(x$params$zeta)[-1]))
    }
      }
      
      if (!is.null(x$randomX)){
        x$params$Br <- NULL
        x$params$sigmaB <- object$params$sigmaB[lower.tri(object$params$sigmaB, diag = TRUE)]
        df <- df - length(x$params$sigmaB)
      }else if(!is.null(x$params$B)){
        df <- df-length(x$params$B)
      }
      
      
    df <- df/p
    
    if(x$row.eff=="fixed"){
      df <- df+n-1
    }else if(x$row.eff=="random"){
      df<-df+length(x$params$sigma)
      if(!is.null(x$params$rho)) df<-df+length(x$params$rho)
    }
    if(x$family=="ordinal"){
    if(x$zeta.struc=="common"){
      df<-df+length(unlist(x$params$zeta)[-1])
    }    
    }
    if(length(x$params$beta0)==1){
      df <- df+1
    }
    if (!is.null(x$randomX)){
      df <- df + length(x$params$sigmaB)
    }else if(!is.null(x$params$B)){
      df <- df+length(x$params$B)
    }

    
    df <- df + sapply(1:p,function(j)sum(!x$params$theta[j,1:x$num.lv]==0))
    #still add terms for traits..
    return(df)
    })

  ll.list <- sapply(objects_order, function(x){
    LL<--colSums(x$TMBfn$report()$nll)
    
    # if(x$method == "VA"){
    #   if(x$num.lv > 0) LL = LL + n*0.5*x$num.lv
    #   if(x$row.eff == "random") LL = LL + length(x$params$row.params)*0.5
    #   if(x$family=="gaussian") {
    #     LL <- LL - n*p*log(pi)/2
    #   }
    # }
    return(LL)
  }
    )

  D <- 2 * (ll.list[,-1,drop=F] - ll.list[,1:(ncol(df.list)-1),drop=F])
  df.chisq <- (df.list[,-1,drop=F] - df.list[,1:(ncol(df.list) - 1),drop=F])
  Pval <- 1 - pchisq(D, df.chisq)
  Pval <- apply(Pval,2,p.adjust,method=method)
  
  result <- NULL
  for(j in 1:p){
    result <- rbind(result,Species=c(colnames(object$y)[j],rep("",length(objects)-1)),Resid.Df = n - df.list[j,], D = c(0, signif(D[j,])), Df.diff = c(0, df.chisq[j,]), P.value = c("", signif(Pval[j,])) )
  }
    result <- as.data.frame(t(result))
    row.names(result) <- paste("Model", 1:length(objects_order))
    #}
    if (any(df.chisq > 20))
      warning( "This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.\n")
    for (i in 1:length(objects_order)) {
      formchar <- as.character(formulas[[i]])
      if (length(formchar) == 3)
        formchar <- formchar[c(2, 1, 3)]
      cat("Model ", i, ": ", formchar, "\n")
    }
  
  
  data<-rbind(Resid.Df = n - df.list, D = cbind(0,D), Df.diff = cbind(0,df.chisq), P.value =cbind(NA,Pval))
  colnames(data) <- paste("Model.",1:length(objects))
  result2<-cbind(rep(c("Resid.Df","D","Df.diff","P.val"),each=p),Species=rep(colnames(object$y),times=4), data)
  result2[,4]<-as.numeric(result2[,4])
  print(result)
  return(invisible(list(table=result,data=result2)))
  }
}
