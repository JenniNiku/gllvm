gllvm.iter <- function(...){
  args <- list(...)
  
  # if(!(args$family %in% c("poisson","negative.binomial","binomial","tweedie","ZIP", "ZINB", "gaussian", "ordinal", "gamma", "exponential", "beta", "betaH", "orderedBeta", "ZIB","ZNIB")))
    # stop("Selected family not permitted...sorry!")
  if(!(args$Lambda.struc %in% c("unstructured","diagonal","bdNN","UNN", "diagU", "UU")))
    stop("Lambda matrix (covariance of variational distribution for latent variable) not permitted...sorry!")
  
  if (!is.numeric(args$y))
    stop( "y must be numeric. If ordinal data, please convert to numeric with lowest level equal to 1.")
  # if ((family %in% c("ZIP")) && (method %in% c("VA", "EVA"))) #"tweedie", 
  #   stop("family=\"", family, "\" : family not implemented with VA method, change the method to 'LA'")
  if (is.null(rownames(args$y)))
    rownames(args$y) <- paste("Row", 1:nrow(args$y), sep = "")
  if (is.null(colnames(args$y)))
    colnames(args$y) <- paste("Col", 1:ncol(args$y), sep = "")
  if(any(args$family == "ordinal")){
    max.levels <- apply(args$y[,args$family == "ordinal", drop=FALSE],2,function(x) length(min(x,na.rm=TRUE):max(x,na.rm=TRUE)))
    if(any(max.levels == 1)&args$zeta.struc=="species" || all(max.levels == 2)&args$zeta.struc=="species")
      stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead.")
    if(any(!apply(args$y[,args$family == "ordinal", drop=FALSE],2,function(x)all(diff(sort(unique(x)))==1)))&args$zeta.struc=="species"){
      warning("Can't fit ordinal model if there are species with missing classes. Setting 'zeta.struc = `common`'.")
      args$zeta.struc = "common"
    }
    if(!all(min(args$y[,args$family == "ordinal", drop=FALSE],na.rm=TRUE)==apply(args$y[,args$family == "ordinal", drop=FALSE],2,min,na.rm=TRUE))&args$zeta.struc=="species"){
      warning("For ordinal data and zeta.struc=`species` all species must have the same minimum category. Setting 'zeta.struc = `common`'.")
      args$zeta.struc = "common"
    }
    if(any(diff(sort(unique(c(args$y[,args$family == "ordinal", drop=FALSE]))))!=1)&args$zeta.struc=="common")
      stop("Can't fit ordinal model if there are missing response classes. Please reclassify.")
  }
  
  if(any(args$family == "orderedBeta")) {
    if (!(args$method %in% c("VA", "EVA"))) #"tweedie", 
      stop("family=\"", "orderedBeta", "\" : family not implemented with LA method, change the method to 'VA'.")
    
    if((sum(args$y[,args$family == "orderedBeta", drop=FALSE]==1, na.rm = TRUE) + sum(args$y[,args$family == "orderedBeta", drop=FALSE]==0, na.rm = TRUE))==0){
      stop("No zeros or ones in the data, please use 'family = `beta`' instead.")
    }
    if(!all(colSums(args$y[,args$family == "orderedBeta", drop=FALSE]==1, na.rm = TRUE)>0) & !all(colSums(args$y[,args$family == "orderedBeta", drop=FALSE]==0, na.rm = TRUE)>0)){
      warning("Not all species have zeros and ones. Setting 'zeta.struc = `common`'.")
      args$zeta.struc = "common"
    }
  }
  
  
  model <- args$model
  args <- args[-which(names(args)=="model")]

  ### Seeds
  # If number of seeds is less than n.init, sample the seeds randomly, but using the given seed
  if((length(args$seed) >1) & (length(args$seed) < args$n.init)) {
    stop("Seed length doesn't match with the number of initial starts.")
  }
  seed = NULL
  if(!is.null(args$seed) & (length(args$seed) ==1) & (length(args$seed) < args$n.init)) {
    set.seed(args$seed)
    seed <- sample(1:10000, args$n.init)
  }else if(!is.null(args$seed)){
    seed = args$seed
  }
  
  if("seed" %in% names(args))args <- args[!names(args)=="seed"]
  
  # If no seed, new seed is not set for now on  
  # if(is.null(seed) & args$starting.val!="zero"){
  #   seed <- sample(1:10000, args$n.init)
  # }
  
if(args$n.init>1){
  fitFinal <- NULL

n.i.i <- 0;n.i <- 1

while(n.i <= args$n.init && n.i.i<args$n.init.max){
  if(args$n.init > 1 && args$trace)
    cat("Initial run", n.i, "\n")
  
  if(!is.null(seed)) set.seed(seed[n.i])
  
  if(model == "gllvm.TMB"){
  fit <- do.call(gllvm.TMB, args)
  }else if(model == "trait.TMB"){
    fit <- do.call(trait.TMB, args)
  }
  #### Check if model fit succeeded/improved on this iteration n.i
  
  # Gradient check with n.i >2 so we don't get poorly converged models - relatively relaxed tolerance
  if(!is.null(fitFinal$TMBfn)){
    norm.gr1 <- norm(as.matrix(fitFinal$TMBfn$gr()/length(fitFinal$TMBfn$par)))
  }else{
    norm.gr1 <- NaN
  }
  if(!is.null(fit$TMBfn$gr)){
  gr2 <- fit$TMBfn$gr()/length(fit$TMBfn$par)
  norm.gr2 <- norm(as.matrix(gr2))
  max.gr2 <- max(abs(gr2))
  n.i.i <- n.i.i +1
  grad.similar <- isTRUE(all.equal(norm.gr1, norm.gr2, tolerance = 1,   scale = 1))
  grad.better  <- !grad.similar && norm.gr2 < norm.gr1  # meaningfully better gradient (diff > 1)
  logL.better  <- fit$logL > fitFinal$logL
  logL.similar <- isTRUE(all.equal(fit$logL, fitFinal$logL, tolerance = 0.5, scale = 1))
  accept <- (logL.better  && (grad.better || grad.similar)) || # better logL with non-worse gradient
            (logL.similar && grad.better)                      # similar logL but meaningfully better gradient

  if((n.i==1 || ((is.nan(norm.gr1) && !is.nan(norm.gr2)) || !is.nan(norm.gr2) && accept))  && is.finite(fit$logL)){
    n.i.i <- 0
    fitFinal <- fit
    #Store the seed that gave the best results, so that we may reproduce results, even if a seed was not explicitly provided
    fitFinal$seed <- seed[n.i]
    if(args$trace & n.i>1){
      cat("  -> Improved: logL =", fit$logL, "| grad norm =", norm.gr2, "\n")
    }else if(args$trace){
      cat("  -> logL =", fit$logL, "| grad norm =", norm.gr2, "| max grad =", max.gr2, "\n")
    }
  }
  # else if(((is.nan(norm.gr1) && !is.nan(norm.gr2)) || !is.nan(norm.gr2))  && is.finite(fit$logL)){
  #   if(args$trace & n.i>1){
  #     cat("  -> Rejected: logL =", fit$logL, "| grad norm =", norm.gr2, "| max grad =", max.gr2, "\n")
  #   }
  # }

  if(n.i.i>=args$n.init.max){
    n.init <- n.i
    message("n.init.max reached after ", n.i, " iterations.")
    break
  }
  }else{
    norm.gr2 <- NaN
  }
  
n.i <- n.i+1;

}

}else{
  if(!is.null(seed)) set.seed(seed)
  
  if(model == "gllvm.TMB"){
    fitFinal <- do.call(gllvm.TMB, args)
  }else if(model == "trait.TMB"){
    fitFinal <- do.call(trait.TMB, args)
  }
  fitFinal$seed = seed
}

# In case gllvm.TMB or trait.TMB fail entirely
if(is.null(fitFinal)){
  fitFinal$logL <- -Inf
}

return(fitFinal)
}
