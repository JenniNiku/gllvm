gllvm.iter <- function(...){
  args <- list(...)
  
  if(!(args$family %in% c("poisson","negative.binomial","binomial","tweedie","ZIP", "ZINB", "gaussian", "ordinal", "gamma", "exponential", "beta", "betaH", "orderedBeta")))
    stop("Selected family not permitted...sorry!")
  if(!(args$Lambda.struc %in% c("unstructured","diagonal","bdNN","UNN")))
    stop("Lambda matrix (covariance of variational distribution for latent variable) not permitted...sorry!")
  
  if (!is.numeric(args$y))
    stop( "y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1.")
  # if ((family %in% c("ZIP")) && (method %in% c("VA", "EVA"))) #"tweedie", 
  #   stop("family=\"", family, "\" : family not implemented with VA method, change the method to 'LA'")
  if (is.null(rownames(args$y)))
    rownames(args$y) <- paste("Row", 1:nrow(args$y), sep = "")
  if (is.null(colnames(args$y)))
    colnames(args$y) <- paste("Col", 1:ncol(args$y), sep = "")
  if(args$family == "ordinal"){
    max.levels <- apply(args$y,2,function(x) length(min(x):max(x)))
    if(any(max.levels == 1)&args$zeta.struc=="species" || all(max.levels == 2)&args$zeta.struc=="species")
      stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead.")
    if(any(!apply(args$y,2,function(x)all(diff(sort(unique(x)))==1)))&args$zeta.struc=="species"){
      warning("Can't fit ordinal model if there are species with missing classes. Setting 'zeta.struc = `common`'.")
      args$zeta.struc = "common"
    }
    if(!all(min(args$y)==apply(args$y,2,min))&args$zeta.struc=="species"){
      warning("For ordinal data and zeta.struc=`species` all species must have the same minimum category. Setting 'zeta.struc = `common`'.")
      args$zeta.struc = "common"
    }
    if(any(diff(sort(unique(c(args$y))))!=1)&args$zeta.struc=="common")
      stop("Can't fit ordinal model if there are missing response classes. Please reclassify.")
  }
  
  if(args$family == "orderedBeta") {
    if (!(args$method %in% c("VA", "EVA"))) #"tweedie", 
      stop("family=\"", args$family, "\" : family not implemented with LA method, change the method to 'VA'.")
    
    if((sum(args$y==1, na.rm = TRUE) + sum(args$y==0, na.rm = TRUE))==0){
      stop("No zeros or ones in the data, so use 'family = `beta`'.")
    }
    if(!all(colSums(args$y==1, na.rm = TRUE)>0) & !all(colSums(args$y==0, na.rm = TRUE)>0)){
      warning("All species do not have zeros and ones. Setting 'zeta.struc = `common`'.")
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
    cat("Initial run ", n.i, "\n")
  
  if(!is.null(seed)) set.seed(seed[n.i])
  
  if(model == "gllvm.TMB"){
  fit <- do.call(gllvm.TMB, args)
  }else if(model == "trait.TMB"){
    fit <- do.call(trait.TMB, args)
  }
  #### Check if model fit succeeded/improved on this iteration n.i
  
  # Gradient check with n.i >2 so we don't get poorly converged models - relatively relaxed tolerance
  if(!is.null(fitFinal$TMBfn)){
    gr1 <- fitFinal$TMBfn$gr()
    gr1 <- as.matrix(gr1/length(gr1))
    norm.gr1 <- norm(gr1)
  }else{
    gr1 <- NaN
    norm.gr1 <- NaN
  }
  if(!is.null(fit$TMBfn$gr)){
  gr2 <- fit$TMBfn$gr()
  gr2 <- as.matrix(gr2/length(gr2))
  norm.gr2 <- norm(gr2)
  n.i.i <- n.i.i +1
  grad.test1 <- all.equal(norm.gr1, norm.gr2, tolerance = 1, scale = 1)#check if gradients are similar when accepting on log-likelihood
  grad.test2 <- all.equal(norm.gr1, norm.gr2, tolerance = .1, scale = 1)#check if gradient are (sufficiently) different from each other, when accepting on gradient. Slightly more strict for norm(gr2)<norm(gr1)

  if((n.i==1 || ((is.nan(norm.gr1) && !is.nan(norm.gr2)) || !is.nan(norm.gr2) && ((isTRUE(grad.test1) && fit$logL > (fitFinal$logL)) || (!isTRUE(grad.test2) && norm.gr2<norm.gr1))))  && is.finite(fit$logL)){
    n.i.i <- 0
    fitFinal <- fit
    #Store the seed that gave the best results, so that we may reproduce results, even if a seed was not explicitly provided
    fitFinal$seed <- seed[n.i]
  }

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

return(fitFinal)
}
