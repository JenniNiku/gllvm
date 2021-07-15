#' @title Re-fit a GLLVM with different starting values or optimizers
#' @description The algorithm used is sensitive to the initial values, so that a single model fit might not always provide the best parameter estimates. Frequent re-fitting might be required to ensure the best solution has been found. This function refits GLLVM objects with different starting values and optimizers, for convenience, to ensure an optimal solution. 
#' 
#' Specifically, it implements some of the recommendations from Niku et al. 2019 with three different optimizers: 1) "zero", 2) "res" (without jitter), 3) "res3" (with jitter).
#' 
#' @param object an object of class 'gllvm'.
#' @param seed seed to use for model fitting if n.init = 1.
#' @param sd.errors should standard errors be returned as part of the best fitting model? Defaults to \code{TRUE}.
#' @param starting.vals should the model be re-fitted with different starting values? Defaults to \code{TRUE}.
#' @param optimizers should the model be re-fitted with different optimizers? Defaults to \code{TRUE}.
#' @param return.best return only the best fitting object (\code{TRUE})? Or all fitted objects \code{FALSE}? 
#' @param return.table return a table with information of the fitted model objects? Defaults to \code{FALSE}.
#' @param ...	 additional arguments for a \code{\link{gllvm}} function call.
#'
#' @author Bert van der Veen
#' @references 
#' Niku, J., Brooks, W., Herliansyah, R., Hui, F. K. C., Taskinen, S., and Warton,  D. I. (2018). Efficient estimation of generalized linear latent variable models. PLoS One, 14(5):1-20.
#' @seealso  \code{\link{gllvm}}.
#' @examples
#' \dontrun{
#'## Load a dataset from the mvabund package
#' data(antTraits)
#' y <- as.matrix(antTraits$abund)
#'# Fit gllvm model
#' fit <- gllvm(y = y, family = poisson())
#' bestFit <- allFit(fit)
#'}
#'@export
#'@export allFit.gllvm 


allFit.gllvm <- function(object, seed = NULL, starting.vals = TRUE, optimizers = TRUE, return.best = TRUE, return.table = FALSE, sd.errors = TRUE, ...){
  fits <- list()
  opts <- list(...)
  #these can create issues if not taken care of
  object$call$n.init <- 1
  object$call$seed <- NULL
  n.init <- NULL
  if("n.init"%in%names(opts)){
    n.init <- opts$n.init
  }
  
  if(is.null(seed)&!is.null(n.init)){
    seed <- sample(1:10000, n.init)
  }

  if(starting.vals&!optimizers){
    suppressWarnings(
      {
    fits[[1]] <- try(update(object, starting.val = "zero", n.init = 1, seed = NULL, sd.errors = FALSE, ...),silent = TRUE)
    fits[[2]] <- try(update(object, starting.val = "res", seed = seed, sd.errors = FALSE, ...),silent = TRUE)
    fits[[3]] <- try(update(object, starting.val = "res", seed = NULL, sd.errors = FALSE, n.init = 3, jitter.var = 0.4, ...),silent=TRUE)
    fits[[4]] <- try(update(object, starting.val = "random", seed = seed, sd.errors = FALSE, ...),silent=TRUE)
      })
    if(inherits(fits[[1]], "try-error"))warning("Fit with starting.val=='zero' failed. \n")
    if(inherits(fits[[2]], "try-error"))warning("Fit with starting.val=='res' failed. \n")
    if(inherits(fits[[3]], "try-error"))warning("Fit with starting.val=='res' and n.init = 3 failed. \n")
    if(inherits(fits[[3]], "try-error"))warning("Fit with starting.val=='random' failed. \n")
  }
  if(optimizers&!starting.vals){
    suppressWarnings(
      {
    fits[[1]] <- try(update(object, optimizer = "optim", seed = seed, sd.errors = FALSE, ...),silent=TRUE)
    fits[[2]] <- try(update(object, optimizer = "nlminb", seed = seed, sd.errors = FALSE, ...),silent=TRUE)
    fits[[3]] <- try(update(object, optimizer = "nlm", seed = seed, sd.errors = FALSE, ...),silent=TRUE)
      })
    if(inherits(fits[[1]], "try-error"))warning("Fit with optimizer = 'optim' failed. \n")
    if(inherits(fits[[2]], "try-error"))warning("Fit with optimizer = 'nlminb' failed. \n")
    if(inherits(fits[[3]], "try-error"))warning("Fit with optimizer = 'nlm' failed. \n")
  }
  if(optimizers&starting.vals){
    suppressWarnings(
      {
    #optim
    #
      fits[[1]] <- try(update(object, starting.val = "zero", n.init = 1, seed = NULL, optimizer = "optim", sd.errors = FALSE, ...),silent = TRUE)
      fits[[2]] <- try(update(object, starting.val = "res", seed = seed, optimizer = "optim", sd.errors = FALSE, ...),silent = TRUE)
      fits[[3]] <- try(update(object, starting.val = "res", seed = NULL, optimizer = "optim", sd.errors = FALSE, n.init = 3, jitter.var = 0.4, ...),silent=TRUE)
      fits[[4]] <- try(update(object, starting.val = "random", seed = seed, optimizer = "optim", sd.errors = FALSE, ...),silent=TRUE)
     #nlminb
      fits[[5]] <- try(update(object, starting.val = "zero", n.init = 1, seed = NULL, optimizer = "nlminb", sd.errors = FALSE, ...),silent = TRUE)
      fits[[6]] <- try(update(object, starting.val = "res", seed = seed, optimizer = "nlminb", sd.errors = FALSE, ...),silent = TRUE)
      fits[[7]] <- try(update(object, starting.val = "res", seed = NULL, optimizer = "nlminb", sd.errors = FALSE, n.init = 3 , jitter.var = 0.4, ...),silent=TRUE)
      fits[[8]] <- try(update(object, starting.val = "random", seed = seed, optimizer = "nlminb", sd.errors = FALSE, ...),silent=TRUE)
     #nlm
      fits[[9]] <- try(update(object, starting.val = "zero", n.init = 1, seed = NULL, optimizer = "nlm", sd.errors = FALSE, ...),silent = TRUE)
      fits[[10]] <- try(update(object, starting.val = "res", seed = seed, optimizer = "nlm", sd.errors = FALSE, ...),silent = TRUE)
      fits[[11]] <- try(update(object, starting.val = "res", seed = NULL, optimizer = "nlm", sd.errors = FALSE, n.init = 3, jitter.var = 0.4, ...), silent=TRUE)
      fits[[12]] <- try(update(object, starting.val = "random", seed = seed, optimizer = "nlm", sd.errors = FALSE, ...),silent=TRUE)
      })
      #optim
      if(inherits(fits[[1]], "try-error"))warning("Fit with 'zero' and optimizer = 'optim' failed. \n")
      if(inherits(fits[[2]], "try-error"))warning("Fit with 'res' and optimizer = 'optim' failed. \n")
      if(inherits(fits[[3]], "try-error"))warning("Fit with 'res3' and optimizer = 'optim' failed. \n")
      if(inherits(fits[[4]], "try-error"))warning("Fit with 'random' and optimizer = 'optim' failed \n")
      
      #nlminb
      if(inherits(fits[[5]], "try-error"))warning("Fit with 'zero' and optimizer = 'nlminb' failed. \n")
      if(inherits(fits[[6]], "try-error"))warning("Fit with 'res' and optimizer = 'nlminb' failed. \n")
      if(inherits(fits[[7]], "try-error"))warning("Fit with 'res3' and optimizer = 'nlminb' failed. \n")
      if(inherits(fits[[8]], "try-error"))warning("Fit with 'random' and optimizer = 'nlminb' failed \n")
      
      #nlm
      if(inherits(fits[[9]],"try-error"))warning("Fit with 'zero' and optimizer = 'nlm' failed. \n")
      if(inherits(fits[[10]],"try-error"))warning("Fit with 'res' and optimizer = 'nlm' failed. \n")
      if(inherits(fits[[11]],"try-error"))warning("Fit with 'res3' and optimizer = 'nlm' failed. \n")
      if(inherits(fits[[12]],"try-error"))warning("Fit with 'random' and optimizer = 'nlm' failed \n")
  }
  
    LL <- unlist(lapply(fits,function(x){if(!inherits(x,"try-error")){logLik(x)}else{-Inf}}))
    optimizers <- unlist(lapply(fits,function(x){if(!inherits(x,"try-error")){x$optimizer}else{NA}}))
    starting.vals <- unlist(lapply(fits,function(x){if(!inherits(x,"try-error")){x$starting.val}else{NA}}))
    for(i in which(starting.vals=="zero")){
      fits[[i]]$seed <- NA
    }
    n.inits <- unlist(lapply(fits,function(x){if(!inherits(x,"try-error")){x$n.init}else{NA}}))
    if(any(n.inits>1)){
      starting.vals[which(n.inits>1)] <- paste(starting.vals[which(n.inits>1)], ", n.init = ", n.inits[n.inits>1], sep="")  
    }
    seeds <- unlist(lapply(fits,function(x){if(!inherits(x,"try-error")){x$seed}else{NA}}))
    info <- data.frame(fit = 1:length(LL),LL = LL, optimizer = optimizers, starting.val = starting.vals, best.seed = seeds)
    
    #order by best fit
    info  <- info[order(info$LL,decreasing = T),]
    if(return.best){
    cat("Best fit has log-likelihood:", paste(signif(info[1,2]),",", sep=""), "optimizer:", paste(info[1,3], ",", sep=""), "starting value: ", info[1,4])
    }else{
      cat("Best fit", paste("(#",info[1,1], ")", sep="") ,"has log-likelihood:", paste(signif(info[1,2]),",", sep=""), "optimizer:", paste(info[1,3], ",", sep=""), "starting value: ", info[1,4])
    }
    
    if(sd.errors){
      best.fit  <- fits[[which.max(unlist(lapply(fits,function(x){if(!inherits(x,"try-error")){logLik(x)}else{-Inf}})))]]
      SEs <- se.gllvm(object)
      best.fit$Hess<-SEs$Hess 
      best.fit$sd <- SEs$sd
      best.fit$call$sd.errors <- TRUE
      fits[[which.max(unlist(lapply(fits,function(x){if(!inherits(x,"try-error")){logLik(x)}else{-Inf}})))]] <- best.fit
    }
  #only return fit with best LL
  if(return.best&!return.table){
    fits[[which.max(unlist(lapply(fits,function(x){if(!inherits(x,"try-error")){logLik(x)}else{-Inf}})))]]
  }else if(!return.table){
    fits
  }else if(return.table&return.best){
    list(best.fit = fits[[which.max(unlist(lapply(fits,function(x){if(!inherits(x,"try-error")){logLik(x)}else{-Inf}})))]], table = info)
  }else if(return.table&!return.best){
    list(fits = fits[[which.max(unlist(lapply(fits,function(x){if(!inherits(x,"try-error")){logLik(x)}else{-Inf}})))]], table = info)
  }
  
}

#'@export allFit
allFit <- function(object, ...)
{
  UseMethod(generic = "allFit")
}