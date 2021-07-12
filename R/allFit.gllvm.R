#' @title Re-fit a GLLVM with different starting values or optimizers
#' @description The algorithm used is sensitive to the initial values, so that a single model fit might not always provide the best parameter estimates. Frequent re-fitting might be required to ensure the best solution has been found. This function re-fits GLLVM objects with different starting values or optimizers, for convenience, to ensure an optimal solution.
#' 
#' @param object an object of class 'gllvm'
#' @param starting.vals should the model be re-fitted with different starting values? Defaults to \code{TRUE}.
#' @param optimizers should the model be re-fitted with different optimizers? Defaults to \code{TRUE}.
#' @param return.best return only the best fitting object (\code{TRUE})? Or all fitted objects \code{FALSE}? 
#' @param return.table return a table with information of the fitted model objects? Defaults to \code{FALSE}.
#' @param ...	 other arguments for a \code{\link{gllvm}} function call.
#'
#' @author Bert van der Veen
#'@seealso  \code{\link{gllvm}}.
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
#'@export print.summary.gllvm 


allFit.gllvm <- function(object, starting.vals = TRUE, optimizers = TRUE, return.best = TRUE, return.table = FALSE, ...){
  
  fits <- list()
  suppressWarnings(
    {
  if(starting.vals&!optimizers){
    fits[[1]] <- try(update(object, starting.val = "zero", ...),silent = TRUE)
    fits[[2]] <- try(update(object, starting.val = "res", ...),silent = TRUE)
    fits[[3]] <- try(update(object, starting.val = "res", n.init = 3, seed = NULL, ...),silent=TRUE)
    fits[[4]] <- try(update(object, starting.val = "random"),silent=TRUE)
    if(inherits(fits[[1]], "try-error"))warning("Fit with starting.val=='zero' failed. \n")
    if(inherits(fits[[2]], "try-error"))warning("Fit with starting.val=='res' failed. \n")
    if(inherits(fits[[3]], "try-error"))warning("Fit with starting.val=='res' and n.init = 3 failed. \n")
    if(inherits(fits[[3]], "try-error"))warning("Fit with starting.val=='random' failed. \n")
  }
  if(optimizers&!starting.vals){
    fits[[1]] <- try(update(object, optimizer = "optim", ...),silent=TRUE)
    fits[[2]] <- try(update(object, optimizer = "nlminb", ...),silent=TRUE)
    fits[[3]] <- try(update(object, optimizer = "nlm", ...),silent=TRUE)
    if(inherits(fits[[1]], "try-error"))warning("Fit with optimizer = 'optim' failed. \n")
    if(inherits(fits[[2]], "try-error"))warning("Fit with optimizer = 'nlminb' failed. \n")
    if(inherits(fits[[3]], "try-error"))warning("Fit with optimizer = 'nlm' failed. \n")
  }
  if(optimizers&starting.vals){
    #optim
      fits[[1]] <- try(update(object, starting.val = "zero", optimizer = "optim", ...),silent = TRUE)
      fits[[2]] <- try(update(object, starting.val = "res", optimizer = "optim", ...),silent = TRUE)
      fits[[3]] <- try(update(object, starting.val = "res", optimizer = "optim", n.init = 3, seed = NULL, ...),silent=TRUE)
      fits[[4]] <- try(update(object, starting.val = "random", optimizer = "optim", ...),silent=TRUE)
     #nlminb
      fits[[5]] <- try(update(object, starting.val = "zero", optimizer = "nlminb", ...),silent = TRUE)
      fits[[6]] <- try(update(object, starting.val = "res", optimizer = "nlminb", ...),silent = TRUE)
      fits[[7]] <- try(update(object, starting.val = "res", optimizer = "nlminb", n.init = 3, seed = NULL, ...),silent=TRUE)
      fits[[8]] <- try(update(object, starting.val = "random", optimizer = "nlminb", ...),silent=TRUE)
     #nlm
      fits[[9]] <- try(update(object, starting.val = "zero", optimizer = "nlm", ...),silent = TRUE)
      fits[[10]] <- try(update(object, starting.val = "res", optimizer = "nlm", ...),silent = TRUE)
      fits[[11]] <- try(update(object, starting.val = "res", optimizer = "nlm", n.init = 3, seed = NULL),silent=TRUE)
      fits[[12]] <- try(update(object, starting.val = "random", optimizer = "nlm", ...),silent=TRUE)
      
      #optim
      if(inherits(fits[[1]], "try-error"))warning("Fit with starting.val=='zero' and optimizer = 'optim' failed. \n")
      if(inherits(fits[[2]], "try-error"))warning("Fit with starting.val=='res' and optimizer = 'optim' failed. \n")
      if(inherits(fits[[3]], "try-error"))warning("Fit with starting.val=='res', and n.init = 3, and optimizer = 'optim' failed. \n")
      if(inherits(fits[[4]], "try-error"))warning("Fit with starting.val=='random' and optimizer = 'optim' failed \n")
      
      #nlminb
      if(inherits(fits[[5]], "try-error"))warning("Fit with starting.val=='zero' and optimizer = 'nlminb' failed. \n")
      if(inherits(fits[[6]], "try-error"))warning("Fit with starting.val=='res' and optimizer = 'nlminb' failed. \n")
      if(inherits(fits[[7]], "try-error"))warning("Fit with starting.val=='res', and n.init = 3, and optimizer = 'nlminb' failed. \n")
      if(inherits(fits[[8]], "try-error"))warning("Fit with starting.val=='random' and optimizer = 'nlminb' failed \n")
      
      #nlm
      if(inherits(fits[[9]],"try-error"))warning("Fit with starting.val=='zero' and optimizer = 'nlm' failed. \n")
      if(inherits(fits[[10]],"try-error"))warning("Fit with starting.val=='res' and optimizer = 'nlm' failed. \n")
      if(inherits(fits[[11]],"try-error"))warning("Fit with starting.val=='res', and n.init = 3, and optimizer = 'nlm' failed. \n")
      if(inherits(fits[[12]],"try-error"))warning("Fit with starting.val=='random' and optimizer = 'nlm' failed \n")
      
    
  }
  
}
)
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
    cat("Best fit has log-likelihood:", paste(signif(info[1,2]),",", sep=""), "optimizer:", paste(info[1,3], ",", sep=""), "starting value: ", info[1,4])
  
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