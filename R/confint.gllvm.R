#' @title Confidence intervals for model parameters
#' @description Computes confidence intervals for  parameters in a fitted gllvm model.
#'
#' @param object an object of class 'gllvm'.
#' @param level the confidence level. Scalar between 0 and 1.
#' @param parm a specification of which parameters are to be given confidence intervals, a vector of names. Examples of options are "beta0", "Xcoef",theta", "phi". If missing, all parameters are considered.
#' @param ...	not used.
#'
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#' @examples
#' \dontrun{
#'## Load a dataset from the mvabund package
#'data(antTraits, package = "mvabund")
#'y <- as.matrix(antTraits$abund)
#'X <- as.matrix(antTraits$env[,1:2])
#'# Fit gllvm model
#'fit <- gllvm(y = y, X = X, family = poisson())
#'# 95 % confidence intervals for coefficients of X variables
#'confint(fit, level = 0.95, parm = "Xcoef")
#'}
#'@aliases confint confint.gllvm
#'@method confint gllvm
#'@importFrom stats confint
#'
#'@export
#'@export confint.gllvm

confint.gllvm <- function(object, parm=NULL, level = 0.95, ...) {
  if(is.logical(object$sd)) stop("Standard errors for parameters haven't been calculated, so confidence intervals can not be calculated.");
  n <- NROW(object$y)
  p <- NCOL(object$y)
  nX <- 0; if(!is.null(object$X)) nX <- dim(object$X.design)[2]
  nTR <- 0; if(!is.null(object$TR)) nTR <- dim(object$TR)[2]
  num.lv <- object$num.lv
  num.lv.c <- object$num.lv.c
  num.RR <- object$num.RR
  quadratic <- object$quadratic
  if(!is.null(object$lv.X) && is.null(object$lv.X.design))object$lv.X.design <- object$lv.X #for backward compatibility
  if("Intercept"%in%row.names(object$sd$B))object$sd$B<-object$sd$B[-which(row.names(object$sd$B)=="Intercept")]
  alfa <- (1 - level) / 2
  if(object$row.eff == "random") object$params$row.params = NULL
  if(is.null(object$col.eff$col.eff))object$col.eff$col.eff <- FALSE # backward compatibility
  if(object$col.eff$col.eff == "random" | !is.null(object$randomX))object$params$Br <- NULL
  if(object$beta0com){
    object$params$beta0 <- unique(object$params$beta0)
    names(object$params$beta0) <- "Community intercept"
    object$sd$beta0 <- unique(object$sd$beta0)
  }
  if(is.null(parm)){
    if (object$family == "negative.binomial") {
      object$params$phi <- NULL
      object$sd$phi <- NULL
    }
    if (object$family == "ZINB") {
      object$params$ZINB.phi <- NULL
      object$sd$ZINB.phi <- NULL
    }
    
    if (!is.null(object$params$sigmaB)) {
      object$params$sigmaB <- sqrt(diag(object$params$sigmaB))
      if(object$col.eff$col.eff=="random"){
        object$sd$sigmaB <- diag(object$sd$sigmaB)
      }
      object$sd$corrpar <- NULL
    }
    
    parm_all <- c("sigma.lv","theta", "LvXcoef","beta0", "Xcoef", "B", "row.params", "sigma", "sigmaB", "sigmaLvXcoef", "inv.phi", "phi", "ZINB.phi", "ZINB.inv.phi" ,"p","zeta", "rho.sp")
    if(object$randomB!=FALSE){
      object$params$LvXcoef <- NULL
    }
    
    if(object$family=="ZINB" && "inv.phi" %in% parmincl)parmincl[parmincl=="inv.phi"]<-"ZINB.inv.phi"
    
    parmincl <- parm_all[parm_all %in% names(object$params)]
    if("rho.sp"%in%names(object$params[parmincl])){
      object$params$rho.sp <- log(-log(object$params$rho.sp))
    }
    
    cilow <- unlist(object$params[parmincl]) + qnorm(alfa) * unlist(object$sd[parmincl])
    ciup <- unlist(object$params[parmincl]) + qnorm(1 - alfa) * unlist(object$sd[parmincl])
    
    if("rho.sp"%in%names(object$params[parmincl])){
      a <- exp(-exp(cilow[grepl("rho.sp", names(cilow))])); b <- exp(-exp(ciup[grepl("rho.sp", names(ciup))]))
      cilow[grepl("rho.sp", names(cilow))] <- pmin(a,b)
      ciup[grepl("rho.sp", names(ciup))] <- pmax(a,b)
    }
    
    M <- cbind(cilow, ciup)
    
    colnames(M) <- c(paste(alfa * 100, "%"), paste((1 - alfa) * 100, "%"))
    rnames <- names(unlist(object$params[parmincl]))
    
    cal <- 0
    if ((num.lv.c+num.RR) > 0 & num.lv == 0) {
      nr <- rep(1:(num.lv.c+num.RR), each = p)
      nc <- rep(1:p, (num.lv.c+num.RR))
      if(quadratic == FALSE)
        rnames[grepl("theta",rnames)][1:((num.lv.c+num.RR) * p)] <- paste(paste("theta.CLV", nr, sep = ""), nc, sep = ".")
      if(quadratic != FALSE)
        rnames[grepl("theta",rnames)][1:((num.lv.c+num.RR) * p *2)] <- c(paste(paste("theta.CLV", nr, sep = ""), nc, sep = "."), paste(paste("theta.CLV", nr, "^2",
                                                                                                                      sep = ""
        ), nc, sep = "."))
      if(quadratic==FALSE)cal <- cal + (num.lv.c+num.RR) * p
      if(quadratic!=FALSE)cal <- cal + (num.lv.c+num.RR) * p * 2
    }
    if (num.lv > 0 & (num.lv.c+num.RR) ==0) {
      nr <- rep(1:num.lv, each = p)
      nc <- rep(1:p, num.lv)
      if(quadratic == FALSE)
        rnames[grepl("theta",rnames)][1:(num.lv * p)] <- paste(paste("theta.LV", nr, sep = ""), nc, sep = ".")
      if(quadratic != FALSE)
        rnames[grepl("theta",rnames)][1:(num.lv * p *2)] <- c(paste(paste("theta.LV", nr, sep = ""), nc, sep = "."), paste(paste("theta.LV", nr, "^2",
                                                                                                          sep = ""
        ), nc, sep = "."))
      if(quadratic==FALSE)cal <- cal + num.lv * p
      if(quadratic!=FALSE)cal <- cal + num.lv * p * 2
    }
    if ((num.lv.c+num.RR) > 0 & num.lv>0) {
      nr <- rep(1:(num.lv.c+num.RR), each = p)
      nc <- rep(1:p, (num.lv.c+num.RR))
      if(quadratic == FALSE)
        rnames[grepl("theta",rnames)][1:((num.lv.c+num.RR) * p)] <- paste(paste("theta.CLV", nr, sep = ""), nc, sep = ".")
      if(quadratic != FALSE)
        rnames[grepl("theta",rnames)][1:((num.lv.c+num.RR) * p *2)] <- c(paste(paste("theta.CLV", nr, sep = ""), nc, sep = "."), paste(paste("theta.CLV", nr, "^2",
                                                                                                                      sep = ""
        ), nc, sep = "."))
      if(quadratic==FALSE)cal <- cal + (num.lv.c+num.RR) * p
      if(quadratic!=FALSE)cal <- cal + (num.lv.c+num.RR) * p * 2
      
      
      nr <- rep(1:num.lv, each = p)
      nc <- rep(1:p, num.lv)
      if(quadratic == FALSE)
        rnames[grepl("theta",rnames)][((num.lv.c+num.RR)*p+1):((num.lv.c+num.RR)*p+num.lv * p)] <- paste(paste("theta.LV", nr, sep = ""), nc, sep = ".")
      if(quadratic != FALSE)
        rnames[grepl("theta",rnames)][((num.lv.c+num.RR)*p+1):((num.lv.c+num.RR)*p+num.lv * p)] <- c(paste(paste("theta.LV", nr, sep = ""), nc, sep = "."), paste(paste("theta.LV", nr, "^2",
                                                                                                                                                 sep = ""
        ), nc, sep = "."))
      if(quadratic==FALSE)cal <- cal + num.lv * p
      if(quadratic!=FALSE)cal <- cal + num.lv * p * 2
    }
    if((num.lv.c+num.RR)>0&object$randomB==FALSE){
      rnames[-c(1:cal)][1:(ncol(object$lv.X.design)*(num.lv.c+num.RR))] <- paste(colnames(object$lv.X.design),"LV",rep(1:(num.lv.c+num.RR),each=ncol(object$lv.X.design)),sep=".")
      cal<-cal + ncol(object$lv.X.design)*(num.lv.c+num.RR)
    }
    
    if(!object$beta0com){
      rnames[(cal + 1):(cal + p)] <- paste("Intercept",names(object$params$beta0), sep = ".")
      cal <- cal + p
    } else {
      rnames[(cal + 1)] <- paste("Intercept",names(object$params$beta0), sep = ".")
      cal <- cal + 1
    }
    if (!is.null(object$TR) | object$col.eff$col.eff == "random") {
      nr <- names(object$params$B)
      rnames[(cal + 1):(cal + length(nr))] <- nr
      cal <- cal + length(nr)
    }
    
    if (is.null(object$TR) && !is.null(object$X)) {
      cnx <- rep(colnames(object$X.design), each = p)
      rnc <- rep(rownames(object$params$Xcoef), nX)
      newnam <- paste(cnx, rnc, sep = ":")
      rnames[(cal + 1):(cal + nX * p)] <- paste("Xcoef", newnam, sep = ".")
      cal <- cal + nX * p
    }
    if (object$row.eff %in% c("fixed",TRUE)) {
      rnames[(cal + 1):(cal + n)] <- paste("Row.Intercept", 1:n, sep = ".")
      cal <- cal + n
    }
    if (object$row.eff == "random") {
      rnames[(cal + 1):(cal+length(object$sd$sigma))] <- "sigma"
      cal <- cal + length(object$sd$sigma)
      # if(!is.null(object$params$rho)) {
      #   rnames[(cal + 1)] <- "rho"
      #   cal <- cal + length(object$sd$rho)
      # }
    }
    if (!is.null(object$randomX) | object$col.eff$col.eff == 'random') {
      cal <- cal + length(object$params$sigmaB)
    }
    if (object$randomB!=FALSE) {
      cal <- cal + length(object$params$sigmaLvXcoef)
    }
    if(object$family == "negative.binomial"){
      s <- length(unique(object$disp.group))
      rnames[(cal + 1):(s+cal)] <- paste("inv.phi", names(object$params$inv.phi), sep = ".")
    }
    
    if(object$family == "ZINB"){
      s <- length(unique(object$disp.group))
      rnames[(cal + 1):(s+cal)] <- paste("inv.phi", names(object$params$ZINB.inv.phi), sep = ".")
    }
    
    if(object$family == "tweedie"){
      s <-length(unique(object$disp.group))
      rnames[(cal + 1):(s+cal)] <- paste("Dispersion phi", names(object$params$phi), sep = ".")
    }
    if(object$family %in%c("ZIP","ZINB")){
      s <- length(unique(object$disp.group))
      rnames[(cal + 1):(s+cal)] <- paste("p", names(object$params$p), sep = ".")
    }
    if(object$family == "gaussian"){
      s <- length(unique(object$disp.group))
      rnames[(cal + 1):(s+cal)] <- paste("Standard deviations phi", names(object$params$phi), sep = ".")
    }
    if(object$family == "gamma"){
      s <- length(unique(object$disp.group))
      rnames[(cal + 1):(s+cal)] <- paste("Shape phi", names(object$params$phi), sep = ".")
    }
    
    rownames(M) <- rnames
  } else {
    if ("beta0" %in% parm) {
      object$params$Intercept = object$params$beta0
      object$sd$Intercept = object$sd$beta0
      parm[parm=="beta0"] = "Intercept"
    }
    
    cilow <- unlist(object$params[parm]) + qnorm(alfa) * unlist(object$sd[parm])
    ciup <- unlist(object$params[parm]) + qnorm(1 - alfa) * unlist(object$sd[parm])
    
    if("theta"%in%parm){
      
      if(num.lv>0&(num.lv.c+num.RR)==0){
        if(object$quadratic==FALSE){
          names(cilow)[gsub("theta.*","theta",names(unlist(object$sd[parm])))%in%"theta"] <- paste(rep(paste("theta.LV", 1:num.lv, sep = ""),each=p),1:p,sep=".")
          names(ciup)[gsub("theta.*","theta",names(unlist(object$sd[parm])))%in%"theta"] <-   paste(rep(paste("theta.LV", 1:num.lv, sep = ""),each=p),1:p,sep=".")
        }else{
          names(cilow)[gsub("theta.*","theta",names(unlist(object$sd[parm])))%in%"theta"] <-c(paste(rep(paste("theta.LV", 1:num.lv, sep = ""),each=p),1:p,sep="."),paste(rep(paste("theta.LV", 1:num.lv, "^2",sep = ""),each=p),1:p,sep="."))
          names(ciup)[gsub("theta.*","theta",names(unlist(object$sd[parm])))%in%"theta"] <-   c(paste(rep(paste("theta.LV", 1:num.lv, sep = ""),each=p),1:p,sep="."),paste(rep(paste("theta.LV", 1:num.lv, "^2",sep = ""),each=p),1:p,sep="."))
        }
      }else if(num.lv==0&(num.lv.c+num.RR)>0){
        if(object$quadratic==FALSE){
          names(cilow)[gsub("theta.*","theta",names(unlist(object$sd[parm])))%in%"theta"] <- paste(rep(paste("theta.CLV", 1:(num.lv.c+num.RR), sep = ""),each=p),1:p,sep=".")
          names(ciup)[gsub("theta.*","theta",names(unlist(object$sd[parm])))%in%"theta"] <-   paste(rep(paste("theta.CLV", 1:(num.lv.c+num.RR), sep = ""),each=p),1:p,sep=".")
        }else{
          names(cilow)[gsub("theta.*","theta",names(unlist(object$sd[parm])))%in%"theta"] <-c(paste(rep(paste("theta.CLV", 1:(num.lv.c+num.RR), sep = ""),each=p),1:p,sep="."),paste(rep(paste("theta.CLV", 1:(num.lv.c+num.RR), "^2",sep = ""),each=p),1:p,sep="."))
          names(ciup)[gsub("theta.*","theta",names(unlist(object$sd[parm])))%in%"theta"] <-   c(paste(rep(paste("theta.CLV", 1:(num.lv.c+num.RR), sep = ""),each=p),1:p,sep="."),paste(rep(paste("theta.CLV", 1:(num.lv.c+num.RR), "^2",sep = ""),each=p),1:p,sep="."))
        }
      }else if(num.lv>0&(num.lv.c+num.RR)>0){
        if(object$quadratic==FALSE){
          names(cilow)[gsub("theta.*","theta",names(unlist(object$sd[parm])))%in%"theta"] <- c(paste(rep(paste("theta.CLV", 1:(num.lv.c+num.RR), sep = ""),each=p),1:p,sep="."),paste(rep(paste("theta.LV", 1:num.lv, sep = ""),each=p),1:p,sep="."))
          names(ciup)[gsub("theta.*","theta",names(unlist(object$sd[parm])))%in%"theta"] <-   c(paste(rep(paste("theta.CLV", 1:(num.lv.c+num.RR), sep = ""),each=p),1:p,sep="."),paste(rep(paste("theta.LV", 1:num.lv, sep = ""),each=p),1:p,sep="."))
        }else{
          names(cilow)[gsub("theta.*","theta",names(unlist(object$sd[parm])))%in%"theta"] <-c(c(paste(rep(paste("theta.CLV", 1:(num.lv.c+num.RR), sep = ""),each=p),1:p,sep="."),paste(rep(paste("theta.LV", 1:num.lv, sep = ""),each=p),1:p,sep=".")),paste(rep(paste("theta.CLV", 1:(num.lv.c+num.RR), "^2",sep = ""),each=p),1:p,sep="."),paste(rep(paste("theta.LV", 1:num.lv, "^2",sep = ""),each=p),1:p,sep="."))
          names(ciup)[gsub("theta.*","theta",names(unlist(object$sd[parm])))%in%"theta"] <- c(c(paste(rep(paste("theta.CLV", 1:(num.lv.c+num.RR), sep = ""),each=p),1:p,sep="."),paste(rep(paste("theta.LV", 1:num.lv, sep = ""),each=p),1:p,sep=".")),paste(rep(paste("theta.CLV", 1:(num.lv.c+num.RR), "^2",sep = ""),each=p),1:p,sep="."),paste(rep(paste("theta.LV", 1:num.lv, "^2",sep = ""),each=p),1:p,sep="."))
        }
      }
      
      
    }
    if("LvXcoef"%in%parm&object$randomB==FALSE){
      names(cilow)[gsub("LvXcoef.*","LvXcoef",names(unlist(object$sd[parm])))%in%"LvXcoef"] <-paste(rep(colnames(object$lv.X.design),2),"LV",rep(1:(num.lv.c+num.RR),each=ncol(object$lv.X.design)),sep=".")
      names(ciup)[gsub("LvXcoef.*","LvXcoef",names(unlist(object$sd[parm])))%in%"LvXcoef"] <- paste(rep(colnames(object$lv.X.design),2),"LV",rep(1:(num.lv.c+num.RR),each=ncol(object$lv.X.design)),sep=".")
    }
    if("sigmaLvXcoef"%in%parm&object$randomB!=FALSE){
      if(object$randomB=="LV"){
        names(cilow)[gsub("sigmaLvXcoef.*","sigmaLvXcoef",names(unlist(object$sd[parm])))%in%"sigmaLvXcoef"]<-  c(paste("sigmaLvXcoef.LV", 1:(num.RR+num.lv.c), sep=""))
      }else if(object$randomB=="P"){
        names(cilow)[gsub("sigmaLvXcoef.*","sigmaLvXcoef",names(unlist(object$sd[parm])))%in%"sigmaLvXcoef"]<-  c(paste("sigmaLvXcoef.", colnames(object$lv.X.design), sep=""))
      }
      
    }
    if("sigma.lv"%in%parm){
      if(num.lv>0&num.lv.c>0){
        names(cilow)[gsub("sigma.lv.*","sigma.lv",names(unlist(object$sd[parm])))%in%"sigma.lv"]<-  c(paste("sigma.CLV", 1:num.lv.c, sep=""),paste("sigma.LV", 1:num.lv, sep=""))
      }else if(num.lv>0){
        names(cilow)[gsub("sigma.lv.*","sigma.lv",names(unlist(object$sd[parm])))%in%"sigma.lv"]<-  c(paste("sigma.LV", 1:num.lv, sep=""))
      }else if(num.lv.c>0){
        names(cilow)[gsub("sigma.lv.*","sigma.lv",names(unlist(object$sd[parm])))%in%"sigma.lv"]<-  c(paste("sigma.CLV", 1:num.lv.c, sep=""))
      }
    }
    
    M <- cbind(cilow, ciup)
    
    
  }
  return(M)
}
