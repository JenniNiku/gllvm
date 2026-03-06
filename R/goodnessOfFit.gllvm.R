#' @title Goodness of fit measures for a gllvm
#' @description Several goodness-of-fit measure are currently available and can be calculated for a gllvm model fit and predicted values.
#'
#' @param object an object of class 'gllvm', to calculate goodness of a model fit.
#' @param y a response matrix of new observations
#' @param pred predicted values for response matrix y if you want to calculate prediction accuracy for new values. Note that for ordinal model, you need to give the predicted classes.
#' @param measure a goodness-of-fit measure to be calculated. Options are \code{"cor"} (correlation between observed and predicted values), \code{"scor"} (Spearman correlation between observed and predicted values), \code{"RMSE"} (root mean squared error of prediction), \code{"MAE"} (Mean Absolute Error), \code{"MARNE"} (Mean Absolute Range Normalized Error), \code{"TjurR2"} (Tjur's R2 measure, only for binary data), \code{"R2"} (R-squared as the square of the correlation), "AUC", \code{"sR2"} (R-squared as the square of the spearman correlation). Likelihood based pseudo R2 meaures \code{"NagelkerkeR2"}, \code{"McFaddenR2"}, \code{"CoxSnellR2"} can be calculated currently only for training data to measure the model's goodness of fit for full data, not response specific.
#' @param species logical, if \code{TRUE}, goodness-of-fit measures are calculated for each species separately. If FALSE,  goodness-of-fit measures are calculated for all species together.
#'
#' @details
#' goodnessOfFit is used for evaluating the goodness-of-fit of a model or predictions. Available goodness-of-fit measures are correlation, RMSE, MARNE, and R2 measures. Definitions are below.
#' Denote an observed response j (species) at sample i, \eqn{i=1,...,n}, as \eqn{y_{ij}}, and predicted value as \eqn{\hat y_{ij}}.
#' 
#' \deqn{RMSE(\boldsymbol{y_{j}}, \boldsymbol{\hat y_{j}}) =  \sqrt{\frac{1}{n}\Sigma_{i=1}^{n} {(y_{ij} - \hat y_{ij})^2}} }
#' 
#' \deqn{MAE(\boldsymbol{y_{j}}, \boldsymbol{\hat y_{j}}) =  \frac{1}{n}\Sigma_{i=1}^{n} |y_{ij} - \hat y_{ij}| }
#' 
#' \deqn{MARNE(\boldsymbol{y_{j}}, \boldsymbol{\hat y_{j}}) =  \frac{1}{n}\Sigma_{i=1}^{n} \frac{|y_{ij} - \hat y_{ij}|}{max(\boldsymbol{y_{j}}) - min(\boldsymbol{y_{j}})} }
#' 
#' \deqn{Tjur's R2(\boldsymbol{y_{j}}, \boldsymbol{\hat y_{j}}) =  \frac{1}{n_1}\Sigma \hat y_{ij}\boldsymbol{1}_{y=1}(y_{ij}) - \frac{1}{n_0}\Sigma \hat y_{ij}\boldsymbol{1}_{y=0}(y_{ij}) }
#' 
#' 
#' 
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>
#'
#'@seealso \code{\link{gllvm}}, \code{\link{predict.gllvm}}
#' @examples
#' \dontrun{
#'# Fit gllvm model with Poisson family
#'data(microbialdata)
#'X <- microbialdata$Xenv
#'y <- microbialdata$Y[, order(colMeans(microbialdata$Y > 0), 
#'                      decreasing = TRUE)[21:40]]
#'fit <- gllvm(y, X, formula = ~ pH + Phosp, family = poisson())
#'# Calculate metrics
#'goodnessOfFit(object = fit, measure = c("cor", "RMSE"))
#'
#'}
#'@export
goodnessOfFit <- function(object = NULL, y = NULL, pred = NULL, measure = c("cor", "RMSE", "MAE", "MARNE"), species = FALSE){
  mispred <- missing(pred) # for AUC
  
  if(is.null(pred)){
    if(is.null(object)) stop("If 'pred' is not given the model fit for 'object' needs to be given.")
    if(all(object$family == "ordinal")){
      pred <- predict(object, type = "class")
    } else {
      pred <- predict(object, type = "response")
      if(any(object$family == "ordinal")){
        pred <- pred[1,,]
        pred[,object$family == "ordinal"] <- predict(object, type = "class")[,object$family == "ordinal"]
      }
    }
  }
  if(is.null(y)){
    if(is.null(object)) stop("If 'y' is not given the model fit for 'object' need to be given.")
    y <- object$y
    if(any(object$family == "betaH")){
      yH01 <- (y>0)*1; colnames(yH01) <- paste("H01", colnames(y), sep = "_")
      y <- cbind(y, (yH01))
    }
  }else if(!is.matrix(y)){
    try(y <- as.matrix(y))
  }
  n <- NROW(y)
  p <- NCOL(y)
  
  out <- list()
  if("cor" %in% measure) {
    if(species) {
      out$cor <- rep(NA,p)
      for (j in 1:p) {
        out$cor[j] <- cor(na.omit(cbind(y[,j], pred[,j])))[2,1]
      }
    } else {
      out$cor <- cor(na.omit(cbind(unlist(c(y)), unlist(c(pred)))))[2,1]
    }
  }
  if("scor" %in% measure) {
    if(species) {
      out$scor <- rep(NA,p)
      for (j in 1:p) {
        out$scor[j] <- cor(na.omit(cbind(y[,j], pred[,j])), method = "spearman")[2,1]
      }
    } else {
      out$scor <- cor(na.omit(cbind(unlist(c(y)), unlist(c(pred)))), method = "spearman")[2,1]
    }
  }
  if("RMSE" %in% measure) {
    RMSE <- function(y,pred){ sqrt(mean((y- pred)^2, na.rm = TRUE))}
    if(species) {
      for (j in 1:p) {
        out$RMSE[j] <- RMSE(y[,j], pred[,j])
      }
    } else {
      out$RMSE <- RMSE(unlist(c(y)), unlist(c(pred)))
    }
  }

  if("MAE" %in% measure) {
    if(species) {
      for (j in 1:p) {
        out$MAE[j] <- mean(abs(y[,j]- pred[,j]), na.rm = TRUE)
      }
    } else {
      out$MAE <- mean(abs(y- pred), na.rm = TRUE)
    }
  }
  if("MARNE" %in% measure) {
    MARNE <- function(y,pred){ mean(abs(y- pred)/(max(y, na.rm = TRUE) - min(y, na.rm = TRUE)), na.rm = TRUE)}
    if(species) {
      for (j in 1:p) {
        out$MARNE[j] <- MARNE(y[,j], pred[,j])
      }
    } else {
      out$MARNE <- mean(t(abs(y- pred))/(apply(y,2,max, na.rm = TRUE) - apply(y,2,min, na.rm = TRUE)), na.rm = TRUE)
    }
  }
  if("TjurR2" %in% measure) {
    if(all(unique(y) %in% c(0,1,NA))){
      tjurR2 <- function(y,pred){mean(pred[y==1], na.rm = TRUE) - mean(pred[y==0], na.rm = TRUE)}
      if(species) {
        for (j in 1:p) {
          out$TjurR2[j] <- tjurR2(y[,j], pred[,j])
        }
      } else {
        out$TjurR2 <- tjurR2(unlist(c(y)), unlist(c(pred)))
      }
    }
  }
  if("R2" %in% measure) {
    if(species) {
      out$R2 <- rep(NA,p)
      for (j in 1:p) {
        out$R2[j] <- cor(na.omit(cbind(y[,j], pred[,j])))[2,1]
        out$R2[j] <- sign(out$R2[j])*out$R2[j]^2
      }
    } else {
      out$R2 <- cor(na.omit(cbind(unlist(c(y)), unlist(c(pred)))))[2,1]
      out$R2 <- sign(out$R2)*out$R2^2
    }
  }
  if("sR2" %in% measure) {
    if(species) {
      out$sR2 <- rep(NA,p)
      for (j in 1:p) {
        out$sR2[j] <- cor(na.omit(cbind(y[,j], pred[,j])), method = "spearman")[2,1]
        out$sR2[j] <- sign(out$sR2[j])*out$sR2[j]^2
      }
    } else {
      out$sR2 <- cor(na.omit(cbind(unlist(c(y)), unlist(c(pred)))), method = "spearman")[2,1]
      out$sR2 <- sign(out$sR2)*out$sR2^2
    }
  }
  if("AUC" %in% measure){
    if(any(object$family%in%c("gamma","exponential","beta")))warning("AUC only makes sense for families that include absence.")
    if(any(object$family %in% c("ordinal"))){
      # just making sure the minimum is 0
      y[,object$family == "ordinal"] <- object$y[,object$family == "ordinal"]-apply(object$y[,object$family == "ordinal"],2,min) 
    }
    predAUC <- pred
    if(mispred){
      predAUC <- predict(object, type = "response")
    }else if(!mispred && any(object$family == "ordinal") && length(dim(pred))!=3){
      stop("To calculate AUC, 'pred' should  be an array with 3 dimensions when ordinal responses are included.")
    }
    
    predAUC <- gllvm.presence.prob(predAUC, object)
    
    if(species){
      out$AUC <- rep(NA, p)
      for(j in 1:p){
        ranks <- rank(predAUC[,j])
        n_pos <- sum(y[,j] > 0)
        n_neg <- sum(y[,j] == 0)
        sum_ranks_pos <- sum(ranks[y[,j] >0])
        out$AUC[j] <-  (sum_ranks_pos - n_pos*(n_pos + 1)/2) / (n_pos * n_neg) 
      }
    }else{
      ranks <- rank(c(predAUC))
      n_pos <- sum(c(y) > 0)
      n_neg <- sum(c(y) == 0)
      sum_ranks_pos <- sum(ranks[c(y) >0])
      out$AUC <-  (sum_ranks_pos - n_pos*(n_pos + 1)/2) / (n_pos * n_neg) 
    }
  
  }
  if(any(c("NagelkerkeR2", "McFaddenR2", "CoxSnellR2") %in% measure) & !is.null(object)) {
    #Fit null model to calculate 
    nob <- nobs(object)
    if(any(object$family == "ordinal")){
      modelnull<- update(object, formula = ~1, lv.formula = NULL, row.eff = NULL, num.lv=0, num.lv.c=0, num.RR=0, sd.errors = FALSE, starting.val="zero", zeta.struc=object$zeta.struc)
    } else {
      modelnull<- update(object, formula = ~1, lv.formula = NULL, row.eff = NULL, num.lv=0, num.lv.c=0, num.RR=0, sd.errors = FALSE, starting.val="zero")
    }
    logLik_full <- object$logL
    logLik_null <- modelnull$logL
    if("McFaddenR2"%in% measure) {
      out$McFaddenR2 <- 1 - logLik_full/logLik_null
    }
    if("CoxSnellR2"%in% measure) {
      out$CoxSnellR2 <- 1 - exp((2 / nob) * (logLik_null - logLik_full))
    }
    if("NagelkerkeR2"%in% measure) {
      out$NagelkerkeR2 <- (1 - exp((2 /nob) * (logLik_null - logLik_full)))/(1 - exp((2 / nob) * (logLik_null)))
    }
  }
  return(out)
}