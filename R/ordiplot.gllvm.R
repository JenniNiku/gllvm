#' @title Plot latent variables from gllvm model
#' @description Plots latent variables and their corresponding coefficients (biplot).
#'
#' @param object   an object of class 'gllvm'.
#' @param biplot   \code{TRUE} if both latent variables and their coefficients are plotted, \code{FALSE} if only latent variables.
#' @param ind.spp  the number of response variables (usually, species) to include on the biplot. The default is none, or all if \code{biplot = TRUE}.
#' @param alpha    a numeric scalar between 0 and 1 that is used to control the relative scaling of the latent variables and their coefficients, when constructing a biplot.
#' @param main  main title.
#' @param which.lvs indices of two latent variables to be plotted if number of the latent variables is more than 2. A vector with length of two. Defaults to \code{c(1,2)}.
#' @param jitter   if \code{TRUE}, jittering is applied on points.
#' @param jitter.amount   numeric, positive value indicating an amount of jittering for each point, defaults to 0.2 (jitter range).
#' @param s.colors colors for sites
#' @param symbols logical, if \code{TRUE} sites are plotted using symbols, if \code{FALSE} (default) site numbers are used
#' @param cex.spp size of species labels in biplot
#' @param predict.region logical, if \code{TRUE} prediction regions for the predicted latent variables are plotted, defaults to \code{FALSE}.
#' @param level level for prediction regions.
#' @param lty.ellips line type for prediction ellipses. See graphical parameter lty.
#' @param lwd.ellips line width for prediction ellipses. See graphical parameter lwd.
#' @param col.ellips colors for prediction ellipses.
#' @param ...	additional graphical arguments.
#'
#' @details
#' Function constructs a scatter plot of two latent variables, i.e. an ordination plot. If only one latent
#' variable is in the fitted model, latent variables are plotted against their corresponding row indices.
#' The latent variables are labeled using the row index of the response matrix y.
#'
#' Coefficients related to latent variables are plotted in the same figure with the latent
#' variables if \code{biplot = TRUE}. They are labeled using the column names of y. The number
#' of latent variable coefficients to be plotted can be controlled by ind.spp. An argument alpha
#' is used to control the relative scaling of the latent variables and their coefficients.
#' If \code{alpha = 0.5}, the latent variables and their coefficients are on the same scale.
#' For details for constructing a biplot, see Gabriel (1971).
#' 
#' @note 
#' - If error is occurred when using \code{ordiplot()}, try full name of the function \code{ordiplot.gllvm()} as functions named 'ordiplot' might be found in other packages as well.
#' 
#' @references 
#' Gabriel, K. R. (1971). The biplot graphic display of matrices with application to principal component analysis. Biometrika, 58, 453-467.
#' 
#' @author Jenni Niku <jenni.m.e.niku@@jyu.fi>, Francis K.C. Hui
#'
#' @examples
#' #'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'fit <- gllvm(y, family = poisson())
#'# Ordination plot:
#'ordiplot(fit)
#'# Biplot with 10 species
#'ordiplot(fit, biplot = TRUE, ind.spp = 10)
#'
#'@aliases ordiplot ordiplot.gllvm
#'@export
#'@export ordiplot.gllvm
ordiplot.gllvm <- function(object, biplot = FALSE, ind.spp = NULL, alpha = 0.5, main = NULL, which.lvs = c(1, 2), predict.region = FALSE, level =0.95,
                           jitter = FALSE, jitter.amount = 0.2, s.colors = 1, symbols = FALSE, cex.spp = 0.7, lwd.ellips = 0.5, col.ellips = 4, lty.ellips = 1,...) {
  if (any(class(object) != "gllvm"))
    stop("Class of the object isn't 'gllvm'.")
  a <- jitter.amount
  n <- NROW(object$y)
  p <- NCOL(object$y)
  num.lv <- object$num.lv
  if (!is.null(ind.spp)) {
    ind.spp <- min(c(p, ind.spp))
  } else {
    ind.spp <- p
  }
  if (object$num.lv == 0)
    stop("No latent variables to plot.")
  
  if (is.null(rownames(object$params$theta)))
    rownames(object$params$theta) = paste("V", 1:p)
  
  if (object$num.lv == 1) {
    plot(1:n, object$lvs, ylab = "LV1", xlab = "Row index")
  }
  
  if (object$num.lv > 1) {
    do_svd <- svd(object$lvs)
    svd_rotmat_sites <- do_svd$v
    svd_rotmat_species <- do_svd$v
    
    choose.lvs <- object$lvs
    choose.lv.coefs <- object$params$theta
    bothnorms <- sqrt(colSums(choose.lvs^2)) * sqrt(colSums(choose.lv.coefs^2)) 
    ## Standardize both to unit norm then scale using bothnorms. Note alpha = 0.5 so both have same norm. Otherwise "significance" becomes scale dependent
    scaled_cw_sites <- t(t(choose.lvs) / sqrt(colSums(choose.lvs^2)) * (bothnorms^alpha)) 
    scaled_cw_species <- t(t(choose.lv.coefs) / sqrt(colSums(choose.lv.coefs^2)) * (bothnorms^(1-alpha))) 
    
    # equally: lvstr <- object$lvs%*%(diag((bothnorms^0.5)/sqrt(colSums(object$lvs^2)))%*%svd_rotmat_sites)
    choose.lvs <- scaled_cw_sites%*%svd_rotmat_sites
    choose.lv.coefs <- scaled_cw_species%*%svd_rotmat_species
    # equally: thettr <- object$params$theta%*%(diag((bothnorms^(1-0.5))/sqrt(colSums(object$params$theta^2)))%*%svd_rotmat_species)
    B<-(diag((bothnorms^alpha)/sqrt(colSums(object$lvs^2)))%*%svd_rotmat_sites)
    Bt<-(diag((bothnorms^(1-alpha))/sqrt(colSums(object$params$theta^2)))%*%svd_rotmat_species)
    
    
    # testcov <- object$lvs %*% t(object$params$theta)
    # do.svd <- svd(testcov, object$num.lv, object$num.lv)
    # choose.lvs <- do.svd$u * matrix( do.svd$d[1:object$num.lv] ^ alpha,
    #     nrow = n, ncol = object$num.lv, byrow = TRUE )
    # choose.lv.coefs <- do.svd$v * matrix(do.svd$d[1:object$num.lv] ^ (1 - alpha),
    #     nrow = p, ncol = object$num.lv, byrow = TRUE )
    
    
    if (!biplot) {
      plot(choose.lvs[, which.lvs],
           xlab = paste("Latent variable ", which.lvs[1]),
           ylab = paste("Latent variable ", which.lvs[2]),
           main = main , type = "n", ... )
      
      if (predict.region) {
        if(length(col.ellips)!=n){ col.ellips =rep(col.ellips,n)}
        
        if (object$method == "LA") {
          for (i in 1:n) {
            covm <- (t(B)%*%object$prediction.errors$lvs[i,,]%*%B)[which.lvs,which.lvs];
            #covm <- object$prediction.errors$lvs[i,which.lvs,which.lvs];
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=object$num.lv)), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
          }
        } else {
          sdb<-sdA(object)
          object$A<-sdb+object$A
          r=0
          if(object$row.eff=="random") r=1
          
          for (i in 1:n) {
            if(!object$TMB && object$Lambda.struc == "diagonal"){
              covm <- (t(B)%*%diag(object$A[i,1:num.lv+r])%*%B)[which.lvs,which.lvs];
              # covm <- diag(object$A[i,which.lvs+r]);
            } else {
              covm <- (t(B)%*%object$A[i,1:num.lv+r,1:num.lv+r]%*%B)[which.lvs,which.lvs];
              # covm <- object$A[i,which.lvs+r,which.lvs+r];
            }
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=object$num.lv)), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
          }
        }
      }
      
      if (!jitter)
        if (symbols) {
          points(choose.lvs[, which.lvs], col = s.colors, ...)
        } else {
          text(choose.lvs[, which.lvs], label = 1:n, cex = 1.2, col = s.colors)
        }
      if (jitter)
        if (symbols) {
          points(choose.lvs[, which.lvs][, 1] + runif(n,-a,a), choose.lvs[, which.lvs][, 2] + runif(n,-a,a), col =
                   s.colors, ...)
        } else {
          text(
            (choose.lvs[, which.lvs][, 1] + runif(n,-a,a)),
            (choose.lvs[, which.lvs][, 2] + runif(n,-a,a)),
            label = 1:n, cex = 1.2, col = s.colors )
        }
    }
    
    if (biplot) {
      largest.lnorms <- order(apply(object$params$theta ^ 2, 1, sum), decreasing = TRUE)[1:ind.spp]
      
      plot(
        rbind(choose.lvs[, which.lvs], choose.lv.coefs[, which.lvs]),
        xlab = paste("Latent variable ", which.lvs[1]),
        ylab = paste("Latent variable ", which.lvs[2]),
        main = main, type = "n", ... )
      
      if (predict.region) {
        if(length(col.ellips)!=n){ col.ellips =rep(col.ellips,n)}
        
        if (object$method == "LA") {
          for (i in 1:n) {
            #covm <- object$prediction.errors$lvs[i,which.lvs,which.lvs];
            covm <- (t(B)%*%object$prediction.errors$lvs[i,,]%*%B)[which.lvs,which.lvs];
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=object$num.lv)), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
          }
        } else {
          sdb<-sdA(object)
          object$A<-sdb+object$A
          r=0
          if(object$row.eff=="random") r=1
          
          for (i in 1:n) {
            if(!object$TMB && object$Lambda.struc == "diagonal"){
              covm <- (t(B)%*%diag(object$A[i,1:num.lv+r])%*%B)[which.lvs,which.lvs];
              #covm <- diag(object$A[i,which.lvs+r]);
            } else {
              covm <- (t(B)%*%object$A[i,1:num.lv+r,1:num.lv+r]%*%B)[which.lvs,which.lvs];
              # covm <- object$A[i,which.lvs+r,which.lvs+r];
            }
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=object$num.lv)), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
          }
        }
      }
      
      if (!jitter){
        if (symbols) {
          points(choose.lvs[, which.lvs], col = s.colors, ...)
        } else {
          text(choose.lvs[, which.lvs], label = 1:n, cex = 1.2, col = s.colors)
        }
        text(
          matrix(choose.lv.coefs[largest.lnorms, which.lvs], nrow = length(largest.lnorms)),
          label = rownames(object$params$theta)[largest.lnorms],
          col = 4, cex = cex.spp )
      }
      if (jitter){
        if (symbols) {
          points(choose.lvs[, which.lvs[1]] + runif(n,-a,a), (choose.lvs[, which.lvs[2]] + runif(n,-a,a)), col =
                   s.colors, ...)
        } else {
          text(
            (choose.lvs[, which.lvs[1]] + runif(n,-a,a)),
            (choose.lvs[, which.lvs[2]] + runif(n,-a,a)),
            label = 1:n, cex = 1.2, col = s.colors )
        }
        text(
          (matrix(choose.lv.coefs[largest.lnorms, which.lvs], nrow = length(largest.lnorms)) + runif(2*length(largest.lnorms),-a,a)),
          label = rownames(object$params$theta)[largest.lnorms],
          col = 4, cex = cex.spp )
      }
    }
    
    
    
  }
}


#'@export ordiplot
ordiplot <- function(object, ...)
{
  UseMethod(generic = "ordiplot")
}

