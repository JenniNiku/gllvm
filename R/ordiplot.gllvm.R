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
#' @param spp.colors colors for sites, defaults to \code{"blue"}
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
#' Effects of environmental variables in constrained ordination are indicated with arrows. 
#' If standard errors are available in the provided model, the slopes of environmental variables
#' for which the 95% confidence intervals do not include zero are shown as red, while others 
#' are slightly less intensely colored.
#' 
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
#' #'# Extract subset of the microbial data to be used as an example
#'data(microbialdata)
#'y <- microbialdata$Y[, order(colMeans(microbialdata$Y > 0), 
#'                      decreasing = TRUE)[21:40]]
#'fit <- gllvm(y, family = poisson())
#'fit$logL
#'ordiplot(fit, predict.region = TRUE)
#' \dontrun{
#' #'## Load a dataset from the mvabund package
#'data(antTraits)
#'y <- as.matrix(antTraits$abund)
#'fit <- gllvm(y, family = poisson())
#'# Ordination plot:
#'ordiplot(fit)
#'# Biplot with 10 species
#'ordiplot(fit, biplot = TRUE, ind.spp = 10)
#'}
#'@aliases ordiplot ordiplot.gllvm
#'@export
#'@export ordiplot.gllvm
ordiplot.gllvm <- function(object, biplot = FALSE, ind.spp = NULL, alpha = 0.5, main = NULL, which.lvs = c(1, 2), predict.region = FALSE, level =0.95,
                           jitter = FALSE, jitter.amount = 0.2, s.colors = 1, symbols = FALSE, cex.spp = 0.7, spp.colors = "blue", lwd.ellips = 0.5, col.ellips = 4, lty.ellips = 1,...) {
  if (!any(class(object) %in% "gllvm"))
    stop("Class of the object isn't 'gllvm'.")
  
  if(any(class(object)=="gllvm.quadratic")){
    warning("A biplot does not accurately visualize a GLLVM with quadratic response model. \n")
  }
  a <- jitter.amount
  n <- NROW(object$y)
  p <- NCOL(object$y)
  num.lv <- object$num.lv
  num.lv.c <- object$num.lv.c 
  quadratic <- object$quadratic
  if (!is.null(ind.spp)) {
    ind.spp <- min(c(p, ind.spp))
  } else {
    ind.spp <- p
  }
  if(length(spp.colors)==1){
    spp.colors <- rep(spp.colors,p)
  }else if(length(spp.colors)!=p){
    stop("spp.colors needs to be of length p or 1.")
  }
  if ((num.lv+num.lv.c) == 0)
    stop("No latent variables to plot.")
  
  if (is.null(rownames(object$params$theta)))
    rownames(object$params$theta) = paste("V", 1:p)
  
  if ((num.lv+num.lv.c) == 1) {
    if(num.lv==1)plot(1:n, getLV(object), ylab = "LV1", xlab = "Row index")
    if(num.lv.c==1)plot(1:n, getLV(object), ylab = "CLV1", xlab = "Row index")
  }
  
  if ((num.lv+num.lv.c) > 1) {
    do_svd <- svd(getLV(object))
    svd_rotmat_sites <- do_svd$v
    svd_rotmat_species <- do_svd$v
    
    choose.lvs <- getLV(object)
    if(quadratic == FALSE)choose.lv.coefs <- object$params$theta
    if(quadratic != FALSE){
      testcov <- getLV(object) %*% t(object$params$theta[, 1:(num.lv+num.lv.c)]) + 
        getLV(object)^2 %*% t(object$params$theta[, -c(1:(num.lv+num.lv.c))])
      do_svd <- svd(testcov, (num.lv+num.lv.c), (num.lv+num.lv.c))
      choose.lvs <- do_svd$u
      choose.lv.coefs <- do_svd$v
    }
    
    bothnorms <- sqrt(colSums(choose.lvs^2)) * sqrt(colSums(choose.lv.coefs^2)) 
    ## Standardize both to unit norm then scale using bothnorms. Note alpha = 0.5 so both have same norm. Otherwise "significance" becomes scale dependent
    scaled_cw_sites <- t(t(choose.lvs) / sqrt(colSums(choose.lvs^2)) * (bothnorms^alpha)) 
    scaled_cw_species <- t(t(choose.lv.coefs) / sqrt(colSums(choose.lv.coefs^2)) * (bothnorms^(1-alpha))) 
    
    # equally: lvstr <- object$lvs%*%(diag((bothnorms^0.5)/sqrt(colSums(object$lvs^2)))%*%svd_rotmat_sites)
    choose.lvs <- scaled_cw_sites%*%svd_rotmat_sites
    choose.lv.coefs <- scaled_cw_species%*%svd_rotmat_species
    # equally: thettr <- object$params$theta%*%(diag((bothnorms^(1-0.5))/sqrt(colSums(object$params$theta^2)))%*%svd_rotmat_species)
    
    B<-(diag((bothnorms^alpha)/sqrt(colSums(getLV(object)^2)))%*%svd_rotmat_sites)
    if(quadratic==FALSE)Bt<-(diag((bothnorms^(1-alpha))/sqrt(colSums(object$params$theta^2)))%*%svd_rotmat_species)
    
    # testcov <- object$lvs %*% t(object$params$theta)
    # do.svd <- svd(testcov, num.lv, num.lv)
    # choose.lvs <- do.svd$u * matrix( do.svd$d[1:num.lv] ^ alpha,
    #     nrow = n, ncol = num.lv, byrow = TRUE )
    # choose.lv.coefs <- do.svd$v * matrix(do.svd$d[1:num.lv] ^ (1 - alpha),
    #     nrow = p, ncol = num.lv, byrow = TRUE )
    
    
    if (!biplot) {
      plot(choose.lvs[, which.lvs],
           xlab = ifelse(which.lvs[1]<=num.lv.c,paste("Constrained latent variable",(1:num.lv.c)[which.lvs[1]]),paste("Latent variable",(1:num.lv)[which.lvs[1]])), 
           ylab = ifelse(which.lvs[1]<=num.lv.c,paste("Constrained latent variable",(1:num.lv.c)[which.lvs[2]]),paste("Latent variable",(1:num.lv)[which.lvs[2]])),
           main = main , type = "n", ... )
      
      if (predict.region) {
        if(length(col.ellips)!=n){ col.ellips =rep(col.ellips,n)}
        
        if (object$method == "LA") {
          for (i in 1:n) {
            covm <- (t(B)%*%object$prediction.errors$lvs[i,,]%*%B)[which.lvs,which.lvs];
            #covm <- object$prediction.errors$lvs[i,which.lvs,which.lvs];
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=(num.lv+num.lv.c))), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
          }
        } else {
          sdb<-sdA(object)
          object$A<-sdb+object$A
          r=0
          if(object$row.eff=="random") r=1
          
          for (i in 1:n) {
            if(!object$TMB && object$Lambda.struc == "diagonal"){
              covm <- (t(B)%*%diag(object$A[i,1:(num.lv+num.lv.c)+r])%*%B)[which.lvs,which.lvs];
              # covm <- diag(object$A[i,which.lvs+r]);
            } else {
              covm <- (t(B)%*%object$A[i,1:(num.lv+num.lv.c)+r,1:(num.lv+num.lv.c)+r]%*%B)[which.lvs,which.lvs];
              # covm <- object$A[i,which.lvs+r,which.lvs+r];
            }
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=(num.lv+num.lv.c))), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
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
        xlab = ifelse(which.lvs[1]<=num.lv.c,paste("Constrained latent variable",(1:num.lv.c)[which.lvs[1]]),paste("Latent variable",(1:num.lv)[which.lvs[1]])),
        ylab = ifelse(which.lvs[1]<=num.lv.c,paste("Constrained latent variable",(1:num.lv.c)[which.lvs[2]]),paste("Latent variable",(1:num.lv)[which.lvs[2]])),
        main = main, type = "n", ... )
      
      if (predict.region) {
        if(length(col.ellips)!=n){ col.ellips =rep(col.ellips,n)}
        
        if (object$method == "LA") {
          for (i in 1:n) {
            #covm <- object$prediction.errors$lvs[i,which.lvs,which.lvs];
            covm <- (t(B)%*%object$prediction.errors$lvs[i,,]%*%B)[which.lvs,which.lvs];
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=(num.lv+num.lv.c))), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
          }
        } else {
          sdb<-sdA(object)
          object$A<-sdb+object$A
          r=0
          if(object$row.eff=="random") r=1
          
          for (i in 1:n) {
            if(!object$TMB && object$Lambda.struc == "diagonal"){
              covm <- (t(B)%*%diag(object$A[i,1:(num.lv+num.lv.c)+r])%*%B)[which.lvs,which.lvs];
              #covm <- diag(object$A[i,which.lvs+r]);
            } else {
              covm <- (t(B)%*%object$A[i,1:(num.lv+num.lv.c)+r,1:(num.lv+num.lv.c)+r]%*%B)[which.lvs,which.lvs];
              # covm <- object$A[i,which.lvs+r,which.lvs+r];
            }
            ellipse( choose.lvs[i, which.lvs], covM = covm, rad = sqrt(qchisq(level, df=(num.lv+num.lv.c))), col = col.ellips[i], lwd = lwd.ellips, lty = lty.ellips)
          }
        }
      }
      
      if (!jitter){
        if (symbols) {
          points(choose.lvs[, which.lvs], col = s.colors, ...)
        } else {
          text(choose.lvs[, which.lvs], label = 1:n, cex = 1.2, col = s.colors)
        }
        spp.colors <- spp.colors[largest.lnorms][1:ind.spp]
        text(
          matrix(choose.lv.coefs[largest.lnorms, which.lvs], nrow = length(largest.lnorms)),
          label = rownames(object$params$theta)[largest.lnorms],
          col = spp.colors, cex = cex.spp )
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
        spp.colors <- spp.colors[largest.lnorms][1:ind.spp]
        text(
          (matrix(choose.lv.coefs[largest.lnorms, which.lvs], nrow = length(largest.lnorms)) + runif(2*length(largest.lnorms),-a,a)),
          label = rownames(object$params$theta)[largest.lnorms],
          col = spp.colors, cex = cex.spp )
      }
    }

    if(num.lv.c>1&all(which.lvs<=num.lv.c)){
            #LvXcoef <- LvXcoef/ sqrt(colSums(object$lvs[,which.lvs]^2)) * (bothnorms^alpha)
      # LVcor <- t(cor(choose.lvs+object$lv.X%*%t(svd_rotmat_sites%*%t(object$params$LvXcoef[,which.lvs])),spider$x))
      # LVcor<-t(t(LVcor)/ (bothnorms^alpha) *sqrt(colSums(object$lvs[,which.lvs]^2)))
      LVcoef <- object$params$LvXcoef
      if(any(row.names(LVcoef)%in%"(Intercept)")){
        intercept <- LVcoef[row.names(LVcoef)%in%"(Intercept)",]
        LVcoef <- LVcoef[!row.names(LVcoef)%in%"(Intercept)",]
      }
      if(!is.logical(object$sd)){
        cilow <- object$params$LvXcoef+qnorm( (1 - 0.95) / 2)*object$sd$LvXcoef
        ciup <- object$params$LvXcoef+qnorm(1- (1 - 0.95) / 2)*object$sd$LvXcoef
        lty <- rep("solid",ncol(object$lv.X))
        col <- rep("red", ncol(object$lv.X))
        lty[sign(cilow[,1])!=sign(ciup[,1])|sign(cilow[,2])!=sign(ciup[,2])] <- "solid"
        col[sign(cilow[,1])!=sign(ciup[,1])|sign(cilow[,2])!=sign(ciup[,2])] <- hcl(0, 100, 80)#rgb(1,0,0,alpha=0.3)
        
      }else{
        lty <- rep("dashed",ncol(object$lv.X))
        col<-rep("red",ncol(object$lv.X))
      }
        LVcoef <- LVcoef%*%svd_rotmat_sites
        LVcoef <- LVcoef/apply(object$lv.X[,!colnames(object$lv.X)%in%"(Intercept)",drop=F],2,sd)
        marg<-par("usr")
    
        Xlength<-min(dist(c(mean(marg[1:2]),marg[1])),dist(c(mean(marg[1:2]),marg[2])))
        Ylength<-min(dist(c(mean(marg[3:4]),marg[3])),dist(c(mean(marg[3:4]),marg[4])))
        origin<- c(mean(marg[1:2]),mean(marg[3:4]))

        #scale the largest arrow to 80% of the smallest distance from 0 to the edge of the plot
        LVcoef <- t(t(LVcoef)/apply(abs(LVcoef),2,max))*min(Xlength,Ylength)*0.8
        
        for(i in 1:nrow(LVcoef)){
          arrows(x0=origin[1],y0=origin[2],x1=origin[1]+LVcoef[i,1],y1=origin[2]+LVcoef[i,2],col=col[i],length=0.1,lty=lty[i])  
          text(x=origin[1]+LVcoef[i,1],y=origin[2]+LVcoef[i,2],labels = row.names(LVcoef)[i],col=col[i])
        }
      }

    }
}


#'@export ordiplot
ordiplot <- function(object, ...)
{
  UseMethod(generic = "ordiplot")
}

