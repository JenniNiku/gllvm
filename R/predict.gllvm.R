predict.gllvm <- function (object, newX = NULL, newTR = NULL, newLV = NULL, type = "link", 
                           level = 1, ...) 
{
  r0 <- NULL
  newdata <- newX
  p <- ncol(object$y)
  n <- max(nrow(object$y), nrow(newdata), nrow(newLV))
  if (!is.null(newdata)) 
    n <- nrow(newdata)
  if (is.null(newdata) && !is.null(object$X) && !is.null(newLV)) {
    if (nrow(newLV) != nrow(object$y)) {
      stop("Number of rows in newLV must equal to the number of rows in the response matrix, if environmental variables are included in the model and newX is not included.")
    }
  }
  if (is.null(newdata) & !is.null(newLV)) {
    n <- nrow(newLV)
  }
  if (!is.null(object$X)) {
    formula <- formula(terms(object))
  }
  else {
    formula <- NULL
  }
  if (object$row.eff != FALSE) {
    if (length(object$params$row.params) != nrow(object$y)) 
      object$params$row.params = c(object$TMBfn$env$data$dr0 %*% 
                                     object$params$row.params)
  }
  b0 <- object$params$beta0
  eta <- matrix(b0, n, p, byrow = TRUE)
  if (!is.null(newTR)) 
    if (nrow(newTR) != p) 
      stop("Number of rows in newTR must match to the number of responses in the original data matrix.")
  if (is.null(colnames(object$y))) {
    colnames(object$y) <- paste("y", 1:p, sep = "")
  }
  if (!is.null(object$X) && is.null(object$TR)) {
    B <- object$params$Xcoef
    if (is.null(newdata)) {
      X.d <- Xnew <- object$X
    }
    else {
      X.d <- Xnew <- newdata
    }
    if (!is.null(newdata)) {
      if (is.null(object$call$formula)) {
        n1 <- colnames(newdata)
        formula1 <- paste("~ ", n1[1], sep = "")
        if (length(n1) > 1) {
          for (i1 in 2:length(n1)) {
            formula1 <- paste(formula1, n1[i1], sep = "+")
          }
        }
        formula1 <- paste(formula1, "1", sep = "-")
        formula1 <- formula(formula1)
        Xnew <- as.matrix(model.matrix(formula1, data = data.frame(newdata)))
      }
      formula <- formula(object$formula)
      xb <- as.matrix(model.matrix(formula, data = data.frame(Xnew)))
      X.d <- as.matrix(xb[, !(colnames(xb) %in% c("(Intercept)")), 
                          drop = F])
      colnames(X.d) <- colnames(xb)[!(colnames(xb) %in% 
                                        c("(Intercept)"))]
    }
    else {
      X.d <- object$X.design
    }
    eta <- eta + X.d %*% t(B)
  }
  if (!is.null(object$X) && !is.null(object$TR)) {
    B <- object$params$B
    X.d <- object$X.design
    if (!is.null(newdata) && !is.null(newTR)) {
      Xnew <- newdata
      TRnew <- newTR
      y1 <- object$y[sample(1:nrow(object$y), nrow(Xnew), 
                            replace = TRUE), ]
      yX <- reshape(data.frame(cbind(y1, Xnew)), direction = "long", 
                    varying = colnames(data.frame(y1)), v.names = "y", 
                    timevar = "species")
      TR2 <- data.frame(species = 1:p, TRnew)
      yXT <- merge(yX, TR2, by = "species")
      X.d <- as.matrix(model.matrix(formula, data = yXT))
    }
    if (!is.null(newdata) && is.null(newTR)) {
      formula <- formula(object$formula)
      n1 <- colnames(newdata)
      formula1 <- paste("~", n1[1], sep = "")
      if (length(n1) > 1) {
        for (i1 in 2:length(n1)) {
          formula1 <- paste(formula1, n1[i1], sep = "+")
        }
      }
      formula1 <- paste(formula1, "1", sep = "-")
      formula1 <- formula(formula1)
      Xnew <- as.matrix(model.matrix(formula1, data = data.frame(newdata)))
      TRnew <- object$TR
      y1 <- object$y[sample(1:nrow(object$y), nrow(Xnew), 
                            replace = TRUE), ]
      yX <- reshape(data.frame(cbind(y1, Xnew)), direction = "long", 
                    varying = colnames(data.frame(y1)), v.names = "y", 
                    timevar = "species")
      TR2 <- data.frame(species = 1:p, TRnew)
      yXT <- merge(yX, TR2, by = "species")
      X.d <- as.matrix(model.matrix(formula, data = yXT))
    }
    if (is.null(newdata) && !is.null(newTR)) {
      if (nrow(newTR) != p) 
        stop("Number of rows in TRnew must equal to the number of response variables.")
      formula <- formula(object$formula)
      TRnew <- newTR
      n1 <- colnames(newTR)
      formula1 <- paste("~", n1[1], sep = "")
      if (length(n1) > 1) {
        for (i1 in 2:length(n1)) {
          formula1 <- paste(formula1, n1[i1], sep = "+")
        }
      }
      formula1 <- paste(formula1, "1", sep = "-")
      formula1 <- formula(formula1)
      TRnew <- as.matrix(model.matrix(formula1, data = data.frame(newTR)))
      Xnew <- object$X
      y1 <- object$y[sample(1:nrow(object$y), nrow(Xnew), 
                            replace = TRUE), ]
      yX <- reshape(data.frame(cbind(y1, Xnew)), direction = "long", 
                    varying = colnames(data.frame(y1)), v.names = "y", 
                    timevar = "species")
      TR2 <- data.frame(species = 1:p, TRnew)
      yXT <- merge(yX, TR2, by = "species")
      X.d <- as.matrix(model.matrix(formula, data = yXT))
    }
    X.d <- as.matrix(X.d[, colnames(X.d) != "(Intercept)"])
    eta <- eta + matrix(X.d %*% B, n, p)
    if (!is.null(object$randomX)) {
      if (!is.null(newdata)) {
        tr <- try(X.xr <- as.matrix(model.matrix(object$randomX, 
                                                 data = yXT)), silent = TRUE)
        if (inherits(tr, "try-error")) {
          X.xr <- as.matrix(yXT[, colnames(object$Xrandom)])
          colnames(X.xr) <- colnames(object$Xrandom)
        }
        rnam <- colnames(X.xr)[!(colnames(X.xr) %in% 
                                   c("(Intercept)"))]
        xr <- as.matrix(X.xr[1:n, !(colnames(X.xr) %in% 
                                      c("(Intercept)"))])
        if (NCOL(xr) == 1) 
          colnames(xr) <- rnam
      }
      else {
        xr <- object$Xrandom
      }
      eta <- eta + matrix(xr %*% object$params$Br, n, p)
    }
  }
  if (level == 1) {
    if (is.null(newLV) && !is.null(newdata) & (object$num.lv + 
                                               object$num.lv.c) > 0) {
      if (nrow(newdata) != nrow(object$y)) {
        stop("Level 1 predictions cannot be calculated for new X values if new latent variable values are not given. Change to 'level = 0' predictions.")
      }
    }
    if (object$num.lv > 0 | (object$num.lv.c + object$num.RR) > 
        0) {
      if (!is.null(newLV)) {
        if (ncol(newLV) != (object$num.lv + object$num.lv.c)) 
          stop("Number of latent variables in input doesn't equal to the number of latent variables in the model.")
        if (!is.null(newdata)) {
          if (nrow(newLV) != nrow(newdata)) 
            stop("Number of rows in newLV must equal to the number of rows in newX, if newX is included, otherwise same as number of rows in the response matrix.")
        }
        lvs <- t(t(newLV) * object$params$sigma.lv)
      } else {
        if (object$num.RR == 0) {
          lvs <- t(t(object$lvs) * object$params$sigma.lv)
        }
        else {
          if (object$num.lv.c > 0) {
            lvs <- cbind(t(t(object$lvs[, 1:object$num.lv.c]) * 
                             object$params$sigma.lv[1:object$num.lv.c]), 
                         matrix(0, ncol = object$num.RR, nrow = n), 
                         t(t(object$lvs[, -c(1:object$num.lv.c)]) * 
                             object$params$sigma.lv[1:object$num.lv]))
          }
          else if (object$num.lv > 0 & object$num.lv.c == 
                   0) {
            lvs <- cbind(matrix(0, ncol = object$num.RR, 
                                nrow = n), t(t(object$lvs) * object$params$sigma.lv))
          }
          else {
            lvs <- matrix(0, ncol = object$num.RR, nrow = n)
          }
        }
      }
      if ((object$num.lv.c + object$num.RR) > 0 & !is.null(newdata)) {
        lv.X <- model.matrix(object$lv.formula, as.data.frame(newdata))[, 
                                                                        -1, drop = F]
      } else {
        lv.X <- object$lv.X
      }
      theta <- object$params$theta[, 1:(object$num.lv + 
                                          (object$num.lv.c + object$num.RR))]
      eta <- eta + lvs %*% t(theta)
      if ((object$num.lv.c + object$num.RR) > 0) {
        eta <- eta + lv.X %*% object$params$LvXcoef %*% 
          t(theta[, 1:(object$num.lv.c + object$num.RR)])
      }
      if (object$quadratic != FALSE) {
        if(object$num.lv>0){
          theta2 <- object$params$theta[, -c(1:(object$num.lv.c + object$num.RR)), drop = F]
          eta <- eta + lvs[(ncol(lvs)-object$num.lv+1):ncol(lvs)]^2 %*% t(theta2)
        }
        if ((object$num.lv.c + object$num.RR) > 0) {
          theta2C <- abs(theta2[, 1:(object$num.lv.c + 
                                       object$num.RR), drop = F])
          lvs <- lvs[,1:(object$num.lv.c+object$num.RR)] + lv.X%*%object$params$LvXcoef
          lvsT <- t(lvs)
          for (j in 1:p) {
            eta[i, j] <- eta[, j] - lvs%*%diag(theta2C[j, ])%*%lvsT
          }
        }
      }
    }
  }
  if (object$row.eff %in% c("random", "fixed", 
                            "TRUE") && nrow(eta) == length(object$params$row.params) & 
      is.null(r0)) {
    r0 <- object$params$row.params
    eta <- eta + r0
  }
  else if (!is.null(r0)) {
    eta <- eta + r0
  }
  if (object$family %in% c("poisson", "negative.binomial", 
                           "tweedie", "gamma", "exponential")) 
    ilinkfun <- exp
  if (object$family == "binomial" || object$family == 
      "beta") 
    ilinkfun <- binomial(link = object$link)$linkinv
  if (object$family == "ordinal") 
    ilinkfun <- pnorm
  if (object$family == "ZIP") 
    ilinkfun <- function(eta) exp(eta) * (1 - matrix(object$params$phi, 
                                                     n, p, byrow = TRUE))
  if (object$family == "gaussian") 
    ilinkfun <- gaussian()$linkinv
  out <- NULL
  preds <- NULL
  if ("link" %in% type) 
    out <- eta
  if ("response" %in% type) 
    out <- ilinkfun(eta)
  if (is.null(newdata) && is.null(newTR) && is.null(newLV) && 
      "logL" %in% type) 
    out <- object$logL
  if (object$family == "ordinal" && type == "response") {
    if (object$zeta.struc == "species") {
      k.max <- apply(object$params$zeta, 1, function(x) length(x[!is.na(x)])) + 
        1
      preds <- array(NA, dim = c(max(k.max), nrow(eta), 
                                 p), dimnames = list(paste("level", 1:max(k.max), 
                                                           sep = ""), NULL, NULL))
      for (i in 1:n) {
        for (j in 1:p) {
          probK <- NULL
          probK[1] <- pnorm(object$params$zeta[j, 1] - 
                              eta[i, j], log.p = FALSE)
          probK[k.max[j]] <- 1 - pnorm(object$params$zeta[j, 
                                                          k.max[j] - 1] - eta[i, j])
          if (k.max[j] > 2) {
            j.levels <- 2:(k.max[j] - 1)
            for (k in j.levels) {
              probK[k] <- pnorm(object$params$zeta[j, 
                                                   k] - eta[i, j]) - pnorm(object$params$zeta[j, 
                                                                                              k - 1] - eta[i, j])
            }
          }
          preds[, i, j] <- c(probK, rep(NA, max(k.max) - 
                                          k.max[j]))
        }
      }
      out <- preds
    }
    else {
      k.max <- length(object$params$zeta) + 1
      preds <- array(NA, dim = c(k.max, nrow(eta), p), 
                     dimnames = list(paste("level", 1:max(k.max), 
                                           sep = ""), NULL, NULL))
      for (i in 1:n) {
        for (j in 1:p) {
          probK <- NULL
          probK[1] <- pnorm(object$params$zeta[1] - eta[i, 
                                                        j], log.p = FALSE)
          probK[k.max] <- 1 - pnorm(object$params$zeta[k.max - 
                                                         1] - eta[i, j])
          levels <- 2:(k.max - 1)
          for (k in levels) {
            probK[k] <- pnorm(object$params$zeta[k] - 
                                eta[i, j]) - pnorm(object$params$zeta[k - 
                                                                        1] - eta[i, j])
          }
          preds[, i, j] <- c(probK)
        }
      }
      out <- preds
    }
    dimnames(preds)[[3]] <- colnames(object$y)
  }
  try(rownames(out) <- 1:NROW(out), silent = TRUE)
  return(out)
}
