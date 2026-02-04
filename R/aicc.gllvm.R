#' @title Corrected Akaike information criterion and number of observations
#' @description Calculates corrected Akaike information criterion for small sample sizes, and extracts number of observations.
#' 
#'@name AICc
#'
#' @param object an object of class 'gllvm'.
#' @param ... Not used.
#' 
#' @author Jenni Niku, Bert van der Veen
#' 
#'@aliases AICc AICc.gllvm nobs nobs.gllvm
#'
#'@export
#'@export AICc.gllvm
#'@export nobs.gllvm
#'
AICc.gllvm <- function(object, ...){
     if (!missing(...)) {
      lls <- lapply(list(object, ...), logLik.gllvm)
      vals <- sapply(lls, function(el) {
        c(as.numeric(el), attr(el, "df"), attr(el, "nobs") %||% 
            NA_integer_)
      })
      val <- data.frame(df = vals[2L, ], ll = vals[1L, ])
      nos <- na.omit(vals[3L, ])
      if (length(nos) && any(nos != nos[1L])) 
        warning("models are not all fitted to the same number of observations")
    val <- data.frame(df = val$df, AICc = -2 * val$ll + val$df*2+2*val$df*(val$df+1)/(vals[3L,])-val$df-1)
      Call <- match.call()
      Call$k <- NULL
      row.names(val) <- as.character(Call[-1L])
      val
    }
    else {
      lls <- logLik.gllvm(object)
      df <- attr(lls, "df")
      nobs <- attr(lls, "nobs")
      -2 * as.numeric(lls) + 2 *  df+ 2*df*(df+1)/(nobs-df-1)
    }
  }

#'@method AICc gllvm
#'@export AICc
AICc <- function(object, ...)
{
  UseMethod(generic = "AICc")
}

#'@rdname AICc
#'@export
nobs.gllvm <- function(object, ...){
  n <- prod(dim(object$y))
  return(n)
}