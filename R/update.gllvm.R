#' @title Update and Re-fit a gllvm Model Call
#' @description Models of class gllvm have three formula interfaces: formula, lv.formula, and row.eff. Each can include fixed and random effects.
#'
#'
#' @param object an object of class 'gllvm'.
#' @param formula interface for column-specific effects.
#' @param lv.formula interface for latent variable (reduced-rank) effects.
#' @param row.eff interface for row-specific (i.e., same for all columns) effects.
#' @param eval default to TRUE, in which case the reformulated model is fitted directly. if FALSE returns the new model call.
#' @param ...	not used.
#' @details
#' Updates the formula of a gllvm object.
#'
#' @author Bert van der Veen
#' 
#' @aliases update update.gllvm
#' @method update gllvm
#' @importFrom stats update
#'
#' @export update.gllvm
update.gllvm <- function(object, formula = NULL, lv.formula = NULL, row.eff = NULL, eval = TRUE, ...){
  call <- getCall(object)
  match.call(expand.dots = FALSE, call = call)
  
  args <- list(...)
  
  if(!is.null(formula)){
    formula.new <- formula
    formula_ <- call$formula
    if(!is.null(formula_)){
    call$formula <- update.formula(formula_, formula.new)
    }else{
      call$formula <- formula.new
    }
  }
  if(!is.null(lv.formula)){
    lv.formula.new <- args$lv.formula
    lv.formula_ <- call$lv.formula
    
    if(!is.null(lv.formula_)){
      call$lv.formula <- update.formula(lv.formula_, lv.formula.new)
    }else{
      call$lv.formula <- lv.formula.new
    }
  }
  if(!is.null(row.eff) && is.language(row.eff)){
    row.eff.new <- args$row.eff
    row.eff_ <- call$row.eff
    
    if(!is.null(row.eff_)){
      call$row.eff <- update.formula(row.eff_, row.eff.new)
    }else{
      call$row.eff <- row.eff.new
    }
  }else if(!is.null(row.eff)){
    call[["row.eff"]] <- row.eff
  }
  
  # Add or override other arguments from ...
  for (arg_name in names(args)) {
    call[[arg_name]] <- args[[arg_name]]
  }
  
  if(eval){
    eval(call, parent.frame())
  }else{
    call
  }
}
