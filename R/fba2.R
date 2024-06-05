#' @export
fba2 <- function(model) {

  #----------------------------------------------------------------------------#
  # Initializing and defining LP problem                                       #
  #----------------------------------------------------------------------------#
  LPprob <- new(paste0("LPproblem_",COBRAR_SETTINGS("SOLVER")),
                name = paste0("LP_", model@mod_id),
                method = COBRAR_SETTINGS("METHOD"))

  loadLPprob(LPprob,
             nCols = react_num(model),
             nRows = met_num(model)+constraint_num(model),
             mat   = rbind(model@S, model@constraints@coeff),
             ub    = model@uppbnd,
             lb    = model@lowbnd,
             obj   = model@obj_coef,
             rlb   = c(rep(0, met_num(model)),
                       model@constraints@lb),
             rtype = c(rep("E", met_num(model)),
                       model@constraints@rtype),
             lpdir = COBRAR_SETTINGS("OPT_DIRECTION"),
             rub   = c(rep(NA, met_num(model)),
                       model@constraints@ub),
             ctype = NULL
  )

  #----------------------------------------------------------------------------#
  # Optimizing problem                                                         #
  #----------------------------------------------------------------------------#
  lp_ok   <- solveLp(LPprob)
  # lp_stat <- getSolStat(LPprob)

  #----------------------------------------------------------------------------#
  # Retrieve predictions                                                       #
  #----------------------------------------------------------------------------#
  # objRes <- getObjValue(LPprob)
  # lp_fluxes <- getColsPrimal(LPprob)
  #
  # redCosts <- getRedCosts(LPprob)
  #

  return(lp_ok)
}

