# some GLPK specific codes:
cplexPar <- list(
  CPX_MIN = 1,
  CPX_MAX = -1
)

#' @importClassesFrom cobrar LPproblem
#'
#' @exportClass LPproblem_cplex
setClass(Class = "LPproblem_cplex", slots = c(ptr.env = "externalptr",
                                              ptr.mod = "externalptr",
                                              ptr.x   = "externalptr",
                                              ptr.c   = "externalptr",
                                              ptr.obj = "externalptr"),
         contains = "LPproblem"
)

setMethod(f = "initialize",
          signature = "LPproblem_cplex",
          definition = function(.Object,
                                name,
                                method) {

            ptrlist <- initProb(name, COBRAR_SETTINGS("TOLERANCE"))

            .Object@ptr     <- ptrlist$cpx
            .Object@ptr.env <- ptrlist$env
            .Object@ptr.mod <- ptrlist$mod
            .Object@ptr.x   <- ptrlist$x
            .Object@ptr.c   <- ptrlist$c
            .Object@ptr.obj <- ptrlist$obj

            .Object@solver = "cplex"
            .Object@method = method

            return(.Object)
          }
)

#' @export
setMethod("loadLPprob", signature(lp = "LPproblem_cplex"),

          function(lp, nCols, nRows, mat, ub, lb, obj, rlb, rtype, lpdir,
                   rub = NULL, ctype = NULL) {

            # # problem dimensions
            addCols(lp, ncols = nCols)
            addRows(lp, nrows = nRows)

            # optimization direction
            lpdir <- switch(EXPR = lpdir,
                            "max" = cplexPar$CPX_MAX,
                            "min" = cplexPar$CPX_MIN)
            setObjDirection(lp, lpdir = lpdir)


            # populate constraint matrix
            TMPmat <- as(mat, "TsparseMatrix")
            loadMatrix(lp,
                       ne = length(TMPmat@x),
                       ia = TMPmat@i,
                       ja = TMPmat@j,
                       ra = TMPmat@x)

            # column (variable) bounds and objective function
            setColsBndsObjCoefs(lp,
                                j = c(1:nCols)-1,
                                lb = lb,
                                ub = ub,
                                obj_coef = obj)


            # # variable type
            # if (!is.null(ctype)) {
            #   cctype <- sapply(ctype,
            #                    function(x) switch(EXPR = x,
            #                                       "C" = glpkPar$GLP_CV,
            #                                       "I" = glpkPar$GLP_IV,
            #                                       "B" = glpkPar$GLP_BV,
            #                                       glpkPar$GLP_CV))
            #
            #   setColsKind(lp, j = c(1:nCols), kind = cctype)
            # }
            #
            # right hand side
            # Note: This is the steady state condition: Production of internal
            # metabolites should equal the consumption. In over words:
            # row lower bound = row upper bounds = 0
            if (is.null(rub)) {
              crub <- rlb
            }
            else {
              crub <- rub
            }
            stopifnot(length(rlb) == length(crub))
            setRowsBnds(lp,
                        i = c(1:nRows),
                        lb = rlb,
                        ub = crub,
                        type = rtype)


          }
)

setMethod("setObjDirection", signature(lp = "LPproblem_cplex"),
          function(lp, lpdir) {
            setObjDirLP(lp@ptr.obj, lpdir)
          }
)

setGeneric("addCols", function(lp, ...) {
  standardGeneric("addCols")
})
setMethod("addCols", signature(lp = "LPproblem_cplex"),
          function(lp, ncols) {
            addColsLP(lp@ptr.env, lp@ptr.x, as.integer(ncols))
          }
)


setMethod("addRows", signature(lp = "LPproblem_cplex"),
          function(lp, nrows) {
            addRowsLP(lp@ptr.env, lp@ptr.c, as.integer(nrows))
          }
)


setMethod("loadMatrix", signature(lp = "LPproblem_cplex"),
          function(lp, ne, ia, ja, ra) {
            loadMatrixLP(lp@ptr.x, lp@ptr.c,
                         as.integer(ne),
                         as.integer(ia),
                         as.integer(ja),
                         as.numeric(ra))
          }
)


setMethod("setColsBndsObjCoefs", signature(lp = "LPproblem_cplex"),
          function(lp, j, lb, ub, obj_coef) {

            setColsBndsObjCoefsLP(lp@ptr.obj, lp@ptr.x,
                                  as.integer(j),
                                  as.numeric(lb),
                                  as.numeric(ub),
                                  as.numeric(obj_coef))
          }
)


# setMethod("setColsKind", signature(lp = "LPproblem_cplex"),
#           function(lp, j, kind) {
#             setColsKindLP(lp@ptr,
#                           as.integer(j),
#                           as.integer(kind))
#           }
# )


setMethod("setRowsBnds", signature(lp = "LPproblem_cplex"),
          function(lp, i, lb, ub , type) {


            # type <- sapply(type,
            #                function(x) switch(EXPR = x,
            #                                   "F" = glpkPar$GLP_FR,
            #                                   "L" = glpkPar$GLP_LO,
            #                                   "U" = glpkPar$GLP_UP,
            #                                   "D" = glpkPar$GLP_DB,
            #                                   "E" = glpkPar$GLP_FX,
            #                                   glpkPar$GLP_FX))
            #
            # if (is.null(type)) {
            #   Ctype <- as.null(type)
            # }
            # else {
            #   Ctype <- as.integer(type)
            # }
            indE <- which(type == "E")
            ub[indE] <- lb[indE]

            indU <- which(type == "U")
            lb[indU] <- -Inf

            indL <- which(type == "L")
            ub[indL] <- Inf

            indF <- which(type == "F")
            lb[indF] <- -Inf
            ub[indF] <- Inf

            setRowsBndsLP(lp@ptr.c,
                          as.integer(i-1),
                          #Ctype,
                          as.numeric(lb),
                          as.numeric(ub))

          }
)

setMethod("solveLp", signature(lp = "LPproblem_cplex"),
          function(lp) {
            out <- solveCPLEX(lp@ptr, lp@ptr.mod, lp@ptr.x, lp@ptr.c, lp@ptr.obj)

            term <- switch(EXPR = out+1,
                           "Unknown",
                           "Optimal",
                           "Unbounded",
                           "Infeasible",
                           "InfOrUnbd",
                           "OptimalInfeas",
                           "NumBest",
                           "FeasibleRelaxedSum",
                           "OptimalRelaxedSum",
                           "FeasibleRelaxedInf",
                           "OptimalRelaxedInf",
                           "FeasibleRelaxedQuad",
                           "OptimalRelaxedQuad",
                           "AbortRelaxed",
                           "AbortObjLim",
                           "AbortPrimObjLim",
                           "AbortDualObjLim",
                           "AbortItLim",
                           "AbortTimeLim",
                           "AbortDetTimeLim",
                           "AbortUser",
                           "OptimalFaceUnbounded",
                           "OptimalTol",
                           "SolLim",
                           "PopulateSolLim",
                           "NodeLimFeas",
                           "NodeLimInfeas",
                           "FailFeas",
                           "FailInfeas",
                           "MemLimFeas",
                           "MemLimInfeas",
                           "FailFeasNoTree",
                           "FailInfeasNoTree",
                           "ConflictFeasible",
                           "ConflictMinimal",
                           "ConflictAbortContradiction",
                           "ConflictAbortTimeLim",
                           "ConflictAbortDetTimeLim",
                           "ConflictAbortItLim",
                           "ConflictAbortNodeLim",
                           "ConflictAbortObjLim",
                           "ConflictAbortMemLim",
                           "ConflictAbortUser",
                           "Feasible",
                           "OptimalPopulated",
                           "OptimalPopulatedTol",
                           "RelaxationUnbounded",
                           "FirstOrder",
                           "MultiObjOptimal",
                           "MultiObjNonOptimal",
                           "MultiObjInfeasible",
                           "MultiObjUnbounded",
                           "MultiObjInfOrUnbd",
                           "MultiObjStopped")

            if(is.null(out))
              term <- paste("Failed to obtain solution, unknown error code:", out)


            return(list(code= out,
                        term = term))
          }
)

setMethod("getSolStat", signature(lp = "LPproblem_cplex"),
          function(lp) {

            out <- getSolStatLP(lp@ptr)

            # get term
            term <- switch(EXPR = out+1,
                           "Unknown",
                           "Feasible",
                           "Optimal",
                           "Infeasible",
                           "Unbounded",
                           "InfeasibleOrUnbounded",
                           "Error")
            if(is.null(out))
              term <- paste("unknown status code:", out)

            return(list(code = out,
                        term = term))
          }
)

setMethod("getObjValue", signature(lp = "LPproblem_cplex"),
          function(lp) {
            out <- getObjVal(lp@ptr)

            return(out)
          }
)

setMethod("getColsPrimal", signature(lp = "LPproblem_cplex"),
          function(lp) {

            out <- getColsPrimalLP(lp@ptr, lp@ptr.env, lp@ptr.x)

            return(out)
          }
)

setMethod("getRedCosts", signature(lp = "LPproblem_cplex"),
          function(lp) {

            out <- getColsDualLP(lp@ptr, lp@ptr.env, lp@ptr.x)

            return(out)
          }
)

setMethod("addSingleConstraint", signature(lp = "LPproblem_cplex"),
          function(lp, coeffs, lb, ub, type) {

            # add new row to constraint matrix
            addRows(lp, nrows = 1)
            i_newrow <- getNumRowsLP(lp@ptr.c)

            nz <- which(coeffs != 0)

            # add coeffs to new row
            setMatRowLP(lp@ptr.x, lp@ptr.c,
                        as.integer(i_newrow-1),
                        as.integer(length(nz)),
                        as.integer(nz-1),
                        as.numeric(coeffs[nz]))

            # set bounds
            setRowsBnds(lp,
                        i = i_newrow,
                        lb = lb,
                        ub = ub,
                        type = type)


          }
)

setMethod("fvaJob", signature(lp = "LPproblem_cplex"),
          function(lp, ind) {

            fvares <- fvaLP(lp@ptr, lp@ptr.mod, lp@ptr.x, lp@ptr.c, lp@ptr.obj,
                            as.integer(ind-1))

            return(fvares)
          }
)


setMethod("deleteLP", signature(lp = "LPproblem_cplex"),
          function(lp) {
            out <- TRUE

            lpXPtrFinalizer(lp@ptr)
            lpmodXPtrFinalizer(lp@ptr.mod)
            lpobjXPtrFinalizer(lp@ptr.obj)
            lpxXPtrFinalizer(lp@ptr.x)
            lpcXPtrFinalizer(lp@ptr.c)
            lpenvXPtrFinalizer(lp@ptr.env)

            return(out)
          }
)
