

# some GLPK specific codes:
glpkPar <- list(
  CPX_MIN = 1,
  CPX_MAX = -1
)


setClass(Class = "LPproblem_cplex",
         contains = "LPproblem"
)

setMethod(f = "initialize",
          signature = "LPproblem_cplex",
          definition = function(.Object,
                                name,
                                method) {

            .Object@ptr <- cobrarCPLEX::initProb(name)
            .Object@solver = "cplex"
            .Object@method = method

            return(.Object)
          }
)

setMethod("loadLPprob", signature(lp = "LPproblem_cplex"),

          function(lp, nCols, nRows, mat, ub, lb, obj, rlb, rtype, lpdir,
                   rub = NULL, ctype = NULL) {

            # optimization direction
            lpdir <- switch(EXPR = lpdir,
                            "max" = glpkPar$CPX_MAX,
                            "min" = glpkPar$CPX_MIN)
            setObjDirection(lp, lpdir = lpdir)

            # problem dimensions
            addCols(lp, ncols = nCols)
            addRows(lp, nrows = nRows)


            # constraint matrix
            TMPmat <- as(mat, "TsparseMatrix")
            loadMatrix(lp,
                       ne = length(TMPmat@x),
                       ia = TMPmat@i + 1,
                       ja = TMPmat@j + 1,
                       ra = TMPmat@x)

            # column (variable) bounds and objective function
            setColsBndsObjCoefs(lp,
                                j = c(1:nCols),
                                lb = lb,
                                ub = ub,
                                obj_coef = obj)


            # variable type
            if (!is.null(ctype)) {
              cctype <- sapply(ctype,
                               function(x) switch(EXPR = x,
                                                  "C" = glpkPar$GLP_CV,
                                                  "I" = glpkPar$GLP_IV,
                                                  "B" = glpkPar$GLP_BV,
                                                  glpkPar$GLP_CV))

              setColsKind(lp, j = c(1:nCols), kind = cctype)
            }

            # right hand side
            # Note: This is the steady state condition: Production of internal
            # metabolites should equal the consumption. In over words:
            # row lower bound = row upper bounds = 0
            if (is.null(rub)) {
              # The values in rlb will be copied to rub. GLPK ignores rlb and rub,
              # depending on the constraint type (e.g. an upper bound, if the
              # constraint type says, it has a lower bound):
              # Constraint type "L": ignore rub
              # Constraint type "U": ignore rlb
              # Constraint type "E": ignore rub
              # Constraint type "F": ignore rlb and rub

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
            setObjDirLP(lp@ptr, lpdir)
          }
)

setMethod("addCols", signature(lp = "LPproblem_cplex"),
          function(lp, ncols) {
            addColsLP(lp@ptr, as.integer(ncols))
          }
)

setMethod("addRows", signature(lp = "LPproblem_cplex"),
          function(lp, nrows) {
            addRowsLP(lp@ptr, as.integer(nrows))
          }
)

setMethod("loadMatrix", signature(lp = "LPproblem_cplex"),
          function(lp, ne, ia, ja, ra) {
            loadMatrixLP(lp@ptr,
                         as.integer(ne),
                         as.integer(ia),
                         as.integer(ja),
                         as.numeric(ra))
          }
)

setMethod("setColsBndsObjCoefs", signature(lp = "LPproblem_cplex"),
          function(lp, j, lb, ub, obj_coef, type = NULL) {


            if (is.null(type)) {
              Ctype <- as.null(type)
            }
            else {
              Ctype <- as.integer(type)
            }

            setColsBndsObjCoefsLP(lp@ptr,
                                  as.integer(j),
                                  Ctype,
                                  as.numeric(lb),
                                  as.numeric(ub),
                                  as.numeric(obj_coef))
          }
)

setMethod("setColsKind", signature(lp = "LPproblem_cplex"),
          function(lp, j, kind) {
            setColsKindLP(lp@ptr,
                          as.integer(j),
                          as.integer(kind))
          }
)


setMethod("setRowsBnds", signature(lp = "LPproblem_cplex"),
          function(lp, i, lb, ub , type) {


            type <- sapply(type,
                           function(x) switch(EXPR = x,
                                              "F" = glpkPar$GLP_FR,
                                              "L" = glpkPar$GLP_LO,
                                              "U" = glpkPar$GLP_UP,
                                              "D" = glpkPar$GLP_DB,
                                              "E" = glpkPar$GLP_FX,
                                              glpkPar$GLP_FX))

            if (is.null(type)) {
              Ctype <- as.null(type)
            }
            else {
              Ctype <- as.integer(type)
            }

            setRowsBndsLP(lp@ptr,
                          as.integer(i),
                          Ctype,
                          as.numeric(lb),
                          as.numeric(ub))

          }
)

