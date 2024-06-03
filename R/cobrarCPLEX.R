## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib cobrarCPLEX, .registration = TRUE
#' @importFrom methods .valueClassTest as is new
## usethis namespace: end
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Using cplex version ",getCPLEXVersion())
}

.onLoad <- function(lib, pkg) {
  COBRAR_SETTINGS("SOLVER","cplex")
}

