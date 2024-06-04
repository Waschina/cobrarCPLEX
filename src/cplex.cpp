#include <iostream>
#include <ilcplex/ilocplex.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

IloEnv env;

// Get CPLEX version number
// [[Rcpp::export]]
SEXP getCPLEXVersion() {
  IloModel model(env);
  IloCplex cplex(model);
  return Rf_ScalarInteger(cplex.getVersionNumber());
}

// Delete CPLEX problem
void lpXPtrFinalizer(SEXP lp_ptr) {
  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(lp_ptr);
  delete cplex;
}

// Initialize a CPLEX problem
// [[Rcpp::export]]
SEXP initProb(const char* name) {
  IloModel model(env);
  IloCplex *cplex = new IloCplex(model);

  cplex->setParam(IloCplex::Param::Simplex::Display, 0); // suppress output

  SEXP xp = R_MakeExternalPtr(cplex, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(xp, lpXPtrFinalizer);

  return xp;
}

// Set objective direction
// [[Rcpp::export]]
SEXP setObjDirLP(SEXP xp, int dir) {
  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);
  IloObjective obj = cplex->getObjective();

  if (dir == -1) { // Assume -1 for Max, 1 for Min
    obj.setSense(IloObjective::Maximize);
  } else {
    obj.setSense(IloObjective::Minimize);
  }

  return R_NilValue;
}

// Add columns
// [[Rcpp::export]]
SEXP addColsLP(SEXP xp, SEXP ncols) {
  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);
  int numCols = Rf_asInteger(ncols);

  IloModel model = cplex->getModel();
  IloEnv env = model.getEnv();
  IloNumVarArray x(cplex->getModel().getEnv(), IloInt(numCols));

  return Rf_ScalarInteger(numCols);
}

// Add rows
// [[Rcpp::export]]
SEXP addRowsLP(SEXP xp, SEXP nrows) {
  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);

  IloModel model = cplex->getModel();
  IloEnv env = model.getEnv();

  int numRows = Rf_asInteger(nrows);

  IloRangeArray c(env);

  c.setSize(numRows);

  // Update model
  model.add(c);

  return Rf_ScalarInteger(numRows);
}

// Load matrix
// [[Rcpp::export]]
SEXP loadMatrixLP(SEXP xp, SEXP ne, SEXP ia, SEXP ja, SEXP ra) {
  SEXP out = R_NilValue;

  const int *ria = INTEGER(ia);
  const int *rja = INTEGER(ja);
  const double *rra = REAL(ra);
  const int rne = Rf_asInteger(ne);

  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);
  IloModel model = cplex->getModel();
  IloEnv env = model.getEnv();

  IloNumVarArray x(env);
  IloRangeArray c(env);

  for (int i = 0; i < rne; i++) {
    c[i].setLinearCoef(x[i], rra[i]);
  }

  // Update model
  model.add(c);

  return R_NilValue;
}
//
// // Set column bounds and objective coefficients
// // [[Rcpp::export]]
// SEXP setColsBndsObjCoefsLP(SEXP xp, SEXP indices, SEXP lb, SEXP ub, SEXP objCoefs) {
//   IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);
//   IntegerVector inds(indices);
//   NumericVector lbs(lb);
//   NumericVector ubs(ub);
//   NumericVector coefs(objCoefs);
//
//   int n = inds.size();
//   for (int i = 0; i < n; i++) {
//     IloNumVar var = cplex->getVar(inds[i]);
//     var.setBounds(lbs[i], ubs[i]);
//     cplex->getObjective().setLinearCoef(var, coefs[i]);
//   }
//
//   return R_NilValue;
// }
//
// // Set column kinds (integer, binary, continuous)
// // [[Rcpp::export]]
// SEXP setColsKindLP(SEXP xp, SEXP indices, SEXP kinds) {
//   IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);
//   IntegerVector inds(indices);
//   IntegerVector knds(kinds);
//
//   for (int i = 0; i < inds.size(); i++) {
//     IloNumVar var = cplex->getVar(inds[i]);
//     switch (knds[i]) {
//     case 0: var.setType(IloNumVar::Float); break;  // Continuous
//     case 1: var.setType(IloNumVar::Int); break;    // Integer
//     case 2: var.setType(IloNumVar::Bool); break;   // Binary
//     }
//   }
//
//   return R_NilValue;
// }
//
// // Set row bounds
// // [[Rcpp::export]]
// SEXP setRowsBndsLP(SEXP xp, SEXP indices, SEXP lb, SEXP ub) {
//   IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);
//   IntegerVector inds(indices);
//   NumericVector lbs(lb);
//   NumericVector ubs(ub);
//
//   for (int i = 0; i < inds.size(); i++) {
//     IloRange range = cplex->getRange(inds[i]);
//     range.setBounds(lbs[i], ubs[i]);
//   }
//
//   return R_NilValue;
// }
//



