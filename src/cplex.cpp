#include <iostream>
#include <ilcplex/ilocplex.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// Get CPLEX version number
// [[Rcpp::export]]
SEXP getCPLEXVersion() {
  IloEnv env;
  IloModel model(env);
  IloCplex cplex(model);
  return Rf_ScalarInteger(cplex.getVersionNumber());
}

// Delete CPLEX problem
void lpXPtrFinalizer(SEXP ptr) {
  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(ptr);
  delete cplex;
}

// Delete CPLEX environment
void lpenvXPtrFinalizer(SEXP ptr) {
  IloEnv* env = (IloEnv*)R_ExternalPtrAddr(ptr);
  delete env;
}

// Delete CPLEX model
void lpmodXPtrFinalizer(SEXP ptr) {
  IloModel* model = (IloModel*)R_ExternalPtrAddr(ptr);
  delete model;
}

// Delete CPLEX variable array
void lpxXPtrFinalizer(SEXP ptr) {
  IloNumVarArray* x = (IloNumVarArray*)R_ExternalPtrAddr(ptr);
  delete x;
}

// Delete CPLEX constraint array
void lpcXPtrFinalizer(SEXP ptr) {
  IloRangeArray* c = (IloRangeArray*)R_ExternalPtrAddr(ptr);
  delete c;
}

// Delete CPLEX constraint array
void lpobjXPtrFinalizer(SEXP ptr) {
  IloObjective* obj = (IloObjective*)R_ExternalPtrAddr(ptr);
  delete obj;
}

// Initialize a CPLEX problem
// [[Rcpp::export]]
Rcpp::List initProb(const char* name) {
  IloEnv* env       = new IloEnv();
  IloModel* model   = new IloModel(*env);
  IloCplex* cplex   = new IloCplex(*model);
  IloNumVarArray* x = new IloNumVarArray(*env);
  IloRangeArray* c  = new IloRangeArray(*env);
  IloObjective* obj = new IloObjective(*env);

  cplex->setParam(IloCplex::Param::Simplex::Display, 0); // Suppress output

  // cplex object with pointer
  SEXP xp = R_MakeExternalPtr(cplex, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(xp, lpXPtrFinalizer);

  // model object with pointer
  SEXP xpmod = R_MakeExternalPtr(model, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(xpmod, lpmodXPtrFinalizer);

  // environment object with pointer
  SEXP xpenv = R_MakeExternalPtr(env, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(xpenv, lpenvXPtrFinalizer);

  // variable array object with pointer
  SEXP xpx = R_MakeExternalPtr(x, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(xpx, lpxXPtrFinalizer);

  // constraint array object with pointer
  SEXP xpc = R_MakeExternalPtr(c, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(xpc, lpcXPtrFinalizer);

  // objective object with pointer
  SEXP xpobj = R_MakeExternalPtr(obj, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(xpobj, lpobjXPtrFinalizer);

  Rcpp::List ptrs = Rcpp::List::create(Named("cpx") = xp,
                                       Named("env") = xpenv,
                                       Named("mod") = xpmod,
                                       Named("x")   = xpx,
                                       Named("c")   = xpc,
                                       Named("obj") = xpobj);

  return ptrs;
}

// Set objective direction
// [[Rcpp::export]]
SEXP setObjDirLP(SEXP xpobj, int dir) {
  IloObjective* obj = (IloObjective*)R_ExternalPtrAddr(xpobj);

  if (dir == -1) {
    obj->setSense(IloObjective::Maximize);
  } else {
    obj->setSense(IloObjective::Minimize);
  }

  return R_NilValue;
}

// Add columns
// [[Rcpp::export]]
SEXP addColsLP(SEXP xpenv, SEXP xpx, SEXP ncols) {
  IloEnv* env = (IloEnv*)R_ExternalPtrAddr(xpenv);
  IloNumVarArray* x = (IloNumVarArray*)R_ExternalPtrAddr(xpx);
  int numCols = Rf_asInteger(ncols);

  x->add(IloNumVarArray(*env, IloInt(numCols), -IloInfinity, IloInfinity, ILOFLOAT));

  // std::cout << "Nr. of cols: " << x->getSize() << std::endl;

  return Rf_ScalarInteger(numCols);
}

// Add rows
// [[Rcpp::export]]
SEXP addRowsLP(SEXP xpenv, SEXP xpc, SEXP nrows) {
  IloRangeArray* c = (IloRangeArray*)R_ExternalPtrAddr(xpc);
  IloEnv* env = (IloEnv*)R_ExternalPtrAddr(xpenv);
  int numRows = Rf_asInteger(nrows);

  c->add(IloRangeArray(*env, IloInt(numRows),IloNum(-IloInfinity), IloNum(IloInfinity)));

  // std::cout << "Nr. of rows: " << c->getSize() << std::endl;

  return Rf_ScalarInteger(numRows);
}

// Load matrix
// [[Rcpp::export]]
SEXP loadMatrixLP(SEXP xpx, SEXP xpc, SEXP ne, SEXP ia, SEXP ja, SEXP ra) {
  SEXP out = R_NilValue;

  IloRangeArray* c = (IloRangeArray*)R_ExternalPtrAddr(xpc);
  IloNumVarArray* x = (IloNumVarArray*)R_ExternalPtrAddr(xpx);

  std::cout << "Nr. of rows: " << c->getSize() << std::endl;
  std::cout << "Nr. of cols: " << x->getSize() << std::endl;

  const int *ria = INTEGER(ia);
  const int *rja = INTEGER(ja);
  const double *rra = REAL(ra);
  const int rne = Rf_asInteger(ne);

  for (int i = 0; i < rne; i++) {
    (*c)[ria[i]].setLinearCoef((*x)[rja[i]], rra[i]);
  }

  return R_NilValue;
}


// Set column bounds and objective coefficients
// [[Rcpp::export]]
SEXP setColsBndsObjCoefsLP(SEXP xpobj, SEXP xpx, SEXP indices, SEXP lb, SEXP ub, SEXP objCoefs) {
  IloObjective* obj = (IloObjective*)R_ExternalPtrAddr(xpobj);
  IloNumVarArray* x = (IloNumVarArray*)R_ExternalPtrAddr(xpx);

  IntegerVector inds(indices);
  NumericVector lbs(lb);
  NumericVector ubs(ub);
  NumericVector coefs(objCoefs);

  int n = inds.size();
  for (int i = 0; i < n; i++) {
    obj->setLinearCoef((*x)[inds[i]], coefs[i]);
    (*x)[inds[i]].setBounds(lbs[i],ubs[i]);
  }

  return R_NilValue;
}

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


// Set row bounds
// [[Rcpp::export]]
SEXP setRowsBndsLP(SEXP xpc, SEXP indices, SEXP lb, SEXP ub) {
  IloRangeArray* c = (IloRangeArray*)R_ExternalPtrAddr(xpc);
  IntegerVector inds(indices);
  NumericVector lbs(lb);
  NumericVector ubs(ub);

  for (int i = 0; i < inds.size(); i++) {
    (*c)[inds[i]].setBounds(lbs[i], ubs[i]);
  }

  return R_NilValue;
}

// solve LP with cplex
// [[Rcpp::export]]
SEXP solveCPLEX(SEXP xp,SEXP xpmod, SEXP xpx, SEXP xpc, SEXP xpobj) {
  IloObjective* obj = (IloObjective*)R_ExternalPtrAddr(xpobj);
  IloNumVarArray* x = (IloNumVarArray*)R_ExternalPtrAddr(xpx);
  IloRangeArray* c = (IloRangeArray*)R_ExternalPtrAddr(xpc);
  IloModel* model = (IloModel*)R_ExternalPtrAddr(xpmod);
  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);

  model->add(*x);
  model->add(*c);
  model->add(*obj);

  cplex->solve();


  std::cout << "Obj: " << cplex->getObjValue() << std::endl;

  return R_NilValue;
}


