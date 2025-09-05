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
// [[Rcpp::export]]
void lpXPtrFinalizer(SEXP ptr) {
  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(ptr);
  if (cplex) {
    cplex->end();
    delete cplex;
    R_ClearExternalPtr(ptr);
  }
}

// Delete CPLEX environment
// [[Rcpp::export]]
void lpenvXPtrFinalizer(SEXP ptr) {
  IloEnv* env = (IloEnv*)R_ExternalPtrAddr(ptr);
  if (env) {
    env->end();
    delete env;
    R_ClearExternalPtr(ptr);
  }
}

// Delete CPLEX model
// [[Rcpp::export]]
void lpmodXPtrFinalizer(SEXP ptr) {
  IloModel* model = (IloModel*)R_ExternalPtrAddr(ptr);
  if (model) {
    model->end();
    delete model;
    R_ClearExternalPtr(ptr);
  }
}

// Delete CPLEX variable array
// [[Rcpp::export]]
void lpxXPtrFinalizer(SEXP ptr) {
  IloNumVarArray* x = (IloNumVarArray*)R_ExternalPtrAddr(ptr);
  if (x) {
    x->end();
    delete x;
    R_ClearExternalPtr(ptr);
  }
  delete x;
}

// Delete CPLEX constraint array
// [[Rcpp::export]]
void lpcXPtrFinalizer(SEXP ptr) {
  IloRangeArray* c = (IloRangeArray*)R_ExternalPtrAddr(ptr);
  if (c) {
    c->end();
    delete c;
    R_ClearExternalPtr(ptr);
  }
  delete c;
}

// Delete CPLEX objective coeff array
// [[Rcpp::export]]
void lpobjXPtrFinalizer(SEXP ptr) {
  IloObjective* obj = (IloObjective*)R_ExternalPtrAddr(ptr);
  if (obj) {
    obj->end();
    delete obj;
    R_ClearExternalPtr(ptr);
  }
}

// Initialize a CPLEX problem
// [[Rcpp::export]]
Rcpp::List initProb(const char* name, double tol_bnd) {
  IloEnv* env       = new IloEnv();
  IloModel* model   = new IloModel(*env);
  IloCplex* cplex   = new IloCplex(*model);
  IloNumVarArray* x = new IloNumVarArray(*env);
  IloRangeArray* c  = new IloRangeArray(*env);
  IloObjective* obj = new IloObjective(*env);

  // set central parameters (Output control and algorithm)
  cplex->setOut(env->getNullStream());
  cplex->setWarning(env->getNullStream());

  //set feasibility/bound tolerance for simplex
  cplex->setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, tol_bnd);

  // cplex object with pointer
  SEXP xp = R_MakeExternalPtr(cplex, R_NilValue, R_NilValue);
  // R_RegisterCFinalizer(xp, lpXPtrFinalizer);

  // model object with pointer
  SEXP xpmod = R_MakeExternalPtr(model, R_NilValue, R_NilValue);
  // R_RegisterCFinalizer(xpmod, lpmodXPtrFinalizer);

  // environment object with pointer
  SEXP xpenv = R_MakeExternalPtr(env, R_NilValue, R_NilValue);
  // R_RegisterCFinalizer(xpenv, lpenvXPtrFinalizer);

  // variable array object with pointer
  SEXP xpx = R_MakeExternalPtr(x, R_NilValue, R_NilValue);
  // R_RegisterCFinalizer(xpx, lpxXPtrFinalizer);

  // constraint array object with pointer
  SEXP xpc = R_MakeExternalPtr(c, R_NilValue, R_NilValue);
  // R_RegisterCFinalizer(xpc, lpcXPtrFinalizer);

  // objective object with pointer
  SEXP xpobj = R_MakeExternalPtr(obj, R_NilValue, R_NilValue);
  // R_RegisterCFinalizer(xpobj, lpobjXPtrFinalizer);

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

  // std::cout << "Nr. of rows: " << c->getSize() << std::endl;
  // std::cout << "Nr. of cols: " << x->getSize() << std::endl;

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

  int out;

  IloObjective* obj = (IloObjective*)R_ExternalPtrAddr(xpobj);
  IloNumVarArray* x = (IloNumVarArray*)R_ExternalPtrAddr(xpx);
  IloRangeArray* c = (IloRangeArray*)R_ExternalPtrAddr(xpc);
  IloModel* model = (IloModel*)R_ExternalPtrAddr(xpmod);
  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);

  model->add(*x);
  model->add(*c);
  model->add(*obj);

  cplex->solve();
  out = cplex->getCplexStatus();

  return Rf_ScalarInteger(out);
}


// [[Rcpp::export]]
SEXP getSolStatLP(SEXP xp) {

  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);

  SEXP out = R_NilValue;
  int stat = 0;

  stat = cplex->getStatus();

  out = Rf_ScalarInteger(stat);

  return out;
}

// [[Rcpp::export]]
SEXP getObjVal(SEXP xp) {

  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);

  SEXP out = R_NilValue;
  double obj;
  if (cplex->getStatus() == IloAlgorithm::Optimal || cplex->getStatus() == IloAlgorithm::Feasible) {
    obj = cplex->getObjValue();
  } else {
    return DoubleVector(1, DoubleVector::get_na());
  }

  out = Rf_ScalarReal(obj);

  return out;
}

// [[Rcpp::export]]
SEXP getColsPrimalLP(SEXP xp, SEXP xpenv, SEXP xpx) {
  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);
  IloNumVarArray* x = (IloNumVarArray*)R_ExternalPtrAddr(xpx);
  IloEnv* env = (IloEnv*)R_ExternalPtrAddr(xpenv);

  // Check the solver status before attempting to get values
  if (cplex->getStatus() != IloAlgorithm::Optimal && cplex->getStatus() != IloAlgorithm::Feasible) {
    // Return a vector of NAs if no feasible or optimal solution exists
    return DoubleVector(x->getSize(), DoubleVector::get_na());
  }

  IloNumArray vals(*env);
  cplex->getValues(vals, *x);

  std::vector<double> prim;

  for(IloInt i = 0; i < vals.getSize(); i++) {
    prim.push_back(vals[i]);
  }

  return Rcpp::wrap(prim);
}

// [[Rcpp::export]]
SEXP getColsDualLP(SEXP xp, SEXP xpenv, SEXP xpx) {
  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);
  IloNumVarArray* x = (IloNumVarArray*)R_ExternalPtrAddr(xpx);
  IloEnv* env = (IloEnv*)R_ExternalPtrAddr(xpenv);

  // Check the solver status before attempting to get values
  if (cplex->getStatus() != IloAlgorithm::Optimal && cplex->getStatus() != IloAlgorithm::Feasible) {
    // Return a vector of NAs if no feasible or optimal solution exists
    return DoubleVector(x->getSize(), DoubleVector::get_na());
  }

  IloNumArray vals(*env);
  cplex->getReducedCosts(vals, *x);

  std::vector<double> dual;

  for(IloInt i = 0; i < vals.getSize(); i++) {
    dual.push_back(vals[i]);
  }

  return Rcpp::wrap(dual);
}

/* get number of rows */
// [[Rcpp::export]]
SEXP getNumRowsLP(SEXP xpc) {

  IloRangeArray* c = (IloRangeArray*)R_ExternalPtrAddr(xpc);

  SEXP out = R_NilValue;

  out = Rf_ScalarInteger(c->getSize());

  return out;
}

/* set or replace row of constraint matrix */
// [[Rcpp::export]]
SEXP setMatRowLP(SEXP xpx, SEXP xpc, SEXP i, SEXP len, SEXP ind, SEXP val) {
  SEXP out = R_NilValue;

  IloRangeArray* c = (IloRangeArray*)R_ExternalPtrAddr(xpc);
  IloNumVarArray* x = (IloNumVarArray*)R_ExternalPtrAddr(xpx);

  int ri = Rf_asInteger(i);
  int rlen = Rf_asInteger(len);
  IntegerVector rind(ind);
  NumericVector rval(val);

  for (int k = 0; k < rlen; k++) {
    (*c)[ri].setLinearCoef((*x)[rind[k]], rval[k]);
  }

  return out;
}

/* wrapper for FVA */
// [[Rcpp::export]]
Rcpp::DataFrame fvaLP(SEXP xp,SEXP xpmod, SEXP xpx, SEXP xpc, SEXP xpobj, SEXP ind) {

  IloObjective* obj = (IloObjective*)R_ExternalPtrAddr(xpobj);
  IloNumVarArray* x = (IloNumVarArray*)R_ExternalPtrAddr(xpx);
  IloRangeArray* c = (IloRangeArray*)R_ExternalPtrAddr(xpc);
  IloModel* model = (IloModel*)R_ExternalPtrAddr(xpmod);
  IloCplex* cplex = (IloCplex*)R_ExternalPtrAddr(xp);
  model->add(*c);

  // Get the indices as integers
  IntegerVector indices(ind);

  // overwrite current objective function
  for(int i = 0; i < x->getSize(); i++) {
    obj->setLinearCoef((*x)[i], 0);
  }

  // model->add(*x);
  // model->add(*c);
  // model->add(*obj);

  // Create vectors to store results
  NumericVector minVals(indices.size());
  NumericVector maxVals(indices.size());

  // MAX
  obj->setSense(IloObjective::Maximize);
  for(unsigned int i = 0; i < indices.size(); i++) {
    int columnIndex = indices[i];
    obj->setLinearCoef((*x)[columnIndex], 1.0);

    // Maximize the variable
    cplex->solve();

    // Store the maximum value
    maxVals[i] = cplex->getObjValue();

    // Reset the objective coefficient to 0
    obj->setLinearCoef((*x)[columnIndex], 0.0);
  }

  // MIN
  obj->setSense(IloObjective::Minimize);
  for(unsigned int i = 0; i < indices.size(); i++) {
    int columnIndex = indices[i];
    obj->setLinearCoef((*x)[columnIndex], 1.0);

    // Minimize the variable
    cplex->solve();

    // Store the minimum value
    minVals[i] = cplex->getObjValue();

    // Reset the objective coefficient to 0
    obj->setLinearCoef((*x)[columnIndex], 0.0);
  }

  // Create a DataFrame to store the results
  DataFrame result = DataFrame::create(_["min.flux"] = minVals, _["max.flux"] = maxVals);

  return result;
}

