#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <R.h>
#include "poissonbinom.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
  using namespace Eigen; 
  
  SEXP C_poisbinom(SEXP probSEXP) {
    // 1. INPUT: Map R memory to Eigen (No Rcpp::as)
    double* p_prob = REAL(probSEXP);
    int rows = Rf_nrows(probSEXP);
    int cols = Rf_ncols(probSEXP);
    Eigen::Map<Eigen::MatrixXd> prob(p_prob, rows, cols);
    
    // 2. LOGIC:
    Eigen::MatrixXd out = dpoisbinom(prob);
    
    // 3. OUTPUT: Allocate R matrix and copy Eigen data (No Rcpp::wrap)
    // Create a new numeric matrix in R memory
    SEXP resSEXP = PROTECT(Rf_allocMatrix(REALSXP, out.rows(), out.cols()));
    
    // Copy data from Eigen to the new R matrix memory
    // Eigen's .data() returns a pointer; std::copy is safe and fast
    std::copy(out.data(), out.data() + out.size(), REAL(resSEXP));
    
    UNPROTECT(1); // Release the protection on resSEXP
    return resSEXP;
  }
  
  SEXP C_set_omp_threads(SEXP tSEXP) {
    int t = INTEGER(tSEXP)[0];
#ifdef _OPENMP
    omp_set_num_threads(t);
#endif
    return R_NilValue; // nothing to return
  }
  
  SEXP C_get_omp_threads(void) {
    int n;
#ifdef _OPENMP
    n = omp_get_max_threads();
#else
    n = 1;
#endif
    return Rf_ScalarInteger(n);  // wrap int as SEXP
  }
  
  const static R_CallMethodDef R_CallDef[] = {
    TMB_CALLDEFS,
    CALLDEF(C_poisbinom, 1),
    CALLDEF(C_get_omp_threads, 0),
    CALLDEF(C_set_omp_threads, 1),
    {NULL, NULL, 0}
  };
  
  void R_init_gllvm(DllInfo *dll)
  {
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
#ifdef TMB_CCALLABLES
    TMB_CCALLABLES("gllvm");
#endif
  }
  
}

