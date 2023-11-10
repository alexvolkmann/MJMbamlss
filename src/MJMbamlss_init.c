#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define USE_FC_LEN_T

SEXP survint(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP survint_re(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP psi_mat_multiplication(SEXP, SEXP, SEXP);
SEXP psi_vec_multiplication(SEXP, SEXP, SEXP);


static R_CallMethodDef callMethods[] = {
  {"survint", (DL_FUNC) &survint, 8},
  {"survint_re", (DL_FUNC) &survint_re, 6},
  {"psi_mat_multiplication", (DL_FUNC) &psi_mat_multiplication, 3},
  {"psi_vec_multiplication", (DL_FUNC) &psi_vec_multiplication, 3},
  {NULL, NULL, 0}
};

void R_init_sourcetools(DllInfo* info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

