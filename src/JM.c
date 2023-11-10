#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
// #include <omp.h>

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#define USE_FC_LEN_T

#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Complex.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>


/* (1) Helper functions. */
SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt, names;
  PROTECT(elmt = R_NilValue);
  PROTECT(names = getAttrib(list, R_NamesSymbol));

  for(int i = 0; i < length(list); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }

  UNPROTECT(2);

  return elmt;
}

void pvec(SEXP vec)
{
  int i, n = length(vec);
  double *vecptr = REAL(vec);
  for(i = 0; i < n; i++, vecptr++)
    Rprintf(" %g", *vecptr);
  Rprintf("\n");
}

void pmat(SEXP mat)
{
  int i,j,n = nrows(mat), k = ncols(mat);
  Rprintf("   ");
  for(j = 0; j < k; ++j)
    Rprintf("[%d] ", j);
  Rprintf("\n");
  for(i = 0; i < n; ++i) {
    Rprintf("[%d]", i);
    for(j = 0; j < k; ++j)
      Rprintf(" %g", REAL(mat)[i + j * n]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

void merr()
{
  char *m = "stopped";
  error(m);
}


/* (2) Suvival integral function */
SEXP survint(SEXP pred, SEXP pre_fac, SEXP pre_vec, SEXP omega,
  SEXP int_fac, SEXP int_vec, SEXP weights, SEXP survtime)
{
  int nProtected = 0;
  int nw = length(weights);
  int nsubj = length(survtime);
  int predictor = INTEGER(pred)[0];

  int p = ncols(int_vec);
  if(predictor == 2)
    p = ncols(pre_vec);

  SEXP score_int;
  PROTECT(score_int = allocVector(REALSXP, p));
  ++nProtected;

  SEXP hess_int;
  PROTECT(hess_int = allocVector(REALSXP, p * p));
  ++nProtected;

  SEXP hess;
  PROTECT(hess = allocMatrix(REALSXP, p, p));
  ++nProtected;

  SEXP hess_i;
  PROTECT(hess_i = allocMatrix(REALSXP, p, p));
  ++nProtected;

  SEXP score_i;
  PROTECT(score_i = allocVector(REALSXP, p));
  ++nProtected;

  SEXP hess_vec_i;
  PROTECT(hess_vec_i = allocMatrix(REALSXP, p, p));
  ++nProtected;

  SEXP score_vec_i;
  PROTECT(score_vec_i = allocVector(REALSXP, p));
  ++nProtected;

  // Iterators.
  int i, j, ii, jj, nr, nc;

  // Pointers.
  double *weights_ptr = REAL(weights);
  double *omega_ptr = REAL(omega);
  double *int_vec_ptr = REAL(int_vec);
  double *hess_ptr = REAL(hess);
  double *score_i_ptr = REAL(score_i);
  double *hess_i_ptr = REAL(hess_i);
  double *survtime_ptr = REAL(survtime);
  double *pre_fac_ptr = REAL(pre_fac);
  double *score_int_ptr = REAL(score_int);
  double *hess_int_ptr = REAL(hess_int);
  double *pre_vec_ptr = REAL(pre_vec);
  double *hess_vec_i_ptr = REAL(hess_vec_i);
  double *score_vec_i_ptr = REAL(score_vec_i);
  double *int_fac_ptr = REAL(int_fac);

  // Others.
  double tmp = 0.0;

  nr = nrows(int_vec);
  nc = ncols(int_vec);

  // Initialize output. 
  for(ii = 0; ii < p; ii++){
    score_int_ptr[ii] = 0.0;
  }
  for(jj = 0; jj < p*p; jj++){
    hess_int_ptr[jj] = 0.0;
  }

  // Lambda.
  if(predictor == 1) {
    for(i = 0; i < nsubj; i++) {
      for(ii = 0; ii < p; ii++) {
        score_i_ptr[ii] = 0.0;
        for(jj = 0; jj < p; jj++) {
          hess_i_ptr[ii + p * jj] = 0.0;
        }
      }
      for(j = 0; j < nw; j++) {
        tmp = weights_ptr[j] * omega_ptr[i * nw + j];
        for(ii = 0; ii < p; ii++) {
          for(jj = 0; jj < p; jj++) {
            if(ii <= jj) {
              hess_ptr[ii + jj * p] = int_vec_ptr[i * nw + j + ii * nr] * int_vec_ptr[i * nw + j + jj * nr];
              hess_i_ptr[ii + jj * p] += tmp * hess_ptr[ii + jj * p];
            }
            hess_ptr[jj + ii * p] = hess_ptr[ii + jj * p];
            hess_i_ptr[jj + ii * p] = hess_i_ptr[ii + jj * p];
          }
          score_i_ptr[ii] += tmp * int_vec_ptr[i * nw + j + ii * nr];
        }
      }

      tmp = survtime_ptr[i] / 2.0 * pre_fac_ptr[i];

      for(ii = 0; ii < p * p; ii++) {
        if(ii < p)
          score_int_ptr[ii] += tmp * score_i_ptr[ii];
        hess_int_ptr[ii] += tmp * hess_i_ptr[ii];
      }
    }
  }

  // Gamma.
  if(predictor == 2) {
    nr = nrows(pre_vec);
    for(i = 0; i < nsubj; i++) {
      tmp = 0.0;
      for(j = 0; j < nw; j++) {
        tmp += weights_ptr[j] * omega_ptr[i * nw + j];
      }
      tmp = survtime_ptr[i] / 2.0 * tmp * pre_fac_ptr[i];
      j = 0;
      for(ii = 0; ii < p; ii++) {
        score_int_ptr[ii] += tmp * pre_vec_ptr[i + ii * nr];
        for(jj = 0; jj < p; jj++) {
          hess_int_ptr[j] += tmp * pre_vec_ptr[i + ii * nr] * pre_vec_ptr[i + jj * nr];
          j += 1;
        }
      }
    }
  }

  // Long.
  if(predictor == 3) {
    int nmarker = nrows(int_vec) / (nsubj * nw);

    for(i = 0; i < nsubj; i++) {
      for(ii = 0; ii < p; ii++) {
        score_i_ptr[ii] = 0.0;
        for(jj = 0; jj < p; jj++) {
          hess_i_ptr[ii + p * jj] = 0.0;
        }
      }

      for(j = 0; j < nw; j++) {
        tmp = weights_ptr[j] * omega_ptr[i * nw + j];

        for(ii = 0; ii < p; ii++) {
          score_vec_i_ptr[ii] = 0.0;
        }

        for(ii = 0; ii < nmarker; ii++) {
          for(jj = 0; jj < p; jj++) {
            score_vec_i_ptr[jj] += int_fac_ptr[ii * nsubj * nw + i * nw + j] *
              int_vec_ptr[ii * nsubj * nw + i * nw + j + jj * nr];
          }       
        }

        for(ii = 0; ii < p; ii++) {
          for(jj = 0; jj < p; jj++) {
            hess_vec_i_ptr[ii + jj * p] = score_vec_i_ptr[ii] * score_vec_i_ptr[jj];
          }
        }

        for(ii = 0; ii < p; ii++) {
          score_i_ptr[ii] += tmp * score_vec_i_ptr[ii];
          for(jj = 0; jj < p; jj++) {
            hess_i_ptr[ii + jj * p] += tmp * hess_vec_i_ptr[ii + jj * p];
          }
        }
      }

      tmp = survtime_ptr[i] / 2.0 * pre_fac_ptr[i];

      j = 0;
      for(ii = 0; ii < p; ii++) {
        score_int_ptr[ii] += tmp * score_i_ptr[ii];
        for(jj = 0; jj < p; jj++) {
          hess_int_ptr[j] += tmp * hess_i_ptr[ii + jj * p];
          j += 1;
        }
      }
    }
  }

  /* Stuff everything together. */
  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 2));
  ++nProtected;

  SET_VECTOR_ELT(rval, 0, score_int);
  SET_VECTOR_ELT(rval, 1, hess_int);

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 2));
  ++nProtected;

  SET_STRING_ELT(nrval, 0, mkChar("score_int"));
  SET_STRING_ELT(nrval, 1, mkChar("hess_int"));
        
  setAttrib(rval, R_NamesSymbol, nrval);

  UNPROTECT(nProtected);

  return rval;
}

/* (3) Suvival integral function for FPC-Random Effects */
SEXP survint_re(SEXP pre_fac, SEXP omega,
  SEXP int_fac, SEXP int_vec, SEXP weights, SEXP survtime)
{
  int nProtected = 0;
  int nw = length(weights);

  int p = ncols(int_vec);

  SEXP score_int;
  PROTECT(score_int = allocVector(REALSXP, p));
  ++nProtected;

  SEXP hess_int;
  PROTECT(hess_int = allocVector(REALSXP, p));
  ++nProtected;

  SEXP hess_i;
  PROTECT(hess_i = allocVector(REALSXP, 1));
  ++nProtected;

  SEXP score_i;
  PROTECT(score_i = allocVector(REALSXP, 1));
  ++nProtected;

  SEXP hess_vec_i;
  PROTECT(hess_vec_i = allocVector(REALSXP, 1));
  ++nProtected;

  SEXP score_vec_i;
  PROTECT(score_vec_i = allocVector(REALSXP, 1));
  ++nProtected;

  // Iterators.
  int i, j, ii, jj, nr, nc;

  // Pointers.
  double *weights_ptr = REAL(weights);
  double *omega_ptr = REAL(omega);
  double *pre_fac_ptr = REAL(pre_fac);
  double *int_vec_ptr = REAL(int_vec);
  double *int_fac_ptr = REAL(int_fac);
  double *survtime_ptr = REAL(survtime);
  double *score_i_ptr = REAL(score_i);
  double *hess_i_ptr = REAL(hess_i);
  double *score_int_ptr = REAL(score_int);
  double *hess_int_ptr = REAL(hess_int);
  double *hess_vec_i_ptr = REAL(hess_vec_i);
  double *score_vec_i_ptr = REAL(score_vec_i);

  // Others.
  double tmp = 0.0;

  nr = nrows(int_vec);
  nc = ncols(int_vec);

  // Initialize output. 
  for(i = 0; i < p; i++){
    score_int_ptr[i] = 0.0;
  }
  for(i = 0; i < p; i++){
    hess_int_ptr[i] = 0.0;
  }

  // Long.
  int nmarker = nrows(int_vec) / (p * nw);

  for(i = 0; i < p; i++) {
      
    score_i_ptr[0] = 0.0;
    hess_i_ptr[0] = 0.0;
        

    for(j = 0; j < nw; j++) {
      tmp = weights_ptr[j] * omega_ptr[i * nw + j];

      score_vec_i_ptr[0] = 0.0;
        

      for(ii = 0; ii < nmarker; ii++) {
        score_vec_i_ptr[0] += int_fac_ptr[ii * p * nw + i * nw + j] *
          int_vec_ptr[ii * p * nw + i * nw + j + i * nr];      
      }

      hess_vec_i_ptr[0] = score_vec_i_ptr[0] * score_vec_i_ptr[0];


      score_i_ptr[0] += tmp * score_vec_i_ptr[0];
      hess_i_ptr[0] += tmp * hess_vec_i_ptr[0];
    }

    tmp = survtime_ptr[i] / 2.0 * pre_fac_ptr[i];

    score_int_ptr[i] += tmp * score_i_ptr[0];
    hess_int_ptr[i] += tmp * hess_i_ptr[0];
      
  }

  /* Stuff everything together. */
  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 2));
  ++nProtected;

  SET_VECTOR_ELT(rval, 0, score_int);
  SET_VECTOR_ELT(rval, 1, hess_int);

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 2));
  ++nProtected;

  SET_STRING_ELT(nrval, 0, mkChar("score_int"));
  SET_STRING_ELT(nrval, 1, mkChar("hess_int"));
        
  setAttrib(rval, R_NamesSymbol, nrval);

  UNPROTECT(nProtected);

  return rval;
}

/* (4) Crossproduct for Psi-Matrix in Hesse Calculation */
SEXP psi_mat_multiplication(SEXP X, SEXP ni_obs, SEXP diags)
{
  int nProtected = 0;
  int n = nrows(X);
  int p = ncols(X);
  int *ni_obs_ptr = INTEGER(ni_obs);
  double *X_ptr = REAL(X);
  double *diags_ptr = REAL(diags);

  SEXP out;
  PROTECT(out = allocVector(REALSXP, p));
  ++nProtected;

  SEXP out_i;
  PROTECT(out_i = allocVector(REALSXP, 1));
  ++nProtected;

  SEXP col_it;
  PROTECT(col_it = allocVector(INTSXP, 1));
  ++nProtected;


  // Iterators.
  int i, j;

  // Pointers.
  double *out_ptr = REAL(out);
  double *out_i_ptr = REAL(out_i);
  int *col_it_ptr = INTEGER(col_it);

  // initialize column iteration
  col_it_ptr[0] = 0;
  
  // Initialize output and sum over non-zero elements.
  for(i = 0; i < p; i++){
    
    out_i_ptr[0] = 0.0;
    int n_i = ni_obs_ptr[i];
    for(j = 0; j < n_i; j++){
    	
	out_i_ptr[0] += X_ptr[i*n + col_it_ptr[0] + j] * X_ptr[i*n + col_it_ptr[0] + j] * diags_ptr[col_it_ptr[0] + j];
	    
    }
    
    out_ptr[i] = out_i_ptr[0];
    col_it_ptr[0] += n_i;
    
  }

  /* Output. */

  UNPROTECT(nProtected);

  return out;
}


/* (5) Crossproduct for Psi-Matrix in Score Calculation */
SEXP psi_vec_multiplication(SEXP X, SEXP ni_obs, SEXP y)
{
  int nProtected = 0;
  int n = nrows(X);
  int p = ncols(X);
  int *ni_obs_ptr = INTEGER(ni_obs);
  double *X_ptr = REAL(X);
  double *y_ptr = REAL(y);

  SEXP out;
  PROTECT(out = allocVector(REALSXP, p));
  ++nProtected;

  SEXP out_i;
  PROTECT(out_i = allocVector(REALSXP, 1));
  ++nProtected;

  SEXP col_it;
  PROTECT(col_it = allocVector(INTSXP, 1));
  ++nProtected;


  // Iterators.
  int i, j;

  // Pointers.
  double *out_ptr = REAL(out);
  double *out_i_ptr = REAL(out_i);
  int *col_it_ptr = INTEGER(col_it);

  // initialize column iteration
  col_it_ptr[0] = 0;
  
  // Initialize output and sum over non-zero elements.
  for(i = 0; i < p; i++){
    
    out_i_ptr[0] = 0.0;
    int n_i = ni_obs_ptr[i];
    for(j = 0; j < n_i; j++){
    	
	out_i_ptr[0] += X_ptr[i*n + col_it_ptr[0] + j] * y_ptr[col_it_ptr[0] + j];
	    
    }
    
    out_ptr[i] = out_i_ptr[0];
    col_it_ptr[0] += n_i;
    
  }

  /* Output. */

  UNPROTECT(nProtected);

  return out;
}


