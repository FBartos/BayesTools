#include "invgamma/BTInvGammaCore.h"

#include <cmath>

#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Random.h>

namespace {

SEXP coerce_numeric(SEXP x, const char *name)
{
  if(TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP){
    Rf_error("'%s' must be numeric.", name);
  }

  return Rf_coerceVector(x, REALSXP);
}

double scalar_numeric(SEXP x, const char *name)
{
  if(Rf_length(x) != 1){
    Rf_error("'%s' must be a numeric scalar.", name);
  }
  SEXP x_real = PROTECT(coerce_numeric(x, name));
  double out = REAL(x_real)[0];
  UNPROTECT(1);

  return out;
}

int scalar_int(SEXP x, const char *name)
{
  if(Rf_length(x) != 1){
    Rf_error("'%s' must be an integer scalar.", name);
  }
  SEXP x_int = PROTECT(Rf_coerceVector(x, INTSXP));
  int out = INTEGER(x_int)[0];
  UNPROTECT(1);
  if(out < 0){
    Rf_error("'%s' must be non-negative.", name);
  }

  return out;
}

bool scalar_bool(SEXP x, const char *name)
{
  if(Rf_length(x) != 1){
    Rf_error("'%s' must be a logical scalar.", name);
  }
  int out = Rf_asLogical(x);
  if(out == NA_LOGICAL){
    Rf_error("'%s' must not be NA.", name);
  }

  return out == TRUE;
}

}

extern "C" SEXP BayesTools_invgamma_d(SEXP x, SEXP shape, SEXP scale, SEXP log)
{
  SEXP x_real = PROTECT(coerce_numeric(x, "x"));
  double shape_value = scalar_numeric(shape, "shape");
  double scale_value = scalar_numeric(scale, "scale");
  bool log_value = scalar_bool(log, "log");

  R_xlen_t n = XLENGTH(x_real);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *out_ptr = REAL(out);
  double const *x_ptr = REAL(x_real);
  for(R_xlen_t i = 0; i < n; ++i){
    double log_density = bayestools::invgamma::log_density(
      x_ptr[i], shape_value, scale_value
    );
    out_ptr[i] = log_value ? log_density : std::exp(log_density);
  }

  UNPROTECT(2);
  return out;
}

extern "C" SEXP BayesTools_invgamma_p(SEXP q, SEXP shape, SEXP scale,
                                      SEXP lower_tail, SEXP log_p)
{
  SEXP q_real = PROTECT(coerce_numeric(q, "q"));
  double shape_value = scalar_numeric(shape, "shape");
  double scale_value = scalar_numeric(scale, "scale");
  bool lower_tail_value = scalar_bool(lower_tail, "lower.tail");
  bool log_p_value = scalar_bool(log_p, "log.p");

  R_xlen_t n = XLENGTH(q_real);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *out_ptr = REAL(out);
  double const *q_ptr = REAL(q_real);
  for(R_xlen_t i = 0; i < n; ++i){
    out_ptr[i] = bayestools::invgamma::cdf(
      q_ptr[i], shape_value, scale_value, lower_tail_value, log_p_value
    );
  }

  UNPROTECT(2);
  return out;
}

extern "C" SEXP BayesTools_invgamma_q(SEXP p, SEXP shape, SEXP scale,
                                      SEXP lower_tail, SEXP log_p)
{
  SEXP p_real = PROTECT(coerce_numeric(p, "p"));
  double shape_value = scalar_numeric(shape, "shape");
  double scale_value = scalar_numeric(scale, "scale");
  bool lower_tail_value = scalar_bool(lower_tail, "lower.tail");
  bool log_p_value = scalar_bool(log_p, "log.p");

  R_xlen_t n = XLENGTH(p_real);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *out_ptr = REAL(out);
  double const *p_ptr = REAL(p_real);
  for(R_xlen_t i = 0; i < n; ++i){
    out_ptr[i] = bayestools::invgamma::quantile(
      p_ptr[i], shape_value, scale_value, lower_tail_value, log_p_value
    );
  }

  UNPROTECT(2);
  return out;
}

extern "C" SEXP BayesTools_invgamma_r(SEXP n, SEXP shape, SEXP scale)
{
  int n_value = scalar_int(n, "n");
  double shape_value = scalar_numeric(shape, "shape");
  double scale_value = scalar_numeric(scale, "scale");

  SEXP out = PROTECT(Rf_allocVector(REALSXP, n_value));
  double *out_ptr = REAL(out);
  GetRNGstate();
  for(int i = 0; i < n_value; ++i){
    out_ptr[i] = bayestools::invgamma::rng(unif_rand(), shape_value, scale_value);
  }
  PutRNGstate();

  UNPROTECT(1);
  return out;
}
