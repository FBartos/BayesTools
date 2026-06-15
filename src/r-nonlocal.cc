#include "nonlocal/BTNonlocalCore.h"

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

SEXP nonlocal_d(SEXP x, SEXP location, SEXP tau, SEXP order, SEXP log,
                bool invmoment, SEXP df = R_NilValue)
{
  SEXP x_real = PROTECT(coerce_numeric(x, "x"));
  double location_value = scalar_numeric(location, "location");
  double tau_value = scalar_numeric(tau, "tau");
  double order_value = scalar_numeric(order, "order");
  double df_value = invmoment ? scalar_numeric(df, "df") : 0.0;
  bool log_value = scalar_bool(log, "log");

  R_xlen_t n = XLENGTH(x_real);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *out_ptr = REAL(out);
  double const *x_ptr = REAL(x_real);
  for(R_xlen_t i = 0; i < n; ++i){
    double log_density = invmoment ?
      bayestools::nonlocal::invmoment_log_density(
        x_ptr[i], location_value, tau_value, order_value, df_value
      ) :
      bayestools::nonlocal::moment_log_density(
        x_ptr[i], location_value, tau_value, order_value
      );
    out_ptr[i] = log_value ? log_density : std::exp(log_density);
  }

  UNPROTECT(2);
  return out;
}

SEXP nonlocal_p(SEXP q, SEXP location, SEXP tau, SEXP order,
                SEXP lower_tail, SEXP log_p, bool invmoment,
                SEXP df = R_NilValue)
{
  SEXP q_real = PROTECT(coerce_numeric(q, "q"));
  double location_value = scalar_numeric(location, "location");
  double tau_value = scalar_numeric(tau, "tau");
  double order_value = scalar_numeric(order, "order");
  double df_value = invmoment ? scalar_numeric(df, "df") : 0.0;
  bool lower_tail_value = scalar_bool(lower_tail, "lower.tail");
  bool log_p_value = scalar_bool(log_p, "log.p");

  R_xlen_t n = XLENGTH(q_real);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *out_ptr = REAL(out);
  double const *q_ptr = REAL(q_real);
  for(R_xlen_t i = 0; i < n; ++i){
    out_ptr[i] = invmoment ?
      bayestools::nonlocal::invmoment_cdf(
        q_ptr[i], location_value, tau_value, order_value, df_value,
        lower_tail_value, log_p_value
      ) :
      bayestools::nonlocal::moment_cdf(
        q_ptr[i], location_value, tau_value, order_value,
        lower_tail_value, log_p_value
      );
  }

  UNPROTECT(2);
  return out;
}

SEXP nonlocal_q(SEXP p, SEXP location, SEXP tau, SEXP order,
                SEXP lower_tail, SEXP log_p, bool invmoment,
                SEXP df = R_NilValue)
{
  SEXP p_real = PROTECT(coerce_numeric(p, "p"));
  double location_value = scalar_numeric(location, "location");
  double tau_value = scalar_numeric(tau, "tau");
  double order_value = scalar_numeric(order, "order");
  double df_value = invmoment ? scalar_numeric(df, "df") : 0.0;
  bool lower_tail_value = scalar_bool(lower_tail, "lower.tail");
  bool log_p_value = scalar_bool(log_p, "log.p");

  R_xlen_t n = XLENGTH(p_real);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *out_ptr = REAL(out);
  double const *p_ptr = REAL(p_real);
  for(R_xlen_t i = 0; i < n; ++i){
    out_ptr[i] = invmoment ?
      bayestools::nonlocal::invmoment_quantile(
        p_ptr[i], location_value, tau_value, order_value, df_value,
        lower_tail_value, log_p_value
      ) :
      bayestools::nonlocal::moment_quantile(
        p_ptr[i], location_value, tau_value, order_value,
        lower_tail_value, log_p_value
      );
  }

  UNPROTECT(2);
  return out;
}

SEXP nonlocal_r(SEXP n, SEXP location, SEXP tau, SEXP order,
                bool invmoment, SEXP df = R_NilValue)
{
  int n_value = scalar_int(n, "n");
  double location_value = scalar_numeric(location, "location");
  double tau_value = scalar_numeric(tau, "tau");
  double order_value = scalar_numeric(order, "order");
  double df_value = invmoment ? scalar_numeric(df, "df") : 0.0;

  SEXP out = PROTECT(Rf_allocVector(REALSXP, n_value));
  double *out_ptr = REAL(out);
  GetRNGstate();
  for(int i = 0; i < n_value; ++i){
    out_ptr[i] = invmoment ?
      bayestools::nonlocal::invmoment_rng(
        unif_rand(), unif_rand(), location_value, tau_value, order_value, df_value
      ) :
      bayestools::nonlocal::moment_rng(
        unif_rand(), unif_rand(), location_value, tau_value, order_value
      );
  }
  PutRNGstate();

  UNPROTECT(1);
  return out;
}

}

extern "C" SEXP BayesTools_moment_d(SEXP x, SEXP location, SEXP tau,
                                    SEXP order, SEXP log)
{
  return nonlocal_d(x, location, tau, order, log, false);
}

extern "C" SEXP BayesTools_moment_p(SEXP q, SEXP location, SEXP tau,
                                    SEXP order, SEXP lower_tail, SEXP log_p)
{
  return nonlocal_p(q, location, tau, order, lower_tail, log_p, false);
}

extern "C" SEXP BayesTools_moment_q(SEXP p, SEXP location, SEXP tau,
                                    SEXP order, SEXP lower_tail, SEXP log_p)
{
  return nonlocal_q(p, location, tau, order, lower_tail, log_p, false);
}

extern "C" SEXP BayesTools_moment_r(SEXP n, SEXP location, SEXP tau,
                                    SEXP order)
{
  return nonlocal_r(n, location, tau, order, false);
}

extern "C" SEXP BayesTools_invmoment_d(SEXP x, SEXP location, SEXP tau,
                                       SEXP order, SEXP df, SEXP log)
{
  return nonlocal_d(x, location, tau, order, log, true, df);
}

extern "C" SEXP BayesTools_invmoment_p(SEXP q, SEXP location, SEXP tau,
                                       SEXP order, SEXP df,
                                       SEXP lower_tail, SEXP log_p)
{
  return nonlocal_p(q, location, tau, order, lower_tail, log_p, true, df);
}

extern "C" SEXP BayesTools_invmoment_q(SEXP p, SEXP location, SEXP tau,
                                       SEXP order, SEXP df,
                                       SEXP lower_tail, SEXP log_p)
{
  return nonlocal_q(p, location, tau, order, lower_tail, log_p, true, df);
}

extern "C" SEXP BayesTools_invmoment_r(SEXP n, SEXP location, SEXP tau,
                                       SEXP order, SEXP df)
{
  return nonlocal_r(n, location, tau, order, true, df);
}
