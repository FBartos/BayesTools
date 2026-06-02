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

double scalar_real(SEXP x, const char *name)
{
  if(Rf_length(x) != 1){
    Rf_error("'%s' must be a numeric scalar.", name);
  }

  SEXP x_real = PROTECT(coerce_numeric(x, name));
  double value = REAL(x_real)[0];
  UNPROTECT(1);

  return value;
}

int scalar_int(SEXP x, const char *name)
{
  if(Rf_length(x) != 1){
    Rf_error("'%s' must be an integer scalar.", name);
  }

  SEXP x_int = PROTECT(Rf_coerceVector(x, INTSXP));
  int value = INTEGER(x_int)[0];
  UNPROTECT(1);

  return value;
}

bool scalar_bool(SEXP x, const char *name)
{
  if(Rf_length(x) != 1){
    Rf_error("'%s' must be a logical scalar.", name);
  }

  int value = Rf_asLogical(x);
  if(value == NA_LOGICAL){
    Rf_error("'%s' must not be NA.", name);
  }

  return value == TRUE;
}

void check_moment_parameters(double location, double tau, double order)
{
  if(!bayestools::nonlocal::valid_moment_parameters(location, tau, order)){
    Rf_error("'location' must be finite, 'tau' must be positive and finite, and 'order' must be a positive integer.");
  }
}

void check_invmoment_parameters(double location, double tau, double order, double df)
{
  if(!bayestools::nonlocal::valid_invmoment_parameters(location, tau, order, df)){
    Rf_error("'location' must be finite, 'tau' must be positive and finite, 'order' must be a positive integer, and 'df' must be positive and finite.");
  }
}

SEXP nonlocal_vector_4(SEXP x, SEXP location, SEXP tau, SEXP order, SEXP flag,
                       const char *flag_name,
                       double (*fun)(double, double, double, double, bool))
{
  SEXP x_real = PROTECT(coerce_numeric(x, "x"));
  double location_value = scalar_real(location, "location");
  double tau_value = scalar_real(tau, "tau");
  double order_value = scalar_real(order, "order");
  bool flag_value = scalar_bool(flag, flag_name);

  check_moment_parameters(location_value, tau_value, order_value);

  R_xlen_t n = XLENGTH(x_real);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  for(R_xlen_t i = 0; i < n; ++i){
    REAL(out)[i] = fun(REAL(x_real)[i], location_value, tau_value, order_value, flag_value);
  }

  UNPROTECT(2);
  return out;
}

}

extern "C" SEXP BayesTools_moment_d(SEXP x, SEXP location, SEXP tau, SEXP order, SEXP log)
{
  return nonlocal_vector_4(x, location, tau, order, log, "log", bayestools::nonlocal::dmoment);
}

extern "C" SEXP BayesTools_moment_p(SEXP q, SEXP location, SEXP tau, SEXP order, SEXP lower_tail, SEXP log_p)
{
  SEXP q_real = PROTECT(coerce_numeric(q, "q"));
  double location_value = scalar_real(location, "location");
  double tau_value = scalar_real(tau, "tau");
  double order_value = scalar_real(order, "order");
  bool lower_tail_value = scalar_bool(lower_tail, "lower.tail");
  bool log_p_value = scalar_bool(log_p, "log.p");

  check_moment_parameters(location_value, tau_value, order_value);

  R_xlen_t n = XLENGTH(q_real);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  for(R_xlen_t i = 0; i < n; ++i){
    REAL(out)[i] = bayestools::nonlocal::pmoment(
      REAL(q_real)[i], location_value, tau_value, order_value, lower_tail_value, log_p_value
    );
  }

  UNPROTECT(2);
  return out;
}

extern "C" SEXP BayesTools_moment_q(SEXP p, SEXP location, SEXP tau, SEXP order, SEXP lower_tail, SEXP log_p)
{
  SEXP p_real = PROTECT(coerce_numeric(p, "p"));
  double location_value = scalar_real(location, "location");
  double tau_value = scalar_real(tau, "tau");
  double order_value = scalar_real(order, "order");
  bool lower_tail_value = scalar_bool(lower_tail, "lower.tail");
  bool log_p_value = scalar_bool(log_p, "log.p");

  check_moment_parameters(location_value, tau_value, order_value);

  R_xlen_t n = XLENGTH(p_real);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  for(R_xlen_t i = 0; i < n; ++i){
    REAL(out)[i] = bayestools::nonlocal::qmoment(
      REAL(p_real)[i], location_value, tau_value, order_value, lower_tail_value, log_p_value
    );
  }

  UNPROTECT(2);
  return out;
}

extern "C" SEXP BayesTools_moment_r(SEXP n, SEXP location, SEXP tau, SEXP order)
{
  int n_value = scalar_int(n, "n");
  if(n_value < 1){
    Rf_error("'n' must be a positive integer scalar.");
  }

  double location_value = scalar_real(location, "location");
  double tau_value = scalar_real(tau, "tau");
  double order_value = scalar_real(order, "order");
  check_moment_parameters(location_value, tau_value, order_value);

  SEXP out = PROTECT(Rf_allocVector(REALSXP, n_value));
  GetRNGstate();
  for(int i = 0; i < n_value; ++i){
    REAL(out)[i] = bayestools::nonlocal::rmoment(
      unif_rand(), unif_rand(), location_value, tau_value, order_value
    );
  }
  PutRNGstate();

  UNPROTECT(1);
  return out;
}

extern "C" SEXP BayesTools_invmoment_d(SEXP x, SEXP location, SEXP tau, SEXP order, SEXP df, SEXP log)
{
  SEXP x_real = PROTECT(coerce_numeric(x, "x"));
  double location_value = scalar_real(location, "location");
  double tau_value = scalar_real(tau, "tau");
  double order_value = scalar_real(order, "order");
  double df_value = scalar_real(df, "df");
  bool log_value = scalar_bool(log, "log");

  check_invmoment_parameters(location_value, tau_value, order_value, df_value);

  R_xlen_t n = XLENGTH(x_real);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  for(R_xlen_t i = 0; i < n; ++i){
    REAL(out)[i] = bayestools::nonlocal::dinvmoment(
      REAL(x_real)[i], location_value, tau_value, order_value, df_value, log_value
    );
  }

  UNPROTECT(2);
  return out;
}

extern "C" SEXP BayesTools_invmoment_p(SEXP q, SEXP location, SEXP tau, SEXP order, SEXP df, SEXP lower_tail, SEXP log_p)
{
  SEXP q_real = PROTECT(coerce_numeric(q, "q"));
  double location_value = scalar_real(location, "location");
  double tau_value = scalar_real(tau, "tau");
  double order_value = scalar_real(order, "order");
  double df_value = scalar_real(df, "df");
  bool lower_tail_value = scalar_bool(lower_tail, "lower.tail");
  bool log_p_value = scalar_bool(log_p, "log.p");

  check_invmoment_parameters(location_value, tau_value, order_value, df_value);

  R_xlen_t n = XLENGTH(q_real);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  for(R_xlen_t i = 0; i < n; ++i){
    REAL(out)[i] = bayestools::nonlocal::pinvmoment(
      REAL(q_real)[i], location_value, tau_value, order_value, df_value, lower_tail_value, log_p_value
    );
  }

  UNPROTECT(2);
  return out;
}

extern "C" SEXP BayesTools_invmoment_q(SEXP p, SEXP location, SEXP tau, SEXP order, SEXP df, SEXP lower_tail, SEXP log_p)
{
  SEXP p_real = PROTECT(coerce_numeric(p, "p"));
  double location_value = scalar_real(location, "location");
  double tau_value = scalar_real(tau, "tau");
  double order_value = scalar_real(order, "order");
  double df_value = scalar_real(df, "df");
  bool lower_tail_value = scalar_bool(lower_tail, "lower.tail");
  bool log_p_value = scalar_bool(log_p, "log.p");

  check_invmoment_parameters(location_value, tau_value, order_value, df_value);

  R_xlen_t n = XLENGTH(p_real);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  for(R_xlen_t i = 0; i < n; ++i){
    REAL(out)[i] = bayestools::nonlocal::qinvmoment(
      REAL(p_real)[i], location_value, tau_value, order_value, df_value, lower_tail_value, log_p_value
    );
  }

  UNPROTECT(2);
  return out;
}

extern "C" SEXP BayesTools_invmoment_r(SEXP n, SEXP location, SEXP tau, SEXP order, SEXP df)
{
  int n_value = scalar_int(n, "n");
  if(n_value < 1){
    Rf_error("'n' must be a positive integer scalar.");
  }

  double location_value = scalar_real(location, "location");
  double tau_value = scalar_real(tau, "tau");
  double order_value = scalar_real(order, "order");
  double df_value = scalar_real(df, "df");
  check_invmoment_parameters(location_value, tau_value, order_value, df_value);

  SEXP out = PROTECT(Rf_allocVector(REALSXP, n_value));
  GetRNGstate();
  for(int i = 0; i < n_value; ++i){
    REAL(out)[i] = bayestools::nonlocal::rinvmoment(
      unif_rand(), unif_rand(), location_value, tau_value, order_value, df_value
    );
  }
  PutRNGstate();

  UNPROTECT(1);
  return out;
}
