#include "nonlocal/BTNonlocalCore.h"
#include "gamma/BTGammaRange.h"

#include <cmath>
#include <cstring>

#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Random.h>

// Internal proving facts use the identical checker that enables certificates.
// Only the ten meaningful x87 bytes are copied; padding is never exposed.
extern "C" SEXP BayesTools_native_range_environment()
{
  bayestools::gamma::range::environment env = bayestools::gamma::range::current_environment();
  const char *names[] = {"compiled_profile", "documented_cpu", "hypervisor", "control_word", "mxcsr",
    "ln2_raw", "ln2_lower_raw", "ln2_upper_raw", "log2_raw", "threshold_raw", "available", "representation"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 12));
  SEXP labels = PROTECT(Rf_allocVector(STRSXP, 12));
  for(int i = 0; i < 12; ++i) SET_STRING_ELT(labels, i, Rf_mkChar(names[i]));
  SET_VECTOR_ELT(out, 0, Rf_ScalarLogical(env.compiled_profile));
  SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(env.documented_cpu));
  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(env.hypervisor));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(env.compiled_profile ? env.control_word : NA_INTEGER));
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(env.compiled_profile ? env.mxcsr : NA_INTEGER));
  const long double extended[] = {env.ln2, 0xb.17217f7d1cf79a4p-4L, 0xb.17217f7d1cf79b4p-4L};
  const double doubles[] = {std::log(2.0), bayestools::gamma::range::far_threshold};
  for(int i = 0; i < 5; ++i){
    int size = env.compiled_profile ? (i < 3 ? 10 : 8) : 0;
    SEXP bytes = PROTECT(Rf_allocVector(RAWSXP, size));
    if(size) std::memcpy(RAW(bytes), i < 3 ? static_cast<const void*>(&extended[i]) :
      static_cast<const void*>(&doubles[i - 3]), size);
    SET_VECTOR_ELT(out, i + 5, bytes);
    UNPROTECT(1);
  }
  SET_VECTOR_ELT(out, 10, Rf_ScalarLogical(env.available));
  SET_VECTOR_ELT(out, 11, Rf_mkString(env.compiled_profile ? "x87-extended-little-endian-10" : "unsupported"));
  Rf_setAttrib(out, R_NamesSymbol, labels);
  UNPROTECT(2);
  return out;
}

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
    if(ISNAN(x_ptr[i])){
      out_ptr[i] = x_ptr[i];
      continue;
    }
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
    if(ISNAN(q_ptr[i])){
      out_ptr[i] = q_ptr[i];
      continue;
    }
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
    if(ISNAN(p_ptr[i])){
      out_ptr[i] = p_ptr[i];
      continue;
    }
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
    // Fix the draw order while retaining the existing GCC seeded stream.
    double u_size = unif_rand();
    double u_sign = unif_rand();
    out_ptr[i] = invmoment ?
      bayestools::nonlocal::invmoment_rng(
        u_sign, u_size, location_value, tau_value, order_value, df_value
      ) :
      bayestools::nonlocal::moment_rng(
        u_sign, u_size, location_value, tau_value, order_value
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

extern "C" SEXP BayesTools_nonlocal_log_interval_mass(SEXP lower, SEXP upper,
  SEXP location, SEXP tau, SEXP order, SEXP df, SEXP inverse)
{
  SEXP l = PROTECT(coerce_numeric(lower, "lower"));
  SEXP u = PROTECT(coerce_numeric(upper, "upper"));
  if(XLENGTH(l) != XLENGTH(u)) Rf_error("'lower' and 'upper' must have equal lengths.");
  double location_value = scalar_numeric(location, "location");
  double tau_value = scalar_numeric(tau, "tau");
  double order_value = scalar_numeric(order, "order");
  double df_value = scalar_numeric(df, "df");
  bool inverse_value = scalar_bool(inverse, "invmoment");
  SEXP out = PROTECT(Rf_allocVector(REALSXP, XLENGTH(l)));
  for(R_xlen_t i = 0; i < XLENGTH(l); ++i){
    REAL(out)[i] = ISNAN(REAL(l)[i]) ? REAL(l)[i] : ISNAN(REAL(u)[i]) ? REAL(u)[i] :
      bayestools::nonlocal::log_interval_mass(REAL(l)[i], REAL(u)[i],
        location_value, tau_value, order_value, df_value, inverse_value);
  }
  UNPROTECT(3);
  return out;
}

extern "C" SEXP BayesTools_nonlocal_truncated_quantile(SEXP p, SEXP lower,
  SEXP upper, SEXP location, SEXP tau, SEXP order, SEXP df, SEXP inverse)
{
  SEXP p_real = PROTECT(coerce_numeric(p, "p"));
  double lower_value = scalar_numeric(lower, "lower");
  double upper_value = scalar_numeric(upper, "upper");
  double location_value = scalar_numeric(location, "location");
  double tau_value = scalar_numeric(tau, "tau");
  double order_value = scalar_numeric(order, "order");
  double df_value = scalar_numeric(df, "df");
  bool inverse_value = scalar_bool(inverse, "invmoment");
  double log_sign_masses[2];
  double log_mass = bayestools::nonlocal::log_interval_mass(lower_value,
    upper_value, location_value, tau_value, order_value, df_value, inverse_value,
    log_sign_masses);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, XLENGTH(p_real)));
  for(R_xlen_t i = 0; i < XLENGTH(p_real); ++i){
    REAL(out)[i] = ISNAN(REAL(p_real)[i]) ? REAL(p_real)[i] :
      bayestools::nonlocal::truncated_quantile(REAL(p_real)[i], lower_value,
        upper_value, location_value, tau_value, order_value, df_value, inverse_value,
        log_sign_masses);
  }
  SEXP mass = PROTECT(Rf_ScalarReal(log_mass));
  Rf_setAttrib(out, Rf_install("log_normalizer"), mass);
  UNPROTECT(3);
  return out;
}
