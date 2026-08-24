#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP BayesTools_lkj_cholesky_from_u(SEXP u, SEXP K);
extern SEXP BayesTools_lkj_log_prior_u(SEXP u, SEXP alpha);
extern SEXP BayesTools_lkj_alpha(SEXP K, SEXP eta);
extern SEXP BayesTools_moment_d(SEXP x, SEXP location, SEXP tau, SEXP order, SEXP log);
extern SEXP BayesTools_moment_p(SEXP q, SEXP location, SEXP tau, SEXP order, SEXP lower_tail, SEXP log_p);
extern SEXP BayesTools_moment_q(SEXP p, SEXP location, SEXP tau, SEXP order, SEXP lower_tail, SEXP log_p);
extern SEXP BayesTools_moment_r(SEXP n, SEXP location, SEXP tau, SEXP order);
extern SEXP BayesTools_invmoment_d(SEXP x, SEXP location, SEXP tau, SEXP order, SEXP df, SEXP log);
extern SEXP BayesTools_invmoment_p(SEXP q, SEXP location, SEXP tau, SEXP order, SEXP df, SEXP lower_tail, SEXP log_p);
extern SEXP BayesTools_invmoment_q(SEXP p, SEXP location, SEXP tau, SEXP order, SEXP df, SEXP lower_tail, SEXP log_p);
extern SEXP BayesTools_invmoment_r(SEXP n, SEXP location, SEXP tau, SEXP order, SEXP df);
extern SEXP BayesTools_invgamma_d(SEXP x, SEXP shape, SEXP scale, SEXP log);
extern SEXP BayesTools_invgamma_p(SEXP q, SEXP shape, SEXP scale, SEXP lower_tail, SEXP log_p);
extern SEXP BayesTools_invgamma_q(SEXP p, SEXP shape, SEXP scale, SEXP lower_tail, SEXP log_p);
extern SEXP BayesTools_invgamma_r(SEXP n, SEXP shape, SEXP scale);
extern SEXP BayesTools_structured_cholesky(SEXP rho, SEXP coordinates, SEXP equicorrelation);

static const R_CallMethodDef callMethods[] = {
  {"BayesTools_lkj_cholesky_from_u", (DL_FUNC) &BayesTools_lkj_cholesky_from_u, 2},
  {"BayesTools_lkj_log_prior_u",     (DL_FUNC) &BayesTools_lkj_log_prior_u,     2},
  {"BayesTools_lkj_alpha",           (DL_FUNC) &BayesTools_lkj_alpha,           2},
  {"BayesTools_moment_d",            (DL_FUNC) &BayesTools_moment_d,            5},
  {"BayesTools_moment_p",            (DL_FUNC) &BayesTools_moment_p,            6},
  {"BayesTools_moment_q",            (DL_FUNC) &BayesTools_moment_q,            6},
  {"BayesTools_moment_r",            (DL_FUNC) &BayesTools_moment_r,            4},
  {"BayesTools_invmoment_d",         (DL_FUNC) &BayesTools_invmoment_d,         6},
  {"BayesTools_invmoment_p",         (DL_FUNC) &BayesTools_invmoment_p,         7},
  {"BayesTools_invmoment_q",         (DL_FUNC) &BayesTools_invmoment_q,         7},
  {"BayesTools_invmoment_r",         (DL_FUNC) &BayesTools_invmoment_r,         5},
  {"BayesTools_invgamma_d",          (DL_FUNC) &BayesTools_invgamma_d,          4},
  {"BayesTools_invgamma_p",          (DL_FUNC) &BayesTools_invgamma_p,          5},
  {"BayesTools_invgamma_q",          (DL_FUNC) &BayesTools_invgamma_q,          5},
  {"BayesTools_invgamma_r",          (DL_FUNC) &BayesTools_invgamma_r,          3},
  {"BayesTools_structured_cholesky", (DL_FUNC) &BayesTools_structured_cholesky, 3},
  {NULL, NULL, 0}
};

void R_init_BayesTools(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
