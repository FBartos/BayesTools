#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP BayesTools_lkj_cholesky_from_u(SEXP u, SEXP K);
extern SEXP BayesTools_lkj_corr_from_u(SEXP u, SEXP K);
extern SEXP BayesTools_lkj_log_prior_u(SEXP u, SEXP alpha);
extern SEXP BayesTools_lkj_alpha(SEXP K, SEXP eta);

static const R_CallMethodDef callMethods[] = {
  {"BayesTools_lkj_cholesky_from_u", (DL_FUNC) &BayesTools_lkj_cholesky_from_u, 2},
  {"BayesTools_lkj_corr_from_u",     (DL_FUNC) &BayesTools_lkj_corr_from_u,     2},
  {"BayesTools_lkj_log_prior_u",     (DL_FUNC) &BayesTools_lkj_log_prior_u,     2},
  {"BayesTools_lkj_alpha",           (DL_FUNC) &BayesTools_lkj_alpha,           2},
  {NULL, NULL, 0}
};

void R_init_BayesTools(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

