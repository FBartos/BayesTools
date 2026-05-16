#include "lkj/BTLKJCore.h"

#include <cmath>
#include <vector>

#include <Rinternals.h>
#include <R_ext/Error.h>

namespace {

unsigned int scalar_K(SEXP K)
{
  if(Rf_length(K) != 1){
    Rf_error("'K' must be an integer scalar.");
  }

  SEXP K_int = PROTECT(Rf_coerceVector(K, INTSXP));
  int value = INTEGER(K_int)[0];
  UNPROTECT(1);

  if(value < 1){
    Rf_error("'K' must be a positive integer scalar.");
  }

  return static_cast<unsigned int>(value);
}

double scalar_eta(SEXP eta)
{
  if(Rf_length(eta) != 1){
    Rf_error("'eta' must be a numeric scalar.");
  }

  SEXP eta_real = PROTECT(Rf_coerceVector(eta, REALSXP));
  double value = REAL(eta_real)[0];
  UNPROTECT(1);

  if(!std::isfinite(value) || value <= 0.0){
    Rf_error("'eta' must be positive and finite.");
  }

  return value;
}

bool is_matrix(SEXP x)
{
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  return TYPEOF(dim) == INTSXP && Rf_length(dim) == 2;
}

void check_u_values(double const *u, unsigned int length)
{
  if(!bayestools::lkj::valid_u(u, length)){
    Rf_error("'u' values must be finite and strictly between 0 and 1.");
  }
}

void check_u_values_strided(double const *u, unsigned int length, unsigned int stride)
{
  if(!bayestools::lkj::valid_u_strided(u, length, stride)){
    Rf_error("'u' values must be finite and strictly between 0 and 1.");
  }
}

SEXP coerce_numeric(SEXP x, const char *name)
{
  if(TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP){
    Rf_error("'%s' must be numeric.", name);
  }

  return Rf_coerceVector(x, REALSXP);
}

SEXP lkj_matrix_from_u(SEXP u, SEXP K, bool correlation)
{
  const unsigned int K_value = scalar_K(K);
  const unsigned int n_pairs = bayestools::lkj::n_cpc(K_value);

  SEXP u_real = PROTECT(coerce_numeric(u, "u"));
  const bool matrix_input = is_matrix(u);
  unsigned int n_draws = 1;
  unsigned int u_width = static_cast<unsigned int>(Rf_length(u_real));

  if(matrix_input){
    SEXP dim = Rf_getAttrib(u, R_DimSymbol);
    n_draws = static_cast<unsigned int>(INTEGER(dim)[0]);
    u_width = static_cast<unsigned int>(INTEGER(dim)[1]);
  }

  if(u_width != n_pairs){
    Rf_error("'u' must have K * (K - 1) / 2 columns or elements.");
  }

  double const *u_ptr = REAL(u_real);
  std::vector<double> row_major(K_value * K_value);
  std::vector<double> corr_work;
  if(correlation){
    corr_work.resize(K_value * K_value);
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, n_draws * K_value * K_value));
  double *out_ptr = REAL(out);

  for(unsigned int draw = 0; draw < n_draws; ++draw){
    if(matrix_input){
      double const *draw_u = u_ptr + draw;
      check_u_values_strided(draw_u, n_pairs, n_draws);
      if(correlation){
        bayestools::lkj::fill_cholesky_from_u_strided(&corr_work[0], n_pairs == 0 ? 0 : draw_u, K_value, n_draws);
        bayestools::lkj::fill_corr_from_cholesky(&row_major[0], &corr_work[0], K_value);
      }else{
        bayestools::lkj::fill_cholesky_from_u_strided(&row_major[0], n_pairs == 0 ? 0 : draw_u, K_value, n_draws);
      }
    }else{
      check_u_values(u_ptr, n_pairs);
      if(correlation){
        bayestools::lkj::fill_corr_from_u(&row_major[0], n_pairs == 0 ? 0 : u_ptr, K_value);
      }else{
        bayestools::lkj::fill_cholesky_from_u(&row_major[0], n_pairs == 0 ? 0 : u_ptr, K_value);
      }
    }

    for(unsigned int row = 0; row < K_value; ++row){
      for(unsigned int column = 0; column < K_value; ++column){
        double cell = row_major[bayestools::lkj::flat_index(row + 1, column + 1, K_value)];
        if(matrix_input){
          out_ptr[draw + n_draws * row + n_draws * K_value * column] = cell;
        }else{
          out_ptr[row + K_value * column] = cell;
        }
      }
    }

    if(!matrix_input){
      break;
    }
  }

  if(matrix_input){
    SEXP dim = PROTECT(Rf_allocVector(INTSXP, 3));
    INTEGER(dim)[0] = static_cast<int>(n_draws);
    INTEGER(dim)[1] = static_cast<int>(K_value);
    INTEGER(dim)[2] = static_cast<int>(K_value);
    Rf_setAttrib(out, R_DimSymbol, dim);
    UNPROTECT(1);
  }else{
    SEXP dim = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(dim)[0] = static_cast<int>(K_value);
    INTEGER(dim)[1] = static_cast<int>(K_value);
    Rf_setAttrib(out, R_DimSymbol, dim);
    UNPROTECT(1);
  }

  UNPROTECT(2);
  return out;
}

}

extern "C" SEXP BayesTools_lkj_cholesky_from_u(SEXP u, SEXP K)
{
  return lkj_matrix_from_u(u, K, false);
}

extern "C" SEXP BayesTools_lkj_corr_from_u(SEXP u, SEXP K)
{
  return lkj_matrix_from_u(u, K, true);
}

extern "C" SEXP BayesTools_lkj_log_prior_u(SEXP u, SEXP alpha)
{
  SEXP u_real = PROTECT(coerce_numeric(u, "u"));
  SEXP alpha_real = PROTECT(coerce_numeric(alpha, "alpha"));
  const unsigned int n_pairs = static_cast<unsigned int>(Rf_length(alpha_real));

  if(!bayestools::lkj::valid_alpha(REAL(alpha_real), n_pairs)){
    Rf_error("'alpha' values must be positive and finite.");
  }

  const bool matrix_input = is_matrix(u);
  unsigned int n_draws = 1;
  unsigned int u_width = static_cast<unsigned int>(Rf_length(u_real));
  if(matrix_input){
    SEXP dim = Rf_getAttrib(u, R_DimSymbol);
    n_draws = static_cast<unsigned int>(INTEGER(dim)[0]);
    u_width = static_cast<unsigned int>(INTEGER(dim)[1]);
  }
  if(u_width != n_pairs){
    Rf_error("'u' must have length or column count matching 'alpha'.");
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, n_draws));
  double *out_ptr = REAL(out);
  double const *u_ptr = REAL(u_real);
  double const *alpha_ptr = REAL(alpha_real);

  if(matrix_input){
    for(unsigned int draw = 0; draw < n_draws; ++draw){
      out_ptr[draw] = bayestools::lkj::log_density_u_strided_valid_alpha(
        n_pairs == 0 ? 0 : u_ptr + draw,
        n_pairs == 0 ? 0 : alpha_ptr,
        n_pairs,
        n_draws
      );
    }
  }else{
    out_ptr[0] = bayestools::lkj::log_density_u(
      n_pairs == 0 ? 0 : u_ptr,
      n_pairs == 0 ? 0 : alpha_ptr,
      n_pairs
    );
  }

  UNPROTECT(3);
  return out;
}

extern "C" SEXP BayesTools_lkj_alpha(SEXP K, SEXP eta)
{
  const unsigned int K_value = scalar_K(K);
  const double eta_value = scalar_eta(eta);
  const unsigned int n_pairs = bayestools::lkj::n_cpc(K_value);

  SEXP out = PROTECT(Rf_allocVector(REALSXP, n_pairs));
  bayestools::lkj::fill_alpha(n_pairs == 0 ? 0 : REAL(out), K_value, eta_value);
  UNPROTECT(1);

  return out;
}
