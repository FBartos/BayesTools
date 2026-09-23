#include <algorithm>
#include <cmath>
#include <cstddef>

#include <Rinternals.h>
#include <R_ext/Error.h>

namespace {

SEXP coerce_numeric(SEXP x, const char *name)
{
  if(TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP){
    Rf_error("'%s' must be numeric.", name);
  }

  return Rf_coerceVector(x, REALSXP);
}

bool scalar_logical(SEXP x, const char *name)
{
  if(Rf_length(x) != 1){
    Rf_error("'%s' must be one logical value.", name);
  }

  const int value = Rf_asLogical(x);
  if(value == NA_LOGICAL){
    Rf_error("'%s' must not be missing.", name);
  }

  return value != 0;
}

inline std::size_t array_index(std::size_t draw,
                               std::size_t row,
                               std::size_t column,
                               std::size_t n_draws,
                               std::size_t n_columns)
{
  return draw + n_draws * row + n_draws * n_columns * column;
}

void fill_cs(double *out,
             const double *rho,
             std::size_t n_draws,
             std::size_t n_columns)
{
  for(std::size_t draw = 0; draw < n_draws; ++draw){
    const double value = rho[draw];
    out[array_index(draw, 0, 0, n_draws, n_columns)] = 1.0;

    for(std::size_t row = 1; row < n_columns; ++row){
      out[array_index(draw, row, 0, n_draws, n_columns)] = value;
      for(std::size_t column = 1; column < row; ++column){
        const double numerator = 1.0 - value;
        const double denominator =
          (1.0 + static_cast<double>(column - 1) * value) *
          (1.0 + static_cast<double>(column) * value);
        out[array_index(draw, row, column, n_draws, n_columns)] =
          value * std::sqrt(numerator / denominator);
      }
      const double diagonal = std::sqrt(
        (1.0 - value) *
        (1.0 + static_cast<double>(row) * value) /
        (1.0 + static_cast<double>(row - 1) * value)
      );
      if(!std::isfinite(diagonal) || diagonal <= 0.0){
        Rf_error("Structured CS Cholesky factor has a non-positive or non-finite diagonal.");
      }
      out[array_index(draw, row, row, n_draws, n_columns)] = diagonal;
    }
  }
}

void fill_markov(double *out,
                 const double *rho,
                 const double *coordinates,
                 std::size_t n_draws,
                 std::size_t n_columns)
{
  for(std::size_t draw = 0; draw < n_draws; ++draw){
    const double value = rho[draw];
    out[array_index(draw, 0, 0, n_draws, n_columns)] = 1.0;

    for(std::size_t row = 1; row < n_columns; ++row){
      const double gap = coordinates[row] - coordinates[row - 1];
      if(!std::isfinite(gap) || gap <= 0.0){
        Rf_error("Structured Markov Cholesky coordinates must have positive finite gaps.");
      }

      double log_phi;
      double phi;
      double innovation_variance;
      if(value == 0.0){
        log_phi = R_NegInf;
        phi = 0.0;
        innovation_variance = 1.0;
      }else{
        if(value < 0.0 && gap != std::floor(gap)){
          Rf_error("A negative structured Markov correlation requires integer coordinate gaps.");
        }
        log_phi = gap * std::log(std::fabs(value));
        phi = std::exp(log_phi);
        if(value < 0.0 && std::fmod(gap, 2.0) != 0.0){
          phi = -phi;
        }
        innovation_variance = -std::expm1(2.0 * log_phi);
      }
      if(!std::isfinite(innovation_variance) ||
         innovation_variance <= 0.0){
        Rf_error("Structured Markov Cholesky factor has a non-positive or non-finite innovation variance.");
      }

      for(std::size_t column = 0; column < n_columns; ++column){
        out[array_index(draw, row, column, n_draws, n_columns)] =
          phi * out[array_index(draw, row - 1, column, n_draws, n_columns)];
      }
      out[array_index(draw, row, row, n_draws, n_columns)] =
        std::sqrt(innovation_variance);
    }
  }
}

}

extern "C" SEXP BayesTools_structured_cholesky(SEXP rho,
                                                SEXP coordinates,
                                                SEXP equicorrelation)
{
  SEXP rho_real = PROTECT(coerce_numeric(rho, "rho"));
  SEXP coordinates_real = PROTECT(coerce_numeric(coordinates, "coordinates"));
  const bool is_cs = scalar_logical(equicorrelation, "equicorrelation");
  const std::size_t n_draws = static_cast<std::size_t>(Rf_length(rho_real));
  const std::size_t n_columns =
    static_cast<std::size_t>(Rf_length(coordinates_real));

  if(n_draws < 1){
    Rf_error("'rho' must contain at least one draw.");
  }
  if(n_columns < 1){
    Rf_error("'coordinates' must contain at least one value.");
  }
  const double *rho_ptr = REAL(rho_real);
  const double *coordinates_ptr = REAL(coordinates_real);
  for(std::size_t draw = 0; draw < n_draws; ++draw){
    if(!std::isfinite(rho_ptr[draw])){
      Rf_error("'rho' must contain only finite values.");
    }
  }
  for(std::size_t column = 0; column < n_columns; ++column){
    if(!std::isfinite(coordinates_ptr[column])){
      Rf_error("'coordinates' must contain only finite values.");
    }
  }

  SEXP out = PROTECT(Rf_allocVector(
    REALSXP,
    static_cast<R_xlen_t>(n_draws * n_columns * n_columns)
  ));
  double *out_ptr = REAL(out);
  std::fill(out_ptr, out_ptr + n_draws * n_columns * n_columns, 0.0);
  if(is_cs){
    fill_cs(out_ptr, rho_ptr, n_draws, n_columns);
  }else{
    fill_markov(
      out_ptr,
      rho_ptr,
      coordinates_ptr,
      n_draws,
      n_columns
    );
  }

  SEXP dim = PROTECT(Rf_allocVector(INTSXP, 3));
  INTEGER(dim)[0] = static_cast<int>(n_draws);
  INTEGER(dim)[1] = static_cast<int>(n_columns);
  INTEGER(dim)[2] = static_cast<int>(n_columns);
  Rf_setAttrib(out, R_DimSymbol, dim);

  UNPROTECT(4);
  return out;
}
