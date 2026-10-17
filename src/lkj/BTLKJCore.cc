#include "BTLKJCore.h"

#include <algorithm>
#include <cmath>
#include <JRmath.h>
#include <limits>
#include <vector>

namespace bayestools {
namespace lkj {

unsigned int n_cpc(unsigned int K)
{
  return K * (K - 1) / 2;
}

unsigned int cpc_index(unsigned int i, unsigned int j)
{
  return ((j - 1) * (j - 2)) / 2 + (i - 1);
}

unsigned int flat_index(unsigned int row, unsigned int column, unsigned int K)
{
  return (row - 1) * K + (column - 1);
}

bool valid_u(double const *u, unsigned int length)
{
  return valid_u_strided(u, length, 1);
}

bool valid_u_strided(double const *u, unsigned int length, unsigned int stride)
{
  for(unsigned int i = 0; i < length; ++i){
    double value = u[i * stride];
    if(!std::isfinite(value) || value <= 0.0 || value >= 1.0){
      return false;
    }
  }

  return true;
}

bool valid_alpha(double const *alpha, unsigned int length)
{
  for(unsigned int i = 0; i < length; ++i){
    if(!std::isfinite(alpha[i]) || alpha[i] <= 0.0){
      return false;
    }
  }

  return true;
}

namespace {

// Preserve the LKJ transform on [0, 1] and keep descendant evaluation
// finite while JAGS rejects stochastic proposals outside that support.
double clamp_u(double u)
{
  return std::max(0.0, std::min(1.0, u));
}

// The innovation factor sqrt(1 - cpc^2) of a primitive u in [0, 1], with
// cpc = 2u - 1, equals 2 sqrt(u) sqrt(1 - u). It is computed from u itself:
// the rounded cpc loses it for a small interior u (cpc rounds to -1 and the
// factor to 0, although it is of order sqrt(u) and the prior density of such
// a u is finite), and 1 - u is exact for u >= 1/2.
double innovation_from_u(double u)
{
  u = clamp_u(u);
  return 2.0 * std::sqrt(u) * std::sqrt(1.0 - u);
}

}

double cpc_from_u(double u)
{
  return 2.0 * clamp_u(u) - 1.0;
}

void fill_cholesky_from_u(double *value, double const *u, unsigned int K)
{
  fill_cholesky_from_u_strided(value, u, K, 1);
}

void fill_cholesky_from_u_strided(double *value, double const *u, unsigned int K, unsigned int stride)
{
  std::fill(value, value + K * K, 0.0);
  value[flat_index(1, 1, K)] = 1.0;

  for(unsigned int row = 2; row <= K; ++row){
    double multiplier = 1.0;
    for(unsigned int column = 1; column < row; ++column){
      double u_value = u[cpc_index(column, row) * stride];
      value[flat_index(row, column, K)] = cpc_from_u(u_value) * multiplier;
      multiplier *= innovation_from_u(u_value);
    }
    value[flat_index(row, row, K)] = multiplier;
  }
}

// L is the lower Cholesky factor of a correlation matrix, as built by
// fill_cholesky_from_u(): its rows have unit norm by construction (for every
// primitive, including the clamped extension), so the diagonal is exactly 1.
// The rounded row sum of squares would scatter within 1 +/- a few eps and make
// a constant diagonal monitor look like varying draws.
void fill_corr_from_cholesky(double *value, double const *L, unsigned int K)
{
  for(unsigned int row = 1; row <= K; ++row){
    for(unsigned int column = 1; column < row; ++column){
      double cell = 0.0;
      for(unsigned int m = 1; m <= column; ++m){
        cell += L[flat_index(row, m, K)] * L[flat_index(column, m, K)];
      }
      value[flat_index(row, column, K)] = cell;
      value[flat_index(column, row, K)] = cell;
    }
    value[flat_index(row, row, K)] = 1.0;
  }
}

void fill_corr_from_u(double *value, double const *u, unsigned int K)
{
  std::vector<double> L(K * K);
  fill_cholesky_from_u(&L[0], u, K);
  fill_corr_from_cholesky(value, &L[0], K);
}

double log_density_u_strided_valid_alpha(double const *u, double const *alpha,
                                         unsigned int length, unsigned int stride)
{
  if(!valid_u_strided(u, length, stride)){
    return -std::numeric_limits<double>::infinity();
  }

  double log_density = 0.0;
  for(unsigned int i = 0; i < length; ++i){
    log_density += dbeta(u[i * stride], alpha[i], alpha[i], true);
  }

  return log_density;
}

double log_density_u_strided(double const *u, double const *alpha,
                             unsigned int length, unsigned int stride)
{
  if(!valid_alpha(alpha, length)){
    return -std::numeric_limits<double>::infinity();
  }

  return log_density_u_strided_valid_alpha(u, alpha, length, stride);
}

double log_density_u(double const *u, double const *alpha, unsigned int length)
{
  return log_density_u_strided(u, alpha, length, 1);
}

}
}
