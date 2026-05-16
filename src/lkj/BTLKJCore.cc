#include "BTLKJCore.h"

#include <algorithm>
#include <cmath>
#include <JRmath.h>
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

double cpc_from_u(double u)
{
  return 2.0 * u - 1.0;
}

double pair_alpha(unsigned int i, unsigned int K, double eta)
{
  return eta + (static_cast<double>(K - i - 1)) / 2.0;
}

void fill_alpha(double *alpha, unsigned int K, double eta)
{
  unsigned int p = 0;
  for(unsigned int j = 2; j <= K; ++j){
    for(unsigned int i = 1; i < j; ++i){
      alpha[p] = pair_alpha(i, K, eta);
      ++p;
    }
  }
}

double sqrt_one_minus_square(double x)
{
  double value = 1.0 - x * x;
  if(value < 0.0 && value > -1e-14){
    value = 0.0;
  }
  return std::sqrt(value);
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
      double cpc = cpc_from_u(u[cpc_index(column, row) * stride]);
      value[flat_index(row, column, K)] = cpc * multiplier;
      multiplier *= sqrt_one_minus_square(cpc);
    }
    value[flat_index(row, row, K)] = multiplier;
  }
}

void fill_corr_from_cholesky(double *value, double const *L, unsigned int K)
{
  for(unsigned int row = 1; row <= K; ++row){
    for(unsigned int column = 1; column <= row; ++column){
      double cell = 0.0;
      for(unsigned int m = 1; m <= column; ++m){
        cell += L[flat_index(row, m, K)] * L[flat_index(column, m, K)];
      }
      value[flat_index(row, column, K)] = cell;
      if(column != row){
        value[flat_index(column, row, K)] = cell;
      }
    }
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
    return R_NegInf;
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
    return R_NegInf;
  }

  return log_density_u_strided_valid_alpha(u, alpha, length, stride);
}

double log_density_u(double const *u, double const *alpha, unsigned int length)
{
  return log_density_u_strided(u, alpha, length, 1);
}

}
}
