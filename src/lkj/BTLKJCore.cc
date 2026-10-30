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

// 'position' holds a permutation of 1, ..., K (as doubles, the JAGS argument type).
bool valid_position(double const *position, unsigned int K)
{
  std::vector<bool> seen(K, false);
  for(unsigned int i = 0; i < K; ++i){
    const double value = position[i];
    if(!std::isfinite(value) || value < 1.0 || value > static_cast<double>(K) ||
       std::floor(value) != value){
      return false;
    }
    const unsigned int index = static_cast<unsigned int>(value) - 1;
    if(seen[index]){
      return false;
    }
    seen[index] = true;
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

void fill_cholesky_from_u_ordered(double *value, double const *u, unsigned int K, double const *position)
{
  fill_cholesky_from_u_ordered_strided(value, u, K, 1, position);
}

// The lower Cholesky factor, in design order, of the correlation matrix whose
// vine-order factor fill_cholesky_from_u() builds, where design column d is
// vine variable position[d]. The rows of the vine factor are taken in design
// order and re-triangularized by Givens rotations of their columns (an LQ
// decomposition): the rotations preserve the row inner products, so L L' is the
// reordered correlation matrix up to rounding, and the correlation matrix itself
// is never formed (its Cholesky factorization would lose the small positive
// diagonal that a primitive near 0 or 1 leaves). An entry that is already zero
// is not rotated, so the identity order returns the vine factor unchanged. A
// rank-deficient factor (the clamped extension of an out-of-support proposal)
// stays finite. In another order a diagonal entry is a conditional SD given the
// preceding design columns, which the rotations resolve to absolute rounding
// accuracy only: a conditional SD below rounding of the unit-norm rows (a
// numerically singular correlation matrix) can come out as 0, whereas the vine
// order keeps the relative accuracy of its products of innovations.
void fill_cholesky_from_u_ordered_strided(double *value, double const *u, unsigned int K, unsigned int stride,
                                          double const *position)
{
  std::vector<double> vine(K * K);
  fill_cholesky_from_u_strided(&vine[0], u, K, stride);

  for(unsigned int row = 1; row <= K; ++row){
    const unsigned int source = static_cast<unsigned int>(position[row - 1]);
    for(unsigned int column = 1; column <= K; ++column){
      value[flat_index(row, column, K)] = vine[flat_index(source, column, K)];
    }
  }

  for(unsigned int row = 1; row <= K; ++row){
    for(unsigned int column = K; column > row; --column){
      const double b = value[flat_index(row, column, K)];
      if(b == 0.0){
        continue;
      }
      // The rotation is computed from the pivot pair scaled to a maximum of 1:
      // c and s from a subnormal radius would carry its rounding (c^2 + s^2
      // far from 1) into the rotated rows.
      const double a = value[flat_index(row, row, K)];
      const double scale = std::max(std::fabs(a), std::fabs(b));
      const double scaled_r = std::hypot(a / scale, b / scale);
      const double c = (a / scale) / scaled_r;
      const double s = (b / scale) / scaled_r;
      const double r = scale * scaled_r;
      for(unsigned int other = row; other <= K; ++other){
        const double x = value[flat_index(other, row, K)];
        const double y = value[flat_index(other, column, K)];
        value[flat_index(other, row, K)] = c * x + s * y;
        value[flat_index(other, column, K)] = c * y - s * x;
      }
      value[flat_index(row, row, K)] = r;
      value[flat_index(row, column, K)] = 0.0;
    }
    if(value[flat_index(row, row, K)] < 0.0){
      for(unsigned int other = row; other <= K; ++other){
        value[flat_index(other, row, K)] = -value[flat_index(other, row, K)];
      }
    }
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
