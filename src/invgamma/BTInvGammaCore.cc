#include "BTInvGammaCore.h"
#include "../gamma/BTGammaCore.h"
#include "../gamma/BTGammaRange.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>

namespace bayestools {
namespace invgamma {

namespace {

const double infinity = std::numeric_limits<double>::infinity();

double quiet_nan()
{
  return std::numeric_limits<double>::quiet_NaN();
}

double radial_log(double x, double scale, double *represented_r = 0, double *coordinate_error = 0)
{
  if(represented_r) *represented_r = quiet_nan();
  if(coordinate_error) *coordinate_error = quiet_nan();
  if(x <= 0.0) return infinity;
  if(x == infinity) return -infinity;
  double r = scale / x;
  if(std::isfinite(r) && r >= DBL_MIN){
    if(represented_r) *represented_r = r;
    return std::log(r);
  }
  double h = std::log(scale) - std::log(x);
  if(coordinate_error && h > gamma::range::far_threshold &&
     gamma::range::current_environment().available){
    *coordinate_error = gamma::range::invgamma_error(h);
  }
  return h;
}

bool candidate_available(double x, double shape, double scale, double lower, double upper)
{
  return std::isfinite(x) && x >= lower && x <= upper &&
    std::isfinite(log_density(x, shape, scale));
}

}

bool valid_parameters(double shape, double scale)
{
  return std::isfinite(shape) && shape > 0.0 && std::isfinite(scale) && scale > 0.0;
}

bool valid_value(double x)
{
  return std::isfinite(x) && x > 0.0;
}

double mode(double shape, double scale)
{
  return valid_parameters(shape, scale) ? scale / (shape + 1.0) : quiet_nan();
}

double log_density(double x, double shape, double scale)
{
  if(!valid_parameters(shape, scale)) return -infinity;
  if(std::isnan(x)) return quiet_nan();
  if(x <= 0.0 || !std::isfinite(x)) return -infinity;
  double represented_r;
  double coordinate_error;
  double log_r = radial_log(x, scale, &represented_r, &coordinate_error);
  return gamma::log_prefix(shape, log_r, represented_r, coordinate_error) - std::log(x);
}

double cdf(double q, double shape, double scale, bool lower_tail, bool log_p)
{
  if(!valid_parameters(shape, scale) || std::isnan(q)) return quiet_nan();
  if(q <= 0.0) return log_p ? (lower_tail ? -infinity : 0.0) : (lower_tail ? 0.0 : 1.0);
  if(q == infinity) return log_p ? (lower_tail ? 0.0 : -infinity) : (lower_tail ? 1.0 : 0.0);
  double represented_r;
  double coordinate_error;
  double log_r = radial_log(q, scale, &represented_r, &coordinate_error);
  gamma::tails out = gamma::probabilities(shape, log_r, represented_r, coordinate_error);
  return log_p ? (lower_tail ? out.log_upper : out.log_lower) :
    (lower_tail ? out.upper : out.lower);
}

double quantile(double p, double shape, double scale, bool lower_tail, bool log_p)
{
  if(!valid_parameters(shape, scale) || std::isnan(p) ||
     (log_p ? p > 0.0 : p < 0.0 || p > 1.0)) return quiet_nan();
  bool zero = log_p ? p == -infinity : p == 0.0;
  bool unit = log_p ? p == 0.0 : p == 1.0;
  if(zero) return lower_tail ? 0.0 : infinity;
  if(unit) return lower_tail ? infinity : 0.0;
  double represented_root;
  double root = gamma::log_quantile(log_p ? p : std::log(p), shape, !lower_tail,
                                    &represented_root);
  if(std::isfinite(represented_root)) return scale / represented_root;
  return std::exp(std::log(scale) - root);
}

double rng(double u, double shape, double scale)
{
  u = std::min(1.0 - DBL_EPSILON, std::max(DBL_MIN, u));
  return quantile(u, shape, scale, false, false);
}

double log_interval_mass(double shape, double scale, double lower, double upper)
{
  if(!valid_parameters(shape, scale) || std::isnan(lower) || std::isnan(upper)) return quiet_nan();
  if(upper <= 0.0 || lower >= upper) return -infinity;
  double lower_error, upper_error;
  double lo = radial_log(upper, scale, 0, &lower_error);
  double hi = radial_log(lower, scale, 0, &upper_error);
  return gamma::log_interval_mass(shape, lo, hi, lower_error, upper_error);
}

double truncated_quantile(double p, double shape, double scale, double lower, double upper)
{
  if(!valid_parameters(shape, scale) || !(p >= 0.0 && p <= 1.0 && lower < upper)) return quiet_nan();
  if(p == 0.0) return std::max(lower, 0.0);
  if(p == 1.0) return upper;
  double lower_error, upper_error;
  double lo = radial_log(upper, scale, 0, &lower_error);
  double hi = radial_log(lower, scale, 0, &upper_error);
  double root = gamma::interval_log_quantile(std::log1p(-p), shape,
    lo, hi, lower_error, upper_error);
  double x = std::exp(std::log(scale) - root);
  return std::isfinite(x) && x > 0.0 && x > lower && x < upper ? x : quiet_nan();
}

double typical_value(double shape, double scale, double lower, double upper)
{
  double value = mode(shape, scale);
  if(candidate_available(value, shape, scale, lower, upper)) return value;
  value = truncated_quantile(0.5, shape, scale, lower, upper);
  if(candidate_available(value, shape, scale, lower, upper)) return value;
  if(std::isfinite(lower) && std::isfinite(upper) && lower < upper){
    value = lower + 0.5 * (upper - lower);
    if(candidate_available(value, shape, scale, lower, upper)) return value;
  }
  if(std::isfinite(lower)){
    value = lower + std::max(1.0, std::fabs(lower) * 0.1);
    if(candidate_available(value, shape, scale, lower, upper)) return value;
  }
  if(std::isfinite(upper) && upper > 0.0){
    value = upper * 0.5;
    if(candidate_available(value, shape, scale, lower, upper)) return value;
  }
  return quiet_nan();
}

}
}
