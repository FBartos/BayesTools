#include "BTInvGammaCore.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <Rmath.h>

namespace bayestools {
namespace invgamma {

namespace {

double quiet_nan()
{
  return std::numeric_limits<double>::quiet_NaN();
}

double clamp_unit(double p)
{
  if(p < 0.0){
    return 0.0;
  }
  if(p > 1.0){
    return 1.0;
  }
  return p;
}

double return_probability(double p, bool lower_tail, bool log_p)
{
  p = clamp_unit(p);
  if(!lower_tail){
    p = 1.0 - p;
  }
  return log_p ? std::log(p) : p;
}

bool inside(double x, double lower, double upper)
{
  return x >= lower && x <= upper && std::isfinite(x);
}

}

bool valid_parameters(double shape, double scale)
{
  return std::isfinite(shape) && shape > 0.0 &&
    std::isfinite(scale) && scale > 0.0;
}

bool valid_value(double x)
{
  return std::isfinite(x) && x > 0.0;
}

double mode(double shape, double scale)
{
  if(!valid_parameters(shape, scale)){
    return quiet_nan();
  }
  return scale / (shape + 1.0);
}

double log_density(double x, double shape, double scale)
{
  if(!valid_parameters(shape, scale)){
    return -std::numeric_limits<double>::infinity();
  }
  if(std::isnan(x)){
    return quiet_nan();
  }
  if(x <= 0.0 || !std::isfinite(x)){
    return -std::numeric_limits<double>::infinity();
  }

  return shape * std::log(scale) - lgammafn(shape) -
    (shape + 1.0) * std::log(x) - scale / x;
}

double cdf(double q, double shape, double scale, bool lower_tail, bool log_p)
{
  if(!valid_parameters(shape, scale)){
    return quiet_nan();
  }
  if(std::isnan(q)){
    return quiet_nan();
  }
  if(q <= 0.0){
    return return_probability(0.0, lower_tail, log_p);
  }
  if(q == std::numeric_limits<double>::infinity()){
    return return_probability(1.0, lower_tail, log_p);
  }

  return pgamma(1.0 / q, shape, 1.0 / scale, !lower_tail, log_p);
}

double quantile(double p, double shape, double scale, bool lower_tail, bool log_p)
{
  if(!valid_parameters(shape, scale)){
    return quiet_nan();
  }

  double prob = log_p ? std::exp(p) : p;
  if(std::isnan(prob) || prob < 0.0 || prob > 1.0){
    return quiet_nan();
  }

  if(lower_tail){
    if(prob == 0.0){
      return 0.0;
    }
    if(prob == 1.0){
      return std::numeric_limits<double>::infinity();
    }
    double gamma_q = qgamma(p, shape, 1.0 / scale, false, log_p);
    return 1.0 / gamma_q;
  }

  if(prob == 0.0){
    return std::numeric_limits<double>::infinity();
  }
  if(prob == 1.0){
    return 0.0;
  }
  double gamma_q = qgamma(p, shape, 1.0 / scale, true, log_p);
  return 1.0 / gamma_q;
}

double rng(double u, double shape, double scale)
{
  if(!valid_parameters(shape, scale)){
    return quiet_nan();
  }
  u = std::min(1.0 - DBL_EPSILON, std::max(DBL_MIN, u));
  return 1.0 / qgamma(u, shape, 1.0 / scale, true, false);
}

double typical_value(double shape, double scale, double lower, double upper)
{
  double value = mode(shape, scale);
  if(inside(value, lower, upper)){
    return value;
  }

  double l = std::max(lower, 0.0);
  double cdf_l = cdf(l, shape, scale, true, false);
  double cdf_u = cdf(upper, shape, scale, true, false);
  double mass = cdf_u - cdf_l;
  if(std::isfinite(mass) && mass > 0.0){
    value = quantile(cdf_l + 0.5 * mass, shape, scale, true, false);
    if(inside(value, lower, upper)){
      return value;
    }
  }

  double surv_l = cdf(l, shape, scale, false, false);
  double surv_u = cdf(upper, shape, scale, false, false);
  mass = surv_l - surv_u;
  if(std::isfinite(mass) && mass > 0.0){
    value = quantile(surv_u + 0.5 * mass, shape, scale, false, false);
    if(inside(value, lower, upper)){
      return value;
    }
  }

  if(std::isfinite(lower) && std::isfinite(upper) && lower < upper){
    return lower + 0.5 * (upper - lower);
  }
  if(std::isfinite(lower)){
    return lower + std::max(1.0, std::fabs(lower) * 0.1);
  }
  if(std::isfinite(upper) && upper > 0.0){
    return upper * 0.5;
  }
  return mode(shape, scale);
}

}
}
