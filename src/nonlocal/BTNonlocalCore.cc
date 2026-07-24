#include "BTNonlocalCore.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <Rmath.h>

namespace bayestools {
namespace nonlocal {

namespace {

bool valid_order(double order)
{
  return std::isfinite(order) && order >= 1.0 && std::floor(order) == order;
}

double quiet_nan()
{
  return std::numeric_limits<double>::quiet_NaN();
}

double log_half()
{
  return -std::log(2.0);
}

double log1mexp(double log_p)
{
  if(log_p == 0.0){
    return -std::numeric_limits<double>::infinity();
  }
  if(log_p < -std::log(2.0)){
    return std::log1p(-std::exp(log_p));
  }
  return std::log(-std::expm1(log_p));
}

double clamp_probability(double p)
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
  p = clamp_probability(p);
  if(!lower_tail){
    p = 1.0 - p;
  }
  if(log_p){
    return std::log(p);
  }
  return p;
}

double return_half_probability(double p, bool log_p)
{
  if(log_p){
    return log_half() + p;
  }
  return 0.5 * p;
}

double return_one_minus_half_probability(double p, bool log_p)
{
  if(log_p){
    return log1mexp(log_half() + p);
  }
  return 1.0 - 0.5 * p;
}

bool inside(double x, double lower, double upper)
{
  return x >= lower && x <= upper && std::isfinite(x);
}

}

bool valid_common_parameters(double location, double tau, double order)
{
  return std::isfinite(location) && std::isfinite(tau) && tau > 0.0 &&
    valid_order(order);
}

bool valid_invmoment_parameters(double location, double tau, double order, double df)
{
  return valid_common_parameters(location, tau, order) &&
    std::isfinite(df) && df > 0.0;
}

double moment_mode(double tau, double order)
{
  return std::sqrt(2.0 * order * tau);
}

double invmoment_mode(double tau, double order, double df)
{
  return std::sqrt(tau) * std::pow(2.0 * order / (df + 1.0), 1.0 / (2.0 * order));
}

double moment_log_density(double x, double location, double tau, double order)
{
  if(!valid_common_parameters(location, tau, order)){
    return -std::numeric_limits<double>::infinity();
  }
  double delta = x - location;
  if(!std::isfinite(delta) || delta == 0.0){
    return -std::numeric_limits<double>::infinity();
  }

  double log_double_factorial =
    lgammafn(2.0 * order + 1.0) - order * std::log(2.0) - lgammafn(order + 1.0);

  return 2.0 * order * std::log(std::fabs(delta)) +
    dnorm(delta, 0.0, std::sqrt(tau), true) -
    order * std::log(tau) - log_double_factorial;
}

double invmoment_log_density(double x, double location, double tau, double order, double df)
{
  if(!valid_invmoment_parameters(location, tau, order, df)){
    return -std::numeric_limits<double>::infinity();
  }
  double delta = x - location;
  if(!std::isfinite(delta) || delta == 0.0){
    return -std::numeric_limits<double>::infinity();
  }

  return std::log(order) + (df / 2.0) * std::log(tau) -
    lgammafn(df / (2.0 * order)) -
    (df + 1.0) * std::log(std::fabs(delta)) -
    std::pow(tau / (delta * delta), order);
}

double moment_cdf(double q, double location, double tau, double order,
                  bool lower_tail, bool log_p)
{
  if(!valid_common_parameters(location, tau, order)){
    return quiet_nan();
  }
  if(std::isnan(q)){
    return quiet_nan();
  }
  if(q == -std::numeric_limits<double>::infinity()){
    return return_probability(0.0, lower_tail, log_p);
  }
  if(q == std::numeric_limits<double>::infinity()){
    return return_probability(1.0, lower_tail, log_p);
  }

  double delta = q - location;
  double z = delta * delta / tau;
  if(delta < 0.0){
    if(lower_tail){
      return return_half_probability(
        pchisq(z, 2.0 * order + 1.0, false, log_p),
        log_p
      );
    }
    return return_one_minus_half_probability(
      pchisq(z, 2.0 * order + 1.0, false, log_p),
      log_p
    );
  }
  if(delta > 0.0){
    if(lower_tail){
      return return_one_minus_half_probability(
        pchisq(z, 2.0 * order + 1.0, false, log_p),
        log_p
      );
    }
    return return_half_probability(
      pchisq(z, 2.0 * order + 1.0, false, log_p),
      log_p
    );
  }

  return return_probability(0.5, lower_tail, log_p);
}

double invmoment_cdf(double q, double location, double tau, double order,
                     double df, bool lower_tail, bool log_p)
{
  if(!valid_invmoment_parameters(location, tau, order, df)){
    return quiet_nan();
  }
  if(std::isnan(q)){
    return quiet_nan();
  }
  if(q == -std::numeric_limits<double>::infinity()){
    return return_probability(0.0, lower_tail, log_p);
  }
  if(q == std::numeric_limits<double>::infinity()){
    return return_probability(1.0, lower_tail, log_p);
  }

  double delta = q - location;
  double s = std::pow(tau / (delta * delta), order);
  if(delta < 0.0){
    if(lower_tail){
      return return_half_probability(
        pgamma(s, df / (2.0 * order), 1.0, true, log_p),
        log_p
      );
    }
    return return_one_minus_half_probability(
      pgamma(s, df / (2.0 * order), 1.0, true, log_p),
      log_p
    );
  }
  if(delta > 0.0){
    if(lower_tail){
      return return_one_minus_half_probability(
        pgamma(s, df / (2.0 * order), 1.0, true, log_p),
        log_p
      );
    }
    return return_half_probability(
      pgamma(s, df / (2.0 * order), 1.0, true, log_p),
      log_p
    );
  }

  return return_probability(0.5, lower_tail, log_p);
}

double moment_quantile(double p, double location, double tau, double order,
                       bool lower_tail, bool log_p)
{
  if(!valid_common_parameters(location, tau, order)){
    return quiet_nan();
  }
  if(std::isnan(p) || (log_p ? p > 0.0 : p < 0.0 || p > 1.0)){
    return quiet_nan();
  }

  bool zero_probability = log_p ?
    p == -std::numeric_limits<double>::infinity() :
    p == 0.0;
  bool unit_probability = log_p ? p == 0.0 : p == 1.0;
  double half_probability = log_p ? log_half() : 0.5;

  if(!lower_tail){
    if(zero_probability){
      return std::numeric_limits<double>::infinity();
    }
    if(unit_probability){
      return -std::numeric_limits<double>::infinity();
    }
    if(p == half_probability){
      return location;
    }
    if(p < half_probability){
      double p_tail = log_p ? p + std::log(2.0) : 2.0 * p;
      double z = qchisq(p_tail, 2.0 * order + 1.0, false, log_p);
      return location + std::sqrt(tau * z);
    }

    double p_tail = log_p ?
      log1mexp(p) + std::log(2.0) :
      2.0 * (1.0 - p);
    double z = qchisq(p_tail, 2.0 * order + 1.0, false, log_p);
    return location - std::sqrt(tau * z);
  }

  if(zero_probability){
    return -std::numeric_limits<double>::infinity();
  }
  if(unit_probability){
    return std::numeric_limits<double>::infinity();
  }
  if(p == half_probability){
    return location;
  }

  if(p < half_probability){
    double p_tail = log_p ? p + std::log(2.0) : 2.0 * p;
    double z = qchisq(p_tail, 2.0 * order + 1.0, false, log_p);
    return location - std::sqrt(tau * z);
  }

  double p_tail = log_p ?
    log1mexp(p) + std::log(2.0) :
    2.0 * (1.0 - p);
  double z = qchisq(p_tail, 2.0 * order + 1.0, false, log_p);
  return location + std::sqrt(tau * z);
}

double invmoment_quantile(double p, double location, double tau, double order,
                          double df, bool lower_tail, bool log_p)
{
  if(!valid_invmoment_parameters(location, tau, order, df)){
    return quiet_nan();
  }
  if(std::isnan(p) || (log_p ? p > 0.0 : p < 0.0 || p > 1.0)){
    return quiet_nan();
  }

  double shape = df / (2.0 * order);
  bool zero_probability = log_p ?
    p == -std::numeric_limits<double>::infinity() :
    p == 0.0;
  bool unit_probability = log_p ? p == 0.0 : p == 1.0;
  double half_probability = log_p ? log_half() : 0.5;

  if(!lower_tail){
    if(zero_probability){
      return std::numeric_limits<double>::infinity();
    }
    if(unit_probability){
      return -std::numeric_limits<double>::infinity();
    }
    if(p == half_probability){
      return location;
    }
    if(p < half_probability){
      double p_tail = log_p ? p + std::log(2.0) : 2.0 * p;
      double s = qgamma(p_tail, shape, 1.0, true, log_p);
      return location + std::sqrt(tau / std::pow(s, 1.0 / order));
    }

    double p_tail = log_p ?
      log1mexp(p) + std::log(2.0) :
      2.0 * (1.0 - p);
    double s = qgamma(p_tail, shape, 1.0, true, log_p);
    return location - std::sqrt(tau / std::pow(s, 1.0 / order));
  }

  if(zero_probability){
    return -std::numeric_limits<double>::infinity();
  }
  if(unit_probability){
    return std::numeric_limits<double>::infinity();
  }
  if(p == half_probability){
    return location;
  }

  if(p < half_probability){
    double p_tail = log_p ? p + std::log(2.0) : 2.0 * p;
    double s = qgamma(p_tail, shape, 1.0, true, log_p);
    return location - std::sqrt(tau / std::pow(s, 1.0 / order));
  }

  double p_tail = log_p ?
    log1mexp(p) + std::log(2.0) :
    2.0 * (1.0 - p);
  double s = qgamma(p_tail, shape, 1.0, true, log_p);
  return location + std::sqrt(tau / std::pow(s, 1.0 / order));
}

double moment_rng(double u_sign, double u_size, double location, double tau,
                  double order)
{
  u_size = std::min(1.0 - DBL_EPSILON, std::max(DBL_MIN, u_size));
  double sign = u_sign < 0.5 ? -1.0 : 1.0;
  double z = qchisq(u_size, 2.0 * order + 1.0, true, false);
  return location + sign * std::sqrt(tau * z);
}

double invmoment_rng(double u_sign, double u_size, double location, double tau,
                     double order, double df)
{
  u_size = std::min(1.0 - DBL_EPSILON, std::max(DBL_MIN, u_size));
  double sign = u_sign < 0.5 ? -1.0 : 1.0;
  double s = qgamma(u_size, df / (2.0 * order), 1.0, true, false);
  return location + sign * std::sqrt(tau / std::pow(s, 1.0 / order));
}

double typical_value(double location, double mode_abs, double lower,
                     double upper)
{
  if(std::isfinite(mode_abs) && mode_abs > 0.0){
    double upper_mode = location + mode_abs;
    if(inside(upper_mode, lower, upper)){
      return upper_mode;
    }
    double lower_mode = location - mode_abs;
    if(inside(lower_mode, lower, upper)){
      return lower_mode;
    }
  }

  if(std::isfinite(lower) && std::isfinite(upper) && lower < upper){
    double width = upper - lower;
    double candidate = location <= (lower + upper) / 2.0 ?
      lower + 0.75 * width : lower + 0.25 * width;
    if(candidate == location){
      candidate = lower + 0.25 * width;
    }
    return candidate;
  }

  double step = std::isfinite(mode_abs) && mode_abs > 0.0 ?
    mode_abs : std::max(1.0, std::fabs(location) * 0.1);
  if(std::isfinite(lower)){
    return lower + step;
  }
  if(std::isfinite(upper)){
    return upper - step;
  }
  return location + step;
}

}
}
