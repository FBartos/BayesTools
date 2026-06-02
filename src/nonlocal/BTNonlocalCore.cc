#include "BTNonlocalCore.h"

#include <cmath>
#include <limits>
#include <JRmath.h>

namespace bayestools {
namespace nonlocal {

namespace {

double maybe_log(double p, bool log_p)
{
  return log_p ? std::log(p) : p;
}

double log_two()
{
  return std::log(2.0);
}

double log_half()
{
  return -std::log(2.0);
}

double half_times_tail(double tail, bool log_tail, bool log_p)
{
  if(log_tail){
    return log_p ? log_half() + tail : 0.5 * std::exp(tail);
  }

  return log_p ? log_half() + std::log(tail) : 0.5 * tail;
}

double half_plus_half(double p, bool log_p)
{
  return log_p ? std::log1p(p) + log_half() : 0.5 + 0.5 * p;
}

double open_unit(double u)
{
  if(u <= 0.0){
    return std::nextafter(0.0, 1.0);
  }
  if(u >= 1.0){
    return std::nextafter(1.0, 0.0);
  }
  return u;
}

double nan_value()
{
  return std::numeric_limits<double>::quiet_NaN();
}

bool inside(double x, double lower, double upper)
{
  return std::isfinite(x) && x >= lower && x <= upper;
}

double nudge_from_location(double x, double location, double lower, double upper)
{
  if(x != location){
    return x;
  }

  double right = std::nextafter(location, upper);
  if(inside(right, lower, upper) && right != location){
    return right;
  }

  double left = std::nextafter(location, lower);
  if(inside(left, lower, upper) && left != location){
    return left;
  }

  return x;
}

}

bool valid_order(double order)
{
  return std::isfinite(order) && order >= 1.0 && std::fabs(order - std::round(order)) < 1e-8;
}

bool valid_moment_parameters(double location, double tau, double order)
{
  return std::isfinite(location) && std::isfinite(tau) && tau > 0.0 && valid_order(order);
}

bool valid_invmoment_parameters(double location, double tau, double order, double df)
{
  return valid_moment_parameters(location, tau, order) && std::isfinite(df) && df > 0.0;
}

double moment_mode(double tau, double order)
{
  return std::sqrt(2.0 * order * tau);
}

double invmoment_mode(double tau, double order, double df)
{
  return std::sqrt(tau) * std::pow(2.0 * order / (df + 1.0), 1.0 / (2.0 * order));
}

double typical_value(double location, double mode, double const *lower, double const *upper)
{
  double lower_value = lower == nullptr ? R_NegInf : *lower;
  double upper_value = upper == nullptr ? R_PosInf : *upper;

  double right_mode = location + mode;
  if(inside(right_mode, lower_value, upper_value)){
    return right_mode;
  }

  double left_mode = location - mode;
  if(inside(left_mode, lower_value, upper_value)){
    return left_mode;
  }

  double candidate;
  if(std::isfinite(lower_value) && std::isfinite(upper_value)){
    candidate = lower_value / 2.0 + upper_value / 2.0;
  }else if(std::isfinite(lower_value)){
    candidate = lower_value + mode;
    if(!std::isfinite(candidate)){
      candidate = std::nextafter(lower_value, R_PosInf);
    }
  }else if(std::isfinite(upper_value)){
    candidate = upper_value - mode;
    if(!std::isfinite(candidate)){
      candidate = std::nextafter(upper_value, R_NegInf);
    }
  }else{
    candidate = right_mode;
  }

  return nudge_from_location(candidate, location, lower_value, upper_value);
}

double dmoment(double x, double location, double tau, double order, bool give_log)
{
  if(std::isnan(x)){
    return x;
  }
  if(!valid_moment_parameters(location, tau, order)){
    return give_log ? R_NegInf : 0.0;
  }

  double delta = x - location;
  if(!std::isfinite(delta) || delta == 0.0){
    return give_log ? R_NegInf : 0.0;
  }

  double log_double_factorial = lgammafn(2.0 * order + 1.0) -
    order * std::log(2.0) - lgammafn(order + 1.0);
  double log_density = 2.0 * order * std::log(std::fabs(delta)) +
    dnorm(delta, 0.0, std::sqrt(tau), true) -
    order * std::log(tau) -
    log_double_factorial;

  return give_log ? log_density : std::exp(log_density);
}

double pmoment(double q, double location, double tau, double order, bool lower_tail, bool log_p)
{
  if(std::isnan(q)){
    return q;
  }
  if(!valid_moment_parameters(location, tau, order)){
    return nan_value();
  }

  double delta = q - location;
  if(delta == 0.0){
    return maybe_log(0.5, log_p);
  }

  double z = delta * delta / tau;
  double df = 2.0 * order + 1.0;
  if(delta < 0.0){
    if(lower_tail){
      return half_times_tail(pchisq(z, df, false, log_p), log_p, log_p);
    }
    return half_plus_half(pchisq(z, df, true, false), log_p);
  }

  if(lower_tail){
    return half_plus_half(pchisq(z, df, true, false), log_p);
  }
  return half_times_tail(pchisq(z, df, false, log_p), log_p, log_p);
}

double qmoment(double p, double location, double tau, double order, bool lower_tail, bool log_p)
{
  if(std::isnan(p)){
    return p;
  }
  if(!valid_moment_parameters(location, tau, order)){
    return nan_value();
  }

  double df = 2.0 * order + 1.0;
  double z;
  double direction;

  if(log_p){
    if(p > 0.0){
      return nan_value();
    }
    if(std::isinf(p) && p < 0.0){
      return lower_tail ? R_NegInf : R_PosInf;
    }
    if(p == 0.0){
      return lower_tail ? R_PosInf : R_NegInf;
    }
    if(p == log_half()){
      return location;
    }

    if(p < log_half()){
      direction = lower_tail ? -1.0 : 1.0;
      z = qchisq(log_two() + p, df, false, true);
    }else{
      direction = lower_tail ? 1.0 : -1.0;
      z = qchisq(log_two() + std::log(-std::expm1(p)), df, false, true);
    }

    return location + direction * std::sqrt(tau * z);
  }

  if(p < 0.0 || p > 1.0){
    return nan_value();
  }
  if(p <= 0.0){
    return lower_tail ? R_NegInf : R_PosInf;
  }
  if(p >= 1.0){
    return lower_tail ? R_PosInf : R_NegInf;
  }
  if(p == 0.5){
    return location;
  }

  if(p < 0.5){
    direction = lower_tail ? -1.0 : 1.0;
    z = qchisq(2.0 * p, df, false, false);
  }else{
    direction = lower_tail ? 1.0 : -1.0;
    z = qchisq(2.0 * (1.0 - p), df, false, false);
  }

  return location + direction * std::sqrt(tau * z);
}

double rmoment(double sign_u, double magnitude_u, double location, double tau, double order)
{
  double sign = sign_u < 0.5 ? -1.0 : 1.0;
  double y = qchisq(open_unit(magnitude_u), 2.0 * order + 1.0, true, false);
  return location + sign * std::sqrt(tau * y);
}

double dinvmoment(double x, double location, double tau, double order, double df, bool give_log)
{
  if(std::isnan(x)){
    return x;
  }
  if(!valid_invmoment_parameters(location, tau, order, df)){
    return give_log ? R_NegInf : 0.0;
  }

  double delta = x - location;
  if(!std::isfinite(delta) || delta == 0.0){
    return give_log ? R_NegInf : 0.0;
  }

  double log_density = std::log(order) +
    (df / 2.0) * std::log(tau) -
    lgammafn(df / (2.0 * order)) -
    (df + 1.0) * std::log(std::fabs(delta)) -
    std::pow(tau / (delta * delta), order);

  return give_log ? log_density : std::exp(log_density);
}

double pinvmoment(double q, double location, double tau, double order, double df, bool lower_tail, bool log_p)
{
  if(std::isnan(q)){
    return q;
  }
  if(!valid_invmoment_parameters(location, tau, order, df)){
    return nan_value();
  }

  double delta = q - location;
  if(delta == 0.0){
    return maybe_log(0.5, log_p);
  }

  double u = std::pow(tau / (delta * delta), order);
  double shape = df / (2.0 * order);
  if(delta < 0.0){
    if(lower_tail){
      return half_times_tail(pgamma(u, shape, 1.0, true, log_p), log_p, log_p);
    }
    return half_plus_half(pgamma(u, shape, 1.0, false, false), log_p);
  }

  if(lower_tail){
    return half_plus_half(pgamma(u, shape, 1.0, false, false), log_p);
  }
  return half_times_tail(pgamma(u, shape, 1.0, true, log_p), log_p, log_p);
}

double qinvmoment(double p, double location, double tau, double order, double df, bool lower_tail, bool log_p)
{
  if(std::isnan(p)){
    return p;
  }
  if(!valid_invmoment_parameters(location, tau, order, df)){
    return nan_value();
  }

  double shape = df / (2.0 * order);
  double u;
  double direction;

  if(log_p){
    if(p > 0.0){
      return nan_value();
    }
    if(std::isinf(p) && p < 0.0){
      return lower_tail ? R_NegInf : R_PosInf;
    }
    if(p == 0.0){
      return lower_tail ? R_PosInf : R_NegInf;
    }
    if(p == log_half()){
      return location;
    }

    if(p < log_half()){
      direction = lower_tail ? -1.0 : 1.0;
      u = qgamma(log_two() + p, shape, 1.0, true, true);
    }else{
      direction = lower_tail ? 1.0 : -1.0;
      u = qgamma(log_two() + std::log(-std::expm1(p)), shape, 1.0, true, true);
    }

    return location + direction * std::sqrt(tau) * std::pow(u, -1.0 / (2.0 * order));
  }

  if(p < 0.0 || p > 1.0){
    return nan_value();
  }
  if(p <= 0.0){
    return lower_tail ? R_NegInf : R_PosInf;
  }
  if(p >= 1.0){
    return lower_tail ? R_PosInf : R_NegInf;
  }
  if(p == 0.5){
    return location;
  }

  if(p < 0.5){
    direction = lower_tail ? -1.0 : 1.0;
    u = qgamma(2.0 * p, shape, 1.0, true, false);
  }else{
    direction = lower_tail ? 1.0 : -1.0;
    u = qgamma(2.0 * (1.0 - p), shape, 1.0, true, false);
  }

  return location + direction * std::sqrt(tau) * std::pow(u, -1.0 / (2.0 * order));
}

double rinvmoment(double sign_u, double magnitude_u, double location, double tau, double order, double df)
{
  double sign = sign_u < 0.5 ? -1.0 : 1.0;
  double u = qgamma(open_unit(magnitude_u), df / (2.0 * order), 1.0, true, false);
  return location + sign * std::sqrt(tau) * std::pow(u, -1.0 / (2.0 * order));
}

}
}
