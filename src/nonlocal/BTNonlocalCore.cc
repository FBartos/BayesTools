#include "BTNonlocalCore.h"
#include "../gamma/BTGammaCore.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>

namespace bayestools {
namespace nonlocal {

namespace {

const double infinity = std::numeric_limits<double>::infinity();
const double log_two = std::log(2.0);

double quiet_nan()
{
  return std::numeric_limits<double>::quiet_NaN();
}

bool valid_order(double order)
{
  return std::isfinite(order) && order >= 1.0 && std::floor(order) == order;
}

double shape_value(double order, double df, bool inverse)
{
  return inverse ? (df / order) * 0.5 : order + 0.5;
}

double log_distance(double x, double location)
{
  double distance = std::fabs(x - location);
  if(std::isfinite(distance)) return std::log(distance);
  return gamma::logaddexp(std::log(std::fabs(x)), std::log(std::fabs(location)));
}

double radial_log(double x, double location, double tau, double order, bool inverse,
                  double *represented_r = 0)
{
  if(represented_r) *represented_r = quiet_nan();
  if(x == location) return inverse ? infinity : -infinity;
  if(!std::isfinite(x)) return inverse ? -infinity : infinity;
  double distance = std::fabs(x - location);
  double first = inverse ? tau / distance : distance / tau;
  double second = inverse ? first / distance : distance / 2.0;
  double ratio = inverse ? second : first * second;
  if(std::isfinite(distance) && first >= DBL_MIN && std::isfinite(first) &&
     second >= DBL_MIN && std::isfinite(second) && ratio >= DBL_MIN && std::isfinite(ratio)){
    if(!inverse){
      if(represented_r) *represented_r = ratio;
      return std::log(ratio);
    }
    double power = std::pow(ratio, order);
    if(power >= DBL_MIN && std::isfinite(power)){
      if(represented_r) *represented_r = power;
      return std::log(power);
    }
  }
  double ld = log_distance(x, location);
  double out = inverse ? order * (std::log(tau) - 2.0 * ld) :
    2.0 * ld - log_two - std::log(tau);
  return std::isfinite(out) ? out : quiet_nan();
}

double reconstruct(double log_r, double location, double tau, double order,
                   bool inverse, bool negative, double represented_root = -1.0)
{
  if(std::isnan(log_r)) return quiet_nan();
  double log_d = inverse ? (std::log(tau) - log_r / order) / 2.0 :
    (log_two + std::log(tau) + log_r) / 2.0;
  double d = std::exp(log_d);
  if(std::isfinite(represented_root) && represented_root >= DBL_MIN){
    double factor = inverse ? std::pow(represented_root, 0.5 / order) :
      std::sqrt(2.0) * std::sqrt(represented_root);
    if(std::isfinite(factor) && factor >= DBL_MIN){
      d = inverse ? std::sqrt(tau) / factor : std::sqrt(tau) * factor;
    }
  }
  double out = negative ? location - d : location + d;
  return out == location ? quiet_nan() : out;
}

double invmoment_reconstruct(double log_r, double represented_root, double log_lower_p,
                            double leading_gamma_ratio, double location, double tau,
                            double order, double df, bool negative)
{
  if(!std::isfinite(leading_gamma_ratio)){
    return reconstruct(log_r, location, tau, order, true, negative, represented_root);
  }
  // Only the Gamma core's certified leading small-root branch supplies g.
  // Retain the declared df/order and lower log probability: log_r can overflow
  // before its division by order even when the final magnitude is ordinary.
  double log_tau_half = std::log(tau) / 2.0;
  double positive_ratio = -log_lower_p / df;
  double correction = -(leading_gamma_ratio / order) / 2.0;
  if(positive_ratio == infinity && std::isfinite(log_lower_p) && log_lower_p < 0.0 &&
     std::isfinite(df) && df > 0.0 && leading_gamma_ratio <= 0.0 &&
     std::isfinite(correction) && correction >= 0.0 &&
     log_tau_half >= std::log(std::numeric_limits<double>::denorm_min()) / 2.0 &&
     log_tau_half <= std::log(DBL_MAX) / 2.0){
    // The nonnegative correction and bounded finite-tau term cannot offset
    // this positive quotient overflow; the final natural magnitude overflows.
    return negative ? -infinity : infinity;
  }
  double log_d = log_tau_half + positive_ratio + correction;
  if(!std::isfinite(log_d)) return quiet_nan();
  double d = std::exp(log_d);
  double out = negative ? location - d : location + d;
  return out == location ? quiet_nan() : out;
}

double density(double x, double location, double tau, double order, double df, bool inverse)
{
  if(!(inverse ? valid_invmoment_parameters(location, tau, order, df) :
       valid_common_parameters(location, tau, order))) return -infinity;
  if(std::isnan(x)) return quiet_nan();
  if(!std::isfinite(x) || x == location) return -infinity;
  double represented_r;
  double log_r = radial_log(x, location, tau, order, inverse, &represented_r);
  double prefix = gamma::log_prefix(shape_value(order, df, inverse), log_r, represented_r);
  return (inverse ? std::log(order) : 0.0) + prefix - log_distance(x, location);
}

double probability(double q, double location, double tau, double order, double df,
                   bool lower_tail, bool log_p, bool inverse)
{
  if(!(inverse ? valid_invmoment_parameters(location, tau, order, df) :
       valid_common_parameters(location, tau, order)) || std::isnan(q)) return quiet_nan();
  if(!std::isfinite(q)){
    bool zero = (q < 0.0) == lower_tail;
    return log_p ? (zero ? -infinity : 0.0) : (zero ? 0.0 : 1.0);
  }
  if(q == location) return log_p ? -log_two : 0.5;
  double represented_r;
  double log_r = radial_log(q, location, tau, order, inverse, &represented_r);
  gamma::tails tails = gamma::probabilities(shape_value(order, df, inverse), log_r, represented_r);
  double small_log = (inverse ? tails.log_lower : tails.log_upper) - log_two;
  bool small = (q < location) == lower_tail;
  if(log_p) return small ? small_log : gamma::log1mexp(small_log);
  double small_probability = 0.5 * (inverse ? tails.lower : tails.upper);
  return small ? small_probability : 1.0 - small_probability;
}

double quantile(double p, double location, double tau, double order, double df,
                bool lower_tail, bool log_p, bool inverse)
{
  if(!(inverse ? valid_invmoment_parameters(location, tau, order, df) :
       valid_common_parameters(location, tau, order)) || std::isnan(p) ||
     (log_p ? p > 0.0 : p < 0.0 || p > 1.0)) return quiet_nan();
  bool zero = log_p ? p == -infinity : p == 0.0;
  bool unit = log_p ? p == 0.0 : p == 1.0;
  if(zero) return lower_tail ? -infinity : infinity;
  if(unit) return lower_tail ? infinity : -infinity;
  double half = log_p ? -log_two : 0.5;
  if(p == half) return location;
  bool small = p < half;
  bool negative = small == lower_tail;
  double log_tail = (small ? (log_p ? p : std::log(p)) :
    (log_p ? gamma::log1mexp(p) : std::log1p(-p))) + log_two;
  double represented_root;
  double leading_gamma_ratio;
  double log_r = gamma::log_quantile(log_tail, shape_value(order, df, inverse), inverse,
                                    &represented_root, &leading_gamma_ratio);
  if(inverse){
    return invmoment_reconstruct(log_r, represented_root, log_tail, leading_gamma_ratio,
      location, tau, order, df, negative);
  }
  return reconstruct(log_r, location, tau, order, inverse, negative, represented_root);
}

double piece_mass(double lower, double upper, double location, double tau,
                  double order, double df, bool inverse)
{
  if(!(lower < upper)) return -infinity;
  double r1 = radial_log(lower, location, tau, order, inverse);
  double r2 = radial_log(upper, location, tau, order, inverse);
  if(std::isnan(r1) || std::isnan(r2)) return quiet_nan();
  return gamma::log_interval_mass(shape_value(order, df, inverse),
    std::min(r1, r2), std::max(r1, r2)) - log_two;
}

bool candidate_available(double value, double lower, double upper, double location,
                         double tau, double order, double df, bool inverse)
{
  return std::isfinite(value) && value >= lower && value <= upper &&
    std::isfinite(density(value, location, tau, order, df, inverse));
}

bool exact_symmetric_bounds(double lower, double upper, double location)
{
  if(!(std::isfinite(lower) && std::isfinite(upper))) return lower == -infinity && upper == infinity;
  double left = lower * 0.5;
  double right = upper * 0.5;
  double target = location;
  if(left * 2.0 != lower || right * 2.0 != upper){
    left = lower;
    right = upper;
    target = location * 2.0;
  }
  // An error-free two-sum certifies the supplied binary64 midpoint. Rounded
  // distance equality can otherwise falsely turn a nearby interior p into .5.
  double sum = left + right;
  if(!std::isfinite(sum) || !std::isfinite(target)) return false;
  double right_part = sum - left;
  double error = (left - (sum - right_part)) + (right - right_part);
  return sum == target && error == 0.0;
}

}

bool valid_common_parameters(double location, double tau, double order)
{
  return std::isfinite(location) && std::isfinite(tau) && tau > 0.0 && valid_order(order);
}

bool valid_invmoment_parameters(double location, double tau, double order, double df)
{
  return valid_common_parameters(location, tau, order) && std::isfinite(df) && df > 0.0;
}

double moment_mode(double tau, double order)
{
  double scaled_tau = 2.0 * order * tau;
  if(std::isfinite(scaled_tau) && scaled_tau >= DBL_MIN) return std::sqrt(scaled_tau);
  return std::exp((log_two + std::log(order) + std::log(tau)) / 2.0);
}

double invmoment_mode(double tau, double order, double df)
{
  double ratio = 2.0 * order / (df + 1.0);
  if(std::isfinite(ratio) && ratio >= DBL_MIN && std::isfinite(2.0 * order)){
    return std::sqrt(tau) * std::pow(ratio, 1.0 / (2.0 * order));
  }
  return std::exp(std::log(tau) / 2.0 +
    (log_two + std::log(order) - std::log(df + 1.0)) / (2.0 * order));
}

double moment_log_density(double x, double location, double tau, double order)
{
  return density(x, location, tau, order, 0.0, false);
}

double invmoment_log_density(double x, double location, double tau, double order, double df)
{
  return density(x, location, tau, order, df, true);
}

double moment_cdf(double q, double location, double tau, double order, bool lower_tail, bool log_p)
{
  return probability(q, location, tau, order, 0.0, lower_tail, log_p, false);
}

double invmoment_cdf(double q, double location, double tau, double order, double df,
                     bool lower_tail, bool log_p)
{
  return probability(q, location, tau, order, df, lower_tail, log_p, true);
}

double moment_quantile(double p, double location, double tau, double order, bool lower_tail, bool log_p)
{
  return quantile(p, location, tau, order, 0.0, lower_tail, log_p, false);
}

double invmoment_quantile(double p, double location, double tau, double order, double df,
                          bool lower_tail, bool log_p)
{
  return quantile(p, location, tau, order, df, lower_tail, log_p, true);
}

double moment_rng(double u_sign, double u_size, double location, double tau, double order)
{
  u_size = std::min(1.0 - DBL_EPSILON, std::max(DBL_MIN, u_size));
  double represented_root;
  double log_r = gamma::log_quantile(std::log(u_size), order + 0.5, true, &represented_root);
  return reconstruct(log_r, location, tau, order, false, u_sign < 0.5, represented_root);
}

double invmoment_rng(double u_sign, double u_size, double location, double tau, double order, double df)
{
  u_size = std::min(1.0 - DBL_EPSILON, std::max(DBL_MIN, u_size));
  double represented_root;
  double leading_gamma_ratio;
  double log_lower_p = std::log(u_size);
  double log_r = gamma::log_quantile(log_lower_p, (df / order) * 0.5, true,
                                    &represented_root, &leading_gamma_ratio);
  return invmoment_reconstruct(log_r, represented_root, log_lower_p, leading_gamma_ratio,
    location, tau, order, df, u_sign < 0.5);
}

double log_interval_mass(double lower, double upper, double location, double tau,
                         double order, double df, bool inverse, double *log_sign_masses)
{
  if(log_sign_masses){
    log_sign_masses[0] = quiet_nan();
    log_sign_masses[1] = quiet_nan();
  }
  if(std::isnan(lower) || std::isnan(upper) ||
     !(inverse ? valid_invmoment_parameters(location, tau, order, df) :
       valid_common_parameters(location, tau, order))) return quiet_nan();
  if(lower >= upper) return -infinity;
  double negative = piece_mass(lower, std::min(upper, location), location, tau, order, df, inverse);
  double positive = piece_mass(std::max(lower, location), upper, location, tau, order, df, inverse);
  if(log_sign_masses){
    log_sign_masses[0] = negative;
    log_sign_masses[1] = positive;
  }
  return gamma::logaddexp(negative, positive);
}

double truncated_quantile(double p, double lower, double upper, double location,
                          double tau, double order, double df, bool inverse,
                          double const *log_sign_masses)
{
  if(!(inverse ? valid_invmoment_parameters(location, tau, order, df) :
       valid_common_parameters(location, tau, order)) ||
     !(p >= 0.0 && p <= 1.0 && lower < upper)) return quiet_nan();
  if(p == 0.0) return lower;
  if(p == 1.0) return upper;
  double ln = log_sign_masses ? log_sign_masses[0] :
    piece_mass(lower, std::min(upper, location), location, tau, order, df, inverse);
  double lp = log_sign_masses ? log_sign_masses[1] :
    piece_mass(std::max(lower, location), upper, location, tau, order, df, inverse);
  double total = gamma::logaddexp(ln, lp);
  if(!std::isfinite(total)) return quiet_nan();
  // Symmetric declared bounds certify the exact sign boundary; log equality alone does not.
  if(p == 0.5 && lower < location && upper > location &&
     exact_symmetric_bounds(lower, upper, location)) return location;
  double log_u = std::log(p);
  bool negative = log_u < ln - total;
  double relative = negative ? log_u - (ln - total) :
    gamma::logdiffexp(log_u, ln - total) - (lp - total);
  double left = negative ? lower : std::max(lower, location);
  double right = negative ? std::min(upper, location) : upper;
  double rleft = radial_log(left, location, tau, order, inverse);
  double rright = radial_log(right, location, tau, order, inverse);
  bool descending = negative != inverse;
  double radial_probability = descending ? gamma::log1mexp(relative) : relative;
  double root = gamma::interval_log_quantile(radial_probability, shape_value(order, df, inverse),
    std::min(rleft, rright), std::max(rleft, rright));
  double out = reconstruct(root, location, tau, order, inverse, negative);
  return std::isfinite(out) && out > lower && out < upper ? out : quiet_nan();
}

double typical_value(double location, double mode_abs, double lower, double upper,
                     double tau, double order, double df, bool inverse)
{
  double value = location + mode_abs;
  if(candidate_available(value, lower, upper, location, tau, order, df, inverse)) return value;
  value = location - mode_abs;
  if(candidate_available(value, lower, upper, location, tau, order, df, inverse)) return value;
  if(std::isfinite(lower) && std::isfinite(upper) && lower < upper){
    double width = upper - lower;
    value = location <= (lower + upper) / 2.0 ? lower + 0.75 * width : lower + 0.25 * width;
    if(value == location) value = lower + 0.25 * width;
    if(candidate_available(value, lower, upper, location, tau, order, df, inverse)) return value;
  }
  double step = std::isfinite(mode_abs) && mode_abs > 0.0 ? mode_abs :
    std::max(1.0, std::fabs(location) * 0.1);
  value = std::isfinite(lower) ? lower + step : std::isfinite(upper) ? upper - step : location + step;
  if(candidate_available(value, lower, upper, location, tau, order, df, inverse)) return value;
  value = truncated_quantile(0.5, lower, upper, location, tau, order, df, inverse);
  return candidate_available(value, lower, upper, location, tau, order, df, inverse) ? value : quiet_nan();
}

}
}
