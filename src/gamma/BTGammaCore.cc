#include "BTGammaCore.h"
#include "BTGammaRange.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <JRmath.h>

namespace bayestools {
namespace gamma {

namespace {

const double neg_inf = -std::numeric_limits<double>::infinity();
const double pos_inf = std::numeric_limits<double>::infinity();
const double log_half = -std::log(2.0);

double unavailable()
{
  return std::numeric_limits<double>::quiet_NaN();
}

bool valid_log_probability(double x)
{
  return !std::isnan(x) && x <= 0.0;
}

tails missing_tails()
{
  tails out = {unavailable(), unavailable(), unavailable(), unavailable()};
  return out;
}

tails from_lower(double log_p)
{
  tails out = {log_p, log1mexp(log_p), std::exp(log_p), -std::expm1(log_p)};
  return out;
}

}

double log1mexp(double x)
{
  if(!valid_log_probability(x)) return unavailable();
  return x < log_half ? std::log1p(-std::exp(x)) : std::log(-std::expm1(x));
}

double logaddexp(double x, double y)
{
  if(std::isnan(x) || std::isnan(y)) return unavailable();
  if(x == neg_inf) return y;
  if(y == neg_inf) return x;
  double hi = std::max(x, y);
  return hi + std::log1p(std::exp(std::min(x, y) - hi));
}

double logdiffexp(double x, double y)
{
  if(std::isnan(x) || std::isnan(y) || x < y) return unavailable();
  if(y == neg_inf) return x;
  if(x == y) return neg_inf;
  return x + log1mexp(y - x);
}

double log_gamma_ratio(double a)
{
  if(!(std::isfinite(a) && a > 0.0)) return unavailable();
  if(a > 0.5){
    double value = lgammafn(a + 1.0) / a;
    return std::isfinite(value) ? value : unavailable();
  }
  // NIST 5.7.3, factored about 1+a. Retained binary64 coefficients (172).
  const double coefficients[] = {
    0x1.4a34cc4a60fa6p-2, -0x1.13e001a557607p-4,
    0x1.51322ac7d8483p-6, -0x1.e404fc218f5f2p-8,
    0x1.7add6eadb6c30p-9, -0x1.38ac5c2bf8e08p-10,
    0x1.0b36af86396e9p-11, -0x1.d3fd4c76d2fc8p-13,
    0x1.a127b0f17d65ap-14, -0x1.78de5bd7c81efp-15,
    0x1.580dcee66eb02p-16, -0x1.3cbc963ce2243p-17,
    0x1.2597a39f34aacp-18, -0x1.11b2eb7679541p-19,
    0x1.0064cdeb22f0fp-20, -0x1.e2600d93cfd2fp-22,
    0x1.c76bbb3f07a4dp-23, -0x1.af5a6cbbf8a97p-24,
    0x1.99b93c2070b0fp-25, -0x1.862c734df3eacp-26,
    0x1.7469daccfadcdp-27, -0x1.6434a8447aeadp-28,
    0x1.555a877ffd2c3p-29, -0x1.47b1679258d0ep-30,
    0x1.3b15d2b2fc10cp-31, -0x1.2f69a9fabe3e0p-32,
    0x1.24932a337434cp-33, -0x1.1a7c26ec2523cp-34,
    0x1.11116e693ed98p-35, -0x1.08424cbc543d8p-36,
    0x1.000026e3f644fp-37
  };
  double polynomial = coefficients[30];
  for(int i = 29; i >= 0; --i) polynomial = coefficients[i] + a * polynomial;
  return (-std::log1p(a) / a + 0x1.b0ee6072093cep-2) + a * polynomial;
}

tails probabilities(double a, double log_r, double represented_r, double coordinate_error)
{
  if(!(std::isfinite(a) && a > 0.0) || std::isnan(log_r)) return missing_tails();
  if(log_r == neg_inf) return from_lower(neg_inf);
  if(log_r == pos_inf) return from_lower(0.0);
  if(range::certified(a, log_r, coordinate_error)) return {0.0, neg_inf, 1.0, 0.0};
  double r = represented_r >= DBL_MIN && std::isfinite(represented_r) ? represented_r : std::exp(log_r);
  if(a <= 0.5 && log_r <= 0.0){
    double term = -r / (a + 1.0);
    double sum = 0.0;
    for(int n = 1; n <= 64; ++n){
      sum += term;
      if(n < 64) term = term * (-r) / (n + 1.0) * (a + n) / (a + n + 1.0);
    }
    double u = a * sum;
    double h = std::fabs(u) < std::sqrt(DBL_EPSILON) ?
      sum * (1.0 - u / 2.0 + u * u / 3.0 - u * u * u / 4.0) :
      std::log1p(u) / a;
    double kappa = -log_r + log_gamma_ratio(a) - h;
    if(!(std::isfinite(kappa) && kappa > 0.0)) return missing_tails();
    double t = a * kappa;
    tails out = from_lower(-t);
    out.log_upper = t < std::sqrt(DBL_EPSILON) ?
      ((std::log(a) + std::log(kappa)) - t / 2.0) + t * t / 24.0 :
      log1mexp(-t);
    return out;
  }
  if(log_r < std::log(DBL_MIN)){
    double g = log_gamma_ratio(a);
    double lp = a * (log_r - g);
    return std::isfinite(g) && valid_log_probability(lp) ? from_lower(lp) : missing_tails();
  }
  if(!std::isfinite(r)) return missing_tails();
  // JRmath loses finite logarithmic tails for subnormal shapes in this domain.
  if(a < DBL_MIN){
    tails out = missing_tails();
    double bound = (a - 1.0) * log_r - r + std::log(a) - a * log_gamma_ratio(a);
    if(bound < std::log(std::numeric_limits<double>::denorm_min()) - std::log(2.0)){
      out.lower = 1.0;
      out.upper = 0.0;
    }
    return out;
  }
  tails out = {pgamma(r, a, 1.0, true, true), pgamma(r, a, 1.0, false, true),
               unavailable(), unavailable()};
  if(!valid_log_probability(out.log_lower) || !valid_log_probability(out.log_upper) ||
     !std::isfinite(out.log_lower) || !std::isfinite(out.log_upper)) return missing_tails();
  out.lower = out.log_lower < log_half ? std::exp(out.log_lower) : -std::expm1(out.log_upper);
  out.upper = out.log_upper < log_half ? std::exp(out.log_upper) : -std::expm1(out.log_lower);
  return out;
}

double log_prefix(double a, double log_r, double represented_r, double coordinate_error)
{
  if(!(std::isfinite(a) && a > 0.0 && std::isfinite(log_r))) return unavailable();
  if(range::certified(a, log_r, coordinate_error)) return neg_inf;
  double r = represented_r >= DBL_MIN && std::isfinite(represented_r) ? represented_r : std::exp(log_r);
  if(!std::isfinite(r)) return unavailable();
  double out = a <= 0.5 ? std::log(a) + a * (log_r - log_gamma_ratio(a)) - r :
    r >= DBL_MIN ? dgamma(r, a, 1.0, true) + log_r :
    a * log_r - r - lgammafn(a);
  return std::isfinite(out) ? out : unavailable();
}

double log_quantile(double log_p, double a, bool lower_tail, double *represented_root,
                    double *leading_gamma_ratio)
{
  if(represented_root) *represented_root = unavailable();
  if(leading_gamma_ratio) *leading_gamma_ratio = unavailable();
  if(!(std::isfinite(a) && a > 0.0) || !valid_log_probability(log_p)) return unavailable();
  if(log_p == neg_inf) return lower_tail ? neg_inf : pos_inf;
  if(log_p == 0.0) return lower_tail ? pos_inf : neg_inf;
  double root = qgamma(log_p, a, 1.0, lower_tail, true);
  if(std::isfinite(root) && root >= DBL_MIN){
    if(represented_root) *represented_root = root;
    return std::log(root);
  }
  if(root == 0.0 || (root > 0.0 && root < DBL_MIN)){
    double lp = lower_tail ? log_p : log1mexp(log_p);
    double g = log_gamma_ratio(a);
    double leading = lp / a + g;
    if(std::isfinite(lp) && lp < 0.0 && std::isfinite(g) && leading <= -36.0){
      if(leading_gamma_ratio) *leading_gamma_ratio = g;
      return leading;
    }
  }
  return unavailable();
}

double log_interval_mass(double a, double lo, double hi, double lower_error, double upper_error)
{
  if(std::isnan(lo) || std::isnan(hi) || lo > hi) return unavailable();
  if(lo == hi) return neg_inf;
  tails lower = probabilities(a, lo, -1.0, lower_error);
  tails upper = probabilities(a, hi, -1.0, upper_error);
  double from_p = logdiffexp(upper.log_lower, lower.log_lower);
  double from_q = logdiffexp(lower.log_upper, upper.log_upper);
  double out = upper.log_lower <= log_half ? from_p :
    lower.log_upper <= log_half ? from_q :
    log1mexp(logaddexp(lower.log_lower, upper.log_upper));
  if(std::isfinite(out)) return out;
  if(std::isfinite(from_p)) return from_p;
  if(std::isfinite(from_q)) return from_q;
  return unavailable();
}

double interval_log_quantile(double log_p, double a, double lo, double hi,
                             double lower_error, double upper_error)
{
  if(!valid_log_probability(log_p) || !(lo < hi)) return unavailable();
  if(log_p == neg_inf) return lo;
  if(log_p == 0.0) return hi;
  double mass = log_interval_mass(a, lo, hi, lower_error, upper_error);
  if(!std::isfinite(mass)) return unavailable();
  tails lower = probabilities(a, lo, -1.0, lower_error);
  tails upper = probabilities(a, hi, -1.0, upper_error);
  double lp = logaddexp(lower.log_lower, log_p + mass);
  double lq = logaddexp(upper.log_upper, log1mexp(log_p) + mass);
  double out = !std::isnan(lp) && (std::isnan(lq) || lp <= lq) ?
    log_quantile(lp, a, true) : log_quantile(lq, a, false);
  return out > lo && out < hi ? out : unavailable();
}

}
}
