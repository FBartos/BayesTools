#ifndef BTGAMMACORE_H_
#define BTGAMMACORE_H_

namespace bayestools {
namespace gamma {

struct tails {
  double log_lower;
  double log_upper;
  double lower;
  double upper;
};

double logaddexp(double x, double y);
double logdiffexp(double x, double y);
double log1mexp(double x);
double log_gamma_ratio(double shape);
tails probabilities(double shape, double log_r, double represented_r = -1.0);
double log_prefix(double shape, double log_r, double represented_r = -1.0);
double log_quantile(double log_p, double shape, bool lower_tail,
                    double *represented_root = 0, double *leading_gamma_ratio = 0);
double log_interval_mass(double shape, double log_lower, double log_upper);
double interval_log_quantile(double log_p, double shape,
                             double log_lower, double log_upper);

}
}

#endif
