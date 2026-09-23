#ifndef BTINVGAMMACORE_H_
#define BTINVGAMMACORE_H_

namespace bayestools {
namespace invgamma {

bool valid_parameters(double shape, double scale);
bool valid_value(double x);

double mode(double shape, double scale);
double log_density(double x, double shape, double scale);
double cdf(double q, double shape, double scale, bool lower_tail, bool log_p);
double quantile(double p, double shape, double scale, bool lower_tail, bool log_p);
double rng(double u, double shape, double scale);
double typical_value(double shape, double scale, double lower, double upper);

}
}

#endif
