#ifndef BTNONLOCALCORE_H_
#define BTNONLOCALCORE_H_

namespace bayestools {
namespace nonlocal {

bool valid_common_parameters(double location, double tau, double order);
bool valid_invmoment_parameters(double location, double tau, double order, double df);

double moment_mode(double tau, double order);
double invmoment_mode(double tau, double order, double df);

double moment_log_density(double x, double location, double tau, double order);
double invmoment_log_density(double x, double location, double tau, double order, double df);

double moment_cdf(double q, double location, double tau, double order,
                  bool lower_tail, bool log_p);
double invmoment_cdf(double q, double location, double tau, double order,
                     double df, bool lower_tail, bool log_p);

double moment_quantile(double p, double location, double tau, double order,
                       bool lower_tail, bool log_p);
double invmoment_quantile(double p, double location, double tau, double order,
                          double df, bool lower_tail, bool log_p);

double moment_rng(double u_sign, double u_size, double location, double tau,
                  double order);
double invmoment_rng(double u_sign, double u_size, double location, double tau,
                     double order, double df);

double typical_value(double location, double mode_abs, double lower,
                     double upper);

}
}

#endif
