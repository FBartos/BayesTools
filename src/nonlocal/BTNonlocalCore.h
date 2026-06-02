#ifndef BTNONLOCALCORE_H_
#define BTNONLOCALCORE_H_

namespace bayestools {
namespace nonlocal {

bool valid_order(double order);
bool valid_moment_parameters(double location, double tau, double order);
bool valid_invmoment_parameters(double location, double tau, double order, double df);

double moment_mode(double tau, double order);
double invmoment_mode(double tau, double order, double df);
double typical_value(double location, double mode, double const *lower, double const *upper);

double dmoment(double x, double location, double tau, double order, bool give_log);
double pmoment(double q, double location, double tau, double order, bool lower_tail, bool log_p);
double qmoment(double p, double location, double tau, double order, bool lower_tail, bool log_p);
double rmoment(double sign_u, double magnitude_u, double location, double tau, double order);

double dinvmoment(double x, double location, double tau, double order, double df, bool give_log);
double pinvmoment(double q, double location, double tau, double order, double df, bool lower_tail, bool log_p);
double qinvmoment(double p, double location, double tau, double order, double df, bool lower_tail, bool log_p);
double rinvmoment(double sign_u, double magnitude_u, double location, double tau, double order, double df);

}
}

#endif
