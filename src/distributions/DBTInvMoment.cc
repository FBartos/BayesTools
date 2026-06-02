#include "DBTInvMoment.h"
#include "../nonlocal/BTNonlocalCore.h"

#include <rng/RNG.h>

namespace jags {
  namespace BayesTools {

    DBTInvMoment::DBTInvMoment() : RScalarDist("dbt_invmoment", 4, DIST_UNBOUNDED, false)
    {
    }

    double DBTInvMoment::d(double x, PDFType type,
                           std::vector<double const *> const &par,
                           bool give_log) const
    {
      return bayestools::nonlocal::dinvmoment(x, par[0][0], par[1][0], par[2][0], par[3][0], give_log);
    }

    double DBTInvMoment::p(double x, std::vector<double const *> const &par,
                           bool lower, bool give_log) const
    {
      return bayestools::nonlocal::pinvmoment(x, par[0][0], par[1][0], par[2][0], par[3][0], lower, give_log);
    }

    double DBTInvMoment::q(double p, std::vector<double const *> const &par,
                           bool lower, bool log_p) const
    {
      return bayestools::nonlocal::qinvmoment(p, par[0][0], par[1][0], par[2][0], par[3][0], lower, log_p);
    }

    double DBTInvMoment::r(std::vector<double const *> const &par, RNG *rng) const
    {
      return bayestools::nonlocal::rinvmoment(rng->uniform(), rng->uniform(), par[0][0], par[1][0], par[2][0], par[3][0]);
    }

    double DBTInvMoment::typicalValue(std::vector<double const *> const &par,
                                      double const *lower, double const *upper) const
    {
      double mode = bayestools::nonlocal::invmoment_mode(par[1][0], par[2][0], par[3][0]);
      return bayestools::nonlocal::typical_value(par[0][0], mode, lower, upper);
    }

    bool DBTInvMoment::checkParameterValue(std::vector<double const *> const &par) const
    {
      return bayestools::nonlocal::valid_invmoment_parameters(par[0][0], par[1][0], par[2][0], par[3][0]);
    }

    bool DBTInvMoment::isLocationParameter(unsigned int index) const
    {
      return index == 0;
    }
  }
}
