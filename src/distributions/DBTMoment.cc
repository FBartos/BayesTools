#include "DBTMoment.h"
#include "../nonlocal/BTNonlocalCore.h"

#include <rng/RNG.h>

namespace jags {
  namespace BayesTools {

    DBTMoment::DBTMoment() : RScalarDist("dbt_moment", 3, DIST_UNBOUNDED, false)
    {
    }

    double DBTMoment::d(double x, PDFType type,
                        std::vector<double const *> const &par,
                        bool give_log) const
    {
      return bayestools::nonlocal::dmoment(x, par[0][0], par[1][0], par[2][0], give_log);
    }

    double DBTMoment::p(double x, std::vector<double const *> const &par,
                        bool lower, bool give_log) const
    {
      return bayestools::nonlocal::pmoment(x, par[0][0], par[1][0], par[2][0], lower, give_log);
    }

    double DBTMoment::q(double p, std::vector<double const *> const &par,
                        bool lower, bool log_p) const
    {
      return bayestools::nonlocal::qmoment(p, par[0][0], par[1][0], par[2][0], lower, log_p);
    }

    double DBTMoment::r(std::vector<double const *> const &par, RNG *rng) const
    {
      return bayestools::nonlocal::rmoment(rng->uniform(), rng->uniform(), par[0][0], par[1][0], par[2][0]);
    }

    double DBTMoment::typicalValue(std::vector<double const *> const &par,
                                   double const *lower, double const *upper) const
    {
      double mode = bayestools::nonlocal::moment_mode(par[1][0], par[2][0]);
      return bayestools::nonlocal::typical_value(par[0][0], mode, lower, upper);
    }

    bool DBTMoment::checkParameterValue(std::vector<double const *> const &par) const
    {
      return bayestools::nonlocal::valid_moment_parameters(par[0][0], par[1][0], par[2][0]);
    }

    bool DBTMoment::isLocationParameter(unsigned int index) const
    {
      return index == 0;
    }
  }
}
