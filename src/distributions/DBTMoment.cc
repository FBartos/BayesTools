#include "DBTMoment.h"
#include "../nonlocal/BTNonlocalCore.h"

#include <cmath>
#include <rng/RNG.h>
#include <util/nainf.h>

namespace jags {
  namespace BayesTools {

    DBTMoment::DBTMoment() : ScalarDist("dbt_moment", 3, DIST_UNBOUNDED)
    {
    }

    bool DBTMoment::checkParameterValue(std::vector<double const *> const &par) const
    {
      return bayestools::nonlocal::valid_common_parameters(*par[0], *par[1], *par[2]);
    }

    bool DBTMoment::checkParameterDiscrete(std::vector<bool> const &mask) const
    {
      return mask[2];
    }

    bool DBTMoment::isDiscreteValued(std::vector<bool> const &mask) const
    {
      return false;
    }

    double DBTMoment::logDensity(double x, PDFType type,
                                 std::vector<double const *> const &par,
                                 double const *lower, double const *upper) const
    {
      double out = bayestools::nonlocal::moment_log_density(x, *par[0], *par[1], *par[2]);
      if(!std::isfinite(out)){
        return JAGS_NEGINF;
      }
      if(lower || upper){
        double l = lower ? *lower : JAGS_NEGINF;
        double u = upper ? *upper : JAGS_POSINF;
        if(x < l || x > u){
          return JAGS_NEGINF;
        }

        double mass = bayestools::nonlocal::moment_cdf(
          u, *par[0], *par[1], *par[2], true, false
        ) - bayestools::nonlocal::moment_cdf(
          l, *par[0], *par[1], *par[2], true, false
        );
        if(!(std::isfinite(mass) && mass > 0.0)){
          mass = bayestools::nonlocal::moment_cdf(
            l, *par[0], *par[1], *par[2], false, false
          ) - bayestools::nonlocal::moment_cdf(
            u, *par[0], *par[1], *par[2], false, false
          );
        }
        if(!(std::isfinite(mass) && mass > 0.0)){
          return JAGS_NEGINF;
        }
        out -= std::log(mass);
      }
      return std::isfinite(out) ? out : JAGS_NEGINF;
    }

    double DBTMoment::randomSample(std::vector<double const *> const &par,
                                   double const *lower, double const *upper,
                                   RNG *rng) const
    {
      if(lower || upper){
        double l = lower ? *lower : JAGS_NEGINF;
        double u = upper ? *upper : JAGS_POSINF;
        double cdf_l = bayestools::nonlocal::moment_cdf(
          l, *par[0], *par[1], *par[2], true, false
        );
        double cdf_u = bayestools::nonlocal::moment_cdf(
          u, *par[0], *par[1], *par[2], true, false
        );
        double mass = cdf_u - cdf_l;
        if(std::isfinite(mass) && mass > 0.0){
          return bayestools::nonlocal::moment_quantile(
            cdf_l + rng->uniform() * mass,
            *par[0], *par[1], *par[2], true, false
          );
        }

        double surv_l = bayestools::nonlocal::moment_cdf(
          l, *par[0], *par[1], *par[2], false, false
        );
        double surv_u = bayestools::nonlocal::moment_cdf(
          u, *par[0], *par[1], *par[2], false, false
        );
        mass = surv_l - surv_u;
        if(std::isfinite(mass) && mass > 0.0){
          return bayestools::nonlocal::moment_quantile(
            surv_u + rng->uniform() * mass,
            *par[0], *par[1], *par[2], false, false
          );
        }
      }
      return bayestools::nonlocal::moment_rng(
        rng->uniform(), rng->uniform(), *par[0], *par[1], *par[2]
      );
    }

    double DBTMoment::typicalValue(std::vector<double const *> const &par,
                                   double const *lower, double const *upper) const
    {
      double l = lower ? *lower : JAGS_NEGINF;
      double u = upper ? *upper : JAGS_POSINF;
      double mode = bayestools::nonlocal::moment_mode(*par[1], *par[2]);
      return bayestools::nonlocal::typical_value(*par[0], mode, l, u);
    }

    bool DBTMoment::canBound() const
    {
      return true;
    }
  }
}
