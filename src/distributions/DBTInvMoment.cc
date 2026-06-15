#include "DBTInvMoment.h"
#include "../nonlocal/BTNonlocalCore.h"

#include <cmath>
#include <rng/RNG.h>
#include <util/nainf.h>

namespace jags {
  namespace BayesTools {

    DBTInvMoment::DBTInvMoment() : ScalarDist("dbt_invmoment", 4, DIST_UNBOUNDED)
    {
    }

    bool DBTInvMoment::checkParameterValue(std::vector<double const *> const &par) const
    {
      return bayestools::nonlocal::valid_invmoment_parameters(*par[0], *par[1], *par[2], *par[3]);
    }

    bool DBTInvMoment::checkParameterDiscrete(std::vector<bool> const &mask) const
    {
      return mask[2];
    }

    bool DBTInvMoment::isDiscreteValued(std::vector<bool> const &mask) const
    {
      return false;
    }

    double DBTInvMoment::logDensity(double x, PDFType type,
                                    std::vector<double const *> const &par,
                                    double const *lower, double const *upper) const
    {
      double out = bayestools::nonlocal::invmoment_log_density(
        x, *par[0], *par[1], *par[2], *par[3]
      );
      if(!std::isfinite(out)){
        return JAGS_NEGINF;
      }
      if(lower || upper){
        double l = lower ? *lower : JAGS_NEGINF;
        double u = upper ? *upper : JAGS_POSINF;
        if(x < l || x > u){
          return JAGS_NEGINF;
        }

        double mass = bayestools::nonlocal::invmoment_cdf(
          u, *par[0], *par[1], *par[2], *par[3], true, false
        ) - bayestools::nonlocal::invmoment_cdf(
          l, *par[0], *par[1], *par[2], *par[3], true, false
        );
        if(!(std::isfinite(mass) && mass > 0.0)){
          mass = bayestools::nonlocal::invmoment_cdf(
            l, *par[0], *par[1], *par[2], *par[3], false, false
          ) - bayestools::nonlocal::invmoment_cdf(
            u, *par[0], *par[1], *par[2], *par[3], false, false
          );
        }
        if(!(std::isfinite(mass) && mass > 0.0)){
          return JAGS_NEGINF;
        }
        out -= std::log(mass);
      }
      return std::isfinite(out) ? out : JAGS_NEGINF;
    }

    double DBTInvMoment::randomSample(std::vector<double const *> const &par,
                                      double const *lower, double const *upper,
                                      RNG *rng) const
    {
      if(lower || upper){
        double l = lower ? *lower : JAGS_NEGINF;
        double u = upper ? *upper : JAGS_POSINF;
        double cdf_l = bayestools::nonlocal::invmoment_cdf(
          l, *par[0], *par[1], *par[2], *par[3], true, false
        );
        double cdf_u = bayestools::nonlocal::invmoment_cdf(
          u, *par[0], *par[1], *par[2], *par[3], true, false
        );
        double mass = cdf_u - cdf_l;
        if(std::isfinite(mass) && mass > 0.0){
          return bayestools::nonlocal::invmoment_quantile(
            cdf_l + rng->uniform() * mass,
            *par[0], *par[1], *par[2], *par[3], true, false
          );
        }

        double surv_l = bayestools::nonlocal::invmoment_cdf(
          l, *par[0], *par[1], *par[2], *par[3], false, false
        );
        double surv_u = bayestools::nonlocal::invmoment_cdf(
          u, *par[0], *par[1], *par[2], *par[3], false, false
        );
        mass = surv_l - surv_u;
        if(std::isfinite(mass) && mass > 0.0){
          return bayestools::nonlocal::invmoment_quantile(
            surv_u + rng->uniform() * mass,
            *par[0], *par[1], *par[2], *par[3], false, false
          );
        }
      }
      return bayestools::nonlocal::invmoment_rng(
        rng->uniform(), rng->uniform(), *par[0], *par[1], *par[2], *par[3]
      );
    }

    double DBTInvMoment::typicalValue(std::vector<double const *> const &par,
                                      double const *lower, double const *upper) const
    {
      double l = lower ? *lower : JAGS_NEGINF;
      double u = upper ? *upper : JAGS_POSINF;
      double mode = bayestools::nonlocal::invmoment_mode(*par[1], *par[2], *par[3]);
      return bayestools::nonlocal::typical_value(*par[0], mode, l, u);
    }

    bool DBTInvMoment::canBound() const
    {
      return true;
    }
  }
}
