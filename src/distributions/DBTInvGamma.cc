#include "DBTInvGamma.h"
#include "../invgamma/BTInvGammaCore.h"

#include <cmath>
#include <rng/RNG.h>
#include <util/nainf.h>

namespace jags {
  namespace BayesTools {

    DBTInvGamma::DBTInvGamma() : ScalarDist("dbt_invgamma", 2, DIST_POSITIVE)
    {
    }

    bool DBTInvGamma::checkParameterValue(std::vector<double const *> const &par) const
    {
      return bayestools::invgamma::valid_parameters(*par[0], *par[1]);
    }

    bool DBTInvGamma::checkParameterDiscrete(std::vector<bool> const &mask) const
    {
      return true;
    }

    bool DBTInvGamma::isDiscreteValued(std::vector<bool> const &mask) const
    {
      return false;
    }

    double DBTInvGamma::logDensity(double x, PDFType type,
                                   std::vector<double const *> const &par,
                                   double const *lower, double const *upper) const
    {
      double out = bayestools::invgamma::log_density(x, *par[0], *par[1]);
      if(!std::isfinite(out)){
        return JAGS_NEGINF;
      }
      if(lower || upper){
        double l = lower ? *lower : 0.0;
        double u = upper ? *upper : JAGS_POSINF;
        if(x < l || x > u){
          return JAGS_NEGINF;
        }

        double mass = bayestools::invgamma::cdf(
          u, *par[0], *par[1], true, false
        ) - bayestools::invgamma::cdf(
          l, *par[0], *par[1], true, false
        );
        if(!(std::isfinite(mass) && mass > 0.0)){
          mass = bayestools::invgamma::cdf(
            l, *par[0], *par[1], false, false
          ) - bayestools::invgamma::cdf(
            u, *par[0], *par[1], false, false
          );
        }
        if(!(std::isfinite(mass) && mass > 0.0)){
          return JAGS_NEGINF;
        }
        out -= std::log(mass);
      }
      return std::isfinite(out) ? out : JAGS_NEGINF;
    }

    double DBTInvGamma::randomSample(std::vector<double const *> const &par,
                                     double const *lower, double const *upper,
                                     RNG *rng) const
    {
      if(lower || upper){
        double l = lower ? *lower : 0.0;
        double u = upper ? *upper : JAGS_POSINF;
        double cdf_l = bayestools::invgamma::cdf(
          l, *par[0], *par[1], true, false
        );
        double cdf_u = bayestools::invgamma::cdf(
          u, *par[0], *par[1], true, false
        );
        double mass = cdf_u - cdf_l;
        if(std::isfinite(mass) && mass > 0.0){
          return bayestools::invgamma::quantile(
            cdf_l + rng->uniform() * mass,
            *par[0], *par[1], true, false
          );
        }

        double surv_l = bayestools::invgamma::cdf(
          l, *par[0], *par[1], false, false
        );
        double surv_u = bayestools::invgamma::cdf(
          u, *par[0], *par[1], false, false
        );
        mass = surv_l - surv_u;
        if(std::isfinite(mass) && mass > 0.0){
          return bayestools::invgamma::quantile(
            surv_u + rng->uniform() * mass,
            *par[0], *par[1], false, false
          );
        }
      }
      return bayestools::invgamma::rng(rng->uniform(), *par[0], *par[1]);
    }

    double DBTInvGamma::typicalValue(std::vector<double const *> const &par,
                                     double const *lower, double const *upper) const
    {
      double l = lower ? *lower : 0.0;
      double u = upper ? *upper : JAGS_POSINF;
      return bayestools::invgamma::typical_value(*par[0], *par[1], l, u);
    }

    bool DBTInvGamma::canBound() const
    {
      return true;
    }
  }
}
