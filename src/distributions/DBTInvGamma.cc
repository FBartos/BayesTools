#include "DBTInvGamma.h"
#include "../invgamma/BTInvGammaCore.h"

#include <cmath>
#include <rng/RNG.h>
#include <util/nainf.h>
#include <module/ModuleError.h>

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
      double l = lower ? *lower : 0.0;
      double u = upper ? *upper : JAGS_POSINF;
      if(x < l || x > u || x <= 0.0 || !std::isfinite(x)) return JAGS_NEGINF;
      double log_mass = lower || upper ? bayestools::invgamma::log_interval_mass(*par[0], *par[1], l, u) : 0.0;
      if(!std::isfinite(log_mass)){
        throwDistError(this, "Prior normalization is numerically unavailable.");
      }
      double out = bayestools::invgamma::log_density(x, *par[0], *par[1]);
      if(!std::isfinite(out) || !std::isfinite(out - log_mass)){
        throwDistError(this, "Prior log density is numerically unavailable.");
      }
      return out - log_mass;
    }

    double DBTInvGamma::randomSample(std::vector<double const *> const &par,
                                     double const *lower, double const *upper,
                                     RNG *rng) const
    {
      double l = lower ? *lower : 0.0;
      double u = upper ? *upper : JAGS_POSINF;
      double out;
      if(lower || upper){
        if(!std::isfinite(bayestools::invgamma::log_interval_mass(*par[0], *par[1], l, u))){
          throwDistError(this, "Prior normalization is numerically unavailable.");
        }
        out = bayestools::invgamma::truncated_quantile(rng->uniform(), *par[0], *par[1], l, u);
      }else{
        out = bayestools::invgamma::rng(rng->uniform(), *par[0], *par[1]);
      }
      if(!std::isfinite(out) || out <= 0.0 || out < l || out > u){
        throwDistError(this, "Prior sampling is numerically unavailable.");
      }
      return out;
    }

    double DBTInvGamma::typicalValue(std::vector<double const *> const &par,
                                     double const *lower, double const *upper) const
    {
      double l = lower ? *lower : 0.0;
      double u = upper ? *upper : JAGS_POSINF;
      if((lower || upper) && !std::isfinite(bayestools::invgamma::log_interval_mass(*par[0], *par[1], l, u))){
        throwDistError(this, "Prior normalization is numerically unavailable.");
      }
      double out = bayestools::invgamma::typical_value(*par[0], *par[1], l, u);
      if(!std::isfinite(out)){
        throwDistError(this, "Prior initialization is numerically unavailable.");
      }
      return out;
    }

    bool DBTInvGamma::canBound() const
    {
      return true;
    }
  }
}
