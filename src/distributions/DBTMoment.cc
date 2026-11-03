#include "DBTMoment.h"
#include "../nonlocal/BTNonlocalCore.h"

#include <cmath>
#include <rng/RNG.h>
#include <util/nainf.h>
#include <module/ModuleError.h>

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
      double l = lower ? *lower : JAGS_NEGINF;
      double u = upper ? *upper : JAGS_POSINF;
      if(x < l || x > u || x == *par[0] || !std::isfinite(x)) return JAGS_NEGINF;
      double log_mass = lower || upper ? bayestools::nonlocal::log_interval_mass(l, u, *par[0], *par[1], *par[2], 0.0, false) : 0.0;
      if(!std::isfinite(log_mass)){
        throwDistError(this, "Prior normalization is numerically unavailable.");
      }
      double out = bayestools::nonlocal::moment_log_density(x, *par[0], *par[1], *par[2]);
      if(!std::isfinite(out) || !std::isfinite(out - log_mass)){
        throwDistError(this, "Prior log density is numerically unavailable.");
      }
      return out - log_mass;
    }

    double DBTMoment::randomSample(std::vector<double const *> const &par,
                                     double const *lower, double const *upper,
                                     RNG *rng) const
    {
      double l = lower ? *lower : JAGS_NEGINF;
      double u = upper ? *upper : JAGS_POSINF;
      double out;
      if(lower || upper){
        if(!std::isfinite(bayestools::nonlocal::log_interval_mass(l, u, *par[0], *par[1], *par[2], 0.0, false))){
          throwDistError(this, "Prior normalization is numerically unavailable.");
        }
        out = bayestools::nonlocal::truncated_quantile(rng->uniform(), l, u, *par[0], *par[1], *par[2], 0.0, false);
      }else{
        double u_size = rng->uniform();
        double u_sign = rng->uniform();
        out = bayestools::nonlocal::moment_rng(u_sign, u_size, *par[0], *par[1], *par[2]);
      }
      if(!std::isfinite(out) || out == *par[0] || out < l || out > u){
        throwDistError(this, "Prior sampling is numerically unavailable.");
      }
      return out;
    }

    double DBTMoment::typicalValue(std::vector<double const *> const &par,
                                     double const *lower, double const *upper) const
    {
      double l = lower ? *lower : JAGS_NEGINF;
      double u = upper ? *upper : JAGS_POSINF;
      if((lower || upper) && !std::isfinite(bayestools::nonlocal::log_interval_mass(l, u, *par[0], *par[1], *par[2], 0.0, false))){
        throwDistError(this, "Prior normalization is numerically unavailable.");
      }
      double out = bayestools::nonlocal::typical_value(*par[0], bayestools::nonlocal::moment_mode(*par[1], *par[2]), l, u, *par[1], *par[2], 0.0, false);
      if(!std::isfinite(out)){
        throwDistError(this, "Prior initialization is numerically unavailable.");
      }
      return out;
    }

    bool DBTMoment::canBound() const
    {
      return true;
    }
  }
}
