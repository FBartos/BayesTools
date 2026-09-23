#include "DBTLKJCPC.h"
#include "../lkj/BTLKJCore.h"

#include <cmath>
#include <JRmath.h>
#include <rng/RNG.h>
#include <util/nainf.h>

namespace jags {
  namespace BayesTools {

    DBTLKJCPC::DBTLKJCPC() : VectorDist("dbt_lkj_cpc", 1)
    {
    }

    bool DBTLKJCPC::checkParameterLength(std::vector<unsigned int> const &len) const
    {
      return len[0] > 0;
    }

    bool DBTLKJCPC::checkParameterValue(std::vector<double const *> const &par,
                                        std::vector<unsigned int> const &len) const
    {
      return bayestools::lkj::valid_alpha(par[0], len[0]);
    }

    bool DBTLKJCPC::checkParameterDiscrete(std::vector<bool> const &mask) const
    {
      return true;
    }

    double DBTLKJCPC::logDensity(double const *x, unsigned int length, PDFType type,
                                 std::vector<double const *> const &par,
                                 std::vector<unsigned int> const &len,
                                 double const *lower, double const *upper) const
    {
      if(length != len[0]){
        return JAGS_NEGINF;
      }

      double log_density = bayestools::lkj::log_density_u(x, par[0], length);
      if(!std::isfinite(log_density)){
        return JAGS_NEGINF;
      }

      return log_density;
    }

    void DBTLKJCPC::randomSample(double *x, unsigned int length,
                                 std::vector<double const *> const &par,
                                 std::vector<unsigned int> const &len,
                                 double const *lower, double const *upper,
                                 RNG *rng) const
    {
      for(unsigned int i = 0; i < length; ++i){
        x[i] = qbeta(rng->uniform(), par[0][i], par[0][i], true, false);
      }
    }

    void DBTLKJCPC::typicalValue(double *x, unsigned int length,
                                 std::vector<double const *> const &par,
                                 std::vector<unsigned int> const &len,
                                 double const *lower, double const *upper) const
    {
      for(unsigned int i = 0; i < length; ++i){
        x[i] = 0.5;
      }
    }

    void DBTLKJCPC::support(double *lower, double *upper, unsigned int length,
                            std::vector<double const *> const &par,
                            std::vector<unsigned int> const &len) const
    {
      for(unsigned int i = 0; i < length; ++i){
        lower[i] = 0.0;
        upper[i] = 1.0;
      }
    }

    unsigned int DBTLKJCPC::length(std::vector<unsigned int> const &len) const
    {
      return len[0];
    }

    bool DBTLKJCPC::isSupportFixed(std::vector<bool> const &fixmask) const
    {
      return true;
    }

    bool DBTLKJCPC::isDiscreteValued(std::vector<bool> const &mask) const
    {
      return false;
    }
  }
}
