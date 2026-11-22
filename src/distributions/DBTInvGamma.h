#ifndef DBTINVGAMMA_H_
#define DBTINVGAMMA_H_

#include <distribution/ScalarDist.h>

namespace jags {
  namespace BayesTools {

    class DBTInvGamma : public ScalarDist
    {
    public:
      DBTInvGamma();
      double logDensity(double x, PDFType type,
                        std::vector<double const *> const &parameters,
                        double const *lower, double const *upper) const override;
      double randomSample(std::vector<double const *> const &parameters,
                          double const *lower, double const *upper,
                          RNG *rng) const override;
      bool checkParameterValue(std::vector<double const *> const &parameters) const override;
      bool checkParameterDiscrete(std::vector<bool> const &mask) const override;
      bool isDiscreteValued(std::vector<bool> const &mask) const override;
      bool canBound() const override;
    };
  }
}

#endif
