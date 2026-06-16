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
                        double const *lower, double const *upper) const;
      double randomSample(std::vector<double const *> const &parameters,
                          double const *lower, double const *upper,
                          RNG *rng) const;
      double typicalValue(std::vector<double const *> const &parameters,
                          double const *lower, double const *upper) const;
      bool checkParameterValue(std::vector<double const *> const &parameters) const;
      bool checkParameterDiscrete(std::vector<bool> const &mask) const;
      bool isDiscreteValued(std::vector<bool> const &mask) const;
      bool canBound() const;
    };
  }
}

#endif
