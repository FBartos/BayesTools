#ifndef DBTINVMOMENT_H_
#define DBTINVMOMENT_H_

#include <distribution/ScalarDist.h>

namespace jags {
  namespace BayesTools {

    class DBTInvMoment : public ScalarDist
    {
    public:
      DBTInvMoment();
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
