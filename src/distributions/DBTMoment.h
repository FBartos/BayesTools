#ifndef DBTMOMENT_H_
#define DBTMOMENT_H_

#include <distribution/RScalarDist.h>

namespace jags {
  namespace BayesTools {

    class DBTMoment : public RScalarDist
    {
    public:
      DBTMoment();

      double d(double x, PDFType type, std::vector<double const *> const &parameters,
               bool give_log) const;
      double p(double x, std::vector<double const *> const &parameters,
               bool lower, bool give_log) const;
      double q(double p, std::vector<double const *> const &parameters,
               bool lower, bool log_p) const;
      double r(std::vector<double const *> const &parameters, RNG *rng) const;
      double typicalValue(std::vector<double const *> const &parameters,
                          double const *lower, double const *upper) const;
      bool checkParameterValue(std::vector<double const *> const &parameters) const;
      bool isLocationParameter(unsigned int index) const;
    };
  }
}

#endif
