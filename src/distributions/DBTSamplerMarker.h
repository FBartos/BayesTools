#ifndef DBTSAMPLERMARKER_H_
#define DBTSAMPLERMARKER_H_

#include <distribution/ScalarDist.h>

#include <string>

namespace jags {
  namespace BayesTools {

    // Sampler declarations are observed marker nodes with a degenerate
    // distribution at zero:
    //   BayesTools_gate_<k>  ~ dbt_gate(<node>, <active>)
    //   BayesTools_block_<k> ~ dbt_block(<node>, <id>)
    // The log density is 0 at the value 0 and minus infinity elsewhere and
    // never reads the parameters, so a marker adds nothing to any full
    // conditional, deviance, or log-likelihood. The sampler factories in
    // src/samplers/ read the parameters as declarations.
    class DBTSamplerMarker : public ScalarDist
    {
    public:
      explicit DBTSamplerMarker(std::string const &name);
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
      bool isSupportFixed(std::vector<bool> const &fixmask) const;
      double l(std::vector<double const *> const &parameters) const;
      double u(std::vector<double const *> const &parameters) const;
    };

    // Distribution names shared by the module registration and the factories.
    extern char const *const gate_marker_name;
    extern char const *const block_marker_name;
  }
}

#endif
