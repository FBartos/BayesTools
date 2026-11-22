#include "DBTSamplerMarker.h"

#include <BayesTools/sampler_declarations.h>

#include <util/nainf.h>

namespace jags {
  namespace BayesTools {

    char const *const gate_marker_name = sampler_declarations::gate_marker_name();
    char const *const block_marker_name = sampler_declarations::block_marker_name();

    DBTSamplerMarker::DBTSamplerMarker(std::string const &name)
      : ScalarDist(name, 2, DIST_SPECIAL)
    {
    }

    double DBTSamplerMarker::logDensity(double x, PDFType type,
                                        std::vector<double const *> const &parameters,
                                        double const *lower, double const *upper) const
    {
      return x == 0.0 ? 0.0 : JAGS_NEGINF;
    }

    double DBTSamplerMarker::randomSample(std::vector<double const *> const &parameters,
                                          double const *lower, double const *upper,
                                          RNG *rng) const
    {
      return 0.0;
    }

    bool DBTSamplerMarker::checkParameterValue(std::vector<double const *> const &parameters) const
    {
      return true;
    }

    bool DBTSamplerMarker::checkParameterDiscrete(std::vector<bool> const &mask) const
    {
      return true;
    }

    bool DBTSamplerMarker::isDiscreteValued(std::vector<bool> const &mask) const
    {
      return true;
    }

    bool DBTSamplerMarker::isSupportFixed(std::vector<bool> const &fixmask) const
    {
      return true;
    }

    double DBTSamplerMarker::l(std::vector<double const *> const &parameters) const
    {
      return 0.0;
    }

    double DBTSamplerMarker::u(std::vector<double const *> const &parameters) const
    {
      return 0.0;
    }
  }
}
