#ifndef BTFINITEREADOUT_H_
#define BTFINITEREADOUT_H_

#include <memory>

namespace jags {
  class Module;
  class MonitorFactory;
  class SamplerFactory;

  namespace BayesTools {

    // Owned by the module that registers the two factories with JAGS.
    class FiniteReadout {
      std::unique_ptr<SamplerFactory> _sampler;
      std::unique_ptr<MonitorFactory> _monitor;
    public:
      explicit FiniteReadout(Module &module);
      ~FiniteReadout();
      void prepare();
    };
  }
}

#endif
