#include "BTSamplerDeclarations.h"

#include <BayesTools/sampler_declarations.h>

namespace jags {
  namespace BayesTools {

    bool has_gate_marker(StochasticNode *node)
    {
      return sampler_declarations::has_gate_marker(node);
    }

    Node const *gate_activity(StochasticNode *node)
    {
      return sampler_declarations::gate_activity(node);
    }

    unsigned int block_id(StochasticNode *node)
    {
      return sampler_declarations::block_id(node);
    }
  }
}
