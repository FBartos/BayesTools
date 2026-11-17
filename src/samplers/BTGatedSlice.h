#ifndef BTGATEDSLICE_H_
#define BTGATEDSLICE_H_

#include <sampler/SamplerFactory.h>

namespace jags {
  namespace BayesTools {

    // Sampler "BayesTools::GatedSlice" for a scalar continuous node declared by
    //   BayesTools_gate_<k> ~ dbt_gate(<node>, <active>)
    // While <active> is 0 the node is drawn directly from its own distribution
    // (its likelihood children are constant in it by declaration); otherwise it
    // gets the univariate slice update of base::RealSlicer.
    class GatedSliceFactory : public SamplerFactory
    {
    public:
      std::vector<Sampler *> makeSamplers(std::list<StochasticNode *> const &nodes,
                                          Graph const &graph) const;
      std::string name() const;
    };

    // Process-lifetime instance; the module never deletes it.
    GatedSliceFactory *gated_slice_factory();
  }
}

#endif
