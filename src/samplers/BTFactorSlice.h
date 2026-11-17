#ifndef BTFACTORSLICE_H_
#define BTFACTORSLICE_H_

#include <sampler/SamplerFactory.h>

namespace jags {
  namespace BayesTools {

    // Sampler "BayesTools::FactorSlice" for a block of scalar continuous nodes
    // with fixed support declared by
    //   BayesTools_block_<k> ~ dbt_block(<node>, <id>)
    // Nodes that share an id form one block (at least two nodes; a node with a
    // dbt_gate marker never joins a block). One update performs a univariate
    // slice update (stepping out with a randomized limit, then shrinkage) along
    // each of a set of orthonormal directions on an unconstrained scale; the
    // directions and widths are learned from the covariance during adaptation
    // and frozen afterwards.
    class FactorSliceFactory : public SamplerFactory
    {
    public:
      std::vector<Sampler *> makeSamplers(std::list<StochasticNode *> const &nodes,
                                          Graph const &graph) const;
      std::string name() const;
    };

    // Process-lifetime instance; the module never deletes it.
    FactorSliceFactory *factor_slice_factory();
  }
}

#endif
