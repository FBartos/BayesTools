#ifndef BTSAMPLERDECLARATIONS_H_
#define BTSAMPLERDECLARATIONS_H_

namespace jags {

  class Node;
  class StochasticNode;

  namespace BayesTools {

    // Reading of the marker declarations (src/distributions/DBTSamplerMarker.h).
    // A marker declares the node it is a direct stochastic child of only when it
    // is observed with value 0 and its first parent is that node itself (pointer
    // identity, so a marker on a deterministic expression declares nothing).

    // Whether some dbt_gate marker names `node` as its first parent, valid or not.
    bool has_gate_marker(StochasticNode *node);

    // The <active> node of the single valid dbt_gate marker of `node`, or null
    // when there is none, more than one, or the <active> node is a deterministic
    // function of `node` (the update would change its own gate).
    Node const *gate_activity(StochasticNode *node);

    // The block id of the single valid dbt_block marker of `node` (a fixed
    // positive integer), or 0 when there is none or more than one.
    unsigned int block_id(StochasticNode *node);
  }
}

#endif
