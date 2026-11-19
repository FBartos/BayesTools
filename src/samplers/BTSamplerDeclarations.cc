#include "BTSamplerDeclarations.h"
#include "../distributions/DBTSamplerMarker.h"

#include <graph/Node.h>
#include <graph/StochasticNode.h>
#include <distribution/Distribution.h>

#include <cmath>
#include <list>
#include <set>
#include <vector>

namespace jags {
  namespace BayesTools {

    namespace {

      // Markers of one kind that name `node` as their first parent.
      std::vector<StochasticNode const *> markers_of(StochasticNode *node, char const *name)
      {
        std::vector<StochasticNode const *> out;
        std::list<StochasticNode *> const *children = node->stochasticChildren();
        for (StochasticNode const *child : *children) {
          if (child->distribution()->name() != name) continue;
          std::vector<Node const *> const &parents = child->parents();
          if (parents.size() < 2 || parents[0] != node) continue;
          out.push_back(child);
        }
        return out;
      }

      // An observed marker holding the value 0 (observed values agree across chains).
      bool observed_zero(StochasticNode const *marker)
      {
        return isObserved(marker) && marker->length() == 1 && marker->value(0)[0] == 0.0;
      }

      // Whether `target` is an ancestor of `start` (or `start` itself) through any
      // path, stochastic nodes included: an <active> node that depends on the
      // gated node, even through one of its stochastic children, would make the
      // closed-gate draw ignore part of the node's full conditional.
      bool depends_on(Node const *start, Node const *target)
      {
        std::vector<Node const *> stack(1, start);
        std::set<Node const *> visited;
        while (!stack.empty()) {
          Node const *current = stack.back();
          stack.pop_back();
          if (current == target) return true;
          if (!visited.insert(current).second) continue;
          for (Node const *parent : current->parents()) stack.push_back(parent);
        }
        return false;
      }
    }

    bool has_gate_marker(StochasticNode *node)
    {
      return !markers_of(node, gate_marker_name).empty();
    }

    Node const *gate_activity(StochasticNode *node)
    {
      std::vector<StochasticNode const *> markers = markers_of(node, gate_marker_name);
      if (markers.size() != 1 || !observed_zero(markers[0])) return nullptr;
      Node const *active = markers[0]->parents()[1];
      if (active->length() != 1) return nullptr;
      if (depends_on(active, node)) return nullptr;
      return active;
    }

    unsigned int block_id(StochasticNode *node)
    {
      std::vector<StochasticNode const *> markers = markers_of(node, block_marker_name);
      if (markers.size() != 1 || !observed_zero(markers[0])) return 0;
      Node const *id = markers[0]->parents()[1];
      if (!id->isFixed() || id->length() != 1) return 0;
      const double value = id->value(0)[0];
      if (!std::isfinite(value) || value < 1.0 || value != std::floor(value) || value > 1.0e9) return 0;
      return static_cast<unsigned int>(value);
    }
  }
}
