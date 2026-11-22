#ifndef BAYESTOOLS_SAMPLER_DECLARATIONS_H
#define BAYESTOOLS_SAMPLER_DECLARATIONS_H

#include <distribution/Distribution.h>
#include <graph/Node.h>
#include <graph/StochasticNode.h>

#include <cmath>
#include <set>
#include <vector>

namespace jags {
namespace BayesTools {
namespace sampler_declarations {

constexpr char const *gate_marker_name() { return "dbt_gate"; }
constexpr char const *block_marker_name() { return "dbt_block"; }

// Markers of one kind that name `node` as their first parent.
inline std::vector<StochasticNode const *> markers_of(StochasticNode *node, char const *name)
{
  std::vector<StochasticNode const *> markers;
  for (StochasticNode const *child : *node->stochasticChildren()) {
    if (child->distribution()->name() != name) continue;
    std::vector<Node const *> const &parents = child->parents();
    if (parents.size() < 2 || parents[0] != node) continue;
    markers.push_back(child);
  }
  return markers;
}

// An observed scalar marker holding zero. Observed values agree across chains.
inline bool observed_zero(StochasticNode const *marker)
{
  return isObserved(marker) && marker->length() == 1 && marker->value(0)[0] == 0.0;
}

// Whether `target` is `start` or an ancestor through any path, including
// through stochastic nodes.
inline bool depends_on(Node const *start, Node const *target)
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

// A gate marker reserves its first parent even if the declaration is invalid.
inline bool has_gate_marker(StochasticNode *node)
{
  return !markers_of(node, gate_marker_name()).empty();
}

inline Node const *gate_activity(StochasticNode *node)
{
  std::vector<StochasticNode const *> markers = markers_of(node, gate_marker_name());
  if (markers.size() != 1 || !observed_zero(markers[0])) return nullptr;
  Node const *active = markers[0]->parents()[1];
  if (active->length() != 1) return nullptr;
  // A gate that depends on the gated node would change its full conditional.
  if (depends_on(active, node)) return nullptr;
  return active;
}

// The id of the single valid block marker, or 0 when there is no such marker.
inline unsigned int block_id(StochasticNode *node)
{
  std::vector<StochasticNode const *> markers = markers_of(node, block_marker_name());
  if (markers.size() != 1 || !observed_zero(markers[0])) return 0;
  Node const *id = markers[0]->parents()[1];
  if (!id->isFixed() || id->length() != 1) return 0;
  const double value = id->value(0)[0];
  if (!std::isfinite(value) || value < 1.0 || value != std::floor(value) || value > 1.0e9) return 0;
  return static_cast<unsigned int>(value);
}

// The corrected scalar sampler defers to either kind of marker child regardless
// of the edge to the node or the validity of the declaration.
inline bool has_sampler_marker_child(StochasticNode *node)
{
  for (StochasticNode const *child : *node->stochasticChildren()) {
    auto const &name = child->distribution()->name();
    if (name == block_marker_name() || name == gate_marker_name()) return true;
  }
  return false;
}

}
}
}

#endif
