#include "BTGatedSlice.h"
#include "BTSamplerDeclarations.h"

#include <graph/DeterministicNode.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <sampler/MutableSampler.h>
#include <sampler/Sampler.h>
#include <sampler/SingletonGraphView.h>
#include <sampler/Slicer.h>

#include <cmath>
#include <memory>
#include <vector>

namespace jags {
  namespace BayesTools {

    namespace {

      const char *const gated_slice_name = "BayesTools::GatedSlice";

      // Same kind of value: equal, or both NaN (two infinities of one sign are equal).
      bool same_log_density(double before, double after)
      {
        return before == after || (std::isnan(before) && std::isnan(after));
      }

      class GatedSlice : public Slicer
      {
        SingletonGraphView const *_graph;
        const unsigned int _chain;
        Node const *_active;
        bool _adapting;
        bool _adapted;  // an update was made in adaptive mode

        // Summed log density of the stochastic children of the node, the marker
        // included (the likelihood part of the graph view's full conditional).
        double children_log_density() const
        {
          double out = 0.0;
          for (StochasticNode const *child : _graph->stochasticChildren()) {
            out += child->logDensity(_chain, PDF_LIKELIHOOD);
          }
          return out;
        }

        bool descendants_finite() const
        {
          for (DeterministicNode const *child : _graph->deterministicChildren()) {
            double const *values = child->value(_chain);
            for (unsigned int i = 0; i < child->length(); ++i) {
              if (!std::isfinite(values[i])) return false;
            }
          }
          return true;
        }

        // Gate closed: by declaration every likelihood child is constant in the
        // node, so its full conditional is its own distribution given its parents
        // and a direct draw is an exact Gibbs update. A draw that makes the node or
        // a deterministic descendant non-finite is outside the region where the
        // model is defined and is rejected (independence Metropolis with the
        // node's distribution as proposal), keeping the previous value.
        void closed_update(RNG *rng)
        {
          StochasticNode *node = _graph->node();
          const double previous = value();
          const double before = _adapting ? children_log_density() : 0.0;
          node->randomSample(rng, _chain);
          const double drawn = node->value(_chain)[0];
          _graph->setValue(&drawn, 1, _chain);
          if (!std::isfinite(drawn) || !descendants_finite()) {
            _graph->setValue(&previous, 1, _chain);
          }
          // Diagnostic for declaration errors during adaptation only; the
          // declaration step in R establishes the condition.
          if (_adapting && !same_log_density(before, children_log_density())) {
            throwNodeError(node, "gated node changes the likelihood while its gate is closed");
          }
        }

      public:
        GatedSlice(SingletonGraphView const *graph, unsigned int chain, Node const *active)
          : Slicer(1.0, 10), _graph(graph), _chain(chain), _active(active), _adapting(true),
            _adapted(false)
        {
          graph->checkFinite(chain);
        }

        double value() const
        {
          return _graph->node()->value(_chain)[0];
        }

        void setValue(double value)
        {
          _graph->setValue(&value, 1, _chain);
        }

        void getLimits(double *lower, double *upper) const
        {
          _graph->node()->support(lower, upper, 1, _chain);
        }

        double logDensity() const
        {
          return _graph->logFullConditional(_chain);
        }

        void update(RNG *rng)
        {
          if (_adapting) _adapted = true;
          if (_active->value(_chain)[0] == 0.0) {
            closed_update(rng);
            return;
          }
          // Gate open: the update of base::RealSlicer; the inherited width adapts
          // only on these updates.
          if (!updateStep(rng)) {
            switch (state()) {
            case SLICER_POSINF:
              throwNodeError(_graph->node(), "Slicer stuck at value with infinite density");
              break;
            case SLICER_NEGINF:
              throwNodeError(_graph->node(), "Current value is inconsistent with data");
              break;
            case SLICER_OK:
              break;
            }
          }
        }

        void adaptOff()
        {
          _adapting = false;
          Slicer::adaptOff();
        }

        // Incomplete until one update (open or closed) in adaptive mode, complete
        // afterwards and unchanged by adaptOff(): runjags' serial method runs the
        // requested adaptation phase only when some sampler still needs adapting.
        // A rarely open gate leaves the width poorly adapted, which costs
        // efficiency only; correctness never depends on adaptation.
        bool checkAdaptation() const
        {
          return _adapted;
        }
      };
    }

    std::vector<Sampler *> GatedSliceFactory::makeSamplers(
        std::list<StochasticNode *> const &nodes, Graph const &graph) const
    {
      std::vector<Sampler *> samplers;
      const std::string sampler_name(gated_slice_name);
      try {
        for (StochasticNode *node : nodes) {
          if (node->length() != 1 || node->isDiscreteValued() || node->df() == 0 ||
              isObserved(node)) continue;
          Node const *active = gate_activity(node);
          if (active == nullptr) continue;
          std::unique_ptr<SingletonGraphView> view(new SingletonGraphView(node, graph));
          std::vector<std::unique_ptr<MutableSampleMethod> > owned_methods;
          std::vector<MutableSampleMethod *> methods;
          for (unsigned int chain = 0; chain < node->nchain(); ++chain) {
            owned_methods.emplace_back(new GatedSlice(view.get(), chain, active));
            methods.push_back(owned_methods.back().get());
          }
          // The Sampler base adopts the view; the methods are released only after
          // the sampler owns them.
          std::unique_ptr<Sampler> sampler(new MutableSampler(view.release(), methods, sampler_name));
          for (std::unique_ptr<MutableSampleMethod> &method : owned_methods) method.release();
          samplers.push_back(sampler.get());
          sampler.release();
        }
      }
      catch (...) {
        for (Sampler *sampler : samplers) delete sampler;
        throw;
      }
      return samplers;
    }

    std::string GatedSliceFactory::name() const
    {
      return gated_slice_name;
    }

    GatedSliceFactory *gated_slice_factory()
    {
      static GatedSliceFactory factory;
      return &factory;
    }
  }
}
