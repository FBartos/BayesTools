#include "BTFactorSlice.h"
#include "BTSamplerDeclarations.h"

#include <BayesTools/factor_slice_kernel.h>

#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <sampler/GraphView.h>
#include <sampler/MutableSampleMethod.h>
#include <sampler/MutableSampler.h>
#include <sampler/Sampler.h>

#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace jags {
  namespace BayesTools {

    namespace {

      // The transforms, the univariate slice update and the adaptation of the
      // axes are the kernel in inst/include/BayesTools/factor_slice_kernel.h,
      // shared with RoBMA's block sampler; this class is its JAGS side.
      namespace kernel = ::BayesTools::factor_slice;

      const char *const factor_slice_name = "BayesTools::FactorSlice";
      const double negative_infinity = -std::numeric_limits<double>::infinity();

      class FactorSlice : public MutableSampleMethod
      {
        GraphView const *_graph;
        const unsigned int _chain;
        const unsigned int _d;
        bool _adapting;
        unsigned int _sampling_updates;  // updates after adaptation
        std::vector<kernel::Coordinate> _coordinates;
        kernel::Axes _axes;
        std::vector<double> _x;  // scratch
        std::exception_ptr _deferred_error;  // a failure of the check at adaptOff()

        // log pi(x | rest) at the current node values: every member's full log
        // density (a member may be a stochastic child of another member) plus the
        // likelihood terms of the other stochastic children of the block. A NaN
        // term, a NaN sum (opposite infinities) or a +Inf sum is a failure to
        // calculate the density and raises a node error, as the stock slice
        // sampler does through GraphView::logFullConditional(): a NaN is no
        // evidence of zero density, so rejecting it could silently truncate the
        // target. A NaN term is reported on its node; a NaN or +Inf sum on the
        // first member, because a module cannot name nodes (JAGS resolves node
        // names from the model's symbol table when it reports the error).
        double log_full_conditional(char const *where) const
        {
          double out = 0.0;
          for (StochasticNode const *node : _graph->nodes()) {
            const double term = node->logDensity(_chain, PDF_FULL);
            if (std::isnan(term)) density_failure(node, where, "the log density term of this node is NaN");
            out += term;
          }
          for (StochasticNode const *child : _graph->stochasticChildren()) {
            const double term = child->logDensity(_chain, PDF_LIKELIHOOD);
            if (std::isnan(term)) density_failure(child, where, "the log density term of this node is NaN");
            out += term;
          }
          if (std::isnan(out) || out > std::numeric_limits<double>::max()) {
            density_failure(_graph->nodes()[0], where, std::isnan(out) ?
              "the summed log density of the block whose first member is this node is NaN" :
              "the summed log density of the block whose first member is this node is +Inf");
          }
          return out;
        }

        static void density_failure(Node const *node, char const *where, char const *cause)
        {
          throwNodeError(node, std::string(
            "Factor slice sampler: failure to calculate the log density of the block at ") + where + ": " + cause);
        }

        // log pi_u of a proposal u; sets the node values to g^-1(u) when it is
        // representable. A proposal outside the support or with log density -Inf
        // has log pi_u = -Inf (rejected; nothing is clamped); a NaN or +Inf
        // density stops the fit (log_full_conditional()). An error JAGS raises
        // while setValue() recomputes the deterministic descendants (invalid
        // parent values) is not caught either: it is a model error and stops the
        // fit as it does with the stock samplers.
        double log_target(std::vector<double> const &u)
        {
          for (unsigned int i = 0; i < _d; ++i) {
            if (!std::isfinite(u[i]) || !_coordinates[i].to_x(u[i], &_x[i])) return negative_infinity;
          }
          _graph->setValue(_x, _chain);
          const double out = log_full_conditional("a proposal") + kernel::log_jacobian(_coordinates, u);
          return out;
        }

        // One univariate slice update along direction k. The node values hold
        // g^-1(u) on entry and on exit; the log density at the current state is
        // recomputed for every direction.
        void slice_along(unsigned int k, std::vector<double> &u, RNG *rng)
        {
          std::vector<double> x_start(_d);
          _graph->getValue(x_start, _chain);
          const double lp0 = log_full_conditional("the current state") + kernel::log_jacobian(_coordinates, u);
          if (!std::isfinite(lp0)) {
            throwNodeError(_graph->nodes()[0],
              "Factor slice sampler: current value of the block has no finite density");
          }
          const kernel::SliceResult result = kernel::slice_along(
            u, _axes.direction(k), _axes.width(k), lp0,
            [this](std::vector<double> const &point) { return log_target(point); },
            [rng]() { return rng->uniform(); },
            [rng]() { return rng->exponential(); });
          // Floating-point exhaustion of the interval: the state stays.
          if (!result.accepted) _graph->setValue(x_start, _chain);
        }

        static std::vector<kernel::Coordinate> coordinates(GraphView const *graph, unsigned int chain)
        {
          std::vector<kernel::Coordinate> out;
          for (StochasticNode *node : graph->nodes()) {
            // The support is fixed (checked by the factory), so it is read once.
            double lower, upper;
            node->support(&lower, &upper, 1, chain);
            out.push_back(kernel::Coordinate(lower, upper));
          }
          return out;
        }

        // Stops the fit when a transformed member's coordinate cannot resolve its
        // posterior at the current state (kernel::resolvable(): the coordinate
        // reaches a value only to a relative precision of about |u| eps of its
        // distance from the anchoring bound, so with bounds far from where the
        // posterior concentrates the sampler would move on a grid coarser than
        // the posterior, or not at all, and silently sample a wrong
        // distribution). A representability limit of the sampler's scale, not a
        // model error: the block must be declared without the member. Runs at
        // construction, at every refresh of the axes, at the first adaptOff()
        // and every kernel::resolvability_interval updates after adaptation (a
        // member can drift to where its coordinate no longer resolves the
        // target); it draws no random numbers and restores the node values,
        // so it changes no draw. A member on a bound is left to update().
        void check_resolution(char const *site)
        {
          std::vector<StochasticNode *> const &nodes = _graph->nodes();
          std::vector<double> x(_d);
          _graph->getValue(x, _chain);
          std::vector<double> trial(x);
          for (unsigned int i = 0; i < _d; ++i) {
            kernel::Coordinate const &coordinate = _coordinates[i];
            if (coordinate.transform() == kernel::Transform::identity ||
                !(x[i] > coordinate.lower() && x[i] < coordinate.upper())) continue;
            bool ok;
            try {
              ok = kernel::resolvable(coordinate, x[i], [&](double value) {
                trial[i] = value;
                _graph->setValue(trial, _chain);
                return log_full_conditional("a resolvability check");
              });
            }
            catch (...) {
              _graph->setValue(x, _chain);
              throw;
            }
            trial[i] = x[i];
            _graph->setValue(x, _chain);
            if (!ok) {
              throwNodeError(nodes[i], std::string(
                "Factor slice sampler: the support bounds of this node are too far from where its posterior "
                "concentrates for the sampler's scale (found ") + site + "); declare the block without this node");
            }
          }
        }

      public:
        FactorSlice(GraphView const *graph, unsigned int chain)
          : _graph(graph), _chain(chain), _d(graph->length()), _adapting(true), _sampling_updates(0),
            _coordinates(coordinates(graph, chain)), _axes(_d), _x(_d)
        {
          graph->checkFinite(chain);
          check_resolution("at the start of the fit");
        }

        void update(RNG *rng)
        {
          if (_deferred_error) std::rethrow_exception(_deferred_error);
          std::vector<StochasticNode *> const &nodes = _graph->nodes();
          _graph->getValue(_x, _chain);
          std::vector<double> u(_d);
          for (unsigned int i = 0; i < _d; ++i) {
            u[i] = _coordinates[i].to_u(_x[i]);
            if (!std::isfinite(u[i]) || !(_x[i] > _coordinates[i].lower() && _x[i] < _coordinates[i].upper())) {
              throwNodeError(nodes[i],
                "Factor slice sampler: current value lies on a bound of its support");
            }
          }
          const std::vector<double> u_start(u);
          for (unsigned int k = 0; k < _d; ++k) slice_along(k, u, rng);
          if (_adapting) {
            if (_axes.adapt(u_start, u)) check_resolution("at a refresh of the axes during adaptation");
          }
          else if (++_sampling_updates % kernel::resolvability_interval == 0) {
            check_resolution("during sampling");
          }
        }

        bool isAdaptive() const
        {
          return true;
        }

        // JAGS calls adaptOff() repeatedly; the check runs at the first call. JAGS
        // reports an exception thrown here only as a message (rjags' adapt()
        // returns normally and the model is deleted), so a failure is raised at
        // the next update() instead, before any retained draw.
        void adaptOff()
        {
          if (!_adapting) return;
          _adapting = false;
          try {
            check_resolution("at the end of adaptation");
          }
          catch (...) {
            _deferred_error = std::current_exception();
          }
        }

        // Correctness never depends on adaptation, but runjags' serial method runs
        // the requested adaptation phase only when some sampler reports that it
        // still needs adapting: incomplete until one adaptive update, complete
        // afterwards (unchanged by adaptOff()).
        bool checkAdaptation() const
        {
          return _axes.updates() > 0;
        }
      };
    }

    std::vector<Sampler *> FactorSliceFactory::makeSamplers(
        std::list<StochasticNode *> const &nodes, Graph const &graph) const
    {
      std::map<unsigned int, std::vector<StochasticNode *> > blocks;
      for (StochasticNode *node : nodes) {
        if (node->length() != 1 || node->isDiscreteValued() || node->df() == 0 ||
            isObserved(node) || !isSupportFixed(node)) continue;
        // A gate declaration takes precedence over a block declaration.
        if (has_gate_marker(node)) continue;
        const unsigned int id = block_id(node);
        if (id == 0) continue;
        double lower, upper;
        node->support(&lower, &upper, 1, 0);
        if (!(lower < upper)) continue;
        blocks[id].push_back(node);
      }

      std::vector<Sampler *> samplers;
      const std::string sampler_name(factor_slice_name);
      try {
        for (auto const &block : blocks) {
          if (block.second.size() < 2) continue;
          std::unique_ptr<GraphView> view(new GraphView(block.second, graph, true));
          std::vector<std::unique_ptr<MutableSampleMethod> > owned_methods;
          std::vector<MutableSampleMethod *> methods;
          for (unsigned int chain = 0; chain < block.second[0]->nchain(); ++chain) {
            owned_methods.emplace_back(new FactorSlice(view.get(), chain));
            methods.push_back(owned_methods.back().get());
          }
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

    std::string FactorSliceFactory::name() const
    {
      return factor_slice_name;
    }

    FactorSliceFactory *factor_slice_factory()
    {
      static FactorSliceFactory factory;
      return &factory;
    }
  }
}
