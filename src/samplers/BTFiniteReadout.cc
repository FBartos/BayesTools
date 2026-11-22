#include "BTFiniteReadout.h"

#include <distribution/Distribution.h>
#include <function/ScalarFunction.h>
#include <graph/Graph.h>
#include <graph/LogicalNode.h>
#include <graph/StochasticNode.h>
#include <model/BUGSModel.h>
#include <model/MonitorFactory.h>
#include <model/MonitorStat.h>
#include <model/NodeArraySubset.h>
#include <model/TraceMonitor.h>
#include <module/Module.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <sampler/ImmutableSampleMethod.h>
#include <sampler/ImmutableSampler.h>
#include <sampler/SamplerFactory.h>
#include <sampler/SingletonGraphView.h>
#include <util/nainf.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <stdexcept>

namespace jags {
  namespace BayesTools {
    namespace {

      constexpr char finite_sampler_name[] = "BayesTools::Finite";
      constexpr char finite_tag_name[] = "bt_finite_tag";
      constexpr char finite_witness_tag_name[] = "bt_finite_witness_tag";
      constexpr char finite_marker_prefix[] = "BayesTools_finite_";
      constexpr char finite_monitor_prefix[] = "BayesTools_finite_readout_";
      constexpr unsigned long finite_max_levels = 100;

      class FiniteTag : public ScalarFunction {
      public:
        FiniteTag(char const *name, unsigned long npar) : ScalarFunction(name, npar) {}
        double evaluate(std::vector<double const *> const &) const override { return 0.0; }
      };

      struct FiniteMarker {
        Node const *tag = nullptr;
        StochasticNode const *witness = nullptr;
      };

      FiniteMarker finite_marker(StochasticNode *target, Graph const &graph)
      {
        static const std::string one_parent = std::string("(") + finite_tag_name + "(x))";
        static const std::string two_parents = std::string("(") + finite_witness_tag_name + "(x,w))";
        FiniteMarker marker;
        for (DeterministicNode *child : *target->deterministicChildren()) {
          if (!dynamic_cast<LogicalNode const *>(child) || graph.contains(child) ||
              !child->stochasticChildren()->empty() || !child->deterministicChildren()->empty()) {
            continue;
          }
          auto const &parents = child->parents();
          if (parents.empty() || parents.front() != target) continue;

          StochasticNode const *witness = nullptr;
          if (parents.size() == 2 && child->deparse({"x", "w"}) == two_parents) {
            witness = dynamic_cast<StochasticNode const *>(parents[1]);
            if (!witness || witness->length() != 1 || witness->isDiscreteValued() ||
                !graph.contains(witness)) return {};
          } else if (parents.size() != 1 || child->deparse({"x"}) != one_parent) {
            continue;
          }
          if (marker.tag) return {};
          marker = {child, witness};
        }
        return marker;
      }

      class FiniteMethod;

      // Both factories see a compiled marker, but JAGS only supplies a model to
      // the monitor factory. Connect them by node identity for the sampler's
      // lifetime. Sampling threads only touch their own rows, outside this map.
      class FiniteOwners {
        std::mutex _mutex;
        std::map<Node const *, FiniteMethod *> _methods;
      public:
        void insert(Node const *marker, FiniteMethod *method)
        {
          std::lock_guard<std::mutex> lock(_mutex);
          if (!_methods.emplace(marker, method).second) {
            throw std::logic_error("Two finite samplers own the same BayesTools marker");
          }
        }

        void erase(Node const *marker)
        {
          std::lock_guard<std::mutex> lock(_mutex);
          _methods.erase(marker);
        }

        FiniteMethod *find(Node const *marker)
        {
          std::lock_guard<std::mutex> lock(_mutex);
          auto method = _methods.find(marker);
          return method == _methods.end() ? nullptr : method->second;
        }
      };

      class FiniteMethod : public ImmutableSampleMethod {
        SingletonGraphView const *_view;
        double const _lower;
        unsigned long const _size;
        FiniteMarker const _marker;
        std::shared_ptr<FiniteOwners> const _owners;
        mutable std::vector<std::vector<double>> _rows;
      public:
        FiniteMethod(SingletonGraphView const *view, double lower, unsigned long size,
                     FiniteMarker marker, std::shared_ptr<FiniteOwners> const &owners)
          : _view(view), _lower(lower), _size(size), _marker(marker), _owners(owners)
        {
          _owners->insert(_marker.tag, this);
        }

        ~FiniteMethod() override { _owners->erase(_marker.tag); }

        StochasticNode const *target() const { return _view->node(); }
        unsigned long width() const { return 2 * _size + 2 + (_marker.witness != nullptr); }
        bool tracing() const { return !_rows.empty(); }

        void start_trace()
        {
          _rows.assign(target()->nchain(), std::vector<double>(width()));
        }

        void stop_trace() { _rows.clear(); }
        std::vector<double> const &latest(unsigned int chain) const { return _rows[chain]; }

        void update(unsigned int chain, RNG *rng) const override
        {
          double *row = tracing() ? _rows[chain].data() : nullptr;
          if (row && _marker.witness) row[2 * _size + 2] = _marker.witness->value(chain)[0];

          // Preserve JAGS 5 base::Finite's enumeration, arithmetic order and
          // one uniform draw. Record mass before it is overwritten by the CDF.
          std::vector<double> likelihood(_size);
          double maximum = JAGS_NEGINF;
          for (unsigned long i = 0; i < _size; ++i) {
            double value = _lower + i;
            _view->setValue(&value, 1, chain);
            likelihood[i] = _view->logFullConditional(chain);
            if (likelihood[i] > maximum) maximum = likelihood[i];
          }

          double total = 0;
          for (unsigned long i = 0; i < _size; ++i) {
            double shifted = likelihood[i] - maximum;
            double mass = std::exp(shifted);
            total += mass;
            likelihood[i] = total;
            if (row) {
              row[i] = mass;
              row[_size + i] = shifted;
            }
          }
          if (!std::isfinite(total)) throwNodeError(target(), "Cannot normalize density");

          double draw = rng->uniform() * total;
          long index = std::upper_bound(likelihood.begin(), likelihood.end(), draw) - likelihood.begin();
          double value = _lower + index;
          _view->setValue(&value, 1, chain);

          if (row) {
            double log_total = std::log(total);
            for (unsigned long i = 0; i < _size; ++i) {
              row[i] /= total;
              row[_size + i] -= log_total;
            }
            row[2 * _size] = value;
            row[2 * _size + 1] = _lower;
          }
        }
      };

      class FiniteSamplerFactory : public SamplerFactory {
        std::shared_ptr<FiniteOwners> const _owners;
      public:
        explicit FiniteSamplerFactory(std::shared_ptr<FiniteOwners> const &owners) : _owners(owners) {}
        std::string name() const override { return finite_sampler_name; }

        std::vector<Sampler *> makeSamplers(std::list<StochasticNode *> const &nodes,
                                            Graph const &graph) const override
        {
          for (StochasticNode *node : nodes) {
            if (node->length() != 1 || !node->isDiscreteValued() || !node->fullRank()) continue;
            auto const &distribution = node->distribution()->name();
            if (distribution != "dcat" && distribution != "dbern") continue;
            FiniteMarker marker = finite_marker(node, graph);
            if (!marker.tag || !isSupportFixed(node)) continue;

            double lower, upper;
            node->support(&lower, &upper, 1, 0);
            double expected_lower = distribution == "dbern" ? 0.0 : 1.0;
            if (lower != expected_lower || !std::isfinite(upper) || upper != std::floor(upper) ||
                upper < lower || upper - lower >= finite_max_levels) continue;

            std::vector<Sampler *> samplers(1);
            auto view = std::make_unique<SingletonGraphView>(node, graph);
            auto method = std::make_unique<FiniteMethod>(
              view.get(), lower, static_cast<unsigned long>(upper - lower + 1), marker, _owners
            );
            std::string sampler_name = name();
            // The Sampler base owns the view as soon as its constructor runs.
            samplers[0] = new ImmutableSampler(view.release(), method.get(), sampler_name);
            method.release();
            return samplers;
          }
          return {};
        }
      };

      class FiniteStat : public MonitorStat {
        FiniteMethod &_method;
      public:
        FiniteStat(std::string const &name, FiniteMethod &method)
          : MonitorStat(method.width()), _method(method)
        {
          std::vector<std::string> labels;
          for (unsigned long col = 1; col <= method.width(); ++col) {
            labels.push_back(name + "[" + std::to_string(col) + "]");
          }
          setNames(labels);
          _method.start_trace();
        }

        // BUGSModel deletes its monitors before Model deletes its samplers.
        ~FiniteStat() override { _method.stop_trace(); }
        std::vector<double> value(unsigned int chain) const override { return _method.latest(chain); }
      };

      class FiniteMonitorFactory : public MonitorFactory {
        std::shared_ptr<FiniteOwners> const _owners;
      public:
        explicit FiniteMonitorFactory(std::shared_ptr<FiniteOwners> const &owners) : _owners(owners) {}
        std::string name() const override { return "BayesTools::FiniteReadout"; }

        Monitor *getMonitor(std::string const &name, Range const &range, BUGSModel *model,
                            std::string const &stat, std::string const &summary,
                            std::string &msg) override
        {
          if (name.compare(0, sizeof(finite_monitor_prefix) - 1, finite_monitor_prefix) != 0) return nullptr;
          if (stat != "value" || summary != "trace" || !isNULL(range) || model->isAdapting()) {
            msg = "A BayesTools finite readout requires a complete value trace after adaptation";
            return nullptr;
          }
          auto id = name.substr(sizeof(finite_monitor_prefix) - 1);
          auto marker = id.empty() || id.find_first_not_of("0123456789") != std::string::npos
            ? nullptr : model->symtab().getVariable(finite_marker_prefix + id);
          if (marker) {
            NodeArraySubset subset(marker, range);
            if (subset.length() == 1 && !subset.missing()[0]) {
              auto method = _owners->find(subset.allnodes()[0]);
              if (method && !method->tracing()) {
                std::vector<Node const *> nodes{method->target()};
                return new TraceMonitor(nodes, new FiniteStat(name, *method));
              }
            }
          }

          msg = "No available BayesTools::Finite sampler owns the marker for '" + name + "'";
          return nullptr;
        }
      };
    }

    FiniteReadout::FiniteReadout(Module &module)
    {
      auto owners = std::make_shared<FiniteOwners>();
      _sampler = std::make_unique<FiniteSamplerFactory>(owners);
      _monitor = std::make_unique<FiniteMonitorFactory>(owners);
      module.insert(new FiniteTag(finite_tag_name, 1));
      module.insert(new FiniteTag(finite_witness_tag_name, 2));
      module.insert(_sampler.get());
      module.insert(_monitor.get());
    }

    FiniteReadout::~FiniteReadout() = default;
  }
}
