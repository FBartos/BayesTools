#include "BTFactorSlice.h"
#include "BTSamplerDeclarations.h"

#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <sampler/GraphView.h>
#include <sampler/MutableSampleMethod.h>
#include <sampler/MutableSampler.h>
#include <sampler/Sampler.h>

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <vector>

namespace jags {
  namespace BayesTools {

    namespace {

      const char *const factor_slice_name = "BayesTools::FactorSlice";

      // Randomized stepping-out limit (Neal 2003, section 4.1).
      const unsigned int stepping_out_limit = 10;
      // Adaptation: the axes are first refreshed after this many adaptive
      // updates and then every refresh_interval updates.
      const unsigned int first_refresh = 100;
      const unsigned int refresh_interval = 50;
      const double width_multiplier = 2.5;
      const double negative_infinity = -std::numeric_limits<double>::infinity();

      // Coordinate-wise bijection u = g(x) to the real line.
      enum class Transform { identity, lower, interval, upper };

      class FactorSlice : public MutableSampleMethod
      {
        GraphView const *_graph;
        const unsigned int _chain;
        const unsigned int _d;
        bool _adapting;
        bool _mean_initialized;
        unsigned int _n;
        std::vector<Transform> _transform;
        std::vector<double> _lower, _upper;
        std::vector<double> _mean, _cov;  // running moments of u (row-major covariance)
        std::vector<double> _dirs;        // column k is direction k (row-major d x d)
        std::vector<double> _width;
        std::vector<double> _x, _point;   // scratch

        double to_u(unsigned int i, double x) const
        {
          switch (_transform[i]) {
          case Transform::lower: return std::log(x - _lower[i]);
          case Transform::interval: return std::log(x - _lower[i]) - std::log(_upper[i] - x);
          case Transform::upper: return std::log(_upper[i] - x);
          default: return x;
          }
        }

        // Back-transform; false when the coordinate is not finite or does not lie
        // strictly inside the support (such a proposal has zero density).
        bool to_x(unsigned int i, double u, double *x) const
        {
          double out;
          switch (_transform[i]) {
          case Transform::lower:
            out = _lower[i] + std::exp(u);
            break;
          case Transform::interval:
            // Evaluate from the nearer bound for accuracy in both tails.
            out = u >= 0.0 ?
              _upper[i] - (_upper[i] - _lower[i]) / (1.0 + std::exp(u)) :
              _lower[i] + (_upper[i] - _lower[i]) / (1.0 + std::exp(-u));
            break;
          case Transform::upper:
            out = _upper[i] - std::exp(u);
            break;
          default:
            out = u;
          }
          *x = out;
          return std::isfinite(out) && out > _lower[i] && out < _upper[i];
        }

        double log_jacobian(unsigned int i, double u) const
        {
          switch (_transform[i]) {
          case Transform::lower:
          case Transform::upper:
            return u;
          case Transform::interval:
            return -std::fabs(u) - 2.0 * std::log1p(std::exp(-std::fabs(u))) +
              std::log(_upper[i] - _lower[i]);
          default:
            return 0.0;
          }
        }

        double log_jacobian(std::vector<double> const &u) const
        {
          double out = 0.0;
          for (unsigned int i = 0; i < _d; ++i) out += log_jacobian(i, u[i]);
          return out;
        }

        // log pi(x | rest) at the current node values: every member's full log
        // density (a member may be a stochastic child of another member) plus the
        // likelihood terms of the other stochastic children of the block. It is
        // summed here because GraphView::logFullConditional() throws on a NaN;
        // a NaN term returns NaN.
        double log_full_conditional() const
        {
          double out = 0.0;
          for (StochasticNode const *node : _graph->nodes()) {
            const double term = node->logDensity(_chain, PDF_FULL);
            if (std::isnan(term)) return std::numeric_limits<double>::quiet_NaN();
            out += term;
          }
          for (StochasticNode const *child : _graph->stochasticChildren()) {
            const double term = child->logDensity(_chain, PDF_LIKELIHOOD);
            if (std::isnan(term)) return std::numeric_limits<double>::quiet_NaN();
            out += term;
          }
          return out;
        }

        // log pi_u of a proposal u; sets the node values to g^-1(u) when it is
        // representable. A proposal outside the support, or whose log density has
        // a NaN term, a NaN sum (opposite infinities) or a +Inf sum, has
        // log pi_u = -Inf (rejected; nothing is clamped). An error JAGS raises
        // while setValue() recomputes the deterministic descendants (invalid
        // parent values) is not caught: it is a model error and stops the fit as
        // it does with the stock samplers.
        double log_target(std::vector<double> const &u)
        {
          for (unsigned int i = 0; i < _d; ++i) {
            if (!std::isfinite(u[i]) || !to_x(i, u[i], &_x[i])) return negative_infinity;
          }
          _graph->setValue(_x, _chain);
          const double out = log_full_conditional() + log_jacobian(u);
          return std::isnan(out) || out > std::numeric_limits<double>::max() ? negative_infinity : out;
        }

        double log_target_along(std::vector<double> const &u, unsigned int k, double t)
        {
          for (unsigned int i = 0; i < _d; ++i) _point[i] = u[i] + t * _dirs[i * _d + k];
          return log_target(_point);
        }

        // One univariate slice update of t -> pi_u(u + t e_k) (Neal 2003,
        // stepping out with randomized limit and shrinkage). The node values hold
        // g^-1(u) on entry and on exit.
        void slice_along(unsigned int k, std::vector<double> &u, RNG *rng)
        {
          std::vector<double> x_start(_d);
          _graph->getValue(x_start, _chain);
          const double lp0 = log_full_conditional() + log_jacobian(u);
          if (!std::isfinite(lp0)) {
            throwNodeError(_graph->nodes()[0],
              "Factor slice sampler: current value of the block has no finite density");
          }
          const double logy = lp0 - rng->exponential();
          const double w = _width[k];
          double left = -w * rng->uniform();
          double right = left + w;
          int steps_left = static_cast<int>(std::floor(stepping_out_limit * rng->uniform()));
          int steps_right = static_cast<int>(stepping_out_limit) - 1 - steps_left;
          while (steps_left > 0 && log_target_along(u, k, left) > logy) {
            left -= w;
            --steps_left;
          }
          while (steps_right > 0 && log_target_along(u, k, right) > logy) {
            right += w;
            --steps_right;
          }
          for (;;) {
            // Floating-point exhaustion of the interval: the state stays.
            const double middle = left + 0.5 * (right - left);
            if (!(middle > left && middle < right)) break;
            const double t = left + (right - left) * rng->uniform();
            if (log_target_along(u, k, t) > logy) {
              u = _point;  // the node values already hold the accepted point
              return;
            }
            if (t < 0.0) left = t;
            else right = t;
          }
          _graph->setValue(x_start, _chain);
        }

        // Directions become the eigenvectors of the covariance (cyclic Jacobi) and
        // widths 2.5 standard deviations along them; the previous axes stay when
        // an eigenvalue is not positive and finite.
        void refresh_axes()
        {
          std::vector<double> a(_cov), v(_d * _d, 0.0);
          for (unsigned int i = 0; i < _d; ++i) v[i * _d + i] = 1.0;
          for (unsigned int sweep = 0; sweep < 100; ++sweep) {
            double off = 0.0, diagonal = 0.0;
            for (unsigned int p = 0; p < _d; ++p) {
              diagonal += a[p * _d + p] * a[p * _d + p];
              for (unsigned int q = p + 1; q < _d; ++q) off += a[p * _d + q] * a[p * _d + q];
            }
            if (!(off > 1e-30 * diagonal)) break;
            for (unsigned int p = 0; p < _d; ++p) {
              for (unsigned int q = p + 1; q < _d; ++q) {
                const double apq = a[p * _d + q];
                if (apq == 0.0) continue;
                const double theta = (a[q * _d + q] - a[p * _d + p]) / (2.0 * apq);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                  (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;
                for (unsigned int k = 0; k < _d; ++k) {
                  const double akp = a[k * _d + p], akq = a[k * _d + q];
                  a[k * _d + p] = c * akp - s * akq;
                  a[k * _d + q] = s * akp + c * akq;
                }
                for (unsigned int k = 0; k < _d; ++k) {
                  const double apk = a[p * _d + k], aqk = a[q * _d + k];
                  a[p * _d + k] = c * apk - s * aqk;
                  a[q * _d + k] = s * apk + c * aqk;
                }
                for (unsigned int k = 0; k < _d; ++k) {
                  const double vkp = v[k * _d + p], vkq = v[k * _d + q];
                  v[k * _d + p] = c * vkp - s * vkq;
                  v[k * _d + q] = s * vkp + c * vkq;
                }
              }
            }
          }
          for (unsigned int k = 0; k < _d; ++k) {
            const double lambda = a[k * _d + k];
            if (!(lambda > 0.0) || !std::isfinite(lambda)) return;
          }
          _dirs = v;
          for (unsigned int k = 0; k < _d; ++k) {
            _width[k] = width_multiplier * std::sqrt(a[k * _d + k]);
          }
        }

        void adapt(std::vector<double> const &u_start, std::vector<double> const &u)
        {
          if (!_mean_initialized) {
            _mean = u_start;
            _mean_initialized = true;
          }
          ++_n;
          const double weight = std::pow(static_cast<double>(_n) + 1.0, -0.7);
          std::vector<double> delta(_d);
          for (unsigned int i = 0; i < _d; ++i) {
            delta[i] = u[i] - _mean[i];
            _mean[i] += weight * delta[i];
          }
          // (1 - w) C + w (1 - w) delta delta', symmetric by construction.
          for (unsigned int i = 0; i < _d; ++i) {
            for (unsigned int j = 0; j <= i; ++j) {
              const double value = _cov[i * _d + j] +
                weight * ((1.0 - weight) * delta[i] * delta[j] - _cov[i * _d + j]);
              _cov[i * _d + j] = value;
              _cov[j * _d + i] = value;
            }
          }
          if (_n >= first_refresh && (_n - first_refresh) % refresh_interval == 0) refresh_axes();
        }

      public:
        FactorSlice(GraphView const *graph, unsigned int chain)
          : _graph(graph), _chain(chain), _d(graph->length()), _adapting(true),
            _mean_initialized(false), _n(0), _transform(_d, Transform::identity),
            _lower(_d), _upper(_d), _mean(_d, 0.0), _cov(_d * _d, 0.0),
            _dirs(_d * _d, 0.0), _width(_d, 1.0), _x(_d), _point(_d)
        {
          std::vector<StochasticNode *> const &nodes = graph->nodes();
          for (unsigned int i = 0; i < _d; ++i) {
            // The support is fixed (checked by the factory), so it is read once.
            nodes[i]->support(&_lower[i], &_upper[i], 1, chain);
            const bool has_lower = std::isfinite(_lower[i]);
            const bool has_upper = std::isfinite(_upper[i]);
            _transform[i] = has_lower && has_upper ? Transform::interval :
              has_lower ? Transform::lower : has_upper ? Transform::upper : Transform::identity;
            _dirs[i * _d + i] = 1.0;
          }
          graph->checkFinite(chain);
        }

        void update(RNG *rng)
        {
          std::vector<StochasticNode *> const &nodes = _graph->nodes();
          _graph->getValue(_x, _chain);
          std::vector<double> u(_d);
          for (unsigned int i = 0; i < _d; ++i) {
            u[i] = to_u(i, _x[i]);
            if (!std::isfinite(u[i]) || !(_x[i] > _lower[i] && _x[i] < _upper[i])) {
              throwNodeError(nodes[i],
                "Factor slice sampler: current value lies on a bound of its support");
            }
          }
          const std::vector<double> u_start(u);
          for (unsigned int k = 0; k < _d; ++k) slice_along(k, u, rng);
          if (_adapting) adapt(u_start, u);
        }

        bool isAdaptive() const
        {
          return true;
        }

        void adaptOff()
        {
          _adapting = false;
        }

        // Correctness never depends on adaptation, but runjags' serial method runs
        // the requested adaptation phase only when some sampler reports that it
        // still needs adapting: incomplete until one adaptive update, complete
        // afterwards (unchanged by adaptOff()).
        bool checkAdaptation() const
        {
          return _n > 0;
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
