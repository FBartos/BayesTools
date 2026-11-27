// Factor slice kernel of the BayesTools JAGS module (BayesTools::FactorSlice),
// header-only so that other modules (RoBMA, through LinkingTo: BayesTools) use
// the same algorithm. Standard C++ only: no JAGS or R types. The numerical
// behaviour is pinned by FactorSlice: its seeded draws must stay bit-identical
// when this file changes (see the package's JAGS guidance). A module that links
// to BayesTools picks up a change of this file only when it is rebuilt.
//
// Contract of the pieces (namespace BayesTools::factor_slice):
//
// - Coordinate: the bijection u = g(x) of one scalar member from its support
//   (lower, upper), lower < upper, to the real line: identity (no finite
//   bound), log(x - lower), log(x - lower) - log(upper - x) (two finite
//   bounds), log(upper - x). to_u(x) is the forward transform; to_x(u, &x) the
//   inverse, false when the result is not finite or not strictly inside the
//   support; log_jacobian(u) is log |dx/du|; log_jacobian(coordinates, u) sums
//   it over a block in member order. A caller refuses a member with two
//   finite bounds whose difference upper - lower is not finite (FactorSlice
//   stops the fit at construction): such a member cannot be put on this scale.
//   Precision limit: the sampler holds a bounded member through its distance
//   from the anchoring bound (for an interval, the nearer bound; near the
//   centre, the span), so it resolves the member only relative to that
//   distance: to about 16 significant digits of it near the centre of an
//   interval and to about 14 to 15 for a one-sided bound or close to a bound
//   (fewer as |log distance| grows), where a double resolves a value to 16
//   digits of its magnitude. A posterior much narrower than that distance is
//   sampled on a coarse grid, with wrong moments or a chain that can freeze.
//   Measured on development (a standard normal, 10,000 draws): symmetric
//   bounds at 1e13 give 1,969 distinct values with correct moments, at 1e14
//   296 with moments within Monte Carlo error, at 1e15 34 and an SD 15% too
//   small, at 1e16 3 to 5 and an SD up to 60% too large (some chains freeze),
//   at 1e17 a frozen chain; one bound at a distance of 1e12 gives 1,416
//   distinct values, 1e13 188, 1e14 11, and 1e15 a frozen chain. In a block
//   of two normal coefficients with unit posterior scale (4 x 5,000 draws,
//   3 seeds) mean, SD and effective sample size are unchanged up to a
//   distance of 1e14; symmetric bounds at 1e15 give an SD 13% too large and
//   15% fewer effective draws and at 1e16 frozen chains; one bound freezes
//   the chains at 1e15. The posterior scale is not known when the sampler
//   starts, so a caller refuses a member whose starting value lies farther
//   than max_anchor_distance (1e10) from its anchoring bound
//   (anchor_distance()): the moments then hold for posteriors down to about
//   1e-4 wide. This is a limit of the sampler, in the units of the node, not
//   of the model: for such a member, remove a bound that excludes no
//   posterior mass, rescale the parameter, or declare the block without the
//   node. Not detected: a posterior narrower than about 1e-14 of the
//   distance within the limit, and a chain that starts within the limit and
//   concentrates far beyond it. An additive constant so large that the log
//   densities lose their differences breaks the slice arithmetic itself.
//   Far tail: beyond exp_overflow_threshold the back-transform of an interval
//   reaches values at subnormal distances from a bound, as a one-sided member
//   always could. A proposal there is evaluated, so a model whose density is
//   not finite at such a value (a precision 1 / x that overflows, say) stops
//   the fit with the caller's non-finite density error.
// - slice_along(u, direction, width, log_density, log_density_at, uniform,
//   exponential): one univariate slice update (Neal 2003: slice level, random
//   placement of the initial interval, randomized division of a stepping-out
//   budget of 10 widths, stepping out, shrinkage) of t -> pi(u + t direction).
//   Preconditions: log_density is log pi(u) and finite (a NaN exhausts the
//   interval silently, -Inf accepts any finite point), width is positive and
//   finite, and direction.size() == u.size(). log_density_at(point) sets the
//   caller's model to the point (in u) when it is in the support and returns
//   log pi there (-Inf outside the support or when a coordinate does not
//   back-transform, possibly without setting the model); uniform() returns a
//   uniform(0, 1) variate and exponential() an exponential(1) variate. Draws:
//   one exponential, then one uniform for the placement, one for the division,
//   and one per shrinkage step. On acceptance u holds the new point, the model
//   was last set to it, and the returned log density can be reused as the
//   current value of the next direction. When the interval is exhausted in
//   floating point (a self-transition of this direction), u is unchanged, the
//   model state is unspecified, and the caller restores its direction-entry
//   snapshot before chaining another direction.
// - Axes: the directions (columns of a d x d matrix) and widths of a block,
//   identity axes and width 1 at the start. adapt(u_start, u) feeds the state
//   after one update (u_start, the state before it, sets the running mean on
//   the first call only): running mean and covariance with weight
//   (n + 1)^-0.7, axes refreshed after 100 updates and then every 50 (cyclic
//   Jacobi eigenvectors, width 2.5 sqrt(eigenvalue); the old axes stay when an
//   eigenvalue is not positive and finite). Nothing changes between adapt()
//   calls.

#ifndef BAYESTOOLS_FACTOR_SLICE_KERNEL_H
#define BAYESTOOLS_FACTOR_SLICE_KERNEL_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace BayesTools {
  namespace factor_slice {

    // Randomized stepping-out limit (Neal 2003, section 4.1).
    constexpr unsigned int stepping_out_limit = 10;
    // Adaptation: the axes are first refreshed after this many adaptive updates
    // and then every refresh_interval updates.
    constexpr unsigned int first_refresh = 100;
    constexpr unsigned int refresh_interval = 50;
    constexpr double width_multiplier = 2.5;

    // log(DBL_MAX) rounded down: the largest |u| whose exp(|u|) is finite.
    // Within it the back-transforms are the original quotients; beyond it they
    // returned the bound (every such proposal was rejected), and the distance
    // from the nearer bound is computed on the log scale instead.
    constexpr double exp_overflow_threshold = 709.782712893383973096;
    constexpr double log_two = 0.693147180559945309417232121458176568;
    // Largest distance of a bounded member's starting value from its anchoring
    // bound at which a caller starts the sampler (the precision limit above).
    constexpr double max_anchor_distance = 1e10;

    // log(b - a) for b > a; log(b / 2 - a / 2) + log(2) only when b - a
    // overflows (a and b of opposite signs near the largest double).
    inline double log_difference(double a, double b)
    {
      const double difference = b - a;
      return std::isfinite(difference) ? std::log(difference) : std::log(0.5 * b - 0.5 * a) + log_two;
    }

    enum class Transform { identity, lower, interval, upper };

    // to_x() rejects a result that is not finite or not strictly inside the
    // support. At the representability limit it can also reject a value that
    // is representable inside: a subnormal distance from the bound is rounded
    // before it is added to the bound, and the one-sided overflow route
    // (exp(u) beyond the largest double, bounds within a factor of two of it)
    // halves first.
    class Coordinate
    {
      Transform _transform;
      double _lower, _upper;
      double _log_span;  // log(upper - lower) with two finite bounds, else 0

    public:
      Coordinate(double lower, double upper)
        : _transform(Transform::identity), _lower(lower), _upper(upper), _log_span(0.0)
      {
        const bool has_lower = std::isfinite(lower);
        const bool has_upper = std::isfinite(upper);
        _transform = has_lower && has_upper ? Transform::interval :
          has_lower ? Transform::lower : has_upper ? Transform::upper : Transform::identity;
        if (_transform == Transform::interval) _log_span = log_difference(lower, upper);
      }

      Transform transform() const { return _transform; }
      double lower() const { return _lower; }
      double upper() const { return _upper; }
      double log_span() const { return _log_span; }

      // Distance of x from the bound that anchors it: the bound of a one-sided
      // member, the nearer bound of an interval, 0 without a finite bound.
      double anchor_distance(double x) const
      {
        switch (_transform) {
        case Transform::lower: return x - _lower;
        case Transform::upper: return _upper - x;
        case Transform::interval: return x - _lower < _upper - x ? x - _lower : _upper - x;
        default: return 0.0;
        }
      }

      double to_u(double x) const
      {
        switch (_transform) {
        case Transform::lower: return log_difference(_lower, x);
        case Transform::interval: return log_difference(_lower, x) - log_difference(x, _upper);
        case Transform::upper: return log_difference(x, _upper);
        default: return x;
        }
      }

      bool to_x(double u, double *x) const
      {
        double out;
        switch (_transform) {
        case Transform::lower:
          out = u <= exp_overflow_threshold ? _lower + std::exp(u) :
            2.0 * (0.5 * _lower + std::exp(u - log_two));
          break;
        case Transform::interval:
          if (std::fabs(u) <= exp_overflow_threshold) {
            const double span = _upper - _lower;
            if (std::isfinite(span)) {
              // the quotient (upper - lower) / (1 + exp(|u|)) from the nearer bound
              out = u >= 0.0 ?
                _upper - span / (1.0 + std::exp(u)) :
                _lower + span / (1.0 + std::exp(-u));
            } else {
              // the same with half the span when the span overflows
              const double half_span = 0.5 * _upper - 0.5 * _lower;
              out = u >= 0.0 ?
                2.0 * (0.5 * _upper - half_span / (1.0 + std::exp(u))) :
                2.0 * (0.5 * _lower + half_span / (1.0 + std::exp(-u)));
            }
          } else {
            // exp(|u|) overflows: the distance from the nearer bound on the log scale
            out = u >= 0.0 ?
              _upper - std::exp(_log_span - u - std::log1p(std::exp(-u))) :
              _lower + std::exp(_log_span + u - std::log1p(std::exp(u)));
          }
          break;
        case Transform::upper:
          out = u <= exp_overflow_threshold ? _upper - std::exp(u) :
            2.0 * (0.5 * _upper - std::exp(u - log_two));
          break;
        default:
          out = u;
        }
        *x = out;
        return std::isfinite(out) && out > _lower && out < _upper;
      }

      double log_jacobian(double u) const
      {
        switch (_transform) {
        case Transform::lower:
        case Transform::upper:
          return u;
        case Transform::interval:
          return -std::fabs(u) - 2.0 * std::log1p(std::exp(-std::fabs(u))) + _log_span;
        default:
          return 0.0;
        }
      }
    };

    inline double log_jacobian(std::vector<Coordinate> const &coordinates, std::vector<double> const &u)
    {
      double out = 0.0;
      for (std::size_t i = 0; i < coordinates.size(); ++i) out += coordinates[i].log_jacobian(u[i]);
      return out;
    }

    struct SliceResult
    {
      bool accepted;       // false: the interval was exhausted in floating point
      double log_density;  // log pi at the accepted point, else the entry value
    };

    template <class LogDensityAt, class Uniform, class Exponential>
    inline SliceResult slice_along(std::vector<double> &u, std::vector<double> const &direction, double width,
                                   double log_density, LogDensityAt &&log_density_at, Uniform &&uniform,
                                   Exponential &&exponential)
    {
      const std::size_t d = u.size();
      std::vector<double> point(d);
      auto along = [&](double t) -> double {
        for (std::size_t i = 0; i < d; ++i) point[i] = u[i] + t * direction[i];
        return log_density_at(static_cast<std::vector<double> const &>(point));
      };
      const double logy = log_density - exponential();
      double left = -width * uniform();
      double right = left + width;
      int steps_left = static_cast<int>(std::floor(stepping_out_limit * uniform()));
      int steps_right = static_cast<int>(stepping_out_limit) - 1 - steps_left;
      while (steps_left > 0 && along(left) > logy) {
        left -= width;
        --steps_left;
      }
      while (steps_right > 0 && along(right) > logy) {
        right += width;
        --steps_right;
      }
      for (;;) {
        // Floating-point exhaustion of the interval: the caller restores its state.
        const double middle = left + 0.5 * (right - left);
        if (!(middle > left && middle < right)) break;
        const double t = left + (right - left) * uniform();
        const double value = along(t);
        if (value > logy) {
          u = point;
          SliceResult accepted = { true, value };
          return accepted;
        }
        if (t < 0.0) left = t;
        else right = t;
      }
      SliceResult exhausted = { false, log_density };
      return exhausted;
    }

    class Axes
    {
      std::size_t _d;
      unsigned int _n;
      bool _mean_initialized;
      std::vector<double> _mean, _cov;  // running moments of u (row-major covariance)
      std::vector<double> _dirs;        // column k is direction k (row-major d x d)
      std::vector<double> _width;

      // Directions become the eigenvectors of the covariance (cyclic Jacobi) and
      // widths 2.5 standard deviations along them; the previous axes stay when
      // an eigenvalue is not positive and finite.
      void refresh()
      {
        const std::size_t d = _d;
        std::vector<double> a(_cov), v(d * d, 0.0);
        for (std::size_t i = 0; i < d; ++i) v[i * d + i] = 1.0;
        for (unsigned int sweep = 0; sweep < 100; ++sweep) {
          double off = 0.0, diagonal = 0.0;
          for (std::size_t p = 0; p < d; ++p) {
            diagonal += a[p * d + p] * a[p * d + p];
            for (std::size_t q = p + 1; q < d; ++q) off += a[p * d + q] * a[p * d + q];
          }
          if (!(off > 1e-30 * diagonal)) break;
          for (std::size_t p = 0; p < d; ++p) {
            for (std::size_t q = p + 1; q < d; ++q) {
              const double apq = a[p * d + q];
              if (apq == 0.0) continue;
              const double theta = (a[q * d + q] - a[p * d + p]) / (2.0 * apq);
              const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
              const double c = 1.0 / std::sqrt(t * t + 1.0);
              const double s = t * c;
              for (std::size_t k = 0; k < d; ++k) {
                const double akp = a[k * d + p], akq = a[k * d + q];
                a[k * d + p] = c * akp - s * akq;
                a[k * d + q] = s * akp + c * akq;
              }
              for (std::size_t k = 0; k < d; ++k) {
                const double apk = a[p * d + k], aqk = a[q * d + k];
                a[p * d + k] = c * apk - s * aqk;
                a[q * d + k] = s * apk + c * aqk;
              }
              for (std::size_t k = 0; k < d; ++k) {
                const double vkp = v[k * d + p], vkq = v[k * d + q];
                v[k * d + p] = c * vkp - s * vkq;
                v[k * d + q] = s * vkp + c * vkq;
              }
            }
          }
        }
        for (std::size_t k = 0; k < d; ++k) {
          const double lambda = a[k * d + k];
          if (!(lambda > 0.0) || !std::isfinite(lambda)) return;
        }
        _dirs = v;
        for (std::size_t k = 0; k < d; ++k) {
          _width[k] = width_multiplier * std::sqrt(a[k * d + k]);
        }
      }

    public:
      explicit Axes(std::size_t d)
        : _d(d), _n(0), _mean_initialized(false), _mean(d, 0.0), _cov(d * d, 0.0),
          _dirs(d * d, 0.0), _width(d, 1.0)
      {
        for (std::size_t i = 0; i < d; ++i) _dirs[i * d + i] = 1.0;
      }

      std::size_t dimension() const { return _d; }
      // number of adapt() calls so far
      unsigned int updates() const { return _n; }
      double width(std::size_t k) const { return _width[k]; }
      std::vector<double> direction(std::size_t k) const
      {
        std::vector<double> out(_d);
        for (std::size_t i = 0; i < _d; ++i) out[i] = _dirs[i * _d + k];
        return out;
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
        for (std::size_t i = 0; i < _d; ++i) {
          delta[i] = u[i] - _mean[i];
          _mean[i] += weight * delta[i];
        }
        // (1 - w) C + w (1 - w) delta delta', symmetric by construction.
        for (std::size_t i = 0; i < _d; ++i) {
          for (std::size_t j = 0; j <= i; ++j) {
            const double value = _cov[i * _d + j] +
              weight * ((1.0 - weight) * delta[i] * delta[j] - _cov[i * _d + j]);
            _cov[i * _d + j] = value;
            _cov[j * _d + i] = value;
          }
        }
        if (_n >= first_refresh && (_n - first_refresh) % refresh_interval == 0) refresh();
      }
    };

  }
}

#endif
