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
//   it over a block in member order. A value is reached to a relative
//   precision of about |u| eps of its distance from the anchoring bound (to_u
//   rounds its log terms), so a member whose posterior is narrower than that
//   cannot be sampled on this scale: resolvable() detects it.
// - resolvable(coordinate, x, log_density_at_x): whether the doubles of u
//   around to_u(x) resolve the target at the member's current native value x.
//   log_density_at_x(value) sets the member to the native value (the other
//   members unchanged) and returns the log full conditional there; the caller
//   restores the member afterwards. No random numbers are drawn.
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
//   floating point, u is unchanged, the model state is unspecified, and the
//   caller restores its direction-entry snapshot before chaining another
//   direction.
// - Axes: the directions (columns of a d x d matrix) and widths of a block,
//   identity axes and width 1 at the start. adapt(u_start, u) feeds the state
//   after one update (u_start, the state before it, sets the running mean on
//   the first call only) and returns true when the axes were due for a
//   refresh: running mean and covariance with weight (n + 1)^-0.7, axes
//   refreshed after 100 updates and then every 50 (cyclic Jacobi eigenvectors,
//   width 2.5 sqrt(eigenvalue); the old axes stay when an eigenvalue is not
//   positive and finite). Nothing changes between adapt() calls.

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
    // Largest change of the log full conditional allowed between the current
    // value and the values of the neighbouring coordinates (resolvable()).
    constexpr double resolvability_tolerance = 0.01;

    // log(b - a) for b > a; log(b / 2 - a / 2) + log(2) only when b - a
    // overflows (a and b of opposite signs near the largest double).
    inline double log_difference(double a, double b)
    {
      const double difference = b - a;
      return std::isfinite(difference) ? std::log(difference) : std::log(0.5 * b - 0.5 * a) + log_two;
    }

    enum class Transform { identity, lower, interval, upper };

    // Values that are not representable strictly inside the support are
    // rejected by to_x(): a distance from the bound below half the spacing of
    // the doubles there (also a subnormal distance rounded before it is added to
    // the bound), or a value beyond the largest double (also a one-sided bound
    // within a factor of two of the largest double, whose overflow route halves
    // first). These are representability limits; a member whose posterior the
    // coordinate cannot resolve stops the fit through resolvable().
    class Coordinate
    {
      Transform _transform;
      double _lower, _upper;
      double _log_span;  // log(upper - lower) with two finite bounds, else 0

      static double spacing(double value)
      {
        const double magnitude = std::fabs(value);
        return std::nextafter(magnitude, std::numeric_limits<double>::infinity()) - magnitude;
      }

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

      double to_u(double x) const
      {
        switch (_transform) {
        case Transform::lower: return log_difference(_lower, x);
        case Transform::interval: return log_difference(_lower, x) - log_difference(x, _upper);
        case Transform::upper: return log_difference(x, _upper);
        default: return x;
        }
      }

      // Absolute rounding error of to_u(x), about the spacing of the doubles at
      // its largest log term; 0 for the identity.
      double u_resolution(double x) const
      {
        switch (_transform) {
        case Transform::lower: return spacing(log_difference(_lower, x));
        case Transform::interval: {
          const double from_lower = spacing(log_difference(_lower, x));
          const double from_upper = spacing(log_difference(x, _upper));
          return from_lower > from_upper ? from_lower : from_upper;
        }
        case Transform::upper: return spacing(log_difference(x, _upper));
        default: return 0.0;
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

    // The coordinate resolves the target at the member's native value x when
    // the value survives the round trip x -> u -> x and the log full
    // conditional at the round-tripped value, and at the values of the
    // neighbouring coordinates (the adjacent doubles of u and u -/+ the rounding
    // error of to_u), differs from its value at x by at most
    // resolvability_tolerance: the points the sampler can reach are then closer
    // than the scale on which the target changes. Otherwise the member's support
    // bounds are too far from where its posterior concentrates for this
    // coordinate (a representability limit, e.g. a standard normal truncated to
    // (-5e14, 5e14) or bounded below by -1e16). An identity coordinate, or a
    // state whose log full conditional is not finite (the sampler's own checks
    // report that), is not judged. A current value whose u is not finite or that
    // does not survive the round trip is not resolvable.
    template <class LogDensityAtX>
    inline bool resolvable(Coordinate const &coordinate, double x, LogDensityAtX &&log_density_at_x)
    {
      if (coordinate.transform() == Transform::identity) return true;
      const double reference = log_density_at_x(x);
      if (!std::isfinite(reference)) return true;
      const double u = coordinate.to_u(x);
      double value;
      if (!std::isfinite(u) || !coordinate.to_x(u, &value)) return false;
      if (!(std::fabs(log_density_at_x(value) - reference) <= resolvability_tolerance)) return false;
      const double infinity = std::numeric_limits<double>::infinity();
      const double step = coordinate.u_resolution(x);
      const double neighbours[4] = { std::nextafter(u, -infinity), std::nextafter(u, infinity), u - step, u + step };
      for (double neighbour : neighbours) {
        if (coordinate.to_x(neighbour, &value) &&
            !(std::fabs(log_density_at_x(value) - reference) <= resolvability_tolerance)) return false;
      }
      return true;
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

      bool adapt(std::vector<double> const &u_start, std::vector<double> const &u)
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
        if (_n >= first_refresh && (_n - first_refresh) % refresh_interval == 0) {
          refresh();
          return true;
        }
        return false;
      }
    };

  }
}

#endif
