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
//   it over a block in member order. Precision limit: the sampler reaches a
//   bounded member to about 16 significant digits of its distance from the
//   anchoring bound (for an interval, of the smaller distance, and of the span
//   near the centre), as a double holds a value to 16 digits of its magnitude.
//   A posterior with features finer than that grid cannot be sampled on this
//   scale; two local, best-effort diagnostics stop its common forms (a slice
//   interval that collapses on collapse_limit consecutive updates, and
//   resolution()), not a guarantee.
// - resolution(coordinate, x, log_density_at_x) and resolvable(): a local
//   diagnostic of the grid of reachable values at the member's current native
//   value x (strictly inside the support; a value on a bound is the caller's
//   own error): the grid step relative to the local scale of the density
//   (curvature) and the round trip of x. resolvable() is true also when
//   nothing could be judged. log_density_at_x(value) sets the member to the
//   native value (the other members unchanged) and returns the block's full
//   log conditional there; the caller restores the member afterwards, also
//   when the callback throws. No random numbers are drawn. A JAGS sampler that
//   runs it at adaptOff() must defer a failure to its next update: JAGS
//   deletes the model when adaptOff() throws. resolvability_due() is the
//   schedule after adaptation.
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
//   direction. With a finite log density at the current point an exhausted
//   interval means that no reachable point near it lies in the slice; on
//   collapse_limit consecutive updates of the same direction the direction is
//   frozen, and the sampler stops the fit.
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
    // resolution(): the grid step h of the reachable values at x is too coarse
    // when |f''| h^2 (f the log full conditional) exceeds this, i.e. when h
    // exceeds about 0.05 of the local scale 1 / sqrt(|f''|) of the density.
    // Calibrated on a standard normal: bounds +-1e13 (h about 0.02) pass,
    // +-1e14 (h about 0.12) stop. Near the limit a probability can be off in
    // the third decimal (with a lower bound at -2e12, P(|x| <= 1) was measured
    // as 0.6857 +- 0.0010 against 0.6827; the check passes it).
    constexpr double curvature_tolerance = 2.5e-3;
    // resolution(): the current value must survive the round trip x -> u -> x
    // to within this many grid steps.
    constexpr double round_trip_steps = 4.0;
    // resolution(): rounding allowance of each log density, in units of the
    // spacing of the doubles at the largest of the three |f|.
    constexpr double log_density_rounding_ulps = 4.0;
    // After adaptation, the samplers repeat the check after 50, 100 and 200
    // updates and then every resolvability_interval updates
    // (resolvability_due()): a member that starts where its coordinate resolves
    // the target (next to a distant bound) can drift to where it does not. One
    // check costs a few full-conditional evaluations per member; every 500
    // updates that is negligible even for RoBMA's corrected block, which makes
    // only two full evaluations per update (every 50 would cost it about a
    // quarter of its time). A chain that freezes earlier stops through its
    // consecutive collapsed slice intervals (collapse_limit).
    // A sampler stops the fit when the same direction's slice interval
    // collapses (slice_along() reports an exhausted interval) on this many
    // consecutive updates; a single collapse is a self-transition of that
    // direction. Consecutive, because isolated collapses occur in valid fits
    // near the precision limit (3 to 5 in about 20,000 direction-updates for a
    // standard normal on (-1e13, 1e13)), while a frozen chain collapses on
    // every update.
    constexpr unsigned int collapse_limit = 10;
    constexpr unsigned int resolvability_early_checks[3] = { 50, 100, 200 };
    constexpr unsigned int resolvability_interval = 500;

    // Whether the check is due after this many updates after adaptation.
    inline bool resolvability_due(unsigned int updates)
    {
      for (unsigned int early : resolvability_early_checks) {
        if (updates == early) return true;
      }
      return updates > 0 && updates % resolvability_interval == 0;
    }

    // log(b - a) for b > a; log(b / 2 - a / 2) + log(2) only when b - a
    // overflows (a and b of opposite signs near the largest double).
    inline double log_difference(double a, double b)
    {
      const double difference = b - a;
      return std::isfinite(difference) ? std::log(difference) : std::log(0.5 * b - 0.5 * a) + log_two;
    }

    enum class Transform { identity, lower, interval, upper };

    // Spacing of the doubles at |value|.
    inline double spacing(double value)
    {
      const double magnitude = std::fabs(value);
      return std::nextafter(magnitude, std::numeric_limits<double>::infinity()) - magnitude;
    }

    // to_x() rejects a result that is not finite or not strictly inside the
    // support. At the representability limit it can also reject a value that
    // is representable inside: a subnormal distance from the bound is rounded
    // before it is added to the bound, and the one-sided overflow route
    // (exp(u) beyond the largest double, bounds within a factor of two of it)
    // halves first. A state that then cannot move stops the sampler through its
    // consecutive collapsed slice intervals (collapse_limit).
    class Coordinate
    {
      Transform _transform;
      double _lower, _upper;
      double _log_span;  // log(upper - lower) with two finite bounds, else 0

      // Absolute rounding error of log(b - a) for b > a: b - a is known to about
      // eps max(|a|, |b|) relative to the distance (the operands are doubles of
      // that magnitude, which cancellation exposes), and the log adds half the
      // spacing of its result.
      static double log_distance_error(double a, double b)
      {
        double distance = b - a;
        double magnitude = std::fmax(std::fabs(a), std::fabs(b));
        if (!std::isfinite(distance)) {
          distance = 0.5 * b - 0.5 * a;
          magnitude *= 0.5;
        }
        return std::numeric_limits<double>::epsilon() * magnitude / distance +
          0.5 * spacing(log_difference(a, b));
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

      // Absolute rounding error of to_u(x): the sum of the errors of its log
      // terms (log_distance_error()), not zero at the centre of a symmetric
      // interval, where both terms are near zero; 0 for the identity.
      double u_resolution(double x) const
      {
        switch (_transform) {
        case Transform::lower: return log_distance_error(_lower, x);
        case Transform::interval: return log_distance_error(_lower, x) + log_distance_error(x, _upper);
        case Transform::upper: return log_distance_error(x, _upper);
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

    struct Resolution
    {
      bool judged;        // false: nothing could be judged (see resolution())
      double curvature;   // |f''| h^2 at x, h the larger grid step around x
      double allowance;   // its rounding allowance from the rounded log densities
      double round_trip;  // |to_x(to_u(x)) - x| / h (infinite when it fails)
    };

    // A local diagnostic of the grid of values the sampler can reach around the
    // member's current native value x, which must lie strictly inside the
    // support. With the nearest distinct back-transforms below and above x of
    // u -/+ delta (delta = the larger of the spacing of u and the rounding error
    // of to_u, doubled until the back-transform moves), the grid steps h1, h2
    // and the log full conditional f at the three points, it measures the
    // second divided difference f'' and reports |f''| max(h1, h2)^2 (the grid
    // step relative to the local scale of the density, squared) and how far x
    // moves in the round trip x -> u -> x. It judges curvature, not the first
    // difference, so a state far in a tail of a sharp posterior is not
    // mistaken for a coarse grid; a step of the density (a discontinuity) next
    // to x reads as curvature, so a member started exactly on a step of its
    // likelihood stops (a false alarm). It sees only these three points: a feature that
    // lies between two reachable values and is never visited cannot be seen,
    // so this is best effort, not a guarantee.
    // log_density_at_x(value) sets the member to the native value (the other
    // members unchanged) and returns the block's full log conditional there
    // (every member's density and every child's likelihood, as the sampler
    // targets it); the caller restores the member afterwards, also when the
    // callback throws. No random numbers are drawn. Not judged (judged false):
    // an identity coordinate, a non-finite f at x (the sampler's own checks
    // report it), a neighbour that is not representable inside the support,
    // and a non-finite f at a neighbour (an edge of the posterior's support).
    // A current value whose u is not finite or that does not back-transform is
    // judged with an infinite round trip. Not handled: an additive constant so
    // large that the log densities lose their differences breaks the slice
    // arithmetic itself (as on development); the allowance then also exceeds
    // the tolerance and the curvature cannot be judged.
    template <class LogDensityAtX>
    inline Resolution resolution(Coordinate const &coordinate, double x, LogDensityAtX &&log_density_at_x)
    {
      Resolution out = { false, 0.0, 0.0, 0.0 };
      if (coordinate.transform() == Transform::identity) return out;
      const double f = log_density_at_x(x);
      if (!std::isfinite(f)) return out;
      const double u = coordinate.to_u(x);
      double back;
      if (!std::isfinite(u) || !coordinate.to_x(u, &back)) {
        out.judged = true;
        out.round_trip = std::numeric_limits<double>::infinity();
        return out;
      }
      const double delta = std::fmax(spacing(u), coordinate.u_resolution(x));
      // the nearest distinct back-transform on each side of u
      double neighbour[2];
      for (int side = 0; side < 2; ++side) {
        const double sign = side == 0 ? -1.0 : 1.0;
        bool found = false;
        double step = delta;
        for (int doubling = 0; doubling < 64 && !found; ++doubling, step *= 2.0) {
          double value;
          if (!coordinate.to_x(u + sign * step, &value)) return out;
          if (value != x) {
            neighbour[side] = value;
            found = true;
          }
        }
        if (!found) return out;
      }
      // the three points in increasing order with their log densities
      double point[3] = { neighbour[0], x, neighbour[1] };
      double value[3] = { log_density_at_x(neighbour[0]), f, log_density_at_x(neighbour[1]) };
      if (!std::isfinite(value[0]) || !std::isfinite(value[2])) return out;
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2 - i; ++j) {
          if (point[j] > point[j + 1]) {
            const double p = point[j], v = value[j];
            point[j] = point[j + 1];
            value[j] = value[j + 1];
            point[j + 1] = p;
            value[j + 1] = v;
          }
        }
      }
      const double h1 = point[1] - point[0], h2 = point[2] - point[1];
      if (!(h1 > 0.0) || !(h2 > 0.0)) return out;
      const double h = std::fmax(h1, h2);
      // |f''| h^2 = 2 h / (h1 + h2) |(f2 - f1) h / h2 - (f1 - f0) h / h1|,
      // arranged so that no h^2 overflows; each f carries a rounding error of
      // log_density_rounding_ulps spacings of the largest |f|
      const double scale = 2.0 * h / (h1 + h2);
      const double rounding = log_density_rounding_ulps *
        spacing(std::fmax(std::fabs(value[0]), std::fmax(std::fabs(value[1]), std::fabs(value[2]))));
      out.judged = true;
      out.curvature = scale * std::fabs((value[2] - value[1]) * (h / h2) - (value[1] - value[0]) * (h / h1));
      out.allowance = scale * 2.0 * rounding * (h / h2 + h / h1);
      out.round_trip = std::fabs(back - x) / h;
      return out;
    }

    // resolution() within its tolerances (curvature beyond the rounding
    // allowance, round trip within round_trip_steps), or not judged.
    template <class LogDensityAtX>
    inline bool resolvable(Coordinate const &coordinate, double x, LogDensityAtX &&log_density_at_x)
    {
      const Resolution r = resolution(coordinate, x, log_density_at_x);
      return !r.judged || (r.curvature <= curvature_tolerance + r.allowance && r.round_trip <= round_trip_steps);
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
