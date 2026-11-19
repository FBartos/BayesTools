// A toy sampler on inst/include/BayesTools/factor_slice_kernel.h for
// test-factor-slice-kernel.R. The header comes first and the rest is the
// standard library, so compiling this file shows that the header is
// self-contained; the entry point uses the .C convention (no R headers).
#include <BayesTools/factor_slice_kernel.h>

#include <cmath>
#include <limits>
#include <random>
#include <vector>

namespace {

  namespace kernel = ::BayesTools::factor_slice;

  // log density (up to a constant) of a bivariate normal with means m, standard
  // deviations s and correlation r
  double log_bivariate_normal(std::vector<double> const &x, double const *m, double const *s, double r)
  {
    const double z1 = (x[0] - m[0]) / s[0];
    const double z2 = (x[1] - m[1]) / s[1];
    return -0.5 * (z1 * z1 - 2.0 * r * z1 * z2 + z2 * z2) / (1.0 - r * r);
  }

}

// One chain of 'n_adapt' adaptive and 'n_iter' retained updates, each a sweep
// of factor slice updates along the block's directions as BayesTools'
// FactorSlice sampler makes them; 'draws' (n_iter x 2, column-major) receives
// the retained states and 'exhausted' the number of exhausted intervals.
extern "C" void bt_factor_slice_kernel_toy(int const *n_adapt, int const *n_iter, int const *seed,
                                           double const *mean, double const *sd, double const *correlation,
                                           double *draws, int *exhausted)
{
  const double infinity = std::numeric_limits<double>::infinity();
  const std::size_t d = 2;
  std::vector<kernel::Coordinate> coordinates(d, kernel::Coordinate(-infinity, infinity));
  kernel::Axes axes(d);
  std::mt19937_64 generator(static_cast<unsigned long long>(*seed));
  std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);
  std::exponential_distribution<double> exponential_distribution(1.0);
  auto uniform = [&]() { return uniform_distribution(generator); };
  auto exponential = [&]() { return exponential_distribution(generator); };

  std::vector<double> x(mean, mean + d), u(d);
  for (std::size_t i = 0; i < d; ++i) u[i] = coordinates[i].to_u(x[i]);
  auto log_density_at = [&](std::vector<double> const &point) {
    for (std::size_t i = 0; i < d; ++i) {
      if (!std::isfinite(point[i]) || !coordinates[i].to_x(point[i], &x[i])) return -infinity;
    }
    return log_bivariate_normal(x, mean, sd, *correlation) + kernel::log_jacobian(coordinates, point);
  };

  *exhausted = 0;
  for (int iteration = 0; iteration < *n_adapt + *n_iter; ++iteration) {
    const std::vector<double> u_start(u);
    for (std::size_t k = 0; k < d; ++k) {
      const std::vector<double> x_entry(x);
      const double log_density = log_bivariate_normal(x, mean, sd, *correlation) +
        kernel::log_jacobian(coordinates, u);
      const kernel::SliceResult result = kernel::slice_along(
        u, axes.direction(k), axes.width(k), log_density, log_density_at, uniform, exponential);
      if (!result.accepted) {
        x = x_entry;
        ++*exhausted;
      }
    }
    if (iteration < *n_adapt) {
      axes.adapt(u_start, u);
    } else {
      const int row = iteration - *n_adapt;
      for (std::size_t i = 0; i < d; ++i) draws[row + static_cast<int>(i) * *n_iter] = x[i];
    }
  }
}

// resolution() and resolvable() for a member with support (lower, upper) at
// native value x of a normal target N(mean, sd^2) (the log density of the
// other members is constant), and the coordinate's u_resolution(x).
extern "C" void bt_factor_slice_resolution_toy(double const *lower, double const *upper, double const *x,
                                               double const *mean, double const *sd, int *judged, double *curvature,
                                               double *allowance, double *round_trip, int *resolvable,
                                               double *u_resolution)
{
  const kernel::Coordinate coordinate(*lower, *upper);
  auto log_density = [&](double value) {
    const double z = (value - *mean) / *sd;
    return -0.5 * z * z;
  };
  const kernel::Resolution r = kernel::resolution(coordinate, *x, log_density);
  *judged = r.judged ? 1 : 0;
  *curvature = r.curvature;
  *allowance = r.allowance;
  *round_trip = r.round_trip;
  *resolvable = kernel::resolvable(coordinate, *x, log_density) ? 1 : 0;
  *u_resolution = coordinate.u_resolution(*x);
}

// exp() is finite at the overflow threshold of the back-transforms and
// overflows at the next double.
extern "C" void bt_factor_slice_exp_threshold_toy(int *finite_at, int *infinite_above)
{
  const double threshold = kernel::exp_overflow_threshold;
  *finite_at = std::isfinite(std::exp(threshold)) ? 1 : 0;
  *infinite_above = std::isfinite(std::exp(std::nextafter(threshold, std::numeric_limits<double>::infinity()))) ? 0 : 1;
}
