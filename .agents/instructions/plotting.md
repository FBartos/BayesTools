# Plotting

Use this guide for prior plots, posterior/model-averaging plots, JAGS
diagnostics, plot-data helpers, and base/ggplot renderers.

## Architecture

Follow the architecture of the plot family being changed. For plot families
supporting both base graphics and ggplot2, keep statistical computation separate
from rendering:

1. the public method validates and dispatches;
2. a focused helper computes backend-neutral plot data;
3. base and ggplot renderers consume the same computed values.

Do not introduce this layering for a trivial single-backend plot, and do not
force unrelated plot families into one abstraction. Normalize shared defaults
once. Backend-specific styling must not change the represented quantities or
semantic layer order.

Use the established `plot_type = "base"` and `plot_type = "ggplot"` contract
where both backends exist. Preserve each public function's documented return
value; ggplot paths return plot objects, while base paths may return invisible
metadata or `NULL` according to the existing family.

Prior plot dispatch and layers live in `R/priors-plot.R` and
`R/priors-plot-layers.R`. Model-averaging plot families live in
`R/model-averaging-plots*.R`. JAGS diagnostic data and rendering live in
`R/JAGS-diagnostics.R`. Inspect the current family before adding helpers or
source files.

Plotting density grids are approximations for display. They are not evidence
for the structural classifications returned by `prior_density_ordinate()`.

## Verification

Test computed plot data, validation, and dispatch separately, then retain
representative human-reviewed `vdiffr` snapshots for both backends where both
are public. Structural tests and render-only checks do not replace visual
regression.

Use the `visual` profile for pure plots and `visual-fixture` for plots that load
cached JAGS fits. Run `fit` first when the required fitted cache is affected.
Never auto-accept snapshots. Make stochastic plot inputs deterministic and ask
the maintainer to review every intentional visual difference.
