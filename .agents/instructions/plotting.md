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

Posterior plots draw point masses only from declared posterior atoms (the
producers' metadata or `posterior_atom_attribute()`), also when a precomputed
posterior density supplies the continuous curve; samples without an atom
declaration stop with the atom-status message. Never infer point masses
from prior lists, component indices, or draws matching prior spikes.

Mixed continuous-and-point plots keep one probability-axis mapping. Base
overlays reuse `par("usr")` and warn when a later point mass is off-scale.
ggplot overlays reuse the mixed plot's stored secondary-axis mapping and issue
the same clipping warning. Do not invent a fake `usr` for ggplot.

Individual weight-function omega plots index the k-th omega in the order the
summary tables (and RoBMA's) print them: ascending p-value intervals,
reference weight first. Every individual-omega selector (`show_figures`,
`show_parameter`) follows that order.

Prior plot dispatch and layers live in `R/priors-plot.R` and
`R/priors-plot-layers.R`. Model-averaging plot families live in
`R/model-averaging-plots*.R`. JAGS diagnostic data and rendering live in
`R/JAGS-diagnostics.R`. Inspect the current family before adding helpers or
source files.

Ordered prior display suppression uses the shared exact-route predicate: omit
the whole continuous curve only for allocation-induced infinity with no
intrinsically infinite total; keep atoms and unknown classifications. Display
decisions precede ordered density generation through a transient logical mask
on a bound prior copy; skip only hidden curves before their grid or KDE work,
keep original component indices, and sample all primitives in their original
order whenever the existing path samples. Ordinary public density inputs
carry no mask and still evaluate every level. Display
densities missing their atom table recover locations from declared route
nodes and probabilities from exact route ordinates. Apply output transformations
only to recovered locations; preserve already transformed atoms and never scale
their masses by a Jacobian. Unknown recovery raises the typed display-unavailable
condition rather than erasing or guessing atoms.
Direct fixed-total allocation shares check both exact support endpoints. Direct
prior selectors retain reference numbering; multiple ordered ggplot selections
return positional lists with NULL holes for unselected or suppressed figures,
and only a single selected visible figure collapses to a ggplot. Factor overlays omit persisted
zero-design references and carry `factor_level_universe` before curve omission,
so retained styles and labels cannot shift. Ordered overlays match posterior
colors with dashed defaults and one posterior-owned legend; mapped ordered
posterior layers use automatic glyph participation when the legend is enabled.
Empty overlays
skip quietly; all-hidden standalone selections fail with the documented class.

`plot_models()` uses semantic ordered levels and transformed model summaries.
Posterior-only plots skip prior summaries and transforms of unused prior fields.
Ordered prior means use declared
allocation expectations; intervals use structural CDF diagnostics, exact atom
jumps and checked generalized inversion, or a typed model/level limitation.

Plotting density grids are approximations for display. They are not evidence
for the structural classifications returned by `prior_density_ordinate()`.
Keep requested display limits separate from transformation support. Custom
transformations can declare `output_support = c(lower, upper)`; inverse-grid
evaluation must exclude coordinates outside that interval, without clamping
actual values or moving point masses. An infinite inverse at a boundary does
not supply a finite continuous-density ordinate.

## Verification

Test computed plot data, validation, and dispatch separately, then retain
representative human-reviewed `vdiffr` snapshots for both backends where both
are public. Structural tests and render-only checks do not replace visual
regression.

Use the `visual` profile for pure plots and `visual-fixture` for plots that load
cached JAGS fits. Run `fit` first when the required fitted cache is affected.
Make stochastic plot inputs deterministic. Retain changed candidates for
maintainer or explicitly delegated review before accepting a baseline change.
