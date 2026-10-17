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
A declared continuous component (or factor level) with fewer than two
distinct continuous draws has no kernel density: `plot_posterior()` (simple,
factor and omega, base and ggplot) and `plot_marginal()` draw its declared atoms,
omit its curve and warn with `BayesTools_posterior_curve_unavailable` (parent
`BayesTools_plot_condition`; fields `parameter`, `continuous_mass`, `n_draws`,
`n_distinct`; the remedy is a larger `n_samples` in `mix_posteriors()`). Prior
overlays are unaffected. Without declared atoms to draw, the plot still stops with
the posterior-density-unavailable message.
Marginal level extraction checks the original typed atom refusal before
discarding draw metadata when no certified scalar law is available.
A linear transformation of declared atoms (factor level cells, formula rows)
declares `{0, mass 1}` for every row of its persisted design that is
identically zero, also when the source declares no scalar marginals, so such
level cells are drawn as spikes at zero; the design decides, never the draws. A usable
independently certified scalar law takes precedence over a joint-law refusal.
Point layers check transformation images before emitting locations and masses.
Finite-source saturation or nonfinite images stop with
`BayesTools_transformation_image_unavailable` (parent `BayesTools_transformation`,
call NULL; transformation/source_values/images fields identify failing entries).
True-zero positive powers and mathematical infinite support limits stay valid;
continuous-grid handling and probability masses are unchanged.

Mixed continuous-and-point plots keep one probability-axis mapping. Base
overlays reuse `par("usr")` and warn when a later point mass is off-scale.
ggplot overlays reuse the mixed plot's stored secondary-axis mapping and issue
the same clipping warning. Do not invent a fake `usr` for ggplot.

Individual weight-function omega plots index the k-th omega in the order the
summary tables (and RoBMA's) print them: ascending p-value intervals,
reference weight first. Every individual-omega selector (`show_figures`,
`show_parameter`) follows that order.

Individual PET-PEESE prior-list plots (`individual = TRUE` in `plot_prior_list()`,
`lines_prior_list()` and `geom_prior_list()`) draw the density of one
coefficient, with the `prior_none()` models as an atom at zero, and need no
`prior_list_mu`: a supplied one is refused with a message that names it as
used only by the regression bands (`individual = FALSE`). A list with both PET
and PEESE priors is refused: its two priors belong to different coefficients.

Prior plot dispatch and layers live in `R/priors-plot.R` and
`R/priors-plot-layers.R`; `plot.prior_linear_density()` in `R/priors-plot.R`
draws a prior linear density from `.prior_linear_density_to_plot_data()` through
`.plot_prior_list.both()`. Model-averaging plot families live in
`R/model-averaging-plots*.R`. JAGS diagnostic data and rendering live in
`R/JAGS-diagnostics.R`. Inspect the current family before adding helpers or
source files.

The density diagnostics reflect at the support of the whole prior: a mixture
takes the union of the continuous supports of its components (points are
atoms, orthonormal and meandif components the real line; a factor component
only resolves the columns), every column of a bias prior (a mixture or a single
prior) takes the bounds of its own component (the omega weights, alpha, pi_null,
and the union of the truncations of the PET or PEESE components), never those of
another component, and a spike-and-slab masks the excluded draws in every
coefficient column. The density and
autocorrelation diagnostics skip columns the prior declares constant and the
draws confirm (fixed weight bins); a parameter with nothing else to show behaves
as a spike prior.

Ordered plots omit only persisted zero-design reference levels by default.
Retain every non-reference prior curve and genuine zero or point level. An
exact infinite density point is not an atom: omit only nonfinite plotting
ordinates without changing the finite curve, its grid, or its mass. Public
`density()` and `prior_density_ordinate()` retain their numerical and
structural meaning. Unknown and unsupported density routes keep their
existing explicit limitations.
Finite displayed peak heights depend on the evaluation grid and do not
represent an infinite mathematical peak. In `plot_posterior()` and
`plot_marginal()` prior overlays and in `plot_transformed_prior()`, unresolved
raw numerical evaluations warn once per continuous route with
`BayesTools_prior_curve_unavailable` (also `BayesTools_plot_condition`),
distinguishing partially available curves from wholly unavailable ones; keep
every available point and declared atom.
The warning's `unresolved_values` field retains evaluation coordinates on the
source scale, before any plotting transformation. Exact infinities, exact
zero/support-zero densities and intentional transformation limits do not
trigger that warning. Direct `plot.prior()`, `lines.prior()` and `geom_prior()`
filter nonfinite plotting ordinates without that numerical-unavailability
warning; this covers ordered levels, simple priors, simplex marginals and
individual omegas, and the display y range is `c(0, max(finite y))`
(`density()` and its `y_range` are unchanged). A display value of an ordered
level whose structural height is unresolved (NA inside the support) is not such
an ordinate: `density()` leaves it out of the level and warns once with
`BayesTools_prior_curve_unavailable` (fields `unresolved_values`,
`structural_values`, `levels`), and `plot()` passes the warning on. A direct
ordered curve without any finite ordinate gives a classed error with the remedy
`x_seq`; any other direct continuous curve gives the same class with the remedy
`xlim`.
Direct prior selectors retain original level numbering; explicit selectors
can select the reference. `show_figures` of ordered and simplex plots holds
non-zero indices of one sign within the number of figures (validated before any
indexing); `show_parameter` of the `lines()` and `geom_prior()` overlays one
positive index. Multiple ordered ggplot selections return positional
lists with NULL holes for unselected figures, including the default reference,
and one explicit selection returns a ggplot. Factor overlays omit persisted
zero-design references and carry `factor_level_universe`, so styles and labels
cannot shift. Ordered overlays match posterior colors with dashed defaults and
one posterior-owned legend; mapped ordered posterior layers use automatic glyph
participation when the legend is enabled. Mixed-measure displays retain exact
atoms with unchanged probability masses and the active probability-axis mapping.
An ordered mixed-measure level honours `ylim2` and `ylab2` through the same
mapping as the mixed plots of density and point priors and stores it for
overlays; without `ylim2` the masses keep their unit scale.
Whole-factor overlays select attached current-scale per-level prior laws before
reconstructing a law from the original factor sources.
Direct overlay mapping and clipping checks inspect selected components only;
an omitted reference cannot create a clipping warning. Explicit reference and
genuine non-reference atom selections still use the active mapping and warn
when their probabilities exceed it.

`plot_models()` uses semantic ordered levels and transformed model summaries;
a decreasing transformation orders the endpoints of every mapped interval (model,
prior and overall) before labels and rendering.
Posterior-only plots skip prior summaries and transforms of unused prior fields.
Ordered prior means use declared
allocation expectations; intervals use structural CDF diagnostics, exact atom
jumps and checked generalized inversion, or a typed model/level limitation.

Plotting density grids are approximations for display. They are not evidence
for the structural classifications returned by `prior_density_ordinate()`.
Automatic prior-coordinate ranges are source ranges, independent of
`transformation_settings`; validate supplied settings before deriving that
interpretation. Preserve the existing range and pretty expansion, replacing
only an endpoint with an infinite transformed image by the existing marginal
tail quantile (`x_range_quant`, or each prior's default). Prior lists use the
union of participating component quantiles on one common source grid, after
their existing weight and sampling-allocation filters. Do not pretty-expand
the corrected endpoint, move atoms, or change support, weights, Jacobians or
sampling budgets. Automatic transformed range metadata is ordered without
separating knots from ordinates; explicit coordinates retain their meaning.
Nonindividual PET/PEESE standard-error and weightfunction p-value axes keep
their specialized dispatch.
Explicit precomputed posterior plots require a matching valid continuous curve
for every selected continuous leaf, resolved before discarding parent parameter
and condition metadata. Pure atoms need no curve; KDE ignores stored curves.
Individual omega plots consume producer-owned scalar atom/support declarations,
including empty and unknown marginals, and never recover atom masses from prior
shares or allocated draw counts. They follow `density_method` (a precomputed
request resolves the stored curve of the omega with its parent samples and
condition; a pure atom needs none) and apply `transformation` to the posterior
and prior data alike: checked atom maps keep the masses, continuous curves use
the checked grid with the Jacobian. Transformed unsupported scalar prior overlays
raise the existing prior-density limitation.
Analytic prior-list plots (`force_samples = FALSE`) keep every positive-weight
component and atom and do not use `n_samples`; only sampled components that
receive fewer than two draws are dropped. A Bernoulli component of a prior list
enters as two atoms at 0 and 1 (its truncation-normalized pmf times its weight),
which merge with equal atoms of other components. Prior-list plots (including
`plot()` of a prior mixture) refuse ordered priors and their mixtures; a single
bound ordered prior is plotted with `plot()`. They also refuse Dirichlet simplex
priors (one marginal density per coordinate, no list density), and the overlays
`lines_prior_list()` and `geom_prior_list()` refuse them alike instead of drawing
nothing; the factor-coded multivariate priors keep plotting.
Factor spike-and-slab and mixture priors built by the constructors plot like
those attached through `JAGS_formula()`: a point factor component without
levels takes the dimension of its sibling factor components for the
point-declaration check and declares no point when none is known.
Prior mixtures bind the common original source knots after component drops,
with exact x/y/Jacobian correspondence and unchanged component weights and atoms.
Approximate named tanh display grids retain tied finite rounded image knots with
their exact source indices and ordinates, warning with the affected source values.
They require two distinct images and consistent nonzero orientation; whole
collapse, custom/affine collapse and invalid interior knots remain unavailable.
This display convention does not relax atom, draw or inverse-region checks.
The KDE curves of `plot_marginal()` and of simple and factor `plot_posterior()`
use the same checked transformation grid as precomputed densities: knots
without a representable image or ordinate (out-of-support tails, saturation) are
dropped and the plot limits come from the retained curve.
A posterior KDE that has no exact support takes its reflection bounds from the
prior of the plotted coordinate only on the fitted scale; a random-effect SD of
`as_mixed_posteriors(transform_scaled = TRUE)` takes the domain of an SD, `[0, Inf)`,
not the bounds of its fitted-scale prior. The prior context rebuilt for a prior with
`multiply_by` keeps the fitted block designs.
The automatic y range of PET-PEESE posterior plots with priors collapses by its
width (below .01): it is then widened around its values (a collapsed range of
small non-negative values keeps the unit interval), and any wider range, also
an all-negative one, is kept.
Posterior PET/PEESE bands use pointwise HD quantiles; sampled-prior bands retain
empirical quantiles and deterministic prior bands retain their CDF inversion.
PET/PEESE continuous-pair CDFs integrate over location probability, retain raw
panel diagnostics, and share a total200 actual subdivisions per original pair.
Each finite panel uses checked affine coordinates with its original physical
endpoints and width Jacobian, preserving the physical partition and criteria.
Default-range continuous primitive quantiles consume the smaller directly
represented lower or upper tail; truncated components retain their existing
lower-probability route and refuse lost interior coordinates. Positive subnormal
weighted CDF-panel terms remain governed by the CDF error criteria.
Optional partition anchors retain numerical omission provenance; required event
arithmetic, quantile ranges and solver limits refuse without sampled rescue.
Deterministic prior bands evaluate that CDF for all x-values and probabilities
together (`.petpeese_prior_quantile_band()`): each continuous pair integrates
over the location-probability window where the bias CDF varies, on panels
anchored per row at the images of fixed bias quantiles, with the same 1e-7 CDF
criterion; cells not resolved within it, and contexts without a continuous
pair, use the scalar quantile above with its diagnostics and refusals.
The bracket endpoints of that solve (`.petpeese_prior_model_range()`) are seeds,
not reported values: they keep the plain sum of the location and scaled bias
endpoints, also when the bias is below half an ulp of the location, while the
band values keep the strict representability check (`.petpeese_prior_sum_checked()`).
Keep requested display limits separate from transformation support. Custom
transformations can declare `output_support = c(lower, upper)`; inverse-grid
evaluation must exclude coordinates outside that interval, without clamping
actual values or moving point masses. An infinite inverse at a boundary does
not supply a finite continuous-density ordinate.

Autocorrelation keeps the original retained iteration positions and every
requested lag. Unavailable ordinates are NA, with same-length `pair_count` and
`unavailability_reason` attributes on the existing x/y object. Reasons prefer
physical `chain_too_short`, then `no_finite_pairs`, then
`nonfinite_autocorrelation`. Emit one `BayesTools_autocorrelation_unavailable`
warning (parent `BayesTools_plot_condition`) per affected parameter/chain with
the unavailable-lag fields. Finite ordinates retain their original estimator;
zero-anchored finite ranges and finite-only bars are display choices, never
compression, NA-to-zero estimates, confidence bands or a new pair-count cutoff.


## Verification

Test computed plot data, validation, and dispatch separately, then retain
representative human-reviewed `vdiffr` snapshots for both backends where both
are public. Structural tests and render-only checks do not replace visual
regression.

Use the `visual` profile for pure plots and `visual-fixture` for plots that load
cached JAGS fits. Run `fit` first when the required fitted cache is affected.
Make stochastic plot inputs deterministic. Retain changed candidates for
maintainer or explicitly delegated review before accepting a baseline change.
