# Priors and Density Semantics

Use this guide when changing prior constructors, distribution methods, density
provenance, transformations, mixtures, or prior-density and Savage-Dickey
ordinates.

## Prior Contract

`prior()` creates the core `prior` S3 family. Distribution aliases and routing
live in `R/priors.R`; constructors and parameter normalization live in
`R/priors-constructors.R`; joint and marginal methods live in
`R/priors-methods-joint.R` and `R/priors-methods-marginal.R`. Specialized
distribution implementations live in `R/distributions-*.R`.

When adding or changing a distribution, keep these paths synchronized:

1. accepted names and aliases in `prior()`;
2. constructor validation, canonical parameters, support, and truncation;
3. applicable `rng`, `quant`, `cdf`/`ccdf`, `pdf`/`lpdf`, and marginal methods;
4. JAGS prior syntax, initialization, monitoring, and bridge-sampling density
   where the distribution is supported for fitting;
5. printing and plotting where applicable;
6. focused unit and visual tests using the existing prior helpers.

Do not claim a joint or marginal method for a prior family unless its
measure and return shape are defined. Point masses, continuous densities,
mixtures, vectors, and factor priors are not interchangeable merely because a
numeric representation can be produced.

Truncation is part of the distribution definition. Preserve exact support and
one-sided boundary behavior; do not replace valid endpoints with nearby values.

`selection_model()` captures fixed source choices and a deferred publication
column reference. Keep `group = NULL` immediately after `weight_rule`; consumers
resolve publication groups only for `"best"` branches. Product selection is
invariant to publication partitions and must not require or bind a group
column. `group` is a deferred column reference rather than a modelling
instruction, and callers parameterize `weight_rule` while forwarding one
`group` for both rules, so a group supplied alongside the product rule is inert
rather than contradictory: keep the combination constructible and let `print()`
report the group as unused. Integration dependencies come from the integrated
sampling and random-effect covariances, independently of publication groups.
Preserve typed reference capture, branch metadata, and source-conditioning
choices.

Generated JAGS syntax and generated initial values are two renderings of one
node layout, so derive both from `.weightfunction_component_node_names()`
rather than rebuilding names on either side. Initial values belong on the
stochastic node; `omega_target` (`omega`, `omega_component_<id>`) is
deterministic and expanded onto the global cut grid, and initializing it either
does nothing or aborts JAGS with a dimension mismatch. Any check that mocks the
model, the sampler, or the backend cannot establish that the two sides still
agree - keep at least one test that compiles the real model.

## Structural Prior-Density Ordinates

`prior_density_ordinate()` classifies the mathematical behavior of a prior at
one exact finite value. Its stable schema and behavior values are public. Keep
`R/prior-density-ordinate.R`, its roxygen documentation, and
`tests/testthat/test-prior-density-ordinate.R` synchronized.

Classification must come from the prior definition and deterministic
provenance. Never infer `regular`, `zero`, `infinite`, `point_mass`, or
`undefined` from:

- prior or posterior samples;
- KDE output;
- a finite grid, interpolation, or FFT clipping;
- probing `value +/- epsilon`;
- binary64 underflow or overflow alone.

A structurally regular ordinate remains regular if its representable density
underflows. A regular ordinate whose value is unavailable (a quadrature
rejected by its diagnostics, or a boundary limit without a structural value)
has `exact = FALSE`, `log_density = NA` and the failure as `reason`:
`exact = TRUE` always comes with an available regular log density, and
consumers (RoBMA's point-test eligibility) read `exact`. Exact point mass
takes precedence at the requested value; the
continuous behavior remains diagnostic provenance: every `point_mass` result
(and only such a result) carries `provenance$continuous_behavior`, the
behavior (`regular`, `zero`, `infinite`, `undefined` or `unknown`) of the
measure without its point masses at that value, with `log_density` the log
density of that continuous part (not renormalized). It is documented in the
roxygen `@return`, and consumers (RoBMA's point-test eligibility) read it, so
keep it on every point-mass path, including stored atoms of linear densities.
Unsupported transformations or convolutions return `unknown` rather than a
guessed structural class.

Full-precision rule: no ordinate is exact when its value is computed from a
subnormal (nonzero, absolute value below `.Machine$double.xmin`), underflowed
or non-finite intermediate (`.prior_density_full_precision()`,
`.prior_density_affine_full_precision()` for a distance from an exact anchor
and its standardizations). Such an ordinate has `exact = FALSE`, no log
density, and the reason "... not representable at full precision"; its
behavior is `regular` where the route has established the support and
`unknown` otherwise (a linear density then shows its grid estimate). One
place per route kind applies it:
- primitive densities (`.prior_density_ordinate_primitive()`, which every
  scalar chain reaches): a nonzero subnormal value, for every family (the
  distribution functions rescale their argument; a gamma(3, 0.7) log density
  was 2.8e-4 off at 1e-320); 0 keeps its structural class; the rescaled
  argument of a normal value whose log the family takes (`value * rate` of
  a gamma with shape other than 1, `value * sdlog` of a lognormal; a
  gamma(2, 1e-20) log density was 1.1e-5 off at 1e-300); a
  t density below `.Machine$double.xmin`, which `extraDistr::dlst()`
  evaluates before its log (a Cauchy log density was 4.1e-4 off at 1e160),
  while a t density that underflows to 0 keeps `-Inf`;
- the scalar affine map and `lin` output transformations: the distance from
  the offset and the inverse value (an underflowed 0 was classified as the
  source bound, a subnormal one rounded); `tanh`: its inverse value;
  `exp_lin`: its inverse value `exp(...)`; the log source: `e^z`, also when
  it overflows (the density is then not evaluated, rather than asserted
  `-Inf`);
- Jacobians (`.prior_density_ordinate_wrap()`): a regular source whose log
  density is `-Inf` (underflowed, e.g. a t density evaluated before its log)
  under a Jacobian factor above 1, which can make the transformed density
  representable;
- quadrature leaves: the scale product (distance from the offset and its
  standardization), the pure scale mixture (the same, standardized by the
  multiplied SD and mean) and the two-term convolution (the distance and its
  standardizations by both weights); their plotted values there are NA; an
  accepted integral value below `.Machine$double.xmin`
  (`.prior_density_quadrature_ordinate()`: the log image of a half-Cauchy
  product was 1.2e-2 off at z = 370), whose plotted value keeps the
  estimate; the log image of a scale product: `e^z`;
- endpoint matching: the neighbouring double of an endpoint counts as the
  endpoint only when both are at full precision (the smallest subnormal is
  not the bound 0).
Safe without the rule: the `exp` output transformation (`log(y)` of a
positive double is finite and 0 or at least 1.1e-16 in absolute value),
closed-form Gaussian sums, the truncated-normal convolution and Gaussian
convolutions with an additive SD (the value enters only Gaussian kernels of
standardized distances, whose relative change under an absolute rounding of
at most the smallest subnormal is negligible), offsets and support bounds
(exact structural points), atoms, mixtures and row mixtures (weighted sums of
the above), and grid ordinates (never exact). Plotted closed-form densities
at subnormal values are not changed (display only).

The density-provenance implementation is shared across
`R/priors-density-context.R`, `R/priors-linear-density.R`,
`R/priors-linear-density-combinations.R`, `R/priors-density.R`, and
`R/prior-density-ordinate.R`. Extend this path instead of building a parallel
density algebra. Provenance must remain deterministic and compact: do not store
draws, large grids, fitted objects, environments, or closures capturing them.

Conditional-normal quadratures are independent 1-D integrals: each distinct
design row, model or conditional mixture component, and mixture leaf receives
the full evaluation budget (`n_grid`) and its own convergence check. Each such
integral is split at the multiplier's (other term's) support bounds, at its
declared-prior quantiles (1e-6, 1e-3, .02, .25, .5, .75, .98, and their
complements), and at the location peak of the conditional normal and +-1, 3,
10 local SDs around it when the multiplied SD is at most half of the
absolute value of its mean (always for Gaussian convolutions; beyond that
guard the window is not a peak and is not used). Next to a finite bound with
infinite density, a breakpoint much closer to the bound than the next one
makes QUADPACK count the mass near
the bound twice: there the extreme quantile is skipped, and up to the quartile
on that side a breakpoint is kept only if its distance to the bound is at
least 1e-3 of the next kept breakpoint's distance. Gaussian-peak breakpoints
at least one local SD from the bound are exempt, so a narrow peak next to the
bound keeps its breakpoints (closer peak points, i.e. rounding residues, are
still dropped). Breakpoints closer than
1e-9 of their magnitude (at least 1) are merged, so no piece is only a few
ulps wide, but the merge width never exceeds half the peak's local SD, so the
peak breakpoints are never merged. Every piece gets the full budget, and the
acceptance criterion applies to the summed value and error: every piece
converged, a positive value, and a reported absolute error of at most 1e-4 of
the value. The criterion is purely relative (Savage-Dickey uses log
ordinates), so a far-tail value is never accepted on an absolute floor; the
first QUADPACK pass stops at an absolute floor of 1e-12 in total, and a total
that misses the relative criterion is refined once against its own value
(every quadrature piece again with its full budget, and an exactly evaluated
Gaussian-convolution piece of value 0 whose bound 2 * Phi(-10) of its mass is
too large for the total by quadrature). A total that still misses it is
rejected (`exact = FALSE`, refused at point hypotheses); plotted curves may
draw its estimate, which never makes an ordinate exact.
Budgets are never divided among them; the budget only caps how many mixture
leaves are expanded. Known limitations: pure scale mixtures with heavy-tailed
multipliers can stop as non-convergent (mostly with a small multiplied SD); a
piece that QUADPACK flags stops the ordinate even when its value is negligible
against the total (e.g. a far tail piece reported as "probably
divergent"); a quadrature that evaluates to exactly zero is rejected (it
cannot be told apart from a missed peak), which also stops a mixture with
such a component, e.g. a narrow component far from the value; very narrow
Gaussian peaks can stop as non-convergent where QUADPACK reaches the
floating-point resolution: next to a bound with infinite density away from
zero (e.g. an SD of 1e-10 at the upper bound of a beta), or when only a few
hundred ulps wide; and the scale peak of a value far in the
multiplier's heavy tail (beyond its extreme quantiles) is not a breakpoint,
also with a nonzero multiplied mean when the multiplied SD is not small
against it, so its mass can be missed without a convergence failure. Such
ordinates are small at ordinary scales, but the missed fraction does not
depend on the units of the value, and the reported error does not bound it.

Product terms (`multiply_by`) and ordered-prior levels share one route
(R/prior-density-route.R). An ordered level (or allocation subset) is the
ordered total times its allocation share (`.prior_ordered_linear_share()`):
a fixed share scales the total, a Beta share multiplies it. A point
multiplier or a deterministic multiplied part folds into an affine term. A
normal multiplied part with a normal or deterministic additive part is the
conditional-normal route; without an additive normal term it is a pure scale
mixture (additive SD 0), split additionally at the multiplier's zero and at
|value - a_m| / b_s (1/10, 1, 10) on both sides of it. At the deterministic
offset a_m its density is phi(b_m / b_s) E[1 / |s|] / b_s: the multiplier's
declared behavior at zero decides it (a vanishing density gives the finite
value, closed-form inverse moments for the untruncated gamma, inverse-gamma,
lognormal and beta families and a quadrature otherwise; a positive or
infinite density an infinite ordinate). A single non-normal multiplied term
(or a non-normal ordered total) with a deterministic additive part is a
`scale_mixture` quadrature over the multiplier, split also at the images of
the term's quantiles and bounds; its offset is classified from both declared
behaviors at zero (infinite when either is infinite or both are positive;
f_L(0) E[1 / |s|] or f_s(0) E[1 / |L|] when only one vanishes; zero when
both vanish; a jump, where a term bounded at zero with a positive finite
limit meets a two-sided other term whose density vanishes at zero, is not
classified, while a two-sided other term with a positive density at zero
makes both one-sided limits infinite, e.g. the level of a t total with a
Beta(1, b) share). The same leaf with the multiplier map m(s) = sqrt(k s) of
a Beta(a, b) share is the continuous part of allocated random-effect SDs
(below): the integral over the share is split at its quantiles and at the
images s = (d / q)^2 / k of the term's quantiles q, and the mapped share's
density at zero behaves like v^(2a - 1) (zero for a > 1/2, 2 / (sqrt(k)
B(1/2, b)) at a = 1/2, infinite below), with E[1 / sqrt(k S)] =
B(a - 1/2, b) / (sqrt(k) B(a, b)). The route provenance of a scale product
records its support hull and its offset behavior, so named transformations
(the square of an SD) classify it. A value whose distance from the offset,
x - c, or standardized distance (x - c) / w is not representable at full
precision has no scale-product ordinate value and its plotted value is NA
(the full-precision rule above: the integrand's factor argument would be
rounded, a gamma(2, 4) factor was off by 1.4e-4 at 1e-320). The rule is
applied at the leaf, so it holds on every route reaching it (products,
ordered levels, allocated SDs, output transformations and log images). Other products (several products, a
non-normal additive term with a product, several multiplied non-normal
terms) have no structural route; their capped product grid is never used for
heights or probabilities, which are then unavailable.

Multivariate t terms: a vector prior `mt` (and `mcauchy`, its one degree of
freedom case) is X = mu 1 + z / sqrt(w), z ~ N(0, s^2 I) and
w ~ Gamma(nu / 2, nu / 2) (the scale matrix s^2 I, not a covariance, as the
JAGS emitter `.JAGS_prior.vector()` draws it and `mvtnorm` evaluates it with
`sigma = s^2 I`). A linear combination a'X of its coordinates is the
univariate t with location mu sum(a), scale s ||a|| and nu degrees of freedom
(zero weights are dropped, so a' S a > 0; a zero combination is the point 0).
Before routing, each such group is rewritten as that scalar t term under the
prior's name, keeping its `multiply_by` scale
(`.prior_density_route_vector_t_terms()`; the ordinate provenance records
the rewritten groups as `multivariate_t`), so levels and contrasts of
mean-difference and orthonormal `mt` factor priors and any linear target of
them are scalar t terms: exact alone (`scalar_affine`, with the t
distribution function for regions and closed-form plotted densities), and
otherwise routed as below (a Gaussian convolution with normal terms, the
Cauchy sum for one degree of freedom, a two-term convolution with one other
non-normal term, a scale product with a `multiply_by` scale). Mixture and
spike-and-slab priors of `mt` components are expanded first (below). Groups
with non-numeric parameters, a truncation (which vector priors do not
support), a source transformation, or a scale that is not representable keep
the general route.

Sums without a product term: untruncated Cauchy terms (Cauchy or t with one
degree of freedom, no log source) first merge into one Cauchy term with the
summed locations and the summed absolute scales, so Cauchy sums are scalar and
exact. A Gaussian part plus one truncated normal term (a normal with a finite
bound, no log source) is closed form (`truncated_normal_convolution`,
`R/prior-density-truncated-normal.R`): with w T' the weighted truncated term,
N(m_1, s_1) before truncation to [l', u'], and the Gaussian part N(a_m, a_s),
f(x) = phi(x; a_m + m_1, s) [Phi(beta(x)) - Phi(alpha(x))] / [Phi(B) - Phi(A)],
s^2 = s_1^2 + a_s^2, alpha(x), beta(x) the bounds standardized by T' given
the sum (mean m_1 + s_1^2 / s^2 (x - a_m - m_1), SD s_1 a_s / s), A and B by
T' itself; both interval masses are evaluated in log space from lower- or
upper-tail probabilities. Its region probabilities keep the Gaussian-
convolution quadrature (they need the bivariate normal distribution
function). The exp of a log-source term X (weight 1, positive support) plus
a Gaussian part G ~ N(m, s), s > 0 (normal terms and points), e.g. the
unscaled intercept of a log-intercept formula scaling, is the scale product
X W with the lognormal multiplier W = exp(G) ~ lognormal(m, s)
(`.prior_density_route_exp_scale_product()`), with the scale-product
ordinates, region probabilities, plotted values and offset classification at
0 (f_X(0) E[1 / W]); mixture priors of X are expanded first. The log-scale
sum log(X) + G itself (e.g. formula marginal posteriors of log-scale linear
predictors) is the log image of that product (`log_scale_product`):
f_Z(z) = f_Y(e^z) e^z with the product's classification at e^z > 0 (no
offset case) and its quadrature error scaled by e^z (no ordinate value where
e^z is not a normal double, z below about -708.4 or above 709.8), region probabilities
P(Z in (a, b)) = P(Y in (e^a, e^b)), and plotted values through the same
ordinate; its exp is the product itself. Other weights of the log-source
term (X^w W) and a non-Gaussian other term keep the general convolution.
Two simple continuous terms (after merging) that are not both normal
and have no Gaussian part are a two-term `convolution`: the 1-D integral over
the first term (the one with an infinite density at a finite bound when only
one has such a bound) of its density times the other term's density, split at
the first term's bounds and quantiles and at the images of the other term's
quantiles and finite bounds, with the conditional-normal budget and
acceptance criterion; its region probability integrates the first term's
density times the other term's exact region probability. Where finite bounds
of both terms meet at the value, the ordinate is classified from their
exponents p (1 for a positive finite density, above 1 for a vanishing one,
the gamma or beta shape for an infinite one): with e = p_A + p_B - 1 an end of
the support is zero for e > 0 and infinite for e < 0, and an inner meeting
point is infinite for e <= 0; the positive finite limit at an end with e = 0
(e.g. two arcsine terms at 0 and 2) is not classified. Combinations that keep
the numerical grid (`unknown` ordinates): three or more non-normal terms
after merging, a Gaussian part with two or more non-normal terms, log-source
terms other than lognormal ones in a sum (except one log-source term with
a Gaussian part, above), the products listed above, custom output
transformations and `bounded_logit`, the unclassified meeting points, and SD
components and totals (and their variances) of nested variance allocations
with a term of the scale prior times two or more independent allocation
shares.

Plotted linear-combination prior densities (`.prior_linear_density_to_plot_data()`)
evaluate the same route at every plotted value: closed forms vectorized over
the plotting grid, and quadrature leaves by one batched quadrature over all
plotted values (`R/prior-density-quadrature.R`): each value keeps its
ordinate's integrand and breakpoints, all pending intervals of all values are
evaluated in one vectorized call with QUADPACK's qk21 / qk15i rules and error
estimates, and each value is accepted at a relative error of 1e-8. Special
values (offsets, support bounds, meeting points), values that do not converge
within the round cap, and all values of a leaf whose integrand has an
integrable singularity (a term with an infinite density at a finite bound;
bisection without extrapolation converges too slowly there) take their
ordinate. A route with quadrature leaves is plotted on at most 200 equally
spaced values plus the values it must include (atoms, offsets, normal means,
meeting points, finite support bounds with a point 1e-6 of the plotted range
outside each bound: the value at a bound is the density's limit inside the
support, so the point as far inside would repeat it and is added only where
the density at the bound is infinite or unavailable) and the parabola vertex
at each interior local maximum.
Every plotted curve (closed forms, quadrature leaves, and the interpolated
numerical grid, whose bounds are those of the exact support hull of its
provenance) leaves a finite support bound of the continuous density by the
vertical edge of `density.prior()`: the bound's value is repeated with density
0 on its outer side (first point `(lower, 0)` then `(lower, f)`, last points
`(upper, f)` then `(upper, 0)`) when the bound lies in the plotted range (up to
a relative 1e-9 of it), the density there is positive (also infinite) and zero
just outside it. A bound outside the range, one where the density is already 0
and one inside another component's support (an interior jump, drawn by the
display values) get no zero value. The repeat is made before the transformation
step, so a transformed curve maps it with the rules of `density.prior()` (a
bound mapped to -Inf is dropped with it).
A batched value can exist where the per-value ordinate's quadrature is rejected
by QUADPACK's flags (heavy-tailed pure scale mixtures): the plot shows it,
while heights and point hypotheses at that value stop. Only a combination
without a structural route, or a density without recorded provenance,
interpolates its numerical grid.

A row mixture of a density context (`.prior_density_from_context_rows()`,
e.g. the prior density of a formula `marginal_posterior()` level) whose
structure has no point mass (only continuous priors, every row with a
nonzero weight, no constant output map) defers its numerical grid, which
holds one exact display grid per row product component: its `density` field
keeps the builder's arguments and a cache, and the grid is built on the first
read of a grid field (plots without a given range, grid heights and
refinement, transformations of the density). Ordinates, exact heights and
region probabilities evaluate the recorded route and never build it; the
atoms (none), continuous mass (1) and grid size are exact without it.

Heights and region probabilities require recorded provenance (the
`adaptive_evaluation` attribute naming the prior measure): a density grid
without it has no error control and cannot be refined, so it is used for
plots only and heights and probabilities stop; without a continuous part a
density is its exact point masses. `prior_density_has_provenance()` is the
public signal of this (also `FALSE` for a recorded product without a
structural route, whose heights stop as well); the ordinate `method`
`"unsupported_provenance"` is not, since it also marks combinations without a
structural route, which have provenance and refined grid probabilities. Producers therefore build densities
through `.prior_linear_combination_density()` or the context builders:
`parameter_prior_density()` passes a semantic transform as the output
transformation (affine as `lin`, `tanh`, square root and square of a
nonnegative source as `exp_lin`, a bounded logit as a recorded custom map on
the refined grid) and a gated variance proportion as one mixture prior of its
atoms and Beta components (with a model-averaged, mixture scale prior as
well). An allocated random-effect SD is its scale prior T times an
independent multiplier M of the fitted model (shares w ~ Dirichlet(a) over
all K components of the allocation, not renormalized over active ones), and
its density is the `allocation_product` measure
(`R/prior-density-allocation.R`): M is a finite mixture of points (0 when a
gate of the chain is off or no component of a total is active, 1 for a total
with every component active or a gate-only chain) and mapped shares
sqrt(k W) (a component SD: W = w_i ~ Beta(a_i, a_- - a_i), k = K for
mean-variance and 1 for total-variance allocations; a total with the partial
active set A: W ~ Beta(a_A, a_- - a_A)), with the gate-configuration
probabilities. The total of a nested allocation is T times the multiplier
chain of its parents (the gates as points, each parent Dirichlet share as a
Beta share, as for the component SDs) times the own active-set factor of the
total of a non-nested allocation (`conditional` conditions on a positive
total: the parent gates on and the own empty set dropped), so RoBMA's
`study / esid` total, T times a gate, is exact. T's mixture and
spike-and-slab components are expanded; each leaf is an atom, the scaled
scale prior, or the square-root share scale product above, and a variance
applies the square (`exp_lin`) to every continuous leaf and to the atoms, so
ordinates, region probabilities and plotted densities are exact and share one
route. SD components and totals of nested allocations with a term of two or
more independent shares keep a product of density grids without provenance:
a plotting density whose ordinates are `unknown` with the reason recorded
(`provenance_unavailable`), refused for heights and point hypotheses. Fitted
coordinates and factor levels use the prior-density context of the priors owning their
coordinates (`multiply_by` stripped). A pairwise correlation of a K-column
LKJ(eta) block is the LKJ marginal, `2 B - 1` with
B ~ Beta(eta - 1 + K / 2, eta - 1 + K / 2), passed as that Beta prior with
the affine (`lin`) output transformation: on the fitted scale (gated blocks
included, since the LKJ primitives are a priori independent of the SDs and
their gates) and for an original-scale pair whose coefficients are each one
rescaled fitted coefficient (the unscale-matrix rows of the pair have one
nonzero entry each). Original-scale correlations that mix fitted
coefficients (the intercept of centred predictors) combine the correlation
with the SDs and have no density. `parameter_mixed_posterior()` (and
`random_effects_summary_posterior()`, built on it) attaches this density,
the catalog support, and atoms declared from structure, never from draw
values: the atoms on which the per-draw states of the quantity's gates and
point components (`parameter_gate_states()`: allocation gates, the indicator
of a mixture scale prior with a point component, the component indicator
of a mixture or spike-and-slab prior of its coordinates, and the branch of
a selection prior: the constant weights of reference bins, fixed weights and
branches without a selection, 0 for `PET`/`PEESE` and the p-hacking
parameters in branches without them, and the `phack_kind` code) put it, with
their shares of the draws as masses (checked against the density's point
masses when there is a density; `parameter_gate_states()` also gives the
prior masses of point components from the component weights); none when the
density has no point mass or, without a density, when no coordinate of the
quantity can take a point mass (e.g. an original-scale LKJ correlation with
continuous SD priors). A coordinate cannot take one when its prior has no
point component (continuous priors, Dirichlet weights), when it is an LKJ
primitive or a standardized random effect, or when it is a generated
random-effect SD or correlation node (registry families `random_sd`,
`random_rho`, `lkj`) whose registry dependencies (`.bt_deterministic_nodes()`)
all cannot: the SDs of an allocation without inclusion gates whose scale
prior has no point component, so the original-scale correlations of an
allocated `us()` block declare none; a gate, component indicator, or point
component among the dependencies propagates. Correlations are
scale-invariant (cor(cS) = cor(S) for c > 0 multiplying every SD of the
block), so for a correlation the coordinates that enter its block only as a
common factor of every SD count as point-free
(`.bt_parameter_correlation_common_factors()`, from the registry chains of
the block's allocated SDs: the shared scale source and the gates of factors
shared by every chain, e.g. a gate-only allocation split by an
`sd_component` child); where the factor is 0 the draw is undefined and
dropped. Gates and point components acting on only some SDs of the block
still propagate (original-scale correlations can then be -1 or 1 where one
SD is 0). Atoms stay undeclared for unmonitored point states, composites of
coordinates with point masses (e.g. original-scale correlations of a block
whose SDs have their own spike-and-slab priors), and
coordinates whose point structure is not classified (coordinates without a
prior other than LKJ primitives, standardized random effects, and generated
random-effect SD and correlation nodes, e.g. a formula's linear predictor,
which is constant on rows whose design is zero; and selection coordinates,
e.g. the `_indicator`, `_inclusion`, and `_variable` coordinates of mixture
priors).
`parameter_mixed_posterior(conditional = TRUE)` keeps the draws of the
quantity's inclusion event and uses the density restricted to that event (the
gate atom drops out, atoms inside the event are renormalized).

A product component of a combination (a `multiply_by` product, or an ordered
level as its total times its Beta allocation share) with a structural route
is constructed from that route (`.prior_linear_density_route_product()`): its
grid is the route's continuous density on at most 1024 values over the
product range, with the product's exact atoms. The grid values share the
batched quadrature of plotted densities; unlike plotted curves, the values of
a leaf whose integrand has an integrable singularity are batched as well,
accepted at the ordinates' criterion (relative error estimate at most 1e-4),
and only values it does not accept take their ordinate. A leaf with a strong
singularity is not batched: a multiplier, share or factor whose density is
infinite at a finite bound with exponent p below 0.1 (density ~
distance^(p - 1); the shape of a Beta share at that bound or of a gamma term
at 0; an unknown exponent counts as strong), at any bound, including a zero
bound whose singularity the other term cancels. Next to it the bisection's
error estimate is not reliable (accepted grid values of a Dirichlet(0.05)
level share were 4.3e-4 off at an estimate of 1e-4), so all its values take
their ordinate. The capped product grid of the
factors' grids (`.prior_linear_density_product()`), which a heavy-tailed
factor can leave without a finite positive mass, remains only for products
without a structural route (plots only). The route-evaluated grid records
whether it resolves the product (`product_grid_resolution`: the Riemann sum
of the route density on the grid within 10% of the continuous mass). A
heavy-tailed factor (a Cauchy total or multiplier, whose 1e-4 tail range spans
thousands of scales on 1024 values) leaves it unresolved (about 80% low),
while normal, t3, t2, gamma, lognormal and inverse-gamma factors stay within
about 1%. A plotted sum without a structural route that relies on an
unresolved product grid (a prior curve 129-456% off) is omitted with a
`BayesTools_prior_curve_unavailable` warning (also
`BayesTools_plot_condition`); the rest of the plot is drawn, and heights stay
refused as for every product grid. Posterior plots (`plot_posterior()`,
`plot_marginal()`) of draws that declare no prior (`prior_none()`, e.g.
catalog mixed posteriors of quantities without a prior density) and carry no
prior density warn with the same class and draw the posterior alone; draws
with a prior list but without their prior density still stop. The
mixed-measure `density()` of an ordered
prior (a total with a spike) evaluates each level with a structural route on
that route at the display values, never by interpolating its level grid.

Mixture ordinates (model, conditional, and mixture or spike-and-slab terms of
one linear combination) are weighted sums of per-component ordinates, each from
its own exact or regular method, so a numerical grid never spans a density jump
between components. A coordinate of a mixture whose components are vector
(factor) priors, e.g. a model-averaged mean-difference level of a single-fit
product-space prior, is no scalar term: the mixture is expanded, and each
component is routed as its vector prior (normal, the univariate t of a
multivariate t, or an atom for a point component). A Gaussian term plus one
other continuous scalar term
(and two non-normal terms) uses its quadrature (the closed form for a
truncated normal term) wherever it occurs, so the same
combination is never exact in one context and a grid approximation in
another. A mixture
height sums its components' exact or regular heights with those of components
that have none; each of the latter uses its own grid (never spanning another
component's jump), all such grids are refined in lockstep, and the documented
refinement criterion (a change of at most 1e-4 of the height, purely
relative) applies to the mixture height with the components'
weighted absolute changes (no cancellation between components). Single grids,
mixture component grids and grid region probabilities share this one
refinement loop (`.prior_linear_density_refine_grids()`). A refinement-
change criterion is not an error bound: grid-based heights of components with
singular source densities (e.g. gamma shape < 1) remain approximate within it.
Structural behavior is never flagged from grids or attached evaluators (the
former product flags and scalar density evaluators are gone); an `exp_lin`
image of a zero source has the structural boundary limit C / (b exp(a)) of a
source density ~ C x^(b - 1) (beta, gamma, and exponent-one families).

Scalar prior-region probabilities (`R/prior-density-region.R`) follow the
ordinate's classification wherever it has a structural route, for regions
whose relations are linear in the quantity (unions of intervals): point masses
add their exact mass (the exact condition decides strict and inclusive
bounds), scalar priors use their exact distribution function (through a log
source and named monotone output transformations), normal sums their normal
distribution function, two-term convolutions the 1-D integral described
above, and Gaussian convolutions and conditional-normal scale
mixtures the 1-D integral of the other term's density times the Gaussian
probability of the region, with the ordinate's breakpoints (the Gaussian-peak
window, under the same guard and so always for Gaussian convolutions, at every
finite region bound), the full budget per piece, and the same acceptance
criterion on the total; a rejected quadrature stops. Mixture,
spike-and-slab, model and conditional components and distinct design rows are
expanded as for the ordinate and summed with their probabilities; a
combination is evaluated this way only when every component has such a route.
Keep the two paths in step: a combination whose ordinate is structurally
classified has a structural region probability, and one whose ordinate is
`unknown` keeps the grid. A quadrature total may exceed [0, 1] by at most its
absolute error. For a Gaussian convolution G + w T, a piece that lies
entirely outside every window t* +- 10 s / |w| of the finite region bounds is
evaluated exactly as 0 or 1 times its mass under T's declared distribution
function (the Gaussian region probability is constant there to within
2 * Phi(-10) ~= 1.5e-23), so heavy tails and region bounds far beyond T's
quantiles need no quadrature. Known limitation: conditional-normal scale
mixtures integrate every piece, and with a heavy-tailed multiplier (e.g.
half-Cauchy or inverse-gamma(1)) QUADPACK flags the infinite end piece of a
region with an infinite bound ("roundoff error", "probably divergent"), so
such probabilities stop as they did on the grid. A region whose probability
underflows to zero stops like a zero ordinate, and so does a Gaussian
convolution whose region probability comes only from pieces evaluated as 0
(a representable probability below 2 * Phi(-10) of their mass under T, e.g. a
region far below a bounded or nonnegative T).

For the other prior-region probabilities, integrate the continuous grid's
piecewise-linear interpolant up to the exact region boundaries (the comparison
value, or boundaries located by bisection within the cells where a composite
condition changes), normalize by the interpolant's total over the same grid,
and evaluate point masses with the requested strict or inclusive operator. Do
not move a region boundary to nearby grid knots. Keep the existing refinement
gate.

## Savage-Dickey Ordinates

`hypothesis_BF()` applies one exactness rule to every point hypothesis: the
prior ordinate at the null must be `regular` and exactly classified by
`prior_density_ordinate()`; otherwise it stops with a classed error
(`BayesTools_point_mass_at_null`, `BayesTools_infinite_ordinate`,
`BayesTools_zero_ordinate`, `BayesTools_undefined_ordinate`,
`BayesTools_inexact_ordinate`, all also `BayesTools_hypothesis_ordinate`),
which callers match by class, never by message. `prior_ordinate_status()` is
this rule as data (per value: eligibility, class, message, and the continuous
behavior); the stopping check is its stop-at-first-ineligible wrapper, so the
two cannot diverge. Linear expressions of
marginal posteriors use the joint-context density, affine expressions of a
prior object its transformed density; nonlinear expressions of a
deterministic prior are inexact. Only user-supplied prior draws keep a
kernel (or normal) estimate of the prior ordinate, signalled by a warning of
the same inexact class (draw-only inputs have no structural prior density).
`Savage_Dickey_BF()` applies the same rule and classes to its prior ordinate
through the same check: a zero or infinite prior ordinate makes the density
ratio a 0/0 or singular limit that a posterior kernel estimate cannot
estimate, so a scalar call stops, and a level of a list posterior (marginal
inference and its tables) gets an NA Bayes factor with the reason, as a level
fixed at the null.

The posterior ordinate of a Savage-Dickey ratio (`R/marginal-savage-dickey.R`)
is a kernel estimate from the continuous draws: the exact Gaussian kernel sum
at the null with bandwidth `bw.nrd0()` of the (per-component) continuous
draws, reflected at finite exact support bounds, never an evaluation grid,
interpolation or binning (a grid spanning a long-tailed draw range is coarser
than the bandwidth). The raw-draw ordinates of `hypothesis_BF()` use the same
kernel sum without reflection. A kernel sum below the double range (draws
many bandwidths from the null) is 0, and the Bayes factor is then +Inf; no
numerical floor is added. A null outside the prior's exact support has a
zero prior ordinate (classified structurally, never from the range of a
grid) and follows the zero-ordinate rule above.
Everything around the ordinate comes
from declared metadata: posterior atoms, exact support, and for mixtures each
draw's component with the components' exact supports (`components`:
the model of a `mix_posteriors()` ensemble, or the indicator tuple of the
mixture terms entering the quantity for `as_mixed_posteriors()`). Never infer
atoms, supports, or components from the draws.

Draw metadata (supports, atoms, components, undefined draws, prior densities
and their context, precomputed posterior densities and ordinates, formula
flags, linear weights, conditioning) live in one validated attribute,
`bayestools_meta` (`R/draws-metadata.R`). Read and write them only through
`.bt_meta_get()`/`.bt_meta_set()` (the public `posterior_metadata()` for
downstream packages), never as free attributes; a unit lint test enforces it.
Producers of the `condition` field set its logical `averaged` element
(unconditional draws), which consumers read instead of comparing condition
keys with a literal.
Posterior atoms come only from the `atoms` field; a precomputed posterior
density describes only the continuous part and carries no point masses. The
component of a draw has one encoding: `component` indexes the declared
component list, with `component_source` the list it indexes (`"model"`: the
models of a `mix_posteriors()` ensemble,
`"mixture"` or `"spike_and_slab"`: the components of the prior, whose slab and
spike positions come from its `components` attribute); an ordered total's
component is `ordered_total_component`. Draws without a mixture carry no
component and form one component. Metadata build failures propagate (never
caught into missing metadata and a silently different estimator); metadata
are absent only by explicit rules (e.g. no linear support for the mixed
draws of a log-intercept term, or for columns outside the prior-density
context). With log-intercept formula scaling the unscaled
intercept is exp(M[intercept, ] L), L the fitted coefficients with the log of
the fitted intercept, so prior-density contexts accept a combination with the
unscaled intercept only through its log (a `log` source transformation, whose
weights M maps onto the fitted coefficients), and the simple marginal
posterior of the unscaled intercept takes its prior density, support and
component supports from the exp of that log (support algebra maps the log
source's support). Formula marginal posteriors of log-scale linear
predictors, whose rows are linear in L, derive their supports, components and
declared atoms through the same log source (atom locations of a point
intercept enter as their log), also for unscaled (`transform_scaled`)
coefficients; their `joint_prior_transformation` marks that their linear
weights act on the log of the intercept. `Ops`/`Math` group generics, `c()`,
`as.numeric()` and subsetting of draws return plain numerics without
metadata; consumers that need the metadata stop on plain draws, and
producers that transform draws transform their metadata explicitly
(`posterior_transform()`, which `marginal_posterior(transformation = )`
applies to the untransformed marginal posterior: strictly monotone maps only;
supports and atoms mapped, prior densities rebuilt with the transformation as
output transformation or transformed on top of their recorded provenance,
stored posterior densities and ordinates changed by the Jacobian, and the
transformation appended to the label parts, whose relation to the term is
their first `transformation` element). The container of draws records a
fingerprint of the values it describes (length, missing count, and two
sums); reading or setting the metadata of draws whose values changed while
their attributes were kept (`x[] <- `, `x[i] <- `, `pmin()`) stops with class
`BayesTools_stale_metadata`, so a producer that replaces values under an
existing container rebuilds it (`.bt_meta_refresh()`) and then transforms
or removes the fields that no longer apply. The fingerprint is computed in one
native pass per read or update, or by an R evaluator of the same definition
when the native routines are not loaded (the package loads without its DLL
when JAGS cannot be located, and metadata never need JAGS); the package
retains no draws between calls.
Subsetting a list of mixed posteriors with `[` keeps the list's container,
with the prior densities of the kept elements.

- When continuous components have different supports, estimate the ordinate
  per component (boundary-reflected on its own support, zero when that
  support excludes the null, the one-sided limit on a bound) and mix by the
  components' shares of the continuous draws. Shared or unavailable supports
  keep the pooled estimate.
- A null outside the continuous draws gives the finite kernel-tail value with
  a warning that it is not reliable evidence, once per parameter or level;
  never `Inf`.
- A declared posterior point mass at the null leaves the ratio undefined:
  level lists and marginal inference return `NA` with the reason in the
  `"warnings"` attribute and compute the other levels; a scalar call stops
  with `BayesTools_posterior_point_mass_at_null` (also
  `BayesTools_hypothesis_ordinate`).

## Numerical Evidence

Use distribution-library references, analytic identities, or independently
derived results. Cover interior values, exact support boundaries, outside-
support values, invalid parameters, transformations, mixtures, and meaningful
extreme tails. Test identities such as CDF/CCDF complementarity and
quantile/CDF inversion only where they are mathematically defined.

For prior changes, run `unit` first and `visual` when plots change. Add `fit`
and `fixture` when generated JAGS syntax, initialization, monitoring, marginal
likelihoods, or cached fitted objects can change.
