# version 0.3.1
### Fixes
- Converged prior-density quadrature with an unavailable full-precision height
  reports its recorded reason without suggesting a larger sample budget.
- Interpretation keeps expanded rows and their evidence/estimate pairs together
  before subsequent plan items, while preserving explicit order overrides.
- Ensemble and marginal estimate tables preserve probability column labels in
  scientific notation so interpretation can read both interval endpoints.
- Interpretation of fallback credible intervals retains the probabilities and
  coverage of the source endpoints.
- Simplex projections with distinct coefficients retain their exact interval
  support even when those coefficients are numerically close.
- Spike-and-slab posterior atoms include declared point-valued slabs and merge
  a zero slab with the excluded atom using the recorded component shares.
- Linear projections of multivariate normal and t priors preserve tiny and large
  nonzero coefficient scales, including representable scales with an overflowing
  or subnormal norm. Empty coefficient vectors retain a point mass at zero.
- `JAGS_bridgesampling()` preserves an explicit `use_neff` bridge-sampler control,
  including `FALSE`, while retaining `TRUE` as the default.

### Breaking changes
These changes affect code and saved objects written for BayesTools 0.3.0.
This version is released together with RoBMA 4.1.0; RoBMA 4.0.0 relies on the
old behaviour.
- installation and saved objects:
  - requires R >= 4.3.0 (was 4.1.0) and the `reformulas` package, and
    contains compiled code: a BayesTools JAGS module (the distributions
    `dbt_invgamma`, `dbt_moment`, `dbt_invmoment`, and `dbt_lkj_cpc`, and the
    functions `bt_lkj_cholesky()` and `bt_lkj_corr()`) that is built against
    the installed JAGS, so an installation from source needs a compiler and
    the JAGS headers and libraries. The module implements the JAGS 4 module
    interface: BayesTools requires JAGS 4.x (>= 4.3.0, < 5.0.0), and
    installation stops with a message naming the reported version when
    `configure` finds another major version through `JAGS_VERSION`,
    `--with-jags-version`, pkg-config, a versioned JAGS prefix, the
    `jags_version()` probe, or the `JAGS_MAJOR` declared by the selected
    headers. On Windows the build selects the newest installed `JAGS-4.*`
    (also under paths with spaces), ignores other major versions, and rejects
    a `JAGS_ROOT` or `JAGS_VERSION` that points to another one.
  - fits created by BayesTools 0.3.0 or earlier must be refitted: they lack
    the parameter map (now version 9) and the fitted-object contract that
    post-fit functions require. `runjags_estimates_table()` /
    `JAGS_estimates_table()` (also for models without formulas),
    `JAGS_check_convergence()`, `JAGS_diagnostics()`, `JAGS_extend()`,
    `JAGS_bridgesampling()`, `transform_scale_samples()`,
    `JAGS_evaluate_formula()`, `as_mixed_posteriors()`, and `mix_posteriors()`
    stop on such fits with an error of class `BayesTools_refit_required` that
    asks to refit the model with this version of BayesTools (every error that
    asks for a refit has this class; callers match it, not the message).
    `JAGS_check_convergence()` and `mix_posteriors()`
    no longer accept `runjags` objects that were not created by `JAGS_fit()`.
    Mixed posteriors saved from 0.3.0 cannot be passed to
    `marginal_posterior()` (mixed factor posteriors without the factor
    metadata of this version stop with class `BayesTools_refit_required`),
    nor marginal posteriors saved from 0.3.0 to
    `Savage_Dickey_BF()`; rebuild them from refitted models. This includes
    models stored by RoBMA 4.0.0.
  - marginal likelihoods are `BayesTools_marglik` objects.
    `JAGS_bridgesampling()` and `bridgesampling_object()` return them instead
    of bridgesampling `"bridge"` objects, and `models_inference()`,
    `ensemble_inference()`, and `mix_posteriors()` accept only them as
    `marglik`: saved 0.3.0 results and `bridgesampling::bridge_sampler()`
    output are rejected, and bridgesampling methods such as `bf()` do not
    apply to the new objects. Recompute marginal likelihoods with
    `JAGS_bridgesampling()` or wrap a known natural-log value with
    `bridgesampling_object(logml)`.
  - weight-function priors saved from 0.3.0 (alone or inside prior lists,
    mixtures, and bias priors) lack the selection-model specification and are
    rejected by `print()`, `rng()`, `mcdf()`, `mean()`, and `prior_mixture()`;
    recreate them with `prior_weightfunction()`.
  - inverse-gamma priors use the module distribution `dbt_invgamma` instead
    of a gamma prior on `inv_<parameter>`. Syntax from `JAGS_add_priors()`
    therefore needs the BayesTools module (loaded automatically by
    `JAGS_fit()` and in sessions with BayesTools; plain JAGS runs fail with
    "Unknown distribution"), `JAGS_get_inits()` initializes the parameter
    itself instead of `inv_<parameter>`, `inv_<parameter>` is no longer
    monitored, and `rng()` draws from inverse-gamma priors differ from 0.3.0
    under the same seed (the distribution is unchanged). Posterior samples are
    read only from the parameter itself: `JAGS_bridgesampling_posterior()`,
    `JAGS_marglik_parameters()`, `JAGS_marglik_priors()`, and their row
    versions no longer invert the `inv_<parameter>` columns of 0.3.0 fits and
    stop with a missing-parameter error.
- removed functions and arguments:
  - removes the weight-function distribution functions `mdone.sided()`,
    `mdtwo.sided()`, `mpone.sided()`, `mptwo.sided()`, `mqone.sided()`,
    `mqtwo.sided()`, `rone.sided()`, `rtwo.sided()`, and their `_fixed`
    variants (16 functions); use `mpdf()`, `mcdf()`, `mquant()`, and `rng()`
    on a `prior_weightfunction()` prior.
  - removes `JAGS_diagnostics_density()`, `JAGS_diagnostics_trace()`, and
    `JAGS_diagnostics_autocorrelation()`; use `JAGS_diagnostics(type =
    "density")`, `"trace"`, or `"autocorrelation"`.
  - removes `interpret2()` and `interpret_tables()`; use `interpret()` or
    `interpret_records()`.
  - removes the `transform_orthonormal` argument (deprecated since 0.2.14) of
    `runjags_estimates_table()` / `JAGS_estimates_table()`,
    `ensemble_estimates_table()`, and `JAGS_diagnostics()`; use
    `transform_factors = TRUE`. `JAGS_diagnostics()` stops when it receives
    the argument, and the other functions report it as unused.
  - removes the `seed` argument of `JAGS_extend()`. It only set R's random
    seed: the extended chains always continue their own backend random-number
    state and cannot be reseeded. Calls that pass `seed` (including
    `update(fit, sample_extend = )` in RoBMA 4.0.0) stop with "unused
    argument"; drop the argument.
  - new arguments and the removed ones change positional calls:
    `JAGS_fit()` takes `formula_random_prior_list` and
    `formula_random_effects_compile_list` before `chains`;
    `JAGS_bridgesampling()` takes `formula_random_prior_list`,
    `formula_random_effects_compile_list`,
    `formula_random_effects_marginalize_list`, `bridge_context`,
    `bridge_context_node_names`, `repetitions`, and `method` before `maxiter`;
    `runjags_estimates_table()` / `JAGS_estimates_table()` take
    `random_effects_summary`, `simplify_names`, `random_effects_metadata`,
    `remove_random_effects`, `keep_random_effects`,
    `remove_random_structures`, and `keep_random_structures` before
    `remove_diagnostics`; `JAGS_extend()` takes `runtime_setup`,
    `runtime_cache`, and `worker_output` where `seed` was; and the arguments
    after `transform_orthonormal` of `runjags_estimates_table()`,
    `ensemble_estimates_table()`, and `JAGS_diagnostics()` move one position
    forward. Pass these and later arguments by name.
  - `check_bool()` rejects `NA` by default (`allow_NA = FALSE`, was `TRUE`),
    as a missing value cannot be used as a switch. Pass `allow_NA = TRUE` to
    keep the previous behaviour. The `allow_NA` defaults of the other
    `check_*()` helpers are unchanged.
- stricter inputs of released functions:
  - `check_int()` rejects infinite values ("must contain only finite
    values"); 0.3.0 accepted `Inf` and `-Inf` as integers.
  - `JAGS_evaluate_formula()` (and `JAGS_predict_formula()`) evaluates a
    formula only through the formula design that `JAGS_fit()` stores for the
    formula parameter. Posterior samples without it, such as a `coda` `mcmc`
    object of draws, stop with a message naming `JAGS_formula_draws()`, which
    attaches to draws the design that `JAGS_fit()` builds from the same
    formula, data, and priors. 0.3.0 evaluated such samples from the factor
    levels and contrasts of the `prior_list` metadata, and standardized the
    predictors only when the samples carried a `formula_scale` attribute
    (otherwise it used the raw predictors).
  - `transform_scale_samples()` accepts only models fitted with `JAGS_fit()`;
    matrices of posterior samples and plain `runjags` or `mcmc` objects stop.
    Pass the fitted object, whose stored design verifies the transformation.
  - original-scale transforms (`transform_scaled = TRUE` of ensemble and
    marginal tables, `plot_transformed_prior()`, and prior densities with
    `formula_scale`) require the fitted design that the `formula_scale`
    attribute of a fitted model carries. Standardization information without
    it (a `formula_scale` list built by hand) stops with a
    `BayesTools_formula_transform_unavailable` error instead of pairing
    coefficient names, which left the intercept and the factor levels of terms
    such as the level slopes of `~ g + g:x` on the standardized scale; pass
    `attr(fit, "formula_scale")`. `transform_scale_samples()` and
    `transform_prior_samples()` take the fitted structure (design, log
    intercept, random-effect SD structure) from `fit` also when
    `formula_scale` is passed; a list built by hand transformed a log
    intercept as a linear one.
  - `JAGS_marglik_parameters_formula()` requires `formula_design_list` (the
    `formula_design` element of `JAGS_formula()`); the reconstruction from
    formula data alone is removed.
  - `add_bounds` of `JAGS_bridgesampling()` and
    `JAGS_bridgesampling_posterior()` must give `lb` and `ub` as numeric
    vectors named by `add_parameters`, without `NA` and with every lower bound
    below its upper bound.
  - factor priors used outside a formula, e.g., in the `prior_list` of
    `JAGS_fit()`, must carry their factor levels, set with the new
    `prior_factor_levels()`: replace `attr(prior, "levels") <- K` by
    `prior <- prior_factor_levels(prior, K)` (or pass the level names).
    `JAGS_fit()`, `as_mixed_posteriors()`, `mix_posteriors()`, and the
    parameter catalog stop on a factor prior with only a `levels` attribute
    instead of inferring its levels from the number of coefficients, and
    coefficient names come from the factor design instead of choosing, by that
    number, between all level cells and the cells after the reference levels.
  - `JAGS_formula()` (and therefore `JAGS_fit()`) rejects mean-difference and
    orthonormal factor priors, other than a point mass at zero, on a term that
    codes a factor by level indicators because one of its lower-order terms is
    missing, such as `g:x` in `~ g + g:x` or `~ g / x`. BayesTools 0.3.0
    fitted such a prior as independent priors on the per-level coefficients
    (iid N(0, 1) for `prior_factor("mnormal", list(0, 1), contrast =
    "meandif")`) and labelled them as differences from the mean
    (`g:x[dif: L]`). Use `~ g * x` to keep the contrast, or a treatment or
    independent prior for one coefficient per level. Such a term does not use
    the contrast of that factor and records the independent coding, so its
    prior may use another contrast than the factor's other terms; when no term
    codes the factor by its contrast (e.g., `~ g:x + g:z`), its
    indicator-coded terms must agree on the contrast.
  - `JAGS_formula()` (and `JAGS_fit()` formulas) reject a `multiply_by`
    attribute on the intercept prior, including one supplied through
    `__default_continuous`. 0.3.0 accepted it silently: the JAGS model never
    scaled the intercept, while `JAGS_evaluate_formula()`,
    marginal-likelihood reconstruction, and marginal posteriors did, so their
    results did not match the fitted model. Term priors keep `multiply_by`.
  - formula random effects (the `(x || g)` terms supported since 0.2.20) need
    `prior_random()` (`formula_random_prior_list` of `JAGS_fit()`); the
    undocumented 0.3.0 route that read their SD priors from `prior_list`
    entries named `"term|group"` is removed, and formulas with random-effect
    terms but without `prior_random()` stop.
- random numbers:
  - derives each chain's JAGS random-number seed (`.RNG.seed`) in `JAGS_fit()`
    and `JAGS_get_inits()` from `seed` through R's random-number generator,
    `set.seed(seed)` followed by `sample.int(.Machine$integer.max, chains)`,
    instead of `seed + chain`. Chain k + 1 of seed s no longer reuses the
    stream of chain k of seed s + 1, which made prior-only draws of fits with
    adjacent seeds identical, and a chain's seed does not depend on the number
    of chains. The chain seeds therefore depend on R's `RNGkind()` settings.
    Automatic restarts of `JAGS_fit()` after a failed initialization draw
    their seeds from a separate stream of `seed` (the first L'Ecuyer-CMRG
    substream) instead of using `seed + i`, so a restart no longer repeats the
    fit with `seed + i`. Initial values and RNG names are unchanged, but every
    seeded fit, and every result computed from its draws, differs from 0.3.0
    under the same seed. `JAGS_get_inits()` also returns each chain's
    `.RNG.name` and `.RNG.seed` for an empty `prior_list` (0.3.0 returned
    `list()`), so `JAGS_fit(seed = )` is reproducible for models whose priors
    are only in the model syntax. Parallel fits are reproducible for the same
    seed, chains, and cores, but need not equal serial fits.
  - seeded functions no longer reset R's global random-number stream to a
    state determined by `seed`: `JAGS_get_inits()`, `JAGS_fit()` (serial,
    parallel, and autofit), `mix_posteriors()`, `marginal_inference()`, and
    `transform_prior_samples()` restore the caller's `.Random.seed` and
    `RNGkind()`. Seeded results are unchanged.
    Code that relied on the reset, for example an unseeded
    `JAGS_bridgesampling()` right after `JAGS_fit(seed = )`, or unseeded draws
    after a seeded `marginal_inference()`, must pass its own seed. Unseeded
    `JAGS_get_inits()` and `mix_posteriors()` take their seed with one draw
    from the caller's stream.
- inference and convergence results:
  - conditioning on a parameter without an inclusion indicator (a prior that
    is neither spike-and-slab nor a null/alternative mixture), or on an
    unknown label, stops with "The parameter '...' is not a conditional
    parameter." instead of warning and using all draws. This applies to
    `as_mixed_posteriors(conditional = )`,
    `as_marginal_inference(conditional_list = )`, and the conditional
    summaries and plots built on them; remove such parameters from the
    conditioning set.
  - `Savage_Dickey_BF()`, `marginal_inference()`, `as_marginal_inference()`,
    and `marginal_estimates_table()` return a finite Bayes factor when the null
    hypothesis lies outside the posterior draws: the kernel density at the
    null from the Gaussian kernel tails, with a warning that it is not
    reliable evidence. 0.3.0 returned `Inf`. The warning is emitted once per
    parameter or level, labelled by the row label of the summary tables
    (`(mu) x[-1SD]`, `(mu) x:g[-1SD, a]`), and tables list it among their
    warnings.
  - `marginal_inference()`, `as_marginal_inference()`, and
    `Savage_Dickey_BF()` on a list of marginal posteriors return `NA` with the
    reason for a level whose posterior has a declared point mass at the null
    hypothesis (e.g., the reference level of treatment-coded coefficients,
    "fixed at the null hypothesis value"), and compute the other levels; 0.3.0
    returned the invalid density ratio with a warning. Summary tables list the
    reason among their warnings. A scalar `Savage_Dickey_BF()` call stops with
    an error of class `BayesTools_posterior_point_mass_at_null`.
  - `Savage_Dickey_BF()` requires an exactly classified, regular prior
    ordinate at the null hypothesis, the rule of point hypotheses of the new
    `hypothesis_BF()`: an ordinate without an exact structural value stops
    with class `BayesTools_inexact_ordinate` (0.3.0 used a numerical grid
    height), an undefined ordinate with `BayesTools_undefined_ordinate`, and a
    zero or infinite ordinate with `BayesTools_zero_ordinate` or
    `BayesTools_infinite_ordinate` (0.3.0 returned 0 or `Inf`),
    because the density ratio is then a limit that a posterior kernel estimate
    cannot estimate. All have the parent class
    `BayesTools_hypothesis_ordinate`. In list posteriors,
    `marginal_inference()`, `as_marginal_inference()`, and their tables, zero
    and infinite ordinates give `NA` with the reason instead (e.g., a null
    outside the prior support). The "Prior density does not span both sides
    of the null hypothesis" warning is removed, as such a null has a zero
    prior ordinate.
  - prior densities without recorded provenance (e.g., a density grid built
    from prior draws and attached by hand) are used only for plotting: their
    heights, ordinates, and region probabilities stop.
  - `var()` and `sd()` of a truncated normal prior (and of spike-and-slab
    priors with such a slab) are computed in closed form (see Fixes) and stop
    with an error when that value may carry a relative rounding error above
    1e-6: one-sided truncations more than about 33 prior standard deviations
    from the mean, where the numerical integration of 0.3.0 failed as well,
    and very narrow truncation intervals (e.g., 0 to 0.001 SD), for which
    0.3.0 returned a value. The message suggests estimating the variance from
    `rng()` draws, which sample such truncations exactly.
  - `JAGS_check_convergence()` and the autofit of `JAGS_fit()` and
    `JAGS_extend()` classify the fitted parameters by the convergence roles
    stored with the fit (`convergence_role` of `parameter_coordinates()`),
    which are derived from the declared model, never from the draws or name
    suffixes, and share one default selection:
    - declared constants are structural and reported as
      `"structural_constant"`: point priors (now monitored, see below),
      one-state indicators, reference and fixed publication-weight bins
      (including mirrored two-sided bins and the constant bins of composed
      bias priors and bias mixtures), a p-hacking kind shared by every branch,
      the point total of an ordered prior and the coefficients it fixes, unit
      correlation diagonals and Cholesky constants, monitored fully observed
      data, and monitored deterministic nodes whose ancestors in the model
      syntax are all data or constants;
    - a sampled parameter whose draws never change, including a constant
      element of a user `add_parameters` node with a stochastic ancestor, is
      not assessable and fails the check unless `allow_not_assessable = TRUE`
      (0.3.0 counted constant columns as converged); leave such elements out
      with `monitor` or `autofit_control$monitor`;
    - indicators and inclusion probabilities are identified from their priors,
      not from the `_indicator` and `_inclusion` name suffixes; with
      `check_indicators = TRUE`, binary indicators are checked as Bernoulli
      occupancies and categorical indicators separately for every observed
      state;
    - deterministic monitors that formulas generate (correlation matrices,
      Cholesky factors, partial correlations, derived coefficients, latent
      effects, and SDs) are checked only when requested with `monitor`, while
      the user's `add_parameters` are checked as in 0.3.0;
    - `prior_list` is optional; when supplied, it must be the prior list stored
      with the fit (`attr(fit, "prior_list")`), otherwise the call stops.
- labels, tables, and posterior draws:
  - square brackets after a factor term always hold a level label. Contrast
    coefficients that are not a level (mean-difference and orthonormal
    coefficients) are labelled `term{j}`, coefficient j of the contrast coding,
    instead of `term[j]` in `transform_factors = FALSE` summaries and
    diagnostics, in `as_mixed_posteriors()` and `mix_posteriors()` columns and
    the ensemble tables and `plot_models()` built from them; ordinary factor
    priors label theirs `p1{j}`. The parameter catalog selects factor levels by
    label: `mu_g[2]`, `g[2]`, and `(mu) g[2]` mean the level labelled "2",
    never the second fitted coordinate. Every level of every contrast, and of
    ordinary factor priors, is a catalog quantity, and the transformed-summary
    labels `mu_g[dif: 2]`, `(mu) g[dif: 2]`, and `g[dif: 2]` are its aliases.
  - every parameter label is rendered by one renderer from structured label
    parts (see `parameter_labels()`), so a quantity has the same label on every
    route, and labels that names were parsed into change:
    - two-level treatment interactions name their level in every table
      (`(mu) g[hi]:x`, was `(mu) g:x` in model tables);
    - interaction-cell columns of `as_mixed_posteriors()` and
      `mix_posteriors()` are catalog selectors (`mu_g__xXx__h[g=b, h=v]`, was
      `mu_g[b]__xXx__h[v]`), and their ensemble rows are the model-table rows;
    - formula prefixes come from the formula parameter of each prior, not from
      name prefixes: with formula parameters `mu` and `mu_tau`, `mu_tau_z` is
      `(mu_tau) z` (was `(mu) tau_z`), and a predictor `mu_income` of `mu` is
      `(mu) mu_income` in ensemble tables (was `(mu) (mu) income`);
      `format_parameter_names()` replaces a formula prefix only at the start
      of a name (with formula parameter `mu`, `mu_mu_income` is
      `(mu) mu_income` and `tau_mu_x` is unchanged; 0.3.0 gave
      `(mu) (mu) income` and `tau_(mu) x`);
    - marginal names and `marginal_estimates_table()` rows name the levels of
      every factor (`(mu) x:g[-1SD, a]`), and the table's warnings start with
      the row label (was the backend node, `mu_x__xXx__g[-1SD, a]:`); a
      marginal of a parameter without levels is `mu` (was `mu[]`);
    - inference-table rows of mixture components are `(mu) z[narrow]` (was
      `(mu) z [narrow]`), a factor level of a mixture component is
      `(mu) x[b][narrow]` (was `(mu) x[narrow][b]`), and the Bayes factor
      MC-error warnings name the row label (was the backend node, e.g.
      `mu_intercept`);
    - factor plot legends show the level text of each level cell (both levels
      of an interaction cell, `b, v`; without the leading space of transformed
      contrasts; level names with brackets whole).
  - the index of an individual weight-function omega in
    `plot(individual = TRUE, show_figures = )`,
    `lines(individual = TRUE, show_parameter = )`, and
    `geom_prior(individual = TRUE, show_parameter = )` is the k-th omega in the
    order summary tables and RoBMA print them (ascending p-value intervals,
    reference weight first). `plot()` previously counted from the least
    significant interval, and its default `show_figures = -1` now omits the
    reference weight fixed at 1 instead of the least significant weight.
  - `runjags_estimates_table()` / `JAGS_estimates_table()` identify inclusion
    rows from the prior list and the formula metadata instead of the row label,
    and every inclusion row reports only the posterior inclusion probability
    (`Mean`) and the MCMC diagnostics, with `SD` and the quantiles left empty.
    0.3.0 emptied these cells on every row whose label contained "inclusion",
    so a parameter or factor level whose name merely contains that word now
    keeps its summary. Inclusion rows are the spike-and-slab and mixture
    indicators, random-effect inclusion quantities such as
    `(mu) inclusion(sd(x))`, raw variance-allocation gate indicators, and the
    indicators of spike-and-slab totals of ordered priors.
    `format_parameter_names()` likewise treats only the literal "(inclusion)"
    marker as an inclusion row.
  - `JAGS_fit()` monitors point priors (scalar, multivariate, and factor),
    which 0.3.0 left out of the fit (`JAGS_to_monitor()` now returns them):
    their constant values are posterior columns, reported as structural
    constants in convergence checks.
  - posterior draws from `mix_posteriors()`, `as_mixed_posteriors()`,
    `marginal_posterior()`, `random_effects_summary_posterior()`, and
    `parameter_draws()` keep their metadata (supports, atoms, undefined draws,
    prior densities and their context, precomputed posterior densities and
    ordinates, conditioning, level weights, and column quantities) in one
    validated attribute, read and set with the new `posterior_metadata()`.
    The 0.3.0 attributes `models_ind`, `sample_ind`, `prior_density`,
    `prior_density_context`, `conditional`, `formula_parameter`,
    `formula_scale`, `transform_scaled`, and `linear_weights` are no longer set
    or read on draws.
  - arithmetic and mathematical functions (the `Ops` and `Math` group
    generics) of mixed and marginal posterior draws return plain numerics
    without their class and metadata; 0.3.0 kept the attributes, so
    `Savage_Dickey_BF(mp * 2, .5)` silently used the prior density of the
    untransformed draws `mp`. Draws whose values were replaced while
    their attributes were kept (`x[] <- `, `x[i] <- `, `pmin()`, `pmax()`)
    stop with an error of class `BayesTools_stale_metadata` (parent
    `BayesTools_metadata`) wherever their metadata are read or set, detected
    from a fingerprint of the values that the metadata record.
    `Savage_Dickey_BF()`, `marginal_posterior()`, `plot_posterior()`, and the
    other functions that need the metadata stop on plain numeric draws. Both
    errors name `posterior_transform()` (or the `transformation` argument of
    `marginal_posterior()`) as the route for transformed posterior
    distributions.
  - posterior plots take point masses only from the declared posterior atoms
    of the draws (`posterior_atom_attribute()` or the BayesTools functions
    that create draws). `plot_posterior()` and `plot_marginal()` no longer
    infer point masses from the prior list, model indicators, or repeated draw
    values, and draws without a declared atom status stop with the
    atom-status message; constant draws declared continuous are no longer
    shown as point masses.
  - conditional and transformed prior overlays of `plot_posterior()` stop when
    the prior-density context cannot be built instead of plotting the raw,
    possibly unconditional, prior.

### Features
- formula random effects:
  - adds the `prior_random()` interface for formula random effects, with the
    helpers `random_block()`, `random_covariance()`, `random_monitor()`,
    `random_new_levels()`, `random_variance_allocation()`, `allocation_ref()`,
    `random_sd_source()`, and `random_group_covariance()`, and compact print
    methods for these specifications, LKJ priors, parameter sources, random-SD
    sources, and variance-allocation references. `JAGS_fit()` and
    `JAGS_bridgesampling()` take them as `formula_random_prior_list`, and
    `JAGS_formula()` as `prior_random` (with `random_effects_compile`).
    Scalar SD priors with negative support are rejected when the specification
    is constructed (spike-and-slab, mixture, and ordered SD priors component by
    component).
  - parses lme4-like random-effect terms through `reformulas`
    (`random_effects_formula()`): ordinary and independent (`||`) random
    effects, named blocks, nested (`/`, in outer-to-inner order as in metafor)
    and interacted (`:`) grouping factors, factor random slopes, and the
    covariance structures diagonal, shared-SD independent, unstructured,
    compound symmetry, heterogeneous compound symmetry, discrete AR(1),
    heterogeneous AR(1), and continuous-time AR(1).
  - completes the covariance priors of a block after its structure is known:
    unstructured blocks get an LKJ(1) prior on their correlation matrix, set
    with `prior_lkj()` and sampled through the compiled BayesTools JAGS module
    (`JAGS_lkj_corr_cholesky()` writes its syntax; the correlation matrices
    have an exactly unit diagonal), and compound-symmetry, AR, and
    continuous-time AR blocks a uniform prior over their admissible
    correlation interval; explicit scalar correlation priors keep the Fisher-z
    default scale. Factor random slopes have one SD per
    generated coefficient, `sd(g{j})` for mean-difference, orthonormal, and
    ordered bases; scalar SD templates expand to indexed independent priors;
    and random-coefficient blocks reuse the contrast of the fixed factor unless
    `random_block(contrasts = )` overrides it. Factor designs of random slopes
    are built under the memory limit of option
    `BayesTools.random_effects_memory_limit_bytes`.
  - `random_variance_allocation()` splits a total variance across random-effect
    components with Dirichlet weights, optionally with independent Bernoulli
    inclusion gates (also a gate-only allocation of one component). Its
    quantities are the totals `sd_total` / `var_total`, the common scale
    `sd_common` / `var_common`, and the components' `var_prop`, `var_mult`,
    and `sd_mult`; gated totals are summarized on their realized
    model-averaged scale (with the all-off branch at zero), and `var_prop` is
    normalized over the active components and summarized conditional on a
    positive total.
  - block parameterizations `"noncentered"`, `"centered"`, `"mean_centered"`
    (for an eligible scalar random intercept), and `"auto"`, which falls back
    to noncentered for centered-CAR blocks whose SD prior support is not
    bounded away from zero and infinity and for known group covariance with
    several columns (an explicit `"centered"` stops there). The
    parameterization changes neither the priors nor the reported quantities.
  - `random_group_covariance()` gives known covariance or correlation kernels
    across grouping levels. A kernel is converted to a correlation with exact
    symmetry, and kernels that are not positive definite at working precision
    after subsetting to the fitted levels are rejected; the fitted kernel
    multiplier is the block's `sd` / `var`.
  - prediction: `JAGS_predict_formula()` gives fixed, conditional, and
    marginalized predictions with selectable random-effect `blocks` (by
    default the random-effect terms in `formula` for
    `formula_target = "marginal"`), explicit `new_levels` policies, and
    covariance- or sampling-based marginal output. `JAGS_evaluate_formula()`
    gains `formula_target`, `blocks`, `new_levels`, and `fitted_rows`, takes
    `formula`, `data`, and `prior_list` from the fit by default, and evaluates
    fitted random-effect formulas for existing grouping levels when latent
    effects or group-level coefficients were monitored; its default
    `formula_target = NULL` evaluates fits without random effects and stops for
    fits with random effects. Levels declared in a grouping factor but without
    rows in the fitting data are new levels in prediction (they stay fitted
    groups), except in blocks with a known group covariance. Stored
    random-effect term formulas do not keep the environment they were written
    in (saved fits do not carry the caller's objects), and prediction takes
    random-effect predictors from `data` only, stopping when one is missing.
    `random_effects_marginal_vcov()` returns posterior observation-level
    covariance or, with `diagonal_only = TRUE`, variance draws, and
    `random_effects_correlation_draws()` the dense correlation matrices of
    scalar-structured blocks.
  - marginalized random effects: `random_effects_compile()` and
    `random_effects_marginal_variance_factors()` are the contract for
    structurally marginalized blocks (unsupported factor representations
    signal `BayesTools_random_effects_marginal_variance_unavailable`);
    `JAGS_formula_random_marginal_covariance()` compiles the full observation
    covariance for all structures, including row-scale sources, variance
    allocations, and an automatic factor representation;
    `random_effects_marginal_diagonal_factor()` (with `diagonal_support`),
    `random_effects_marginal_factor_diagonal()`,
    `random_effects_marginal_factor_states()` (with an optional `cache` and a
    `contract_id`), `random_effects_marginal_factor_product()`, and
    `random_effects_marginal_factor_vcov()` evaluate exact covariance factors
    across all posterior draws in batches, without dense draw-by-row-by-row
    arrays; `random_effects_marginal_update_plan()` (affine, factor, Markov,
    or unsupported updates, with `source_transform_spec`) and
    `random_effects_marginal_update_grid()` describe exact scalar updates of
    the covariance; and `random_effects_dependency_matrix()`,
    `random_effects_source_roles()`, and `random_effects_level_roles()` report
    the structural row dependencies, source roles, and estimate-level terms.
    Reconstructed covariances are exactly symmetric and use the persisted
    Cholesky factors; invalid or saturated coordinates are rejected rather than
    clamped or repaired.
  - bridge sampling of formula random effects uses the standardized latent
    effects, scalar correlation coordinates, and LKJ primitive coordinates as
    bridge parameters. The `formula_random_effects_marginalize_list` argument
    of `JAGS_bridgesampling()` integrates selected sampled Gaussian blocks
    exactly: their latent coordinates are removed and the likelihood callback
    receives the draw-specific covariance or an exact block factor
    (`bridge_context`, `bridge_context_node_names`, including a nodes-only
    context).
  - summaries: `runjags_estimates_table()` / `JAGS_estimates_table()` gain
    `random_effects_summary`, `simplify_names`, `random_effects_metadata`,
    `remove_random_effects`, `keep_random_effects`,
    `remove_random_structures`, and `keep_random_structures`, and parameter
    filters such as `"random"`, `"random_sd"`, `"random_cor"`,
    `"random_var_prop"`, `"random_var_mult"`, `"random_allocation"`, and
    `"random_sd_mult"`. Random-effect quantities are named
    `(formula) owner: quantity(parameter[level], ...)` (e.g.,
    `(mu) id: sd(x)`, `cor(intercept,x)`); standard tables show the
    prior-facing quantities, ordered by scale, allocation, and correlation,
    full tables every deterministic representation, and raw tables label LKJ
    primitives by their component pair (`(mu) lkj_u(intercept,x | g)`);
    `simplify_names = TRUE` shows a sole random intercept's `sd(intercept)` as
    `sd` (and `parameter_catalog_resolve()` and `hypothesis_parse()` accept
    these simplified aliases with `simplify_names = TRUE`). With
    `transform_scaled = TRUE`, SDs and correlations of standardized random
    slopes are reported on the original scale, and a correlation that is
    undefined in some draws is summarized over its defined draws with a
    footnote giving their share.
    `random_effects_summary_posterior()` returns these quantities as mixed
    posteriors, and `transform_prior_samples()` also returns the random-effect
    monitors that the model defines deterministically (allocation-derived SDs,
    Fisher-z or logit scalar correlations, LKJ Cholesky factors, correlation
    matrices, and partial correlations); latent effects and the auxiliary nodes
    of Dirichlet priors have no prior draws.
  - adds the `RandomEffects` vignette comparing BayesTools formula random
    effects with lme4 and rstanarm.
- fitted metadata, parameter catalog, and labels:
  - every `JAGS_fit()` result stores one versioned parameter map
    (`parameter_map()`) with linked coordinate, quantity, and alias tables.
    `parameter_coordinates()` lists the fitted coordinates with their roles,
    including `convergence_role` (`"sampled"`, `"indicator"`, `"structural"`,
    `"derived"`, or `"auxiliary"`); `parameter_catalog()` lists the public
    quantities with their `label_parts`, exact `support` (from the prior
    provenance, or `NULL`), and `definedness` (`"always"`, `"correlation"`, or
    `"allocation_active"`); the fit also stores a formula name map
    (`JAGS_formula_name_map()`), a fitted-object contract
    (`JAGS_fit_contract()`, `JAGS_validate_fit_contract()`), and the draw
    geometry (`JAGS_draw_geometry()`, `JAGS_materialize_draws()`). Their
    versioned schemas are returned by `parameter_map_schema()`,
    `parameter_coordinates_schema()`, `parameter_catalog_schema()`,
    `JAGS_fit_contract_schema()`, `JAGS_draw_geometry_schema()`, and
    `JAGS_formula_coefficient_transform_schema()`.
    `fit_backend_fingerprint()` identifies the backend code and native
    libraries for cache invalidation. Fits whose metadata are missing or of
    another version stop with class `BayesTools_refit_required`; fits whose
    posterior samples lack coordinates a request needs because of their
    monitoring or sampling settings (a `parameter_draws()` quantity whose
    source coordinates were not monitored, random effects without
    `random_monitor(latent = TRUE)` or `random_monitor(coefficients = TRUE)`
    in `JAGS_evaluate_formula()` and `JAGS_bridgesampling()`, or a
    marginalized random-effect block) stop with its child class
    `BayesTools_refit_monitoring`.
  - `parameter_catalog_resolve()` resolves selectors exactly, with classed
    errors for ambiguous selectors (`BayesTools_parameter_ambiguous`) and for a
    contrast-coefficient selector `term{j}` of a coordinate that is a level
    (`BayesTools_selector_unavailable`, naming the level form, e.g. `g1[10]`
    for `g1{1}` of a treatment factor with levels 5, 10, and 20).
    `parameter_catalog_extend()` adds validated provider quantities (canonical
    names are unique per namespace and component), `parameter_draws()`
    extracts the draws of a selection (also from an already materialized
    posterior matrix) and names the remedy for unavailable quantities,
    `parameter_transform()` with `parameter_transform_forward()`,
    `parameter_transform_inverse()`, and `parameter_transform_jacobian()`
    evaluates one-to-one quantity maps, `JAGS_with_draws()` replaces a fit's
    draws (e.g., by prior draws) while keeping its draw geometry, and
    `parameter_map_cache()` caches map-derived values in a bounded session
    registry, keyed by a caller-supplied token.
  - adds `parameter_labels()`, the renderer of every parameter label (catalog
    selectors, table rows, plot legends, and warnings), from the label parts of
    catalog quantities and of the columns of mixed and marginal posteriors (the
    `quantities` draw metadata, with the fitted coordinates, weights, and
    catalog quantity of each column); `vocabulary` renders random-effect
    quantities under other names. Estimates tables carry a per-row
    `quantities` attribute (row, catalog quantity id, label parts), and the
    rows of ensemble tables carry the quantity ids that the model tables give
    the same values (also transformed factor levels and one-coordinate columns
    such as weight-function weights, `PET`, and `PEESE`); rows whose values
    `transformations` changed have no id.
  - adds the parameter-name helpers `JAGS_indexed_parameter_columns()`,
    `JAGS_indexed_parameter_matrix()`, `JAGS_indexed_parameter_vector()`, and
    `JAGS_regex_escape()`, `JAGS_formula_coordinate_dependencies()`, which
    reports the fitted formula predictors that depend on given coordinates,
    and `JAGS_formula_internal_coordinate_priors()`, the exact priors of
    stochastic formula coordinates outside the fitted `prior_list` (the beta
    primitives of LKJ blocks).
  - adds `JAGS_deterministic_nodes()`, which lists the deterministic nodes
    BayesTools generates in a fitted JAGS model (allocation-derived
    random-effect SDs, scalar and LKJ correlations, publication weights
    `omega`, spike-and-slab and mixture parameters, and formula linear
    predictors) with their family, coordinates, and dependencies;
    `JAGS_evaluate_deterministic()`, which recomputes them from supplied draws;
    and `JAGS_deterministic_evaluator()`, which returns such an evaluator with
    the node registry resolved once. Each node is defined once: the model
    syntax (the formula linear predictor written by `JAGS_formula()` and the
    weights written by `selection_backend_spec()`) is emitted from these
    definitions, and prior draws, the parameter catalog, bridge sampling,
    marginal-likelihood parameters, prediction, marginal posteriors of formula
    parameters, random-effect unscaling, and convergence roles evaluate them.
  - adds `parameter_source()`, which describes an existing JAGS node that
    generated model components may reference, and `random_sd_source()`, which
    marks one as an external random-effect SD source (e.g., a scalar or
    row-shaped `tau`). A `values` function requires `inputs`, the posterior
    coordinates it reads (`character()` for a function of `data` alone): they
    are the source's dependencies in `JAGS_deterministic_nodes()`, and the
    function receives only them; without them it stops with class
    `BayesTools_missing_source_inputs` (also `BayesTools_parameter_source`).
- posterior draws and their metadata:
  - adds `posterior_metadata()` and `posterior_metadata<-`, which read and set
    the metadata of BayesTools posterior draws: `"support"`, `"atoms"`,
    `"undefined_draws"`, `"prior_density"`, `"prior_densities"`,
    `"prior_context"`, `"posterior_density"`, `"posterior_densities"`,
    `"posterior_ordinate"`, `"posterior_ordinates"`, `"condition"` (with the
    logical `averaged`, `TRUE` for unconditional model-averaged draws),
    `"linear_weights"`, `"quantities"`, and `"output_transformations"`. The
    metadata record a fingerprint of the values they describe (see the
    Breaking changes), and subsetting a list of mixed posteriors with `[`
    keeps the list and its metadata, with the prior densities of the kept
    parameters.
  - adds the constructors `posterior_support_attribute()`,
    `posterior_atom_attribute()`, `posterior_density_attribute()`,
    `posterior_ordinate_attribute()`, and `posterior_ordinate_append()`, and
    the predicates `posterior_ordinate_has_value()`,
    `posterior_ordinate_supports_bf()`, and `posterior_atoms_free()` (whether
    draws declare their atoms and have none). Metadata are accepted only as
    these objects, and alias field names (e.g., `estimator` or `height`) are
    not read. A precomputed posterior density describes only the
    continuous part of the posterior; its point masses are the `atoms` of the
    draws, which `plot_posterior()` and `plot_marginal()` with
    `density_method = "precomputed"` draw with the stored density.
  - adds the `density_method` argument (`"KDE"` or `"precomputed"`) to
    `Savage_Dickey_BF()`, `plot_posterior()`, and `plot_marginal()`, and to
    `marginal_inference()` and `as_marginal_inference()` (only `"KDE"`), and
    `posterior_density_method_match()` and
    `posterior_density_method_uses_precomputed()` for packages with their own
    posterior density estimators (such as qCMDE or IWMDE in RoBMA), which map
    their method names to `"precomputed"`.
  - adds `posterior_transform()`, which transforms posterior draws together
    with their metadata: supports, atoms, prior densities, stored posterior
    densities and ordinates (by the Jacobian), component supports, linear
    weights, and label parts, and records the applied transformations in the
    `output_transformations` metadata. Transformations that are not strictly
    monotone and invertible stop with class
    `BayesTools_nonmonotone_transformation`, and draws outside the domain with
    `BayesTools_transformation_domain` (both with parent class
    `BayesTools_transformation`). `marginal_posterior(transformation = )`
    applies it to the untransformed marginal posterior, so stored posterior
    densities and ordinates are transformed rather than dropped, an attached
    prior density is transformed instead of rebuilt from `prior_list`, and the
    quantities keep their labels with the transformation recorded; mixed
    posteriors already transformed with `posterior_transform()` are refused.
  - adds `parameter_mixed_posterior(fit, selection, conditional = FALSE)`: the
    posterior draws of one catalog quantity as a mixed posterior with its
    catalog `support`, its `prior_density` from `parameter_prior_density()`,
    its `condition`, `undefined_draws` (undefined draws are omitted), and
    declared `atoms`, derived from the quantity's structure and never from
    draw values: the point components of mixture and spike-and-slab priors,
    with masses from their indicators; the inclusion-gate atoms of gated
    random-effect SDs, variances, and allocation totals; the value each branch
    of weight-function, publication-bias, and bias-mixture priors gives
    `omega`, `PET`, `PEESE`, and the p-hacking `alpha`, `pi_null`,
    `beta_null`, and `phack_kind`, with masses the shares of the draws of that
    branch; and no atoms for quantities that cannot take a point mass, decided
    from their dependencies in `JAGS_deterministic_nodes()` (e.g.,
    original-scale `us()` correlations of blocks whose SDs have no point
    component or own gate; a gate or point scale shared by every SD of the
    block does not put an atom on these scale-invariant correlations). The
    correlations of blocks whose SDs have their own gates stay undeclared.
    `conditional = TRUE` keeps the draws of the quantity's inclusion event and
    restricts the prior density to it (a total with an ungated component,
    whose event is certain, cannot be conditioned). `parameter_gate_states()`
    returns the per-draw gate and point states and, as `prior_atoms`, the
    prior masses of the atoms of mixture, spike-and-slab, and selection
    priors. `random_effects_summary_posterior()` is built on it, so totals,
    common scales, and gated variance proportions carry their prior measure and
    declared atoms.
  - `plot_posterior(prior = TRUE)` and `plot_marginal(prior = TRUE)` of draws
    that declare no prior (`prior_none()`, such as the
    `parameter_mixed_posterior()` draws of original-scale correlations without
    a prior density) draw the posterior alone with a warning of class
    `BayesTools_prior_curve_unavailable`.
  - adds `JAGS_formula_draws()`, which attaches to posterior or prior draws
    without a fit (e.g., the prior draws of a model that was not fitted) the
    formula design that `JAGS_fit()` builds from the same formula, data,
    priors, and standardization, and from the model data (`model_data`, the
    `data` of `JAGS_fit()`) that `expression()` terms read, so that
    `JAGS_evaluate_formula()` and `JAGS_predict_formula()` evaluate them
    exactly as the fit, on the fitting data or on new data. Draws that lack a
    coefficient of the formula stop with a message naming it.
  - `transform_prior_samples()` also returns the `_indicator` columns of
    mixture priors and of the totals of ordered priors, and the `_indicator`,
    `_inclusion`, and `_variable` columns of spike-and-slab priors, taken from
    the components of their `rng()` draws, so catalog quantities of these
    priors can be evaluated on prior draws; the other columns keep their
    random-number stream.
- prior densities and hypotheses:
  - adds `parameter_prior_density()`, the prior density of a catalog quantity:
    fitted coordinates and one-to-one transforms, factor levels, squared
    nonnegative scales (with integrable singularities at zero),
    variance-allocation marginals, variance proportions, and allocated
    random-effect SDs, variances, and totals (the scale prior times the
    Dirichlet multiplier of the fitted model, for a nested allocation also the
    gates and Dirichlet shares of its parents, with exact atoms at zero for
    gates that are off and continuous parts from the one-dimensional integral
    of the scale prior over each Beta share margin, so their ordinates are
    exact; e.g., the total SD of a gated allocation split into two nested
    blocks is the scale prior times the gate), the SDs and variances of
    gate-only allocations (the scale prior times the gate), the Bernoulli
    prior of `inclusion(...)` quantities
    (support {0, 1}), and the LKJ marginal `2 B - 1`, `B ~ Beta(eta - 1 + K/2,
    eta - 1 + K/2)`, of the correlations of unstructured blocks with K columns
    (also with gated or point-mass SDs, and on the original scale for pairs of
    coefficients that are each one rescaled fitted coefficient, e.g., two
    scaled slopes). Original-scale correlations that mix fitted coefficients
    (the intercept and a slope of a centred predictor) have no prior density,
    and SD components and totals of nested allocations with two or more
    independent allocation shares keep a numerical grid used for plotting,
    with their point hypotheses refused. It returns `NULL` only when
    no fitted prior owns the quantity's coordinates and stops for an owning
    entry that is not a BayesTools prior.
  - adds `prior_density_ordinate()`, the exact-value structural
    classification of prior ordinates (point masses, continuous limits,
    deterministic mixtures, analytic normal combinations, and named
    transformations) with `provenance$continuous_behavior` (the behaviour of
    the continuous part at a point mass). `exact = TRUE` always comes with an
    available `log_density`; ordinates without a value (a rejected quadrature,
    a boundary limit without a structural value) and ordinates computed from
    subnormal, underflowed, or overflowed intermediate values (e.g., a
    `multiply_by` product at a subnormal distance from its offset, or
    log-scale terms beyond about +/-709) are reported with `exact = FALSE` and
    the reason, so point hypotheses there are refused as inexact.
  - linear combinations of multivariate t priors (`"mt"`, and `"mcauchy"`
    with one degree of freedom) are univariate t: the levels, level
    contrasts, and other linear targets of mean-difference and orthonormal
    `prior_factor("mt", ...)` priors have exact prior densities, ordinates,
    and region probabilities (location mu sum(a), scale s ||a||, and the
    prior's degrees of freedom; alone or combined with other terms through
    the exact routes, and as weighted sums in mixtures and spike-and-slab
    priors), so point hypotheses on them are no longer refused with
    `BayesTools_inexact_ordinate`. `prior_density_ordinate()` records the
    reduction in `provenance$multivariate_t` (in each component's provenance
    for mixtures).
  - adds `prior_density_has_provenance()`, which tells whether BayesTools
    evaluates the heights, ordinates, and region probabilities of a prior
    density (`FALSE` for a density grid without its recorded prior measure and
    for a product without a structural route); the ordinate method
    `"unsupported_provenance"` is not that signal, as it also marks
    combinations without a structural route, which have refined grid
    probabilities.
  - adds `prior_ordinate_status()`, the exactness rule of point hypotheses as
    data: one row per value (at least one is required) with `eligible`, the
    `condition` class and `reason` a point hypothesis at that value stops with
    (`NA` when eligible), and the `continuous_behavior` of the prior without
    its point masses.
  - adds `JAGS_formula_coefficient_transform()` and
    `JAGS_formula_prior_density()`: the fitted-to-original coefficient map of a
    scaled formula (with each target's `map_type`, `"identity"`, `"affine"`,
    `"exp_affine"`, or `"unsupported"`, and the image of its map in
    `support`) and the induced original-scale prior density of a target or of a
    weighted combination of targets (`weights`, e.g., a mean-difference level
    of a scaled formula), including structural point values, interactions,
    log-intercept Jacobians, and model-mixture atoms; nonlinear maps stop with
    class `BayesTools_formula_prior_density_unavailable` (reason
    `"nonlinear_map"`). `JAGS_formula_predictor_basis()` returns exact
    observation-level update bases from the fitted design, with a `coordinate`
    field that is `"log"` for a logged formula intercept that is the only
    moving coordinate (its update is additive in the log of the intercept) and
    `"identity"` otherwise; nonlinear or coupled directions are reported as
    `non_affine`.
  - adds `hypothesis_BF()` for expression-based point and region hypotheses on
    posterior and prior draws or marginal posterior objects, including
    comparisons between named factor levels, with `density_method` `"KDE"`,
    `"normal"`, or `"precomputed"`. Explicit comparisons accept a region with
    prior mass one (e.g., `"mu > 0.5 vs mu > 0"` under a half-normal prior);
    implicit statements require a complement with positive prior mass. Point
    hypotheses require an exactly classified, regular prior ordinate at the
    null and otherwise stop with the classes of `Savage_Dickey_BF()` (see the
    Breaking changes); a linear expression of marginal-posterior levels (e.g.,
    `mu[B] - mu[A] = 0`) is tested as a Savage-Dickey ratio of the linear
    combination with the joint prior density and a boundary-reflected
    posterior ordinate estimated per component; nonlinear expressions of a
    deterministic prior are refused; and user-supplied prior draws keep their
    kernel estimate with a `BayesTools_inexact_ordinate` warning. Region
    hypotheses integrate the prior up to the exact region boundaries, exactly
    or by the quadratures of `prior_density_ordinate()` where a structural
    route exists; other densities and regions (e.g., `abs()` conditions) use a
    refined grid, a region with an infinite bound under a heavy-tailed
    multiplier and a region whose probability underflows to zero stop, and
    combinations with very many mixture components keep the grid. Labels show
    numeric literals by their shortest representation that round-trips to the
    exact value (`mu <= 0.2`), several statements on one quantity are labelled
    `mu (1)`, `mu (2)`, a region's `error%(BF)` omits the prior Monte Carlo
    variance when the prior mass is exact, `columns` takes the result column
    names, `prior` accepts mixture and spike-and-slab priors for numeric draws,
    and a seeded call restores the caller's random-number state. References
    to quantities or levels that the posterior does not contain, and a
    `parameter` it does not contain, stop with class
    `BayesTools_parameter_not_found` (also
    `BayesTools_parameter_resolution_error`), as unresolved catalog selectors
    do, also in `hypothesis_linear_target()`.
  - adds a versioned hypothesis syntax tree: `hypothesis_parse()` (with
    parameter-catalog resolution of aliases, unquoted non-syntactic aliases,
    several level references in one hypothesis, colon-separated interaction
    levels such as `factor:moderator[level]`, level names in backticks or the
    catalog's escaped form such as `mu[w%7B1%7D]` for the level "w{1}", and
    symbolic equalities such as `theta = phi` normalized to a difference from
    zero),
    `hypothesis_ast_schema()`, `hypothesis_render()`, `hypothesis_symbols()`,
    `hypothesis_rewrite()`, `hypothesis_resolve()`, and the parser helpers
    `hypothesis_parse_point_reference()`, `hypothesis_parse_level_reference()`,
    and `hypothesis_normalize_level_references()`. `hypothesis_resolve()`
    stops a hypothesis without parameter symbols (e.g. `"1 > 0"`) with class
    `BayesTools_hypothesis_no_parameters`, and a level-qualified symbol whose
    level differs from `component` with class
    `BayesTools_hypothesis_component_mismatch` (both also
    `BayesTools_parameter_resolution_error`).
  - adds `hypothesis_linear_target()`, which compiles a point or simple region
    hypothesis linear in the levels of one parameter (level differences,
    scaled levels, level averages) to a scalar target with its combined
    weights and offset, exact joint prior density, and conditioning. Targets
    that cannot be certified stop with class
    `BayesTools_linear_target_unavailable` (also
    `BayesTools_hypothesis_target`) with `reason` `"posterior_atoms"` (the
    target's prior has a point mass), `"atom_declarations"` (a level lacks its
    posterior-atom declaration), or `"prior_context"` (no valid joint prior
    context).
  - gives refusals classed conditions, so callers match the class instead of
    the message: `BayesTools_point_mass_at_null`,
    `BayesTools_infinite_ordinate`, `BayesTools_zero_ordinate`,
    `BayesTools_undefined_ordinate`, and `BayesTools_inexact_ordinate` (each
    also `BayesTools_hypothesis_ordinate`) for point hypotheses;
    `BayesTools_posterior_point_mass_at_null` (also
    `BayesTools_hypothesis_ordinate`) for a declared posterior point mass at
    the null of a scalar `Savage_Dickey_BF()` call; and
    `BayesTools_missing_monitored_columns` (also `BayesTools_marglik_input`)
    for marginal-likelihood samples that lack monitored coordinates the priors
    read in `JAGS_marglik_parameters()`, `JAGS_marglik_priors()`, and their
    row and formula versions.
- priors and selection models:
  - ordered prior plots and layers quietly omit whole continuous curves whose
    infinite density is introduced by random allocation; atoms and intrinsically
    singular totals remain. Factor posterior/marginal overlays omit declared
    reference rows consistently, retain zero and point levels, preserve the full
    retained level color/label list, and default to matching dashed prior curves.
    Fixed-total allocation shares check both exact endpoints. Multiple ordered
    ggplot selections preserve figure indices with `NULL` holes for omitted
    levels; each posterior legend key uses only its level's curve or point glyph.
    Mixed point-total displays recover declared atoms from exact routes before
    omitting allocation-singular curves, including transformed atom locations
    with unchanged probability mass; unavailable recovery reports
    `BayesTools_ordered_prior_display_unavailable`.
    All-hidden standalone selections report `BayesTools_ordered_prior_display_empty`.
    Certified curve omission is decided before plotting density generation,
    avoiding unused quadrature and KDE evaluations while public `density()`
    continues to evaluate every level.
  - supported scalar ordered literal point totals, point-slab totals and
    mixtures whose total components are all literal points use the existing exact
    mixed-measure density route: intermediate Dirichlet shares have scaled-Beta
    curves and declared atoms, and fixed-share/full-total levels retain their
    exact point masses. `force_samples = FALSE` does not sample, store draws,
    or advance the RNG for these totals. `force_samples = TRUE` attaches the
    original seeded sampled values while retaining these analytic curves and
    atoms. Other discrete totals, including Bernoulli, retain the existing
    sampled fallback. Mixed-measure plot axis labels reflect whether the
    displayed component has continuous density, atoms, or both. Grid-integral
    diagnostics report the actual trapezoid integral without renormalization.
  - `plot_models()` uses ordered semantic level summaries and selected level
    draws. Posterior-only plots skip prior summaries and unused prior transforms;
    ordered prior means use
    declared allocation expectations and intervals use structural generalized
    quantiles, respecting exact atom jumps for signed discrete targets.
    Unavailable prior means/intervals report the model and level with
    `BayesTools_prior_interval_unavailable` and the `prior = FALSE` remedy.
  - adds `prior_ordered()` for ordered-factor priors that separate a scalar
    total effect from fixed or Dirichlet allocations across cumulative level
    increments, with `is.prior.ordered()` and the contrasts
    `contr.ordered_cumulative()` and `contr.ordered_cumulative_levels()`. The
    slices of a spike-and-slab total (an ordered factor in an interaction with
    another factor) share one inclusion indicator per draw, also in prior
    draws; standard and full estimates tables omit the internal gamma
    allocation nodes (`prior_par_eta_*`); and a level contrast with
    Dirichlet(1, ...) shares has an infinite prior density at 0. Estimates
    tables summarize ordered terms like the other contrasts: with
    `transform_factors = TRUE` they show the level effects labelled by their
    level cells (`g[b]`, not `g[dif: b]`), without the zero reference level
    of a cumulative contrast and without the total, which equals the last
    level; with `transform_factors = FALSE`, model and ensemble tables show
    the sampled parameters, the total (`g_ordered_total`) and the shares of
    random allocations (`g_ordered_allocation[<level>]`, the share of the
    total added when reaching the level), on the fitted scale. Ensembles retain
    exact selected source rows: fixed allocations fill sampled union rows and
    absent-model shares are undefined, with their defined fraction footnoted.
    Zero totals do not erase an ordered model's primitive allocation parameters.
  - adds `JAGS_ordered_parameter_spec()` for authoritative ordered source labels,
    primitive coordinates and tensor projections, and
    `JAGS_ordered_density_kernel()` for batched active-total and allocation
    density evaluation, separate from bridge eligibility. Ordered allocation
    and coefficient families share JAGS emission and deterministic replay.
    Fixed numeric arguments round trip through JAGS syntax; prior draws and
    seeded initialization retain their random-number streams. Fits missing the
    persisted literal provenance require refitting.
  - ordered scalar marginal atoms follow declared tensor identities and total
    states, including zero allocations, fixed totals and the full-simplex last
    level. Semantic point-state draws equal their declared locations exactly.
    Marginal declarations survive supported extraction, transforms and row
    subsetting; they do not fabricate joint atoms. Inclusion conditioning uses
    the stored ordered total event, shared across slices. Continuous total
    families remain continuous with expression parameters, and spike totals
    retain their spike atoms when the slab has expression parameters. Expression
    point locations remain structurally unavailable without a certified ancestor
    recipe; scalar consumers report `BayesTools_ordered_expression_unavailable`
    with a supported-prior or fitted-snapshot remedy. Observed constant draws
    never supply such a recipe.
  - `parameter_draws(model_samples = )` for ordered quantities needs the total,
    Gamma allocations and component indicators listed in the ordered accessor's
    `source_coordinates`, as well as the selected increment columns. Supply the
    fitted primitive sources or omit `model_samples` to use complete fitted draws;
    missing supplied sources report `BayesTools_ordered_coordinates_unavailable`.
    Supplied primitives determine resolved ordered semantic values even when
    monitored increment columns are stale; declared transformed views keep their
    existing transformation handling. Whole mixed posteriors, factor levels and
    marginal views use the same primitive projections for every available value,
    including continuous cumulative effects, while preserving undefined rows.
  - conditional ordered mixtures reweight model probabilities by the declared
    inclusion event as well as selecting eligible original draws. Posterior
    event fractions and declared prior event probabilities condition their
    respective model weights; totals without an inclusion event keep all
    states. Independently conditioned ordered parameters are not joint draws.
  - adds Dirichlet simplex priors (`prior("dirichlet", ...)` /
    `prior("simplex", ...)`, `is.prior.simplex()`), with random generation,
    log-density, marginal distribution helpers, JAGS syntax, initialization,
    posterior extraction, and bridge-sampling support.
  - adds moment and inverse-moment nonlocal priors, with R density,
    distribution, quantile, and random-generation functions and JAGS module
    support for prior-only and formula models.
  - adds `prior_factor_levels(prior, levels)`, which gives a factor prior
    (including mixture, spike-and-slab, and ordered factor priors) the complete
    factor metadata that formula factor terms carry: the number and names of
    its levels, its contrast, and the design mapping its coefficients to the
    levels. `levels` is the number of levels or their names.
  - adds `selection_model()` (with `check_selection_model()` and
    `selection_model_spec()`) and `prior_weightfunction(model = )`, which
    specify separate choices for the estimate random effect, the other random
    effects, and the complete sampling error (integrate, condition, integrate
    by default). Conditioning retains the whole sampling-error vector, also in
    univariate models; at most one estimate-level random term is resolved
    through its one-to-one grouping map; publication groups (`group`) apply
    only to `weight_rule = "best"`, and under the default `"product"` rule,
    which is invariant to publication partitions, `print()` reports them as
    unused. Product and best p-value rules of vector selection compile to fixed
    branch metadata without additional sampled parameters, and the
    weight-function prior itself is unchanged.
  - adds the selection helpers `selection_qmc_design()` (deterministic
    shifted-Halton designs with a local random-number stream),
    `selection_event_support()` (structural support of Gaussian selection
    events), `selection_context_validate()`, `selection_context_subset_rows()`,
    `selection_context_subset_observations()`, `selection_row_arg()`,
    `selection_native_kernel_args()`, and `selection_native_static_args()`,
    and `selection_backend_spec(include_init = )`: with `include_init = FALSE`
    the specification has no initial values and its compilation does not
    change the random-number state.
- fitting, runtime, bridge sampling, and model averaging:
  - `JAGS_fit()` gains `jags_modules` and carries the generated
    `add_parameters`, `required_packages`, and `jags_modules` of formulas into
    fitting, extension, convergence checks, and parallel workers;
    `BayesTools_load_JAGS_module()` loads the BayesTools JAGS module.
  - `JAGS_fit()` and `JAGS_extend()` gain `runtime_setup` (a callback that
    receives the chain and process topology, coordinator and worker roles, and
    the start and finish of the run), `runtime_cache` (process-sharded cache
    capture and restoration across saved fits and extensions; a failed capture
    keeps the valid draws, and a dead parallel worker drops the whole
    capture), and `worker_output` (captures worker output without storing a
    machine-specific log path in the fit). Parallel workers
    are checked for the required package versions, loaded R function
    definitions, and native builds before fitting or extending, a connection
    failure is not retried and keeps the original backend error, and cleanup
    reaches every worker after a failed stop. `JAGS_runtime_cluster()` and
    `JAGS_runtime_cluster_stop()` export this cluster handling for packages
    built on BayesTools (workers start with `Rscript --vanilla` and the calling
    session's library paths, receive the given `options` and all current
    `BayesTools.*` options, and load the calling session's builds of
    `packages`), and `JAGS_package_builds()` exports the build fingerprints
    (version, R code hash, native library checksums) these checks compare.
  - `JAGS_check_convergence()` gains `monitor` (request parameters for the
    checks, including derived and auxiliary ones) and `allow_not_assessable`,
    and returns a `diagnostics` attribute that classifies every available
    parameter as `"assessable"`, `"structural_constant"`, `"not_assessable"`,
    `"not_requested"`, or `"not_checked"` (after `fail_fast` stops). With
    `check_indicators = TRUE`, eligible indicators are added to a `monitor`
    selection. The `autofit_control` of `JAGS_fit()` and `JAGS_extend()` takes
    the same `monitor` and `allow_not_assessable`, and unknown `monitor` names
    are rejected before sampling or extension.
  - `JAGS_bridgesampling()` gains `seed` (a seeded call restores the caller's
    random-number state), `repetitions`, `method`, `cores` (forwarded to
    `bridgesampling::bridge_sampler()`), and `nonfinite`
    (`"drop"` aggregates only the finite repetitions with a warning), derives
    the effective sample size from the fitted chains, evaluates fixed-parameter
    models without sampled parameters exactly, and takes formula priors from
    the fitted formula design: a `prior_list` entry for a formula parameter is
    ignored when it specifies the same prior (distribution, parameters,
    truncation, contrast, `multiply_by`, mixture components) and is an error
    when it differs.
  - adds `JAGS_marglik_priors_rows()` and
    `JAGS_marglik_priors_rows_evaluator()` for row-preserving prior-density
    evaluation with reusable compiled evaluators. Ordered priors accept numeric
    matrices and data frames with the same strict stochastic Gamma boundaries
    and missing-coordinate conditions; shared allocations contribute once.
    `JAGS_marglik_parameters()`
    evaluates spike-and-slab, mixture, and publication-bias mixture priors for
    a draw from its active component instead of refusing them (bias mixtures
    return `omega`, `PET`, `PEESE`, and the active branch's p-hacking
    parameters); `JAGS_marglik_parameters_formula()` gains `model_data` and
    `JAGS_marglik_priors_formula()` gains `prior_list`.
  - `compute_inference()`, `inclusion_BF()`, `ensemble_inference()`,
    `models_inference()`, and `mix_posteriors()` gain `on_failure` for missing
    (`NA`) marginal likelihoods of models with positive prior probability:
    `"error"` (the default), `"drop"` (remove the failed models and
    renormalize), or `"zero"` (zero evidence), with a warning and audit
    metadata. `bridgesampling_object(NA)` stores a failed computation in a
    model list.
  - `as_marginal_inference(compute_BF = FALSE)` returns the averaged and
    conditional marginal posteriors without computing inclusion Bayes factors.
- credible bands by Harrell-Davis quantiles:
  - adds `harrell_davis_quantile()`, which estimates quantiles of equally
    weighted draws (a vector, or a matrix with the draws in rows) as a weighted
    average of all order statistics (Harrell and Davis, 1982). Applied to the
    columns of a draws-by-grid matrix it gives smooth pointwise bands: the
    kinks of the empirical quantile along the grid are attenuated at the same
    estimand and, in simulations of smooth curves, about the same accuracy.
    The bands are pointwise, and every draw contributes: below roughly 200
    draws and with heavy tails the estimate can differ much from the empirical
    quantile, and next to the edge of a point mass it blends the point mass
    with the continuous draws. A column with an infinite draw uses the
    empirical quantile; an undefined result stops with an error of class
    `BayesTools_harrell_davis_undefined`.
  - the median and 95% band of posterior PET-PEESE lines
    (`plot_posterior()` with `parameter = "PET"`, `"PEESE"`, or `"PETPEESE"`)
    are Harrell-Davis quantiles of the draws of the line (one call, so the
    median lies inside the band) instead of empirical quantiles; only the line,
    the band, and the automatic y range that follows them change. Prior lines
    computed from sampled draws keep the empirical quantiles.

### Fixes
- prior densities and distribution methods:
  - prior densities of linear combinations (the prior densities of formula
    levels and marginal means from `marginal_posterior()`, of
    `as_mixed_posteriors(transform_scaled = TRUE)` and
    `plot_transformed_prior()`, and the prior ordinates of Savage-Dickey Bayes
    factors) use an exact route wherever one exists instead of a numerical
    grid: point masses and scalar distribution functions; sums of normal terms;
    a normal part plus one other continuous scalar term (a Gaussian-convolution
    quadrature split at the other term's support bounds and quantiles and
    around the Gaussian peak, and a closed form for a truncated normal term);
    `multiply_by` products and the levels of ordered factors (total times
    share) as conditional-normal scale mixtures over the multiplier's support
    (their capped product grids were off by up to 20%); the exp of a
    log-transformed positive term plus normal terms (e.g., the unscaled
    intercept of a `log(intercept)` formula scaling) as the product of the term
    and the lognormal exp of the normal part, and the log-scale sum as the log
    image of that product (for a log-transformed term of weight one with
    normal other terms; other such combinations keep a numerical grid); sums
    of two non-normal terms by a two-term convolution quadrature, and sums of
    Cauchy terms in closed form; and mixtures (model-averaged, conditional,
    spike-and-slab, and mixture priors,
    mixtures of factor priors, and covariate rows) as weighted sums of
    per-component ordinates, each component, design row, and mixture leaf with
    its own evaluation budget and convergence check. Region probabilities and
    plotted curves use the same routes. Heights are zero outside a known prior
    support, and a null at a density jump of one component uses the one-sided
    limit inside that component's support.
  - combinations without an exact route keep a numerical grid that starts from
    a spacing resolving the narrowest component, halves the spacing at every
    refinement step while shrinking the omitted tail, and is accepted only
    when its change, like the reported error of a quadrature, is at most 1e-4
    of the value (values that cannot meet the bound are inexact, and point
    hypotheses there stop with `BayesTools_inexact_ordinate`, while plotted
    curves still draw them); grids that would exceed 2^21 knots and
    mixtures of incompatible scales stop with a clear error instead of
    returning a biased height. Known limitations: scale mixtures with
    heavy-tailed multipliers can stop as non-convergent; a component
    quadrature that evaluates to exactly zero (e.g., a narrow component far
    from the value) stops the ordinate, also of a mixture; very narrow Gaussian
    peaks can stop where the quadrature reaches floating-point resolution; the
    scale peak of a value far in a heavy-tailed multiplier's tail can be missed
    without a convergence failure; and grid heights of components with
    singular source densities (e.g., gamma shape below 1) remain approximate
    within the refinement criterion, which is not an error bound.
  - linear-combination prior densities, `marginal_posterior(prior_samples =
    TRUE)`, and allocation margins accept priors whose density is infinite but
    integrable at a truncation bound (gamma or beta shapes below one): the
    boundary grid cell keeps its exact mass.
  - `exp_lin` transformations use the analytic limit at a zero source value.
    Nonlinear transformations keep every strictly increasing knot, saturating
    ones omit outer knots whose transformed value or density is not
    representable (e.g., `tanh` at +/-1, `exp` under- or overflow) instead of
    returning `Inf` or failing, finite endpoints of transformed supports are
    kept exactly (e.g., for negative powers of truncated positive priors), and
    negative powers of inverse-gamma variables have the correct density at the
    zero boundary (infinite, finite nonzero, or zero according to the shape and
    the power).
  - row-varying prior densities are mixed on the linear-predictor scale and
    then transformed once, instead of building grids of tens of millions of
    knots under `exp` or `tanh`.
  - linear combinations in which a `multiply_by` scale also enters as its own
    term (e.g., `x * sigma + sigma`) stop instead of being convolved as
    independent components.
  - prior densities of monitored formula coefficients are the coefficients'
    own priors: a `multiply_by` scaling, which applies only to the linear
    predictor, is no longer applied in `marginal_posterior(use_formula =
    FALSE)`, the stored densities of `as_mixed_posteriors(transform_scaled =
    TRUE)`, `plot_transformed_prior()`, and conditional prior overlays of
    `plot_posterior()`. Linear-predictor results are unchanged.
  - density grids report the FFT clipping mass in probability units
    (including grid spacing and component mass); a continuous density that
    underflows to an all-zero grid, overflows, or has a collapsed range fails
    explicitly instead of returning a purported normalized density or
    inventing a point mass.
  - distribution methods (`rng()`, `cdf()`, `quant()`, `mquant()`, `mean()`,
    `var()`, `range()`, `density()`, ...) stop with a clear message for prior
    classes they do not support instead of returning a function or the prior
    object.
  - far-tail truncated normal, t, and Cauchy priors sample and invert
    exactly, and `rng()` of factor spike-and-slab priors honours
    `transform_factor_samples = FALSE`.
  - `mean()` and `var()` of truncated normal priors are computed analytically
    instead of by numerical integration, which returned wrong values a few
    prior standard deviations into a tail (at 8 SD, a mean of 7.58 instead of
    8.12 and a variance 284 times too large) and for narrow truncation
    intervals (a variance 12 times too large for 5 to 5.001 SD), and failed
    from about 10 SD. The variance stops where its rounding error may exceed
    1e-6 (see Breaking changes).
  - natural prior-support bounds, structural probability endpoints,
    reference-bin weights, and overlapping bridge bounds are compared by exact
    equality instead of a numerical tolerance, so finite truncations close to
    a support bound stay truncated.
- Savage-Dickey Bayes factors and marginal inference:
  - the posterior ordinate of `Savage_Dickey_BF()`, `marginal_inference()`,
    and `as_marginal_inference()` is the exact Gaussian kernel sum at the null
    (bandwidth `bw.nrd0()` of the continuous draws) instead of a 512-point
    grid estimate, which overestimated the posterior density at the null of
    long-tailed posteriors up to several-fold and so understated their Bayes
    factors; an ordinate below the double range gives an infinite Bayes
    factor. The kernel density is reflected at the exact support bounds when
    the marginal posterior carries exact support metadata, so Bayes factors of
    boundary nulls for bounded parameters differ from the standard kernel
    estimate of 0.3.0. When the continuous components of a marginal posterior
    have different exact supports (e.g., a prior truncated in only some models
    or mixture components), the ordinate is estimated per component, each on
    its own support and mixed by the components' shares of the continuous
    draws, instead of one kernel density smoothed across the support boundary
    (in prior-only checks, where the true Bayes factor at a boundary null is
    1, the estimates are 0.97 for a model mixture and 1.00 for a single-fit
    mixture). Marginal posteriors record each draw's component (the model
    for `mix_posteriors()` ensembles, the indicator tuple of the mixture or
    spike-and-slab terms entering the parameter or level for
    `as_mixed_posteriors()`) and the components' supports.
  - `marginal_estimates_table()` reports the Bayes factor warnings of scalar
    (non-formula) parameters, which it dropped.
  - metadata are read by their exact names: `Savage_Dickey_BF()` on a
    posterior with a prior-density context but no prior density reports the
    missing density instead of misreading the context as a prior.
  - inclusion Bayes factors keep a log-space value (`log_BF` and
    `inclusion_log_BF` attributes) for log-scale output, and the documentation
    of `compute_inference(conditional = TRUE)` states that the inclusion Bayes
    factor stays the unconditional inclusion odds.
- formulas, scaling, and original-scale quantities:
  - `ensemble_estimates_table(transform_scaled = TRUE)` and
    `marginal_estimates_table(transform_scaled = TRUE)` transform mixed
    columns through the fitted design by their fitted coordinates. For
    `~ g + g:x` with standardized `x`, the intercept and the levels of `g`
    were left on the standardized scale, also in `mix_posteriors()` mixtures;
    columns that do not identify their fitted coordinates, coefficients
    outside the design of the passed `formula_scale`, and random-effect SDs
    whose block the passed `formula_scale` does not describe now stop.
  - the original-scale transformation of fixed coefficients is derived from
    the fitted formula design and verified to reproduce the linear predictor
    exactly, in `transform_scale_samples()`, `transform_prior_samples()`,
    `JAGS_estimates_table(transform_scaled = TRUE)`,
    `as_mixed_posteriors(transform_scaled = TRUE)`, and marginal and hypothesis
    prior draws. Nested slopes such as `~ f/x` are transformed correctly, and
    formulas whose centered terms have no original-scale representation (e.g.,
    `~ x + x:f` with standardized `x`) stop with an informative error instead
    of returning wrong coefficients. `JAGS_formula()` and `JAGS_fit()` store
    the fitted design with the formula-scale metadata.
  - `formula_scale` centers scaled predictors also in terms without a free
    intercept (documented): `~ 0 + x` fits a line through the predictor mean,
    with original-scale intercept `-b * mean(x) / sd(x)`, and independent
    random slopes of scaled predictors, `(1 + x || g)`, are independent on the
    centered scale and imply a non-zero original-scale intercept-slope
    correlation.
  - factor point priors (e.g., `prior_factor("spike", ...)`) no longer add a
    spurious unindexed column to `transform_scale_samples()` and
    `transform_prior_samples()` output or a duplicate row to
    `JAGS_estimates_table(transform_scaled = TRUE, remove_spike_0 = FALSE)`.
  - `log(intercept)` formulas: `marginal_posterior()` levels use
    log(intercept) for draws, prior densities, and point masses, matching
    `JAGS_evaluate_formula()` (the flag comes from the fitted formula design,
    and an explicit `"log(intercept)"` formula attribute must agree with it);
    formula marginal posteriors of their linear predictors declare their exact
    support, mixture components, and posterior atoms (derived through the log
    of the intercept), also for `transform_scaled` samples;
    `marginal_posterior()` of the unscaled intercept of `transform_scaled`
    samples (`use_formula = FALSE`) has the prior density, support (0, Inf),
    and component supports of the exp of its log-scale combination; and
    `as_mixed_posteriors(conditional = , transform_scaled = TRUE)` works and
    conditions the intercept prior density on the event.
  - `marginal_posterior()` accepts intercept-only formulas (`formula = ~ 1`)
    and returns the intercept level, and `marginal_posterior()` with a formula,
    `marginal_inference()`, and `as_marginal_inference()` accept predictors
    that enter the formula only through interactions, such as `x` in
    `~ g + g:x`, using the type, levels, and contrast recorded on those
    interaction terms.
  - coefficients of a factor interaction whose other component has no main
    effect, such as `g:x` in `~ g + g:x` (one slope per level of `g`) or `g:h`
    in `~ g + g:h`, are named by their level cells from the term design in
    `as_mixed_posteriors()`, `mix_posteriors()`, summary tables, and
    hypothesis tests on these levels.
  - formulas whose parameters and terms generate the same JAGS node name (for
    example formula `a` with term `b_x` and formula `a_b` with term `x`) stop
    with a message naming them instead of a duplicate `prior_list` name, and
    fixed formulas with dot expansion, `offset()`, inline transformations, or
    other calls stop with a clear message before the data are looked up.
  - `expression()` formula terms are part of the formula parameter that
    `JAGS_bridgesampling()` reconstructs for `log_posterior`; 0.3.0 left them
    out (for `~ x + expression(0.25 * z[i])` it passed the intercept plus the
    `x` term only), so marginal likelihoods of such models were wrong. The
    terms may reference sampled scalar or one-dimensional indexed parameters
    (also individually monitored indexed coordinates); their parsed syntax and
    data and parameter dependencies are stored and replayed draw by draw in
    prediction and marginal-likelihood reconstruction, and `JAGS_fit()`
    rejects dependencies that cannot be replayed (opaque deterministic nodes
    and formula outputs). Scaled formula prediction keeps the raw expression
    inputs and honours explicit no-intercept prediction subsets without
    changing the fitted contrasts.
  - `JAGS_evaluate_formula()` and `marginal_posterior()` evaluate terms with a
    `multiply_by` multiplier with the arithmetic of the JAGS model (the
    multiplier times the coefficient, times the predictor), so values can
    change in the last bit.
  - package-defined factor contrasts are resolved inside the namespace across
    fixed, prediction, and marginal-posterior design matrices, so
    `BayesTools::` calls do not require attaching the package.
- marginal posteriors, model averaging, and posterior atoms:
  - `marginal_posterior()` works for model averages in which a model omits a
    factor term, no longer fails with `prior_samples = FALSE` when optional
    support or point-mass metadata are unavailable, reports no exact support
    for terms whose support is unknown (bias, p-hacking, and derived point
    priors), and stops with a clear message for weight-function, bias, and
    p-hacking posteriors.
  - affine transformations are carried into the joint marginal-prior weights
    and offsets, and unscaled posterior atoms are rebuilt from the joint
    coefficient structure; nonlinear transformed level combinations fail
    explicitly when the joint prior is unavailable instead of using the
    untransformed prior.
  - `as_mixed_posteriors(transform_scaled = TRUE)` rescales the point masses
    of treatment and independent factor coefficients. Under an `OR` condition
    with several labels (e.g., `conditional = c("mu", "omega")`),
    `as_mixed_posteriors()` keeps the bias columns of every branch present in
    the conditioned draws, so `plot_posterior(., "PETPEESE")` no longer stops,
    and a list of bias priors counts as one condition-label branch each.
  - point masses and density-grid values are merged by exact location instead
    of their 15-digit labels, so nearby distinct atoms keep separate masses.
- summary tables:
  - estimates tables blank the MCMC diagnostics of every publication-weight bin
    that its prior declares constant (the reference bin, every fixed weight,
    mirrored bins, and constant bins of bias mixtures), identified from the
    prior instead of rows whose label starts with `omega[0,`; the p-value
    interval rows (`omega[0,0.05]`) are catalog selectors of their bins.
  - estimates tables: `conditional = TRUE` supports factor mixtures with named
    components (rows `par[component][level]`); with `transform_scaled = TRUE`,
    removal of spike-at-zero coefficients is decided on the transformed
    values; requested transformations apply to independent factor
    coefficients; point priors with expression locations no longer fail; and
    log-scale inclusion Bayes factors are formatted from log-space values and
    stay finite beyond the double range. `interpret_records()` labels fallback
    intervals with the probabilities of the columns it uses.
  - ensemble tables print `n_models` as one denominator per row, and
    BayesTools tables keep their metadata and safe row names when subsetted by
    row.
  - Stan draws are extracted per parameter element, so
    `stan_estimates_table()` handles matrix-valued parameters.
- plots:
  - `plot_posterior()` bias prior overlays (full PET-PEESE and weight
    function, and individual PET, PEESE, and omega) follow the samples' full
    condition, including conditions on the effect and `OR` combinations with
    other parameters: each bias branch is weighted by its probability given the
    condition event, weight-function and no-bias branches count as PET = PEESE
    = 0, and the PET-PEESE overlay keeps the weights of mixture or
    spike-and-slab effect priors.
  - `plot()`, `geom_prior()`, and `plot_posterior()` prior overlays weight
    spike-and-slab priors by the inclusion probability instead of 50/50.
  - individual omega posterior plots take point masses and draws from the
    correct models when an ensemble contains duplicate or zero-weight bias
    priors or a conditioned bias mixture, and weight-function posterior plots
    take their bins from the mixed samples' columns, so a model with prior
    weight 0 no longer shifts or drops bins.
  - factor posterior and prior plots label an interaction cell by the levels
    of all its factors (a 3 x 3 treatment interaction showed 2 labels for 4
    curves), factor posterior plots use the declared posterior point masses of
    each level, PET-PEESE prior plots no longer apply the effect-size
    transformation to the standard-error axis when
    `transformation_settings = TRUE`, and the overall diamond of
    `plot_models()` uses the requested factor column.
  - mixed density plots resolve one density-to-probability mapping from
    `ylim` and `ylim2` for the secondary axis and all point masses, reserve
    room for the secondary-axis labels, and take their ranges from the joint
    prior and posterior densities and point masses; later overlays (`lines()`,
    `geom_prior()` with the default `scale_y2 = NULL`) reuse that mapping, and
    an overlay whose point mass falls outside the displayed axis warns instead
    of silently rescaling the plot. `geom_prior()` draws point-mass arrows of
    mixed overlays also without a secondary-axis scale, and weight-function
    step coordinates are shared by the plot data and the renderers.
  - transformed prior densities are evaluated on the displayed plotting range
    (including valid zero boundaries) and within the transformed support,
    while wider display limits and boundary point masses stay available;
    marginal-prior plots no longer exponentiate user-facing axis limits as
    fitted-scale coordinates.
  - prior curves of combinations without an exact route are omitted, with a
    warning of class `BayesTools_prior_curve_unavailable`, when their numerical
    grid cannot resolve a heavy-tailed product term (e.g., a t term plus a
    `multiply_by` product of a Cauchy and a beta term), instead of drawing a
    curve that is far off; the rest of the plot is drawn.
- fitting, convergence, and bridge sampling:
  - `JAGS_extend()` recompiles the model from the stored chain states instead
    of continuing a compiled model left in the session, so extending the same
    object twice, or a saved and reloaded copy, gives identical draws. It
    keeps the fit's warnings through successful extensions, and `JAGS_fit()`
    keeps (and shows during visible fitting) the messages of failed backend
    runs when it succeeds after a restart.
  - `JAGS_bridgesampling()` stops on rank-deficient draws of the bridge
    coordinates (e.g., a deterministic row-shaped SD source given through
    `add_parameters`; describe such sources with `parameter_source(values = )`)
    instead of returning a meaningless estimate, reports the original fitting
    error for a failed `JAGS_fit()` result, resolves the JAGS scalar names of
    complete singleton coordinate arrays, and stops explicitly on missing
    bridge coordinates and undefined prior densities.
  - `JAGS_marglik_parameters()` stops when the free coordinates of an
    independent weight function are missing from `samples`, instead of
    returning `NA` weights.
- selection and bias priors:
  - `selection_backend_spec(names = )` uses the custom node names of
    single-branch specifications in the generated syntax and initial values as
    in the monitors, and validates them as distinct JAGS node names. JAGS
    initial values of weight functions are drawn independently for the free
    omega bins, on their stochastic nodes (`omega_local*` / `log_omega*`)
    rather than the deterministic weight array.
  - p-hacking null masses are computed from tail probabilities, so far-tail
    cut points no longer give negative masses (non-representable masses are
    rejected), and mixed p-hacking source and destination geometry is reported
    as unsupported.

### Performance
- `JAGS_fit()` resolves the JAGS executable once per session instead of twice
  per fit (runjags searched for it with a system command on every call); a
  path set with `runjags::runjags.options(jagspath = )` is used as given. A
  small fit takes about 0.8 s less.
- prior-density curves of plots and the display grids of prior densities
  (e.g., the prior density of `marginal_posterior()`) evaluate the exact
  routes of Gaussian convolutions, `multiply_by` products, ordered-factor
  levels, and sums of two non-normal terms in one batched quadrature
  (relative error 1e-8) on at most 200 equally spaced values plus their
  support bounds, jumps, peaks, and atoms; shares or multipliers with a strong
  singularity at a bound (shape below 0.1) use per-value quadrature. Row-wise
  prior densities of formula `marginal_posterior()` levels build their
  numerical grid only when a plot or a grid-based summary needs it, while
  ordinates, heights, and region probabilities use the exact route.
- `JAGS_bridgesampling()` compiles the linear-predictor plans of formulas and
  the prior evaluators (support and normalization of simple priors, one
  vectorized evaluation of independent factor, PET, and PEESE priors) once per
  call, and skips formula reconstruction for models without formulas.
- two-bin cumulative weight-function priors are sampled through their exact
  Beta marginal, which removes the non-identifiable auxiliary Gamma scale and
  one likelihood-updating JAGS coordinate (also in bridge sampling); seeded
  fits with such priors differ from 0.3.0.
- saved fits are smaller: `JAGS_fit()` and `JAGS_extend()` return fits
  without runjags' compiled rjags model (`method.options$rjags`, closures that
  hold another copy of the model and its data: about 40 KB with 20
  observations and 360 KB with 20,000). It is continued only within the
  fitting call and is not alive after saving; `JAGS_extend()` and
  `runjags::extend.jags()` recompile the model from the stored chain states,
  as they do for a reloaded fit.
- `as_mixed_posteriors(transform_scaled = TRUE)` builds the transformed prior
  densities only for the requested `parameters` and keeps them for the fitted
  object (a bounded per-session cache that is dropped when the parameter map
  of the fit is replaced and never saved with the fit); the formula
  coefficient transforms of one scaling are built once per distinct input.
  Repeated calls with the same fit, priors, conditioning, and
  `n_prior_samples` no longer recompute them (ten scaled parameters: 0.05 s
  per call instead of 0.4 s). The `"prior_densities"` draw metadata hold the
  densities of the requested parameters only.

# version 0.3.0
### Features
- major refactoring and speed-up of unit tests
- adds support for `__default_factor` and `__default_continuous` priors in `JAGS_formula()` - when specified in the `prior_list`, these are used as default priors for factor and continuous predictors that are not explicitly specified
- adds automatic standardization of continuous predictors via `formula_scale` parameter in `JAGS_formula()` and `JAGS_fit()` - improves MCMC sampling efficiency and numerical stability
- adds `transform_scale_samples()` function to transform posterior samples back to original scale after standardization
- adds `transform_prior_samples()` function to generate and transform prior samples using the same matrix transformation as posterior samples - enables correct visualization of priors on the original (unscaled) predictor scale, including proper handling of the intercept which depends on multiple coefficient priors
- adds `transform_scaled` argument to `plot_posterior()` for visualizing prior and posterior distributions on the original (unscaled) scale when using formula-based models with auto-scaling
- adds `exp_lin` transformation type for log-intercept unscaling in density/plotting functions: `exp(a + b * log(x))`
- adds `log(intercept)` formula attribute for specifying models of the form `log(intercept) + sum(beta_i * x_i)` - useful for parameters that must be positive (e.g., standard deviation) while keeping the intercept on the original scale. Set via `attr(formula, "log(intercept)") <- TRUE`. Supported in `JAGS_formula()`, `JAGS_evaluate_formula()`, and marginal likelihood computation
- adds advanced parameter filtering options to `runjags_estimates_table()`:
  - `remove_parameters = TRUE` to remove all non-formula parameters
  - `remove_formulas` to remove all parameters from specific formulas
  - `keep_parameters` to keep only specified parameters
  - `keep_formulas` to keep only parameters from specified formulas
  - when `bias` is specified in `remove_parameters` or `keep_parameters`, the corresponding bias-related parameters (`PET`, `PEESE`, `omega`, `alpha`, `pi_null`, and `phack_kind`) are automatically included based on the bias prior type
- adds `probs` argument to `runjags_estimates_table()` and `runjags_estimates_empty_table()` for custom quantiles (default: `c(0.025, 0.5, 0.975)`)
- adds `effect_direction` argument to `plot_posterior()`, `plot_prior_list()`, `lines_prior_list()`, and `geom_prior_list()` for PET-PEESE regression plots - use `"positive"` (default) for `mu + PET*se + PEESE*se^2` or `"negative"` for `mu - PET*se - PEESE*se^2`
- redesigns `prior_weightfunction()` around a unified `side`, `steps`, and `weights` specification, with `wf_cumulative()`, `wf_fixed()`, and `wf_independent()` constructors for cumulative Dirichlet, fixed, independent, and log-independent weightfunction priors
- adds p-hacking and composed selection-bias priors via `prior_phacking()`, `prior_bias()`, calibration helpers, and `selection_backend_spec()` for compiling active step/p-hacking backend parameters
- adds error % for inclusion BF calculation

### Changes
- changes quantile column names in `runjags_estimates_table()` and `stan_estimates_table()` from `lCI`/`Median`/`uCI` to numeric values (e.g., `0.025`/`0.5`/`0.975`) for consistency with ensemble summary tables
- implied prior distributions for estimated marginal means, unstandardized coefficients, and PET-PEESE no longer require prior samples 
- implied prior distributions for weightfunction weights now use analytical forms for cumulative Dirichlet, fixed, independent, and log-independent priors, including mixture and model-averaged weightfunctions where possible
- independent weightfunction priors now allow non-reference weights above one via non-negative omega-scale priors or unrestricted log-omega priors
- replaces the legacy dot-named weightfunction prior specifications with the unified weightfunction prior API and updates JAGS generation, marginal likelihood computation, posterior extraction, diagnostics, and summary tables to use the new component-local `omega` representation
- composed selection-bias priors and publication-bias mixtures now support prior sampling and explicit unsupported-operation errors for ambiguous scalar prior generics

### Fixes
- reports inclusion Bayes factors as `NA` when the prior assigns probability 0 or 1 to inclusion, while keeping finite-sample bounds for posterior inclusion probabilities of 0 or 1
- fixes incorrect ordering the printed mixture priors
- fixes formula with no intercepts coded as `0` (instead of only `-1`)
- fixes bug in `.is.wholenumber` with NAs and `na.rm = TRUE`
- fixes ggplot prior spike layers for marginal factor plots with density and point components

# version 0.2.23
### Fixes
- `JAGS_diagnostics` functions now correctly handle factor parameters nested within mixture priors

# version 0.2.22
### Fixes
- `plot_posterior()` function with spike and slab priors 

### Changes
- unifies back-end of `prior_mixture()` and `prior_spike_and_slab()` 

# version 0.2.21
### Fixes
- `JAGS_formula()` function now replaces removed missing intercept with 0 (so the model matrix remains unchanged)
-  resetting `silent = FALSE` argument in the `JAGS_fit()` function now fits the model non-silently again 

# version 0.2.20
### Features
- extending prior functions to accept `expression()` instead of a parameter, such objects can be use to create prior distributions that depend on other parameters in JAGS
- extending the formula interface of `JAGS_fit()` function to accept expressions that are appended as literal text to the generated JAGS formula 
- extending the formula interface of `JAGS_fit()` function to handle uncorrelated random effects via `(x||y)` (lme4-like) notation

### Fixes
- `JAGS_estimates_table` not printing formula prefix when only spike and slab priors are supplied 

# version 0.2.19
### Features
- adds `max_extend` option to `autofit_control` argument in `JAGS_fit()` to limit the number of iterations for the model extension
- adds JASP progress bar integration

### Fixes
- `JAGS_diagnostics_density()` plots for mixture distributions
- prior and posterior `plot_posterior()` for simple `as_mixed_posteriors` objects
- `JAGS_evaluate_formula()` for mixture and spike and slab priors
- set Bayes factors based on alternative only prior distributions to NA
- better handling of posterior samples in `.fit_to_posterior()`

# version 0.2.18
### Features
- adding `prior_mixture()` function for creating a mixture of prior distributions
- adding `as_mixed_posteriors()` and `as_marginal_inference()` functions for a single JAGS models (with spike and slab or mixture priors) to enabling tables and figures based on the corresponding output
- adding `interpret2()` function for another way of creating textual summaries without the need of inference and samples objects
- speedup and improvements to the `runjags_estimates_table()` function

### Fixes
- small fixes for expansion of the RoBMA functionality

## version 0.2.17
### Features
- adding informed prior distributions for dichotomous and time to event outcomes based on Cochrane Database of Systematic Reviews to `prior_informed()` function
- adding bridge object convenience function `bridge_object()` (fixes: https://github.com/FBartos/BayesTools/issues/28)
- adding `Na/NaN` tests for `check_` functions (fixes: https://github.com/FBartos/BayesTools/issues/26)

### Fixes
- ability to run more than 4 chains (fixes: https://github.com/FBartos/BayesTools/issues/20)

## version 0.2.16
### Features
- update an existing JAGS fit with `JAGS_extend()` function
- new element of the `autofit_control` argument in `JAGS_fit()`: `"restarts"` allows to restart model initialization up to `restarts` times in case of failure

## version 0.2.15
### Fixes
- fixing repeated print of previous prior distribution in `model_summary_table()` in case of `prior_none()`

## version 0.2.14
### Features
- adding `contrast = "meandif"` to the `prior_factor` function which generates identical prior distributions for difference between the grand mean and each factor level
- adding `contrast = "independent"` to the `prior_factor` function which generates independent identical prior distributions for each factor level
- `remove_column` function for removing columns from `BayesTools_table` objects without breaking the attributes etc...
- adding empty table functions (https://github.com/FBartos/BayesTools/issues/10)
- adding `remove_parameters` argument to `model_summary_table()`
- adding multivariate point distribution functions
- adding `point` prior distribution as option to `prior_factor` with `"meandif"` and `"orthonormal"` contrasts
- adding `marginal_posterior()` function which creates marginal prior and posterior distributions (according to a model formula specification)
- adding `Savage_Dickey_BF()` function to compute density ratio Bayes factors based on `marginal_posterior` objects
- adding `marginal_inference()` function to combine information from `marginal_posterior()` and `Savage_Dickey_BF()`
- adding `marginal_estimates_table()` function to summarize `marginal_inference()` objects
- adding `plot_marginal()` function to visualize `marginal_inference()` objects

### Changes
- `contrast = "meandif"` is now the default setting for `prior_factor` function
- depreciating `transform_orthonormal` argument in favor of more general `transform_factors` argument 
- switching `dummy` contrast/factor attributes to `treatment` for consistency (https://github.com/FBartos/BayesTools/issues/23)

### Fixes
- zero length inputs to `check_bool()`, `check_char()`, `check_real()`, `check_int()`, and `check_list()` do not throw error if `allow_NULL = TRUE`
- properly aggregating identical priors in the plotting function (previously overlying multiple spikes on top of each other when attributes did not match)
- `student-t` allowed as a prior distribution `name`
- fixing factor contrast settings in `JAGS_evaluate_formula`
- fixing spike prior transformations

## version 0.2.13
### Features
- `runjags_estimates_table()` function can now handle factor transformations 
- `plot_posterior` function can now handle factor transformations 
- ability to remove parameters from the `runjags_estimates_table()` function via the `remove_parameters` argument

### Fixes
- inability to deal with constant intercept in marglik formula calculation
- `runjags_estimates_table()` function can now remove factor spike prior distributions
- marginal likelihood calculation for factor prior distributions with spike 
- mixing samples from vector priors of length 1
- same prior distributions not always combined together properly when part of them was generated via the formula interface

## version 0.2.12
### Features
- `stan_estimates_summary()` function
- reducing dependency on runjags/rjags

### Fixes
- dealing with posterior samples from rstan
- dealing with vector posterior samples
- fixing MCMC error of SD calculation for transformed samples (previously reported 100 times lower)

## version 0.2.11
### Features
- adding Bernoulli prior distribution
- adding spike and slab type of prior distributions (without marginal likelihood computations/model-averaging capabilities)
- new vignette comparing Bayes factor computation via marginal likelihood and spike and slab priors

### Fixes
- when a transformation is applied, JAGS summary tables now produce the mean of the transformed variable (previous versions incorrectly returned transformation of the mean) 

### Changes
- runjags_XXX_table functions are now also exported as JAGS_XXX_functions for consistency with the rest of the code

## version 0.2.10
### Features
- trace, density, and autocorrelation diagnostic plots for JAGS models

## version 0.2.9
### Fixes
- dealing with NaNs in inclusion Bayes factors due to overflow with very large marginal likelihoods

## version 0.2.8
### Fixes
- dealing with point prior distributions in `JAGS_marglik_parameters_formula` function
- posterior samples dropping name in `runjags_estimates_table` function
- `ensemble_summary_table` and `ensemble_diagnostics_table` function can create table without model components

## version 0.2.7
### Features
- `JAGS_evaluate_formula` for evaluating formulas based on data and posterior samples (for creating predictions etc)  
- `JAGS_parameter_names` for transforming formula names into the JAGS syntax

## version 0.2.6
### Features
- `plot_models` implementation for factor predictors
- `format_parameter_names` for cleaning parameter names from JAGS
- `mean`, `sd`, and `var` functions now return the corresponding values for differences from the mean for the orthonormal prior distributions

### Fixes
- proper splitting of transformed posterior samples based on orthonormal contrasts in `runjags_summary_table` function (previous version crashed under other than default `fit_JAGS` settings)
- always showing name of the comparison group for treatment contrasts  in `runjags_summary_table` function
- better handling of transformed parameter names in `plot_models` function

## version 0.2.5
### Features
- `add_column` function for extending `BayesTools_table` objects without breaking the attributes etc...
- ability to suppress the formula parameter prefix in `BayesTools_table` functions with with `formula_prefix` argument

### Fixes
- allowing to pass point prior distributions for factor type predictors

## version 0.2.4
### Features
- adding possibility to multiply a (formula) prior parameter by another term (via `multiply_by` attribute passed with the prior)
- t-test example vignette

## version 0.2.3
### Fixes
- fixing error from trying to rename formula parameters in BayesTools tables when multiple parameters were nested within a component

## version 0.2.2
### Fixes
- fixing layering of prior and posterior plots in `plot_posterior` (posterior is now plotted over the prior)

## version 0.2.1
### Fixes
- fixing JAGS code for multivariate-t prior distribution

## version 0.2.0
### Changes
- ensemble inference, summary, and plot functions now extract the prior list from attribute of the fit objects (previously, the prior_list needed to be passed for each model within the model_list as the priors argument

### Features
- adding formula interface for fitting and computing marginal likelihood of JAGS models
- adding factor prior distributions (with treatment and orthonormal contrasts)

## version 0.1.4
### Fixes
- fixing DOIs in the references file
- adds marglik argument `inclusion_BF` to deal with over/underflow (Issue #9)
- better passing of BF names through the `ensemble_inference_table()` (Issue #11)

### Features
- adding logBF and BF01 options to `ensemble_summary_table` (Issue #7)

## version 0.1.3
### Features
- `prior_informed` function for creating informed prior distributions based on the past psychological and medical research

## version 0.1.2
### Fixes
- `prior.plot` can't plot "spike" with `plot_type == "ggplot"` (Issue #6)
- `MCMC error/SD` print names in BayesTools tables (Issue #8)
- `JAGS_bridgesampling_posterior` unable to add a parameter via `add_parameters`

### Features
- `interpret` function for creating textual summaries based on inference and samples objects

## version 0.1.1
### Fixes
- `plot_posterior` fails with only mu & PET samples (Issue #5)
- ordering by "probabilities" does not work in 'plot_models' (Issue #3)
- BF goes to NaN when only a single model is present in 'models_inference' (Issue #2)
- summary tables unit tests unable to deal with numerical precision
- problems with aggregating samples across multiple spikes in `plot_posterior'

### Features
- allow density.prior with range lower == upper  (Issue #4)
- moving rstan towards suggested packages

## version 0.1.0
- published on CRAN

## version 0.0.0.9010
- plotting functions for models

## version 0.0.0.9009
- plotting functions for posterior samples

## version 0.0.0.9008
- plotting functions for mixture of priors

## version 0.0.0.9007
- improvements to prior plotting functions

## version 0.0.0.9006
- ensemble and model summary tables functions

## version 0.0.0.9005
- posterior mixing functions

## version 0.0.0.9004
- model-averaging functions

## version 0.0.0.9003
- JAGS fitting related functions

## version 0.0.0.9002
- JAGS bridgesampling related functions

## version 0.0.0.9001
- JAGS model building related functions

## version 0.0.0.9000
- priors and related methods
