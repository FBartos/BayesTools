# Fit-level hypothesis Bayes factors

Contract of `hypothesis_BF()` for a fitted model, the engine behind it
(`hypothesis_BF_model()`), the provider that describes a model to the engine
(`hypothesis_provider()`), and the registry of density backends. Draw-based
`hypothesis_BF.default()`, the parser and the Bayes-factor algebra are in
[priors.md](priors.md) (Savage-Dickey Ordinates); formula and parameter-map
metadata are in [jags-formula.md](jags-formula.md).

## Files

- `R/hypothesis-BF.R`: `hypothesis_BF()` generic; `.default` (draws, marginal
  posteriors, marginal inference; unchanged signature, unknown arguments stop)
  and `.BayesTools_fit` (calls the engine; unused `...` warn).
- `R/hypothesis-model.R`: entry `hypothesis_BF_model()`, group execution,
  binding, row names, `hypothesis_cache_clear()`.
- `R/hypothesis-model-plan.R`: refusals, plan cache, plan skeleton, scalar,
  linear (levels, combinations), exp-affine routes, formula route
  (`JAGS_formula_coefficient_route()`), linear-weight checks
  (`posterior_formula_measure_check()`, `posterior_linear_weight_space()`,
  `posterior_fitted_linear_weights()`), `hypothesis_stop()`.
- `R/hypothesis-model-random.R`: random-effect routes (SD, variance through its
  SD, variance point against region, prior-draw region route).
- `R/hypothesis-model-marginal.R`: marginal-means plans (provider kind
  `"marginal_means"`).
- `R/hypothesis-model-select.R`: parsing against the provider catalog, parameter
  selection, component vocabulary, rewrite.
- `R/hypothesis-model-ordinates.R`: stored-ordinate bookkeeping, attachment
  through a backend, warning and diagnostics lifting to table rows.
- `R/hypothesis-model-provider.R`: `hypothesis_provider()` generic, the default
  provider of a `BayesTools_fit`, provider validation.
- `R/hypothesis-model-backends.R`: backend constructor and registry.
- `R/hypothesis-model-quantities.R`: `hypothesis_quantities()`.

## Stages of `hypothesis_BF_model()`

1. Arguments are checked; the provider is validated; `density_method` is
   matched (case-insensitive, partial) over `provider$density_methods` or the
   registered backends that are not of type `"stored"`; the backend's
   `control()` normalizes `density_control` (ordinate backends always, other
   backends when a control is supplied).
2. Every statement is parsed with the provider catalog (`simplify_names =
   TRUE`) and planned (resolve, select the entry, rewrite onto the tested
   parameter, route, targets with `prior_ordinate_status()`). The draws of a
   statement are extracted while it is planned, before its point targets
   classify their prior ordinates, and without the public `seed`.
3. Every plan is checked against the requested backend
   (`.hypothesis_model_plan_status()`): plan refusal; region-only statements
   pass; point-target refusals; `source_refusals`; the backend's
   `available(object, plan, control)`; then `method_refusals` (name-keyed, then
   type-keyed). Nothing is evaluated before all plans are checked. Availability
   is computed on every call and never cached.
4. `provider$check_plans(plans)` validates across statements. For marginal
   means the engine first refuses mixed point/region requests and point
   statements on several levels of a model-averaged provider.
5. Statements with the same plan `group` are evaluated together (scalar,
   levels, combination, random, random sd, random sd-region, marginal means).
   Point statements under an ordinate backend: stored ordinates at the requested
   values are dropped unless `reuse()` accepts them, `attach()` computes the
   missing ones, the engine checks them with `posterior_ordinate_supports_bf()`,
   the group is evaluated as `"precomputed"` by `hypothesis_BF.default()`, and
   warnings and diagnostics are lifted. Region statements use the draws.
6. Results are bound in statement order: labels as written, `label (i)` row
   names, `raw_log_BF`-based BF formatting, merged footnotes, row-keyed
   warnings, and the backend's `bind_diagnostics()` for `density_diagnostics`.

## Random numbers

Draws and prior densities are built during planning without the public `seed`
(the fit-level methods never reseed); `seed` reaches `hypothesis_BF.default()`
and the prior-draw route of random-effect region tests only. The plan cache
does not contain `seed`.

## Caching

One cache per model object (`rlang::hash(provider$object)`), bounded (4 objects,
256 entries, 256 MB each), holds structural plans, draws, and random-effect
selections; the key of a plan contains the component, level component, scale,
conditioning, `n_samples`, and the provider metadata. Backends own their caches
(`clear_cache()`); `hypothesis_cache_clear()` drops the plan caches and calls
every registered backend's `clear_cache()`. The cache is lost with the
session and when the package is reloaded.

## Provider (`hypothesis_provider()`)

`hypothesis_provider.BayesTools_fit()` is the default provider; a model class
starts from `hypothesis_provider(<its fit>)` and replaces callbacks with
`utils::modifyList()`. Fields (defaults are those of the fit provider):

| field | meaning |
|---|---|
| `kind` | `"model"`; `"marginal_means"` for stored marginal means |
| `object`, `fit`, `capability_object` | model object (hashed, passed to backends as `request$model`), its `BayesTools_fit`, object backends compute from (`request$object`; default `object`) |
| `build`, `refit_class`, `caller`, `ordinate_label`, `random_point_mass_advice`, `transform_scaled_default` | wording and class hooks: build named in refit requests, classes prefixed to `BayesTools_refit_required`, entry point named in reasons (`"hypothesis_BF()"`), label of the ordinate backends in refusal reasons, sentence ending the random point-mass reason, default scale of `hypothesis_quantities()` |
| `density_methods` | names `density_method` may match, in order (default: registered non-`stored` backends) |
| `components` | `values`, `synonyms`, `default`, `refused` (named reasons), `random`, `advice`; fit provider: `auto`, formula parameters / `model`, `random` |
| `metadata()` | `list(catalog, entries)`; entry columns `quantity_id, parameter, component, term, source, formula_parameter, role, quantity, status, fixed_value, aliases, member_quantity_ids`; alias rows of the extended catalog carry `namespace = component = entry component` |
| `uncovered(resolution, quantity_ids, symbol, metadata)` | refusal list (`reason` plus `type` or complete `class`) for resolved quantities without an entry |
| `sample_parameters(parameter)`, `mixed_posteriors(parameters, conditional, conditional_rule, transform_scaled, n_prior_samples)`, `density_parameter(samples, parameter, sample_parameters)`, `formula_name_map(formula_parameter)` | draws of a fixed target: prior-list names to mix, `as_mixed_posteriors()`, element tested, formula coefficient owners |
| `random` | `select`, `mixed_posterior`, `support`, `density_target`, `sd_pair`, `zero_boundary`, `defined_share`, `defined_footnote`, `defined_draws`, `prior_draws` (signatures in `?hypothesis_provider`) |
| `point_mass_reason()`, `full_ensemble(parameter)`, `combination_hint(plan)`, `check_plans(plans)` | optional hooks (default provider: no point-mass text, mixture rule, hint for factor terms with null components, none) |

Marginal-means providers add `inference`, `model_averaged` (logical),
`select(hypothesis, parameter)` (parameter and aliases of a statement),
`density_method_default()`, `density_control_default(density_method)`,
`source_missing_reason` (used when `capability_object` is `NULL`), and
`metadata()` returning `term_map`. Model-averaged point statements use
`inference$conditional`, others `inference$averaged`.

The default provider builds entries from `parameter_catalog(fit)`: fixed
coefficients (a factor term is one `formula_coefficient_group` entry whose
members are its level cells, with an extension quantity of provider
`"BayesToolsHypothesis"`), scalar `parameter` quantities that own a prior
(not indicators, inclusions, or vector-prior components), and supported
random-effect quantities (`sd`, `var`, totals, `cor`, `var_prop`, multipliers;
component `"random"`, canonical names as row labels). Contrast coefficients and
the other catalog quantities are refused as targets
(`BayesTools_hypothesis_statement` / `BayesTools_hypothesis_target`). Limits of
the default provider: non-formula vector and factor priors have no entries;
`random$zero_boundary()` marks allocation-derived SDs at 0 without naming an
omission alternative; `conditional = TRUE` for random-effect quantities is what
`parameter_mixed_posterior(conditional = TRUE)` supports (allocation gates);
quantities with point components (spike-and-slab SDs and correlations) stop with
its error.

## Scale and conditioning

`transform_scaled = FALSE` (fitted scale) is the default for coefficients and
random-effect SDs, variances, and correlations alike
(`parameter_mixed_posterior(transform_scaled = )` blanks the declared predictor
scales of the fit for draws and prior densities). `conditional = TRUE`
conditions coefficients through `as_mixed_posteriors(conditional = <parameter>,
conditional_rule = "AND")` and random-effect quantities through
`parameter_mixed_posterior(conditional = TRUE)`. Without it, a point null at a
prior point mass stops with `BayesTools_point_mass_at_null`, and a test of a
parameter whose prior has null and alternative components warns once per group
(`full_ensemble()`); pass `conditional = FALSE` explicitly to silence it.

## Density backends

`hypothesis_density_backend(name, type, available, control, attach, reuse,
warnings, diagnostics, bind_diagnostics, rescale_failure, refusal_classes,
clear_cache, method, policy, advertised, package)`; registry in the package
environment (`register_hypothesis_density_backend(backend, overwrite)`,
`unregister_hypothesis_density_backend(name)`,
`hypothesis_density_backends()`). Built-ins: `KDE` and `normal` (type `draws`,
evaluated by `hypothesis_BF.default()`; `normal` is not advertised) and
`precomputed` (type `stored`, draw inputs only). Names must not differ only by
case; built-ins are not replaceable. Types: `draws` (method `KDE`/`normal`),
`ordinate` (needs `attach`), `stored`.

Callbacks: `available(object, plan, control)` returns `NULL` or a refusal
(`reason`, optional `type`, additional `class`, `fields`); `control(control)`
normalizes and validates (default: no settings); `attach(posterior, request,
control, cache)` returns the posterior with one ordinate per requested value in
`posterior_ordinate` metadata (matching the draws' parameter and condition
metadata) or stops with the backend's own condition; `reuse(entry, request,
control)` (default `FALSE`); `warnings(posterior, request)` (records with
`warning, parameter, level, label, value`; default reads `diagnostics$warning`);
`diagnostics(posterior, request)` and `bind_diagnostics(list)`;
`rescale_failure(condition, requested, evaluated, label, density_method)` for
variances evaluated through SDs; `refusal_classes(cause)` for the generic causes
`original_scale`, `conditional_random`, `random_target`; `clear_cache()`.
`request`: `kind` (`scalar`/`level`), `route` (plan kind), `object`, `model`,
`parameter`, `label`, `level`, `values`, `conditional`, `spec` (`primitive`,
`linear` with `weights` and `prior_density`, or provider specs),
`display_transform`, `density_method`, `inference_type`. The engine owns
bookkeeping only; acceptance rules, thresholds, warnings and their text belong
to the backend.

## Conditions

Engine conditions carry BayesTools classes only: `BayesTools_hypothesis_statement`
(restate), `_ambiguous` (parent `_statement`), `_component_mismatch`,
`_parameter_mismatch` (marginal-means providers), `BayesTools_hypothesis_unavailable`
with `_fixed`, `_target`, `_method`, the classes of `prior_ordinate_status()`
(+ `BayesTools_hypothesis_ordinate`), and BayesTools resolution errors after
`BayesTools_hypothesis_statement`. Backend conditions stay the backend's.
`hypothesis_stop()` raises the same refusals for other entry points of a
package.

## Checks

`tests/testthat/test-hypothesis-BF-model.R` (unit, synthetic fits, toy backend,
toy providers) and `test-hypothesis-BF-model-fixture.R` (cached fits:
`fit_hypothesis_hcs`, `fit_hypothesis_car`, `fit_label_random_slope`,
`fit_re_summary_point_components`, `fit_dual_param_regression`,
`fit_formula_treatment`). The one-call result must equal the three-step
extraction exactly for the same seed.
