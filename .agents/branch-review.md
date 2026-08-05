# Deep branch review for `0.3.1-add-random` -> `master`

Review date: 2026-08-05

- Reviewed head: `158930d` (`0.3.1-add-random`)
- Target: local and `origin/master` at `e8533cd`
- Merge base: `7d9d378`
- Scope: 626 changed files, approximately 141,687 additions and 32,778 deletions
- Standard: retain only merge/CI blockers and functional findings backed by a reproduced failure, a violated persisted contract, or an independent mathematical oracle

The previous R03-R11 findings were removed. Their topics (API positioning,
seeding documentation, release-note wording, minimum-R policy, future JAGS
versions, a live-fit warning, package-check notes, and whitespace) are not
included below.

## Merge and CI blockers

### R01 - The branch conflicts with current `master`

**Issue.** A read-only merge simulation reports conflicts in `.Rbuildignore`
and `.gitignore` because `master` changed both files after the branch's merge
base.

**Impact.** The branch cannot merge cleanly. A mechanical resolution could
exclude `tests/` from source packages or discard generated-artifact and cache
rules.

**Suggested change.** Rebase or merge current `master` and combine the intended
rules from both sides. Preserve the branch's inclusion of `tests/`; do not use a
blanket "ours" or "theirs" resolution.

**Decision:** combine rules from both sides

### R02 - All three committed vignette caches are stale

**Issue.** Direct validation reports `ComparisonR.RDS`, `SpikeAndSlab.RDS`, and
`RandomEffects.RDS` as stale against the current source and package versions.
All were produced with BayesTools 0.3.1.7; the branch is 0.3.1.12.

**Impact.** The mandatory cache checks in the vignette and pkgdown workflows
stop before rendering, so the branch is not CI-ready.

**Suggested change.** Regenerate the three caches only after all approved code,
documentation, and version changes are final, then validate and render every
vignette:

```text
Rscript tools/regenerate-precomputed-vignette.R ComparisonR
Rscript tools/regenerate-precomputed-vignette.R SpikeAndSlab
Rscript tools/regenerate-random-effects-vignette.R
```

**Decision:** update cache

## Functional findings

### F01 - Level-qualified hypotheses cannot be resolved to catalog quantities

**Location.** `R/hypothesis-ast.R:238-263` and
`R/JAGS-parameter-catalog.R:396-465`

**Issue.** The AST correctly parses `f[b]` into parameter `f` and level `b`, but
`hypothesis_resolve()` passes only `occurrences$parameter` to
`parameter_catalog_resolve()`. The level is discarded. In addition, the
BayesTools catalog rows for a treatment-coded factor store backend coordinates
such as `mu_f[1]` in `component`, rather than the semantic factor levels needed
to distinguish `f[b]` from `f[c]`.

**Reproduced example.** For `~ 1 + f` with fitted levels `a`, `b`, and `c`, the
AST occurrence is `(parameter = "f", level = "b", symbol = "f[b]")`. Resolving
`f[b] > 0` errors that alias `f` is ambiguous between `mu_f[1]` and `mu_f[2]`.
The same problem occurs when one expression refers to two levels; the public
scalar `component` argument cannot select a different component per occurrence.

**Impact.** The new AST/catalog integration cannot perform its central job for
factor or allocation-level hypotheses. A downstream consumer may have to parse
backend indices again, defeating the semantic catalog contract.

**Suggested change.** Persist the semantic level/component mapping from the
formula design and prior metadata in catalog quantities and aliases. Resolve a
level reference per occurrence using its complete symbol or its exact level as
the component; do not resolve only the root. Add tests for `f[b]`, `f[c]`, and a
single hypothesis containing both levels.

**Decision:** OK

### F02 - Valid JAGS `expression()` terms can fit but cannot be replayed

**Location.** `R/JAGS-formula-helpers.R:136-204`, with callers in
`R/JAGS-bridge-formula-evaluators.R:244` and
`R/JAGS-marglik-formula-parameters.R:293`

**Issue.** `JAGS_formula()` still accepts literal JAGS-scale additions, but
prediction and bridge replay parse the stored text and evaluate it in a base-R
environment containing only data and `i`. JAGS functions and sampled model
parameters are absent, and R is not a general evaluator for JAGS expressions.

**Reproduced example.** `~ expression(step(z[i]))` generates valid model syntax
containing `step(z[i])`. Bridge reconstruction then errors with `could not find
function "step"`. Likewise, a valid expression referring to a sampled scalar,
such as `expression(theta)`, cannot be reconstructed because the posterior
state is not in the evaluation environment.

**Impact.** A model accepted and fitted by `JAGS_fit()` can fail only when
`JAGS_bridgesampling()` or `JAGS_evaluate_formula()` is called. For functions
that exist in both languages but differ semantically, replay could evaluate a
different likelihood from the fitted JAGS target.

**Suggested change.** Define one enforceable replay contract. The smaller safe
change is to accept only data/`i` arithmetic and a documented whitelist whose
R and JAGS semantics are identical, rejecting every other expression during
`JAGS_formula()` construction and directing posterior-dependent terms to
`parameter_source()`. If arbitrary literal JAGS remains supported, compile a
validated expression representation and evaluate it with the reconstructed
posterior/formula state using JAGS-equivalent semantics. Add construction and
replay tests for JAGS-only functions and posterior-dependent expressions.

**Decision:** make expressions fully supported throughout the whole package (jags, marglik, predict, random effects)

### F03 - `JAGS_extend()` relabels stale metadata as the current schema

**Location.** `R/JAGS-fit.R:549-705`

**Issue.** `JAGS_extend()` saves the input fit contract but never validates its
component versions. At the end it preserves the old formula design and other
metadata, then calls `.bt_attach_fit_contract()`, which stamps every component
with the current version.

**Reproduced example.** A fit with both
`fit_contract$formula_design_version = 2L` and an actual formula design at
schema 2 was extended through the no-backend-work `max_time` exit. The returned
fit still contained schema-2 design metadata but claimed
`formula_design_version = 3L`.

**Impact.** Prediction, marginal-likelihood reconstruction, or downstream
consumers can interpret stale serialized metadata under the current schema.
This defeats the new fail-closed fitted-object compatibility contract and can
change reconstructed model semantics.

**Suggested change.** Validate every preserved contract component required by
`JAGS_extend()` before any extension. Reject unsupported metadata with a refit
instruction. Only attach a current contract after the corresponding current
metadata have been validated or rebuilt; never upgrade version numbers without
upgrading the objects. Add a regression test that a stale formula-design
contract fails before the backend is called.

**Decision:** OK

### F04 - Formula coefficient transforms reject ordinary mixed formulas

**Location.** `R/JAGS-formula-coefficient-density.R:209-242`

**Issue.** `.bt_formula_coefficient_sources()` derives its expected fixed
sources from every entry in `design$prior_list`. For a mixed formula that list
also contains generated random-effect SD/covariance priors. It then compares
that full set with registry rows restricted to `role == "fixed_coefficient"`,
so the sets necessarily disagree.

**Reproduced example.** For `~ 1 + random(1 | id, covariance = "diag")`, the
design prior list contains `mu_intercept` and
`mu__xREx__id_intercept`, while the fixed-coefficient registry contains only
`mu_intercept`. `JAGS_formula_coefficient_transform(fit, "mu")` errors with
`Formula coefficient sources ... disagree with the parameter registry`.

**Impact.** The branch's new coefficient-transform and induced-prior-density
APIs are unavailable for ordinary random-effect models whose covariance
parameters have priors, including the standard diagonal random-intercept case,
even though fixed-coefficient unscaling is well-defined and independent of the
random-effect priors.

**Suggested change.** Derive expected sources from authoritative fixed-term
design/name-map metadata or filter prior-owned coordinates by the registry's
`fixed_coefficient` role before checking identity and order. Do not infer fixed
sources from the entire design prior list. Add tests for sampled and
marginalized random-effect formulas, including a scaled fixed predictor.

**Decision:** OK

### F05 - The inverse-gamma negative-power boundary classification is reversed

**Location.** `R/prior-density-ordinate.R:1582-1634`

**Issue.** For an inverse-gamma source and negative power `b`,
`.prior_density_ordinate_exp_lin_boundary()` uses
`-(shape + b)` and applies the positive-exponent rule appropriate to a limit in
the source coordinate. Because `b < 0` reverses the source and output limits,
the zero/infinite decisions are reversed.

**Independent oracle and reproduction.** Let
`X ~ InvGamma(alpha, 1)` and `Y = exp(-2 * log(X)) = X^-2`. Near zero,

```text
f_Y(y) is proportional to y^((alpha - 2) / 2) * exp(-sqrt(y)).
```

Therefore the density is infinite for `alpha < 2`, finite and positive for
`alpha = 2`, and zero for `alpha > 2`. With an immaterial lower truncation at
0.01, `prior_density_ordinate()` reports the opposite: `zero` for shape 1 and
`infinite` for shape 3 (shape 2 is `regular`).

**Impact.** Structural prior behavior at an exact boundary is wrong. Consumers
can accept an invalid regular density ratio or reject a valid one, and the
reported zero/infinite provenance contradicts the actual transformed measure.

**Suggested change.** Classify the output-scale exponent
`-(shape + b) / b` for `b < 0`, or equivalently account explicitly for the
limit reversal. Add analytic tests for inverse-gamma shapes below, equal to,
and above `-b`, plus a numerical interior-density check away from zero.

**Decision:** OK

## Deep-audit coverage and verification

The review traced the changed implementations and their tests across:

- random-effect design generation, covariance construction, centered and
  noncentered reconstruction, new-level prediction, and bridge prior replay;
- formula scaling, expression replay, fitted design identity, registry,
  draw geometry, parameter catalog, and fully fixed models;
- prior-density provenance, transformations, mixtures, point masses,
  Savage-Dickey behavior, and formula coefficient densities; and
- hypothesis parsing, AST round trips, rewriting, catalog resolution, and
  model-averaging alignment.

The earlier branch decision archives (D01-D32 and NF01-NF18) were used only as
an audit map. Every archived item is marked implemented and verified; none was
carried into this report without a fresh current-tree failure.

Existing broad verification on this head:

- unit: 8,852 passed, 0 failed, 0 warnings, 7 profile skips;
- visual: 810 passed, 0 failed, 0 warnings, 11 profile skips;
- fresh live-fit profile: 14,951 passed, 0 failed, 1 data-free backend warning;
- fresh fixture: 10,781 passed, 0 failed, 0 warnings, 1 profile skip;
- fresh visual-fixture: 666 passed, 0 failed, 0 warnings, 26 profile skips; and
- package check without vignette rebuilding: 0 errors, 0 warnings, 1 note.

Targeted independent/reproduction checks were additionally run for F01-F05.
No source, test, fixture, or snapshot changes have been made; this file records
suggested changes for maintainer approval.
