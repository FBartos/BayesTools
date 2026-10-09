# Post-implementation deep review for `0.3.1-add-random` -> `master`

Review date: 2026-08-06

- Reviewed code head: `41b7e59`
- Prior deep-review baseline: `8601948`
- Target and merge base: `e8533cd` (`master`)
- Reviewed implementation delta: 31 files, approximately 2,643 additions and
  216 deletions
- Standard: report only reproduced functional failures, violated persisted
  contracts, independently checkable statistical errors, or concrete
  maintenance pain points encountered while implementing and verifying the
  fixes

## Outcome

All approved PF01-PF04 decisions are implemented. The post-implementation
audit found four additional contract gaps inside the approved scope; each was
reproduced, fixed, regression-tested, and committed separately or folded into
the issue commit before this report. No unresolved inference, fitting,
prediction, bridge-reconstruction, or catalog-resolution defect was reproduced
on the final code head.

PF03 intentionally remains a no-code decision. Publicly constructed fits keep
their existing object-local schema and contract validation; BayesTools does not
add cross-object consistency checks for manually mutated fit metadata.

| Decision or follow-up | Resolution commit(s) |
| --- | --- |
| PF01: transformed finite support endpoints | `3444c32`, `94929d6` |
| PF02: draw-aware formula-expression replay | `97bd538`, `41b7e59` |
| PF03: metadata validation scope | `fc130d7` (decision record only) |
| PF04: named factor term quantities | `b44afca`, `1a39e6b` |
| Exact fixed random-scale bridge rows | `ce61d0c` |
| Named single-source coefficient densities | `2b5affe` |
| Semantic random-summary aliases | `cf37003`, `f31d899` |
| Final vignette cache refresh | `000ddb2` |

## Substantive gaps resolved during the final audit

### PI01 - A finite linear-transform endpoint remained one double outside support

**Issue.** PF01 initially recognized exact forward-transform values and
alternative power forms, but a finite linear endpoint could still differ from
the user's mathematically identical literal by one representable double. For a
truncated source with lower endpoint `0.1` and `Y = 0.1 + 2 * X`, the computed
endpoint is not `identical()` to the literal `0.3`.

**Impact.** `prior_density_ordinate(..., 0.3)` could report an outside-support
zero at a known finite endpoint, affecting structural Savage-Dickey decisions.

**Change.** Provenance-derived endpoints now recognize equality or one-step
representable adjacency using midpoint rounding. The enclosure applies only to
known transformed support endpoints; no global epsilon, clamping, or interior
density tolerance was introduced. Commit: `94929d6`.

### PI02 - Factor-prior wrappers were omitted from named factor quantities

**Issue.** `prior.factor_mixture` and `prior.factor_spike_and_slab` do not
satisfy `is.prior.factor()`. The first PF04 implementation therefore skipped
their persisted factor designs even though ordinary formula handling accepts
both wrappers as factor terms.

**Impact.** Named factor levels resolved for an unwrapped prior but disappeared
when the same factor prior was placed in a mixture or spike-and-slab wrapper.

**Change.** Catalog construction now uses the formula subsystem's complete
factor-prior predicate. Wrapper regressions cover derived orthonormal cells.
This correction is included in `b44afca`.

### PI03 - Interaction cell labels were not injective or always hypothesis-safe

**Issue.** Interaction components were assembled as `term=value` strings
separated by commas. Legal levels can contain the same delimiter text. For
example, the two distinct cells `(f = "a, g=v", g = "u")` and
`(f = "a", g = "v, g=u")` both rendered as `f=a, g=v, g=u`. Empty levels and
levels containing `]` also could not round-trip through level-reference syntax.

**Impact.** Distinct fitted interaction cells became ambiguous during exact
catalog or hypothesis resolution, contradicting PF04's every-cell contract.

**Change.** Catalog components preserve ordinary labels, quote delimiter-
bearing interaction tokens, and percent-escape syntax-sensitive characters.
The representation is injective and round-trips through `hypothesis_parse()`
and `hypothesis_resolve()`. Tests cover the collision, empty levels, closing
brackets, direct catalog resolution, and AST resolution. Commit: `1a39e6b`.

### PI04 - A valid individually monitored indexed expression could fit but not replay

**Issue.** Expression replay required every indexed parameter to contain
contiguous monitored coordinates starting at one. A valid model using
`expression(theta[2])` with only `theta[2]` in `add_parameters` fitted in JAGS,
but `JAGS_evaluate_formula()` rejected the posterior because `theta[1]` was not
monitored.

**Impact.** The PF02 contract depended on unrelated monitor choices. Fitted
prediction and bridge reconstruction could fail even though every coordinate
actually referenced by the expression was available.

**Change.** Indexed draws remain compact and retain their actual indices. One
draw is expanded only while evaluating the expression, so missing unreferenced
coordinates are allowed without allocating a draw-by-maximum-index matrix.
Referencing a genuinely unavailable coordinate still fails through the
non-finite replay check. Unit, bridge-style, and live JAGS regressions cover the
sparse case. Commit: `41b7e59`.

## Findings requiring maintainer approval

None. The final audit did not leave a reproduced correctness issue awaiting a
decision.

## Follow-up recommendations (not current defects)

### S01 - Centralize formula factor-prior capability classification

**Pain point.** Formula code currently repeats variants of
`is.prior.factor()` plus wrapper-class checks. The PF04 wrapper omission was a
direct consequence of one call site using the narrower predicate.

**Impact.** New wrapper types or factor capabilities can work in fitting while
silently disappearing in prediction, scaling, or semantic metadata.

**Suggested change.** Keep one formula-specific factor-capability predicate
near factor metadata construction and use it across formula validation,
prediction, scaling, coefficient transforms, and catalogs. Do not replace
distribution-specific dispatches where wrapper versus component behavior is
intentionally different.

### S02 - Persist an expression access plan, not only dependency roots

**Pain point.** Formula designs persist parsed expressions and dependency
classes, but replay still rediscovers scalar/indexed geometry by scanning each
posterior or bridge row. The sparse-coordinate failure arose at this boundary.

**Impact.** Shape errors are reported later than necessary, and every replay
path must independently interpret coordinate names.

**Suggested change.** Extend the next formula-design schema only when needed
with a compact access plan per parameter dependency: source kind, scalar versus
one-dimensional indexed access, and statically referenced indices where they
can be determined. Keep dynamic `parameter[index_data[i]]` access supported;
do not build a general JAGS interpreter.

### S03 - Reduce fit-cache invalidation and cascading fixture failures

**Pain point.** A test-only edit inside `test-00-model-fits.R` invalidated the
single model-fit completion marker. The fixture lane then emitted 44 downstream
failures for the same stale hash and required a 17-minute 44-second fit-profile
refresh before any fixture assertion could run meaningfully.

**Impact.** Small, isolated live regressions impose near-full cache refresh
cost and produce noisy cascades that obscure the single actionable cause.

**Suggested change.** Split live-fit cache provenance into smaller fixture
groups or hash each cached fit against only its constructor and relevant helper
sources. In the fixture lane, stop after one classed stale-cache precondition
failure instead of repeating it at every consumer.

### S04 - Add batched extraction for derived catalog quantities when consumers need it

**Pain point.** `parameter_draws()` supports an ordinary multi-coordinate
selection but explicitly rejects mixed or multiple derived selections. A
downstream hypothesis comparing several contrast-derived factor cells must
extract each cell separately and realign chains itself.

**Impact.** The current singular API is correct, but repeated extraction adds
consumer complexity precisely for multi-level hypotheses that motivated the
catalog work.

**Suggested change.** If RoBMA integration demonstrates the need, accept a
validated selection containing multiple BayesTools-owned derived quantities
and materialize one aligned column per quantity in selection order. Preserve
the current dependency declarations and chain timing; do not add hypothesis
evaluation logic to the catalog layer.

## Audit coverage

The final review traced and challenged:

- expression parsing, dependency ownership, expression-only data persistence,
  literal JAGS syntax, fixed and random prediction, posterior draw geometry,
  marginal-likelihood reconstruction, bridge evaluators, point parameters,
  full and sparse indexed parameters, and cross-formula rejection;
- factor catalog construction for treatment, independent, orthonormal,
  mean-difference, both ordered encodings, mixtures, spike-and-slab wrappers,
  structural cells, derived cells, interactions, delimiter collisions, and
  hypothesis AST resolution;
- transformed prior-density support and boundary behavior for linear,
  exponential, power, and hyperbolic-tangent transforms; and
- exact fixed random-scale bridge rows, coefficient-density source identity,
  and raw-versus-derived random-summary alias precedence.

Estimated marginal means remain deliberately outside coefficient-level catalog
quantities and continue through `marginal_means()` / `as_marginal_inference()`.
Literal expression indexing remains JAGS syntax: users are responsible for
writing scalar-valid JAGS expressions, while BayesTools guarantees replay for
the documented subset after a model fits.

## Verification

- Focused regression contexts for prior-density ordinates, formula-design
  oracles, and parameter catalogs passed after each corresponding change.
- Clean public live reproduction for `expression(theta[2])` with only
  `theta[2]` monitored passed fitting and fitted prediction.
- Unit profile: 8,953 passed, 0 failed, 0 warnings, 7 profile skips.
- Live-fit profile: 14,960 passed, 0 failed, 1 established warning from the
  deliberately data-free model.
- Fixture profile after the final fit-cache refresh: 10,781 passed, 0 failed,
  0 warnings, 1 profile skip.
- The initial fixture attempt after a test-source edit failed only its stale
  completion-marker precondition; the required fit refresh was performed and
  the fixture lane then passed.
- ComparisonR, SpikeAndSlab, and RandomEffects caches regenerated and validated
  with 4 objects, 5 objects, and 13 models, respectively.
- `R CMD check --no-manual --as-cran`: 0 errors, 0 warnings, 2 environment-only
  notes (CRAN incoming metadata and the worktree's `.git` directory). Package
  installation, compiled code, documentation, examples, tests, and vignette
  rebuilding passed.

The branch remains at development version `0.3.1.14`; `NEWS.md` describes the
final expression, factor-catalog, and transformed-endpoint behavior.
