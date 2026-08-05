# Post-fix deep review for `0.3.1-add-random` -> `master`

Review date: 2026-08-05

- Reviewed head: `8601948`
- Target: local and `origin/master` at `e8533cd`
- Merge base: `e8533cd` (current `master` is merged)
- Scope: 626 changed files, approximately 142,197 additions and 32,762 deletions
- Standard: retain only findings backed by a reproduced failure, a violated
  persisted contract, or a direct semantic contradiction

The approved findings in `.agents/branch-review.md` were implemented before
this review. Each issue has its own commit:

| Original issue | Resolution commit |
| --- | --- |
| R01: merge current `master` | `ab3597d` |
| F01: level-aware hypothesis/catalog resolution | `5254f65` |
| F02: replayable formula-expression contract | `44324d1` |
| F03: pre-extension fit-metadata validation | `7a9b603` |
| F04: mixed-formula coefficient transforms | `842f251` |
| F05: inverse-gamma negative-power boundary | `921418f` |
| R02: release version and vignette caches | `8601948` |

## Findings requiring maintainer approval

### PF01 - Exact finite transformed endpoints can be classified outside support

**Location.** `R/prior-density-ordinate.R:1813-1889`

**Issue.** A finite endpoint is reduced to an ordinary double and later
recovered through an inverse `exp()`/`log()` transformation. Two algebraically
equivalent representations of the same endpoint need not round to the same
double. The classifier then treats the tiny round-trip discrepancy as a real
support violation.

**Reproduced example.** Let a truncated inverse-gamma source have lower bound
`0.01`, and transform it as `Y = exp(-2 * log(X)) = X^-2`. The mathematical
upper support endpoint is `0.01^-2 = 10000`, where the transformed density is
finite. The current implementation reports `behavior = "zero"` and reason
`"The requested value is outside the transformed prior support."` The forward
endpoint is represented as `9999.999999999991`, while inverse-transforming
`10000` produces `0.009999999999999995`.

**Impact.** Exact-boundary density ordinates can be wrong even when the
boundary and transformation are structurally known. This can incorrectly
permit or reject a Savage-Dickey density ratio at a finite transformed
boundary. The F05 regression currently avoids the endpoint with
`(1 - 1e-12)`, so it does not protect this case.

**Suggested change.** Preserve finite source/output endpoint correspondence in
the transformation provenance and classify a recognized endpoint from its
original source endpoint. If interval rounding is needed, constrain it to
provenance-derived endpoint enclosures (for example, adjacent representable
values), not a global epsilon or general support clamping. Add exact lower and
upper endpoint tests for increasing and decreasing transforms.

**Decision:** pending

### PF02 - Accepted unindexed expression columns generate invalid scalar JAGS nodes

**Location.** `R/JAGS-formula-helpers.R:155-205` and
`R/JAGS-formula.R:530-533`

**Issue.** The expression validator accepts a bare data-column symbol anywhere
in the AST. R replay treats that symbol as a row-length vector, but JAGS treats
the same data object as a vector inside each scalar `mu[i]` assignment. The
whitelist therefore establishes syntactic safety without establishing scalar
JAGS shape.

**Reproduced example.** This accepted call:

```r
JAGS_formula(
  ~ expression(log(x)),
  "mu",
  data.frame(x = c(1, 2, 3)),
  list(intercept = prior("normal", list(0, 1)))
)
```

emits:

```text
for(i in 1:N_mu){
  mu[i] = mu_intercept + log(x)
}
```

R replay returns the expected three values, but JAGS compilation fails with
`Cannot insert node into mu[1]. Dimension mismatch`.

**Impact.** `JAGS_formula()` accepts an expression documented as supported, but
the model cannot fit. Other accepted index forms can likewise have different
R and JAGS dimensional or indexing semantics, so the current validator does
not yet guarantee the promised cross-package replay contract.

**Suggested change.** Compile each whitelisted data-only expression to a
validated row vector at formula construction, store it as internal formula
data, and emit only `internal_expression[i]` into JAGS. Retain the parsed AST
for prediction/new-data evaluation and marginal-likelihood reconstruction.
This creates one scalar JAGS path and one R evaluator instead of trying to
prove equivalent shape semantics for two expression languages. Add a live
compile/fit test for bare-column, explicitly indexed, constant, and mixed
random-effect formulas.

**Decision:** pending

### PF03 - `JAGS_extend()` validates metadata objects but not their relationships

**Location.** `R/JAGS-fit.R:558-594` and
`R/JAGS-parameter-catalog.R:944-1018`

**Issue.** The F03 fix validates every stored schema before extension, but the
validators are object-local. In particular, a catalog extraction key only has
to contain a character dependency vector; the dependencies are not checked
against the fit's parameter registry. `JAGS_extend()` then preserves that
catalog and attaches a current fit contract.

**Reproduced example.** Replacing a random-summary catalog dependency with
`"not_a_registry_coordinate"` still passes
`.bt_validate_parameter_catalog()`. A no-work `JAGS_extend()` accepts and
returns the fit with current contract versions. Resolving the preserved
quantity succeeds, but `parameter_draws()` later fails with
`Unknown registry parameter: 'not_a_registry_coordinate'.`

**Impact.** A serialized fit assembled from individually current but mutually
inconsistent metadata can still be relabelled as a valid current fit. The
failure is deferred to a downstream semantic operation, which is the stale-fit
behavior the compatibility contract is intended to prevent.

**Suggested change.** Add one aggregate fit-metadata consistency validator and
call it before extension and other contract-sensitive entry points. At minimum,
require catalog dependencies to be available registry coordinates, ordinary
catalog rows to agree with registry status and canonical names, sampled
registry coordinates to exist in the chains, and formula name maps to agree
with the registry rows owned by each formula parameter. Keep the existing
object-local validators as schema checks.

**Decision:** pending

### PF04 - Level-qualified hypotheses remain partial outside non-reference treatment coefficients

**Location.** `R/JAGS-parameter-catalog.R:399-446`

**Issue.** Semantic factor components are created only for treatment and
independent priors. A treatment-coded reference level has no catalog quantity,
and orthonormal, mean-difference, and ordered/contrast-coded factor levels are
left as backend coordinates. This conflicts with the public catalog
description that fixed factor coefficients use fitted level labels.

**Reproduced examples.** For the fitted two-level treatment fixture, levels are
`A` and `B`: `x_fac2t[B]` resolves to `mu_x_fac2t`, but `x_fac2t[A]` reports no
matching public quantity even though its coefficient is structurally zero. For
the three-level orthonormal fixture, all of `x_fac3o[A]`, `x_fac3o[B]`, and
`x_fac3o[C]` report no matching public quantity; only backend coordinates
`mu_x_fac3o[1]` and `mu_x_fac3o[2]` exist.

**Impact.** The F01 fix supports an important subset of level-qualified
hypotheses but not the full factor-prior family. Downstream code must know the
contrast type and reconstruct level effects itself, bypassing the semantic
catalog and hypothesis AST.

**Suggested change.** Represent every fitted factor level as a catalog
quantity. Use a structural-zero extraction key for a treatment reference
level, direct registry keys where a level is a sampled coordinate, and derived
linear-combination keys based on the persisted contrast/design metadata for
orthonormal, mean-difference, and ordered encodings. If derived level effects
are intentionally out of scope, narrow the documentation and return an
explicit contrast-aware unsupported error instead of a generic no-match error.

**Decision:** pending

## Verification of the implemented fixes

- Focused regression tests for every approved issue passed before its commit.
- Unit profile: 8,881 passed, 0 failed, 0 warnings, 7 profile skips.
- Live-fit profile: 14,951 passed, 0 failed, 1 expected warning from the
  deliberately data-free fully fixed model.
- Fixture profile: 10,781 passed, 0 failed, 0 warnings, 1 profile skip.
- Visual profile: 810 passed, 0 failed, 0 warnings, 11 profile skips.
- Visual-fixture profile: 666 passed, 0 failed, 0 warnings, 26 profile skips.
- All three regenerated vignette caches passed contract validation, and
  `ComparisonR`, `SpikeAndSlab`, and `RandomEffects` rendered from cache without
  fitting.
- CI-equivalent `rcmdcheck --no-manual --as-cran` with the `all` profile:
  0 errors, 0 warnings, 2 environment-only notes. Package installation,
  compiled code, documentation, examples, tests, and vignette rebuilding all
  passed.

No code changes for PF01-PF04 were made during this second review; they await
maintainer approval.
