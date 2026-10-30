# Maintainer decisions and rationale

These decisions were settled with the maintainer on 2026-10-26. Read them when
reviewing the corresponding contracts. They supersede earlier audit proposals
and retained-coordinate reporting exceptions. Fix defects while preserving the
agreed meaning; changing that meaning requires a new maintainer decision, not
an automatic restoration of an older release or test expectation. Record any
later replacement here with its reason. Technical details remain in
[JAGS and formulas](jags-formula.md) and [priors](priors.md).

## One scaling switch for fixed and random summaries

**Decision.** In estimates tables, `transform_scaled` controls the predictor
units of fixed coefficients and random-effect coefficients, SDs, variances and
correlations. `FALSE` means fitted units; `TRUE` requests original units.
`random_effects_summary` controls which random summaries are displayed, not
their units. Allocation parameters retain their declared parameterization.

**Why.** A caller who declares `formula_scale_list = list(mu = list(x = TRUE))`
and selects a display scale expects the whole coefficient block to use that
scale. Previously, a semantic random summary could silently use original units
beside fitted-unit fixed coefficients. The sleepstudy comparison exposed a large
intercept-slope correlation discrepancy that was partly a coordinate mismatch.
Summary depth is not a meaningful reason to change the estimand.

Centering matters as well as units. With `x_s = (x - m) / s`, a random intercept
and slope transform as `b0_original = b0_fitted - (m/s) * b1_fitted` and
`b1_original = b1_fitted/s`. The covariance transforms as `A G A'`; the intercept
SD and its correlation can therefore change markedly. Dividing slope SDs alone
is insufficient. Comparable coordinates do not imply equal estimates across
packages with different priors or posterior sampling.

Do not attach an original-scale catalog quantity ID to a fitted-scale table row.
Catalog quantities retain their own declared semantics; this decision does not
make them depend on table options. Implemented in `6684637`; regression coverage
is `transform_scaled selects fixed and random summary scales at every depth` in
`tests/testthat/test-JAGS-random-effect-scaling.R` and the summary tests of
`tests/testthat/test-random-effects-summary-posterior.R`.

## An unavailable original coordinate is NA in an estimates table

**Decision.** If a requested original-scale random quantity cannot be represented
or identified by the fitted block, keep its table row as `NA` with the specific
reason. Apply this to standard/full/raw displays and returned table samples;
retain identifiable quantities in the same block. The rule is per row (settled
with the maintainer on 2026-10-27 after the pre-release review): for a
correlated block with an SD shared by several columns, only the shared SD (and
its variance) and the rows of unrepresented columns are `NA`; the SD of an
unshared column and the correlations among identified columns keep their
values. A block whose whole design is refused (linearly dependent, or a fallback
that mixes categorical coding) has no identified row, and all its rows are `NA`.
The unrepresented columns are those that a dropped centering contribution would
have to shift: the intercept for `f:x`, and likewise the lower-order slope `z`
for `f:x:z`, and so on for every numeric group of the block.

**Why.** For `diag(1 + f:x | id)` with centered `x`, the transformation can
require factor-specific intercepts absent from the block. Its surviving fitted
intercept remains an effect at the predictor mean, not the requested common
original-scale intercept. Showing that number with a footnote still invites
comparison of different quantities. `NA` makes the boundary visible without
discarding valid slope summaries. This supersedes the earlier retain-and-note
table policy.

This is a reporting policy, not permission to invent lower-order random terms,
alter the fitted model, catch unrelated fixed-effect/input errors, or replace
typed refusals in general inference accessors with NAs. Implemented in
`98c8af4`; the boundary is covered by the table-availability tests of
`tests/testthat/test-JAGS-random-effect-scaling.R` (`tables report unavailable
induced intercepts as NA with a note`, the shifted-slope and shared-leaf
tests) and by the logical-column block tests of
`tests/testthat/test-random-effects-summary-posterior.R`.

## Preserve exact prior laws for pure SD rescaling

**Decision.** When the structural transformation maps an SD to `abs(c)` times
its own fitted SD, preserve the corresponding scaled prior law. Named versus
unnamed index attributes cannot turn that map into a combination of SDs.

**Why.** The prior transformation is exact and is needed by prior-density plots
and prior-ordinate inference. A false composite classification needlessly made
those operations unavailable even though posterior summaries were valid.
Conversely, an intercept SD that truly combines several random coefficients
does not inherit one fitted SD's prior. Do not use numerical smallness, KDEs,
or posterior draws to certify a structural prior law.

The C07 fix (`ad49eac`) corrects classification; it does not introduce generic
vector-owner support into `parameter_prior_density()`. Persisted maps changed
to version 17 because they stored the erroneous classification. Older development
fits must be refitted, not relabelled or resealed. The maintainer authorized the
declared vignette-cache refresh (`1494bb2`); that one-time authorization does not
change the general vignette policy. `named fallback maps preserve pure random-slope
prior rescalings` in `tests/testthat/test-JAGS-random-effect-scaling.R` and the
classification checks of `tests/testthat/test-JAGS-parameter-catalog.R` cover
classification and supported prior-density routes.

## Meandif and orthonormal factor terms are exchangeable units

**Decision** (settled with the maintainer on 2026-10-30). In ID/DIAG/US
random-effect blocks a meandif or orthonormal factor term is one exchangeable
unit under every structure: one shared SD and zero within-unit correlations.
The units are `f`, `f:x` (`x` continuous) and `f:g` (all factors so coded); a
term crossing them with treatment, independent, or ordered factors has one unit
per partner column (`sd(f:h[h2])`). Under US a unit keeps its correlations with
the other columns, with the LKJ(eta) prior restricted to an identity unit block,
`det(R)^(eta - 1)`, implemented by ordering the unit first in the C-vine with
its within-unit primitives fixed at 1/2. Two or more units of three or more
levels in one US block are refused with a classed error rather than made
mutually uncorrelated. Two-level factors (one column) are ordinary LKJ columns
in the design order. Summaries report `sd(f)`, the correlation with each level
`cor(o, f[a])` (also for two-level factors), and the multiple correlation
`cor_total(o, f)`; the coordinates `f{j}` are internal. An SD shared by several
columns (a unit, or the indicator columns of a logical flag's interaction with
a slope) is certified on the original scale whenever its columns are pure
rescalings by one common factor, in diagonal and correlated blocks; the
diagonal equal-variance-weight certificate remains. This revises the earlier
rule that a correlated block never certifies a shared SD.
`random_variance_allocation(target = "block")` counts SD leaves, so a block
that is one unit has one SD component.

**Why.** `contr.orthonormal()` and `contr.meandif()` take an eigenbasis of the
centering matrix, which is not unique for three or more levels; LAPACK rounding
picks it, and it differs across platforms. Fixed effects are unaffected
(their mnormal/mt/mcauchy priors are spherical), but per-column SDs and an LKJ
prior over contrast columns made the random-effect model of the level
deviations depend on that arbitrary basis (CI showed it as platform-dependent
results). An exchangeable unit is the basis-free model the contrast declares.
Silently restricting the model when an implementation is hard (for example,
making two units uncorrelated) would change the statistical target, so such
cases are refused instead. Cross-platform reproducibility of the basis itself
is a separate follow-up.
