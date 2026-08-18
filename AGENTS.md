# AGENTS.md

Guidance for coding assistants working in BayesTools. Be extremely concise.

## Package Overview

BayesTools provides reusable Bayesian infrastructure for prior specification,
JAGS model generation and fitting, marginal likelihoods, model averaging,
posterior summaries, diagnostics, and interpretation.

- Backend: JAGS through `runjags`/`rjags`, with focused native distributions
  and numerical kernels.
- Estimation: ordinary and formula-generated JAGS models, plus bridge sampling
  for marginal likelihoods.
- Consumers: packages such as RoBMA depend on BayesTools' stable metadata and
  numerical contracts.
- System requirement: JAGS 4.3.0 or newer.

## Numerical and Engineering Guardrails

### Numerical Faithfulness

- Correct inference-changing numerical errors; do not pursue arbitrary
  machine-level perfection.
- Implement the statistical model as specified. Never silently alter user
  inputs or computed values with epsilons, clamping, covariance repair, or
  similar heuristics.
- Invalid inputs must fail clearly. Mathematically valid boundary cases must
  retain their defined behavior.
- Prefer established implementations from base R, JAGS, or standard numerical
  libraries.
- Support statistically meaningful extreme cases. Do not promise correctness
  for every representable binary64 input unless the public API requires it.
- Base numerical warnings on structural or provenance information where
  possible. Do not infer prior-density behavior from posterior samples,
  numerical grids, or KDE output.
- Use independent references, analytic identities, or human-verified results
  for numerical tests. Do not test an implementation against the same
  unverified calculation.

### Complexity Budget

- Prefer the smallest boring implementation that solves the observed problem.
- Do not add general exact-arithmetic, arbitrary-precision, or binary64
  emulation subsystems.
- Native C/C++ is justified only by JAGS integration, statistical correctness,
  or a substantial measured performance benefit.
- Before adding native code, abstractions, or source files, verify that existing
  R, JAGS, or package infrastructure cannot solve the problem cleanly.
- Avoid parallel implementations, duplicated numerical paths, speculative
  extensibility, and abstractions serving one simple call site.
- Generic Bayesian functionality belongs in BayesTools. Analysis-specific
  integration belongs in the downstream package that owns that analysis.
- Remove machinery superseded by the current change when it is in scope. Do not
  remove unrelated existing code without maintainer approval.

### Compatibility and Release Notes

- During development of an unreleased patch, do not preserve compatibility with
  earlier iterations of that same unreleased code. Replace inferior architecture
  instead of adding migrations, deprecated aliases, schema adapters, or
  compatibility layers.
- For functionality present in a released version, prefer backward-compatible
  changes.
- If released architecture prevents a clean solution or appears short-sighted,
  ask the maintainer whether to preserve compatibility or make a breaking
  change. Never assume either choice.
- When feature implementation on a branch is complete, increment the
  development version in `DESCRIPTION` and update `NEWS.md` in the same feature
  change.
- Write `NEWS.md` for the final branch state that will be squashed into the
  future release. Do not describe intermediate revisions of unreleased
  features.

### Testing Budget

- Add the smallest high-information regression test proving correctness. Do not
  add tests merely to increase coverage.
- Keep the routine non-refit suite within a 15-minute target on the reference
  machine. Expensive fitting is a separate verification lane.
- Stop an expensive fitting or numerical-certification run approaching two
  hours. Split it into independently runnable pieces and retain cases with a
  high evidence-to-runtime ratio.
- Prefer focused tests during development. Do not repeatedly run the full suite
  without a concrete reason.
- Preserve human-verified `vdiffr` snapshots. Never replace meaningful visual
  comparisons with superficial render-only assertions.
- Do not weaken, skip, regenerate, or loosen an existing regression expectation
  merely to make a changed implementation pass.

### Engineering Behavior

- Before non-trivial work, state consequential assumptions and success
  criteria. Do not narrate obvious assumptions.
- Use a short visible plan for multi-step work; revise it when evidence changes.
- When requirements, code, tests, or documentation conflict, name the conflict
  and seek or record a decision instead of guessing.
- Push back with concrete correctness, complexity, maintenance, or runtime
  costs. Propose the simpler alternative and follow the maintainer's decision.
- Establish a simple correct implementation or reference before optimizing.
  Preserve behavior while optimizing and measure performance claims.
- Distinguish dead code introduced by the current change from unrelated
  existing code. Remove the former; do not remove the latter without approval.
- State uncertainty, verification performed, and verification omitted. Never
  present incomplete evidence as completion.

### Change Discipline

- Keep changes focused and reviewable. Commit related code, tests, and
  documentation together when the maintainer requests a commit.
- Do not refactor unrelated code while fixing a numerical issue.
- For ambiguous statistical or architectural choices, stop and record the
  issue, impact, alternatives, and recommendation. If the choice must remain
  pending, create `.agents/instructions-decisions.md`; remove the entry, and the
  file when empty, after the maintainer decides.
- Before finishing, review changed files for unnecessary abstractions,
  compatibility layers, tests, native code, and unrelated edits.
- Do not run Git operations, including status, diff, add, commit, branch,
  checkout, merge, rebase, reset, tag, push, or pull, unless the maintainer
  explicitly asks for that Git operation in the current turn.
- Keep only durable, cross-task guidance in `AGENTS.md` and
  `.agents/instructions/`. Do not retain completed plans, prompts, transcripts,
  generic engineering boilerplate, scratch files, or resolved decision logs.
- Add a guide only for stable subsystem rules. Preserve any still-valid unique
  rule before deleting an obsolete guide.

## Detailed Instructions

Read only the guide relevant to the files being changed:

- [testing.md](.agents/instructions/testing.md): test profiles, cached fits,
  reference output, coverage, and visual regression.
- [priors.md](.agents/instructions/priors.md): prior implementation and the
  structural prior-density ordinate contract.
- [jags-formula.md](.agents/instructions/jags-formula.md): generated syntax,
  formula scaling, random effects, fitted metadata, and marginal likelihoods.
- [plotting.md](.agents/instructions/plotting.md): plot-data/rendering separation
  and visual verification.
- [vignettes.md](.agents/instructions/vignettes.md): vignette caching,
  regeneration, style, and citations.

Do not load all guides by default. Before changing an instruction, verify every
referenced file and function against the current tree. Create
`.agents/instructions-decisions.md` only while a real maintainer choice remains
unresolved.

## R Code Style

- Read the existing implementation and tests before changing behavior. Match
  local style and keep edits scoped.
- Use `snake_case` for new ordinary functions, arguments, variables, and files.
  Preserve established families such as `JAGS-*.R`; use dots only for S3
  dispatch or an established public API.
- Use `<-` assignment, two-space indentation, and `TRUE`/`FALSE`.
- Do not introduce pipe-heavy code. Name intermediate results when that makes
  control flow clearer.
- Prefer base data structures, `vapply()` for type-stable atomic iteration, and
  explicit loops when they are clearer than forced vectorization.
- Use the validators in `R/tools.R`, including `check_bool()`, `check_char()`,
  `check_int()`, `check_real()`, and `check_list()`. Do not use `stopifnot()` for
  exported or user-facing validation.
- Use concise `stop()`, `warning()`, and `message()` calls. Prefer
  `call. = FALSE` for user-facing conditions where it fits existing style.
- Qualify non-base calls in scripts and examples when clarity matters. Follow
  the established namespace-import style in package source.
- Align adjacent assignments and named arguments when it improves readability;
  do not reformat unrelated code.
- Leave a blank line after an opening brace in function definitions.
- Never call `setwd()` in package code or use absolute local paths, hardcoded
  credentials, or hidden writes to user directories.

## Architecture

- Priors and density algebra: `R/priors*.R`, `R/distributions-*.R`, and
  `R/prior-density-ordinate.R`.
- JAGS runtime, formulas, and bridge sampling: `R/JAGS-*.R`.
- Random-effect contracts: `R/random-effects-*.R` and
  `R/random-group-covariance.R`.
- Marginal inference and model averaging: `R/marginal-*.R` and
  `R/model-averaging*.R`.
- Summaries and validation: `R/summary-tables*.R`, `R/interpret.R`, and
  `R/tools.R`.

### Formula and Random-Effect Semantics

- BayesTools formulas intentionally do not inherit ordinary
  `stats::model.matrix()` intercept/contrast coupling. `prior_factor()` owns a
  fixed factor's contrast family. Removing the intercept represents a
  structural zero intercept while preserving that prior-owned factor basis; it
  must not silently force raw indicator or treatment coding.
- `id()`, `diag()`, and `us()` / `un()` are general random-coefficient
  structures. Their left side is a coefficient formula: `1`, `0`, and `-1`
  control the random intercept, and continuous slopes, factor slopes, and
  interactions are supported. Plain `(expr | group)` defaults to `us()` and
  `||` to `diag()`. `id()` uses one shared SD for independent columns,
  `diag()` uses one SD per independent column, and `us()` estimates one SD per
  column plus an unstructured correlation matrix.
- Random-coefficient factor coding is owned by the block's concrete contrast
  metadata, not by intercept syntax. Existing fixed-factor contrast metadata is
  reused by default when available; `random_block(contrasts = ...)` explicitly
  overrides it for that random block. Thus `us(0 + group | study)` with an
  independent block contrast gives one correlated coefficient per `group`
  level and no random intercept, whereas `0 + group` alone does not force a
  level-indicator basis.
- `cs()` / `hcs()`, `ar1()` / `ar()` / `har()`, and `car()` are
  structure-owned index specifications, not coefficient formulas. `cs()` and
  `hcs()` accept one or more discrete index columns and combine multiple
  columns by their observed interaction. `ar1()` and `har()` accept exactly one
  discrete index column; existing factor levels or sorted unique values define
  its order. These discrete structures accept factor, character,
  numeric/integer, or logical values and persist the resolved levels.
- `car()` accepts exactly one finite numeric/integer coordinate, or an ordered
  factor with numeric level labels, and uses actual coordinate distances.
  Structure-owned index specifications reject explicit `1`, `0`, and `-1` and
  reject `random_block(contrasts = ...)`. `hcs()` has level-specific SDs and one
  common pairwise correlation; it is not equivalent to the unrestricted
  correlation matrix from `us()`.
- Persist basis ownership, concrete index levels, design columns, and public
  labels in formula metadata. Prediction, covariance reconstruction,
  summaries, and downstream packages must consume that metadata rather than
  reconstructing a basis from formula text.
- LKJ primitive coordinates and other covariance-construction nodes are
  internal implementation parameters. Register them as internal through
  `JAGS_parameter_registry()` and expose only semantic SDs, correlations, and
  other declared public summaries.
- Complete omitted correlation priors only after the random-effect structure
  and dimension are resolved: US/UN uses `LKJ(1)`; CS/HCS uses a raw uniform
  prior on `(-1 / (K - 1), 1)`; AR1/HAR uses raw `Uniform(-1, 1)`; and CAR uses
  raw `Uniform(0, 1)`. Explicit scalar priors retain the Fisher-z default scale.
  BayesTools must not invent a generic SD magnitude because that scale belongs
  to the outcome model; direct `JAGS_formula()` use still requires an SD prior,
  SD source, or variance allocation.

`JAGS_parameter_registry()` is the authoritative mapping from posterior columns
to semantic roles, formula terms, fitted scales, and display labels. Downstream
code must use the registry or its accessors instead of parsing JAGS names or
fitted-object internals.

Native JAGS distributions live in `src/distributions/`; shared kernels are in
`src/invgamma/`, `src/lkj/`, and `src/nonlocal/`; registrations are in
`src/BayesTools.cc` and `src/init.c`. Update registration and `Makevars*` files
when adding native sources.

## Documentation and CRAN

- Add or update roxygen2 documentation for exported functions. Use
  `\insertCite{key}{BayesTools}` for package documentation; vignettes use
  Pandoc citations.
- Keep dependencies minimal; prefer base R or existing imports.
- Use `skip_on_cran()` for computationally intensive tests.
- Never regenerate reference output or visual snapshots without maintainer
  review.
- Commit durable agent instructions. Keep `AGENTS.md` and
  `.agents/instructions/` excluded through `.Rbuildignore`; keep CI workflows
  under `.github/workflows/`.
