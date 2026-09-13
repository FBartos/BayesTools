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
- When a product-space Bayes factor has zero or one draw in the rare state,
  require agreement in conclusion, not high numerical precision, if bridge or
  density-based estimates show overwhelming evidence in the same direction.
  Keep reported error percentages and uncertainty visible. This does not relax
  density or integration diagnostics; estimator ESS is not state occupancy.

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
- Refactor and unify unreleased implementations across BayesTools and RoBMA
  within the requested scope when this improves correctness, performance, or
  maintenance. This permission does not extend to breaking released interfaces.
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
- Apply the maintainer's existing authorization without asking again for
  routine in-scope fixes or unreleased refactors. Before substantial work,
  consult `.agents/instructions-decisions.md` when present for pending design
  choices.
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

### Performance Investigations

- Profile suspiciously slow original calls. Prioritize missing optimized paths,
  unused compiled metadata, and repeated construction or calculation of values
  already available before adding new machinery.
- For large cross-package investigations, use several agents for independent
  analysis and review, and coordinate work in BayesTools and RoBMA. Read each
  package's instructions and pending decisions. Parallel agent work does not
  authorize adding parallel computation to scenarios.
- Measure improvements with the same scenario expressions, draws, requested
  grids and sample budgets, seeds, diagnostic criteria, and parallel settings.
  Do not reduce or otherwise alter numerical budgets to manufacture speedups.
  Avoid resource contention during comparative timing runs.
- When defaults are insufficient, explicit sample or integration budget
  increases in the affected scenario calls are authorized, including for
  qCMDE failures. Keep diagnostic criteria unchanged. Do not introduce arbitrary
  hidden sampling changes or silently alter package defaults.
- Record budget increases and their timings separately from matched-workload
  performance evidence. Historical minima or an unchanged artifact name do not
  establish a comparable workload or a fresh measurement.

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
- Persist fitted-parameter metadata as one versioned `parameter_map()` with
  linked `coordinates`, `quantities`, and `aliases` tables. Obtain the concrete
  backend view with `parameter_coordinates()` and the semantic public view with
  `parameter_catalog()`. User selection, display, hypotheses, plotting, and
  density estimation must resolve catalog quantities and obtain draws through
  `parameter_draws()`. Never expose a coordinate as a public alias merely
  because it is monitored.
- LKJ primitive coordinates and other covariance-construction nodes are
  internal implementation parameters. Mark them internal in the map's
  coordinate table and omit them from its semantic quantities; expose only
  declared semantic SDs, correlations, allocations, and other public summaries.
- Canonical random-effect catalog names use
  `(formula) owner: quantity(arguments)`. Downstream packages can explicitly
  request centrally generated simplified names. Simplification removes only a
  sole `intercept` argument (`sd(intercept)` becomes `sd`) and permits omission
  of an owner only when resolution remains unique; non-intercept arguments stay
  explicit. Parentheses contain coefficient or parameter names and square
  brackets contain factor or index levels. Public correlations use `cor`, never
  the backend `rho` coordinate. A known group covariance has a fitted `sd`/`var`
  kernel scale, not an `sd_mult`/`var_mult`. Total-variance
  allocations expose `sd_total`, `var_total`, and `var_prop(...)`;
  mean-variance allocations expose `sd_common`, `var_common`,
  `var_mult(...)`, and `sd_mult(...)`.
- A bare random formula or unnamed one-entry formula list has no redundant
  top-level component prefix. An explicitly named one-entry list retains its
  name. Lists with two or more entries replace missing names with
  `component 1`, `component 2`, and so on. Allocation `name` is a required
  stable backend identifier; `display_name` and `component_names` independently
  own its public semantic labels.
- Complete omitted correlation priors only after the random-effect structure
  and dimension are resolved: US/UN uses `LKJ(1)`; CS/HCS uses a raw uniform
  prior on `(-1 / (K - 1), 1)`; AR1/HAR uses raw `Uniform(-1, 1)`; and CAR uses
  raw `Uniform(0, 1)`. Explicit scalar priors retain the Fisher-z default scale.
  BayesTools must not invent a generic SD magnitude because that scale belongs
  to the outcome model; direct `JAGS_formula()` use still requires an SD prior,
  SD source, or variance allocation.

The stored parameter map is authoritative. Its coordinate view owns concrete
fitted coordinates and provenance; its catalog view owns public quantities and
exact aliases. Code must use these views instead of parsing JAGS names, summary
labels, or fitted-object internals.

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
