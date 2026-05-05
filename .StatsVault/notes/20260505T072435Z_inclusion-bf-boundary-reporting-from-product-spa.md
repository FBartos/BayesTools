# Inclusion BF boundary reporting from product-space samples

Timestamp: 2026-05-05T07:24:35.076359+00:00
Type: research_synthesis

Reviewed local summaries and subagent notes for reporting inclusion Bayes factors from product-space/model-index MCMC when the sampled inclusion indicator is all 0 or all 1. Current BayesTools behavior: `inclusion_BF()` returns the exact odds-ratio boundary (`0` or `Inf`) from posterior model probabilities; `runjags_inference_table(..., BF_diagnostics=TRUE)` computes delta-method `error%(BF)` for interior probabilities and returns `NA` plus a warning when one side has zero posterior visits.

Synthesis:

1. Heck et al. (2019) / MCMCprecision supports diagnostics, not finite boundary correction. Its Markov-model method estimates uncertainty of posterior model probabilities from the full model-index chain by sampling transition matrices and stationary distributions. For inclusion BFs, the principled extension is to fit the full model-index transition matrix, then compute inclusion subset probabilities from stationary-distribution draws. It explicitly does not justify adding pseudocounts to an aggregated all-0/all-1 inclusion bit; never/rarely visited models make BFs unreliable and indicate insufficient product-space mixing or resolution.

2. BEAST/BSSVS and Gámbaro et al. (2025) support an inequality-style reporting convention for saturated posterior inclusion probabilities. When posterior inclusion `p = 1`, replace only for display with the closest resolvable value `(N - 1) / N`, yielding a lower bound `BF > (N - 1) / prior_odds`, where `prior_odds = q / (1 - q)`. A symmetric BayesTools extension for all-0 indicators would use `p = 1 / N`, yielding `BF < 1 / ((N - 1) * prior_odds)`. This should be labelled as a Monte Carlo resolution bound, not a corrected/unbiased BF and not `BFadj`; Gámbaro's BFadj is a domain-specific empirical-prior adjustment for phylogeographic sampling bias.

3. Heck & Davis-Stober / multinomineq uses the encompassing-prior BF framework for inequality-constrained multinomial models. Its stepwise/automatic algorithms use beta uncertainty for estimated constraint-hit proportions and ensure positive hits for components where ratios are needed. This can inspire optional MC interval/bound displays for count-based proportions, but not a finite product-space point BF. If used, smoothing should be labelled as uncertainty regularization and should preferably use ESS-adjusted effective counts because product-space indicators are autocorrelated.

4. Klugkist/Laudy/Hoijtink/Gu informative-hypothesis sources show why an encompassing-prior BF remains finite when all posterior draws satisfy a constraint: `BF = posterior_constraint_mass / prior_constraint_mass`, so plug-in posterior mass 1 gives `1 / prior_constraint_mass`. This finite result is a property of the constrained-vs-encompassing estimand, not a correction for product-space model odds. It should be borrowed for explanation only: product-space all-1/all-0 indicators mean the complement model probability is unresolved by the chain.

Recommended BayesTools reporting direction:

- Keep the internal numeric BF estimator exact: `Inf` for all included and `0` for all excluded.
- Improve table/reporting display by optionally replacing boundary displays with one-sided resolution bounds: `> bound` for all-1 inclusion BFs and `< bound` for all-0 inclusion BFs, with formulas above and with `N` or preferably indicator ESS stated.
- Keep `error%(BF) = NA` at boundaries, but make the warning/footnote visible whenever a boundary BF is printed, not only when BF diagnostics are requested.
- Do not silently use Jeffreys/add-one smoothing as the main BF. If an interval is added, present it as a Monte Carlo uncertainty interval over the indicator probability transformed to BF scale, not as the estimate.
- Longer-term: add a full model-index diagnostic path following MCMCprecision for product-space model indicators, then derive inclusion-BF uncertainty from stationary model-probability draws. This is the principled upgrade for autocorrelated product-space chains.
