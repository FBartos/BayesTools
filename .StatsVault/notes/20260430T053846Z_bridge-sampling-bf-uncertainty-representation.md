# Bridge-sampling BF uncertainty representation

Timestamp: 2026-04-30T05:38:46.601928+00:00
Type: research

Investigated whether bridge-sampling numerical uncertainty for Bayes factors should be represented as relative/percentage error rather than lower/upper BF columns. Materialized bridge-sampling and marginal-likelihood sources: Meng & Wong (1996), Frühwirth-Schnatter (2004), Gronau et al. (2017), Sinharay & Stern (2005), Gelman & Meng (1998), Overstall & Forster (2010), Robert et al. (2009), Llorente et al. (2023), Kass & Raftery (1995). Core rationale: bridge-sampling theory defines estimator accuracy through asymptotic relative mean-squared error / coefficient of variation for ratios of normalizing constants; on the log scale, BF uncertainty is additive and can be propagated from log marginal likelihood errors, with relative BF error approximately sqrt(se_logml_1^2 + se_logml_0^2) when independent. BF-scale intervals are transforms of this quantity and can be skewed or visually over-interpreted as inferential intervals.
