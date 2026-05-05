# BF Monte Carlo error display in BayesTools tables

Timestamp: 2026-04-30T05:46:03.206854+00:00
Type: research
Related papers: `gronau2017jmp_tutorial`, `oberauerndmisc_variance2`, `ardia2012csda_comparative`, `sinharay2005jcgs_empirical`

For BayesTools BF diagnostics, reviewed bridge-sampling and marginal-likelihood error literature. The relevant design point is that numerical uncertainty for marginal likelihood ratios is naturally handled on the log scale, while bridge-sampling accuracy is commonly summarized by relative mean-squared error/coefficient of variation. For inclusion BFs computed from model indicator frequencies, the implementation analog is the delta-method SE(log BF) = SE(p) / (p(1-p)); reporting 100 * SE(log BF) as a relative BF MC error percentage gives a one-column, scale-invariant diagnostic. Caveat: this is numerical/Monte Carlo error only and does not address MCMC non-convergence or poor bridge/proposal overlap.
