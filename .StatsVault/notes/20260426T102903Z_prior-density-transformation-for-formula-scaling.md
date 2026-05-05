# Prior density transformation for formula scaling

Timestamp: 2026-04-26T10:29:03.107441+00:00
Type: research-log

Investigated original-scale prior density overlays for formula_scale regression terms. Current BayesTools path generates standardized prior samples, applies the existing coefficient unscale matrix, and plots a KDE of transformed samples. A deterministic alternative is available for scalar independent coefficient priors: for each original-scale coefficient row beta_j = sum_i M[j,i] theta_i, compute the marginal density as the convolution of scaled standardized prior marginals f_{M[j,i] theta_i}. This handles arbitrary priors with evaluable densities and ranges, including mixtures and point masses by decomposing into continuous and discrete components. Other packages mostly avoid a general arbitrary-prior density transform: rstanarm/brms use adjusted priors or transformed/prior-predictive draws; BayesFactor uses scale-adaptive g-priors; MCMCglmm gelman.prior solves only the multivariate-normal covariance transform case. Added local mockup script dev/prior-unscale-density-mockup.R for a two-predictor interaction case.
