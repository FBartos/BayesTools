skip_if_not_test_profile("unit")

# Inclusion gates and spike-and-slab totals of R-squared budgets.

.R2_gate_data <- function(){
  set.seed(7)
  data.frame(f = factor(rep(c("a", "b", "c"), each = 10)), x = stats::rnorm(30), z = stats::rnorm(30))
}
.R2_gate_budget <- function(inclusion = prior("beta", list(1, 1)), total = prior("beta", list(1.5, 2)),
                            convention = "common", law = "exact", concentration = "inherited"){
  prior_R2(terms = c("x", "z", "f"), total = total, allocation = prior("dirichlet", list(alpha = c(1, 0.5, 2))),
    law = law, convention = convention, concentration = concentration, sigma = 1, contrast = "orthonormal",
    inclusion = inclusion)
}
.R2_gate_bind <- function(budget, formula = ~ x + z + f){
  JAGS_formula(formula, "mu", .R2_gate_data(), list(intercept = prior("normal", list(0, 1)), R2 = budget))$prior_list
}
.R2_spike <- function(p, slab = prior("beta", list(1.5, 2))){
  prior_spike_and_slab(slab, prior_inclusion = prior("spike", list(p)))
}

test_that("prior_R2() accepts spike-and-slab totals and inclusion gates and validates them", {

  budget <- .R2_gate_budget(total = .R2_spike(0.6))
  expect_identical(budget$total_inclusion, list(type = "spike", p = 0.6))
  expect_identical(budget$inclusion$form, "common")
  expect_output(print(budget), "spike:      zero signal in the budget with probability 0.4", fixed = TRUE)
  expect_output(print(budget), "shared rate ~ Beta(1, 1) (beta-binomial model prior)", fixed = TRUE)

  terms <- .R2_gate_budget(inclusion = list(x = prior("spike", list(0.3)), f = prior("beta", list(2, 2))))
  expect_identical(names(terms$inclusion$laws), c("x", "f"))
  expect_equal(terms$inclusion$laws$f$p, 0.5)
  expect_output(print(terms), "x: 0.3, f: Beta(2, 2); other terms always included", fixed = TRUE)

  expect_error(.R2_gate_budget(total = prior_mixture(list(prior("beta", list(1, 1)), prior("spike", list(0))), is_null = c(FALSE, TRUE))),
    "other mixtures are not supported", fixed = TRUE)
  expect_error(.R2_gate_budget(total = prior_spike_and_slab(prior("normal", list(0, 1)))), "'total' must be a Beta prior", fixed = TRUE)
  expect_error(.R2_gate_budget(inclusion = prior("spike", list(1))), "strictly between 0 and 1", fixed = TRUE)
  expect_error(.R2_gate_budget(inclusion = prior("normal", list(0, 1))), "must be a probability prior", fixed = TRUE)
  expect_error(.R2_gate_budget(inclusion = list(w = prior("spike", list(0.5)))), "not in 'terms': 'w'", fixed = TRUE)
  expect_error(.R2_gate_budget(inclusion = list(prior("spike", list(0.5)))), "named by the term labels", fixed = TRUE)
})

test_that("the gate law is exact: configuration probabilities, truncation and marginals", {

  # shared rate Beta(1, 1): beta-binomial, uniform over sizes
  gates <- .R2_gate_bind(.R2_gate_budget())$mu_R2$gates
  grid <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), 3)))
  p <- .bt_R2_gate_configuration_probability(gates, grid)
  expect_equal(sum(p), 1)
  expect_equal(p[rowSums(grid) == 0], 1 / 4)
  expect_equal(unique(round(p[rowSums(grid) == 1], 12)), 1 / 12)

  # with a spike and every term gated the gate prior is conditioned on S nonempty
  spiked <- .R2_gate_bind(.R2_gate_budget(total = .R2_spike(0.6)))$mu_R2$gates
  expect_true(spiked$truncated)
  q <- .bt_R2_gate_configuration_probability(spiked, grid)
  expect_equal(q[rowSums(grid) == 0], 0)
  expect_equal(q[rowSums(grid) > 0], p[rowSums(grid) > 0] / (1 - p[rowSums(grid) == 0]))
  # the sequential conditional probabilities reproduce the truncated law
  sequential <- apply(grid, 1, function(open){
    prob <- 1
    for(t in seq_along(open)){
      g <- .bt_R2_gate_probability(spiked, t, sum(open[seq_len(t - 1L)]))
      prob <- prob * if(open[[t]]) g else 1 - g
    }
    prob
  })
  expect_equal(sequential, q, tolerance = 1e-12)
  marginals <- .bt_R2_gate_marginals(spiked, c("x", "z", "f"))
  expect_equal(marginals$budget, 0.6)
  expect_equal(unname(marginals$open), rep(0.6 * 0.5 / 0.75, 3))

  # per-term gates with an ungated term: no truncation under a spike
  mixed <- .R2_gate_bind(.R2_gate_budget(inclusion = list(x = prior("spike", list(0.3)), f = prior("beta", list(2, 2))),
                                        total = .R2_spike(0.7)))$mu_R2
  expect_false(mixed$gates$truncated)
  expect_identical(unname(mixed$term_gate), c(1L, 0L, 2L))
  expect_equal(unname(.bt_R2_gate_marginals(mixed$gates, mixed$terms)$open), c(0.7 * 0.3, 0.7, 0.7 * 0.5))
})

test_that("the core emits the exact sequential gates, the budget indicator and the gated allocation", {

  priors <- .R2_gate_bind(.R2_gate_budget(total = .R2_spike(0.6)))
  expect_true(is.prior.R2_budget(priors$mu_R2))
  expect_null(priors$mu_R2_allocation)
  syntax <- .JAGS_add_priors.fun(priors)
  expect_match(syntax, "mu_R2_gate[1] ~ dbern(0.6666666666666666)", fixed = TRUE)
  expect_match(syntax, "mu_R2_gate[3] ~ dbern(ifelse(sum(mu_R2_gate[1:2]) > 0, (1 + sum(mu_R2_gate[1:2])) / 4, 1))", fixed = TRUE)
  expect_match(syntax, "mu_R2_indicator ~ dbern(0.6)", fixed = TRUE)
  expect_match(syntax, "mu_R2 <- mu_R2_variable * step(", fixed = TRUE)
  expect_true(all(c("mu_R2", "mu_R2_indicator", "mu_R2_gate", "mu_R2_allocation", "prior_par_eta_mu_R2_allocation") %in%
    JAGS_to_monitor(priors)))

  # linked: term signals g_t eta_t / H; fixed_total: shapes rescaled per configuration
  linked <- .JAGS_add_priors.fun(.R2_gate_bind(.R2_gate_budget(convention = "linked", total = prior("beta", list(3.5, 2)))))
  expect_match(linked, "prior_par_H_mu_R2 ~ dgamma(2, 1)", fixed = TRUE)
  fixed <- .JAGS_add_priors.fun(.R2_gate_bind(.R2_gate_budget(concentration = "fixed_total")))
  expect_match(fixed, "prior_par_eta_mu_R2_allocation[1] ~ dgamma(ifelse(mu_R2_gate[1] > 0, 3.5 / (", fixed = TRUE)
})

test_that("prior draws follow the configuration law and the prior of each model", {

  set.seed(1)
  priors <- .R2_gate_bind(.R2_gate_budget(total = .R2_spike(0.6)))
  core <- priors$mu_R2
  draws <- .bt_R2_core_draws(core, 40000)
  open <- draws[, "mu_R2_indicator"] * draws[, c("mu_R2_gate[1]", "mu_R2_gate[2]", "mu_R2_gate[3]")]
  configurations <- .bt_R2_core_configurations(core)
  expect_equal(mean(rowSums(open) == 0), configurations$null, tolerance = 0.02)
  expect_equal(mean(draws[, "mu_R2"] == 0), configurations$null, tolerance = 0.02)
  expect_equal(unname(colMeans(open)), unname(.bt_R2_gate_marginals(core$gates, core$terms)$open), tolerance = 0.02)

  # given S = {x, f}, the total and the allocation are those of the model x + f
  keep <- open[, 1] == 1 & open[, 2] == 0 & open[, 3] == 1
  expect_equal(stats::ks.test(draws[keep, "mu_R2"], stats::pbeta, 1.5, 2)$p.value > 1e-3, TRUE)
  phi <- draws[keep, "mu_R2_allocation[1]"]
  expect_equal(stats::ks.test(phi, stats::pbeta, 1, 2)$p.value > 1e-3, TRUE)
  expect_equal(draws[keep, "mu_R2_allocation[2]"], rep(0, sum(keep)))

  # linked: Beta(A_S, b) given S
  linked <- .R2_gate_bind(.R2_gate_budget(convention = "linked", total = prior("beta", list(3.5, 2))))$mu_R2
  draws <- .bt_R2_core_draws(linked, 40000)
  keep <- draws[, "mu_R2_gate[1]"] == 1 & draws[, "mu_R2_gate[2]"] == 0 & draws[, "mu_R2_gate[3]"] == 1
  expect_equal(stats::ks.test(draws[keep, "mu_R2"], stats::pbeta, 3, 2)$p.value > 1e-3, TRUE)
})

test_that("gated term priors have an atom at zero in their marginal methods", {

  set.seed(2)
  priors <- .R2_gate_bind(.R2_gate_budget(total = .R2_spike(0.7)))
  prior <- priors$mu_x
  closed <- 1 - 0.7 * 0.5 / 0.75
  expect_equal(mcdf(prior, 0) - mcdf(prior, -1e-12), closed, tolerance = 1e-10)
  expect_equal(mquant(prior, c(0.5 - closed / 4, 0.5 + closed / 4)), c(0, 0))
  draws <- rng(prior, 100000)
  expect_equal(mean(draws == 0), closed, tolerance = 0.01)
  q <- c(-0.8, 0.4)
  expect_equal(mcdf(prior, q), vapply(q, function(v) mean(draws <= v), numeric(1)), tolerance = 0.01)
  expect_equal(var(prior), mean(draws^2), tolerance = 0.03)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(prior))
  expect_silent(plot(priors$mu_R2))
  expect_output(print(priors$mu_R2), "R2 gated budget 'R2': R2 = 0 with probability 0.3, R2 ~ Beta(1.5, 2) in the full model", fixed = TRUE)
  expect_output(print(prior), "allocation 1 of 3.5, included with probability 0.47", fixed = TRUE)
})

test_that("gated budgets fit with inclusion summaries, catalog quantities and refuse bridge sampling", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  set.seed(5)
  data <- .R2_gate_data()
  data$y <- 0.8 * (data$f == "b") + 0.5 * data$x + stats::rnorm(30)
  budget <- .R2_gate_budget(total = .R2_spike(0.5))
  args <- list(
    model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
    data = list(y = data$y, N = nrow(data)), prior_list = list(sigma = prior("invgamma", list(1, 1))),
    formula_list = list(mu = ~ x + z + f), formula_data_list = list(mu = data),
    formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 5)), R2 = budget)))
  fit <- do.call(JAGS_fit, c(args, list(chains = 2, adapt = 200, burnin = 300, sample = 1500, seed = 5, silent = TRUE)))
  post <- as.matrix(.fit_to_posterior(fit))
  open <- post[, "mu_R2_indicator"] * post[, c("mu_R2_gate[1]", "mu_R2_gate[2]", "mu_R2_gate[3]")]

  table <- runjags_estimates_table(fit)
  expect_true(all(c("(mu) x (inclusion)", "(mu) f (inclusion)", "(mu) R2{total} (inclusion)") %in% rownames(table)))
  expect_false(any(grepl("gate|indicator|allocation", rownames(table))))
  expect_equal(table["(mu) x (inclusion)", "Mean"], mean(open[, 1]))
  expect_equal(table["(mu) R2[share: x]", "Mean"], mean((post[, "mu_R2_allocation[1]"])[open[, 1] == 1]))
  conditional <- runjags_estimates_table(fit, conditional = TRUE)
  expect_equal(conditional["(mu) x", "Mean"], mean(post[open[, 1] == 1, "mu_x"]))

  inference <- runjags_inference_table(fit)
  expect_equal(inference["(mu) R2{total}", "prior_prob"], 0.5)
  expect_equal(inference["(mu) x", "prior_prob"], 0.5 * 0.5 / 0.75)
  expect_equal(inference["(mu) x", "post_prob"], mean(open[, 1]))

  catalog <- parameter_catalog(fit)
  total <- parameter_catalog_resolve(catalog, alias = "R2{total}")
  expect_identical(total$quantities$support[[1L]]$points, 0)
  share <- parameter_catalog_resolve(catalog, alias = "R2[share: x]")
  expect_identical(share$quantities$definedness, "allocation_active")
  inclusion <- parameter_catalog_resolve(catalog, alias = "x (inclusion)")
  expect_identical(inclusion$quantities$support[[1L]]$points, c(0, 1))
  expect_equal(as.numeric(as.matrix(parameter_draws(fit, inclusion))), unname(open[, 1]))

  expect_error(hypothesis_BF_model(fit, "R2{total} = 0"), class = "BayesTools_point_mass_at_null")
  region <- hypothesis_BF_model(fit, "R2{total} > 0.1")
  prior_prob <- 0.5 * stats::pbeta(0.1, 1.5, 2, lower.tail = FALSE)
  post_prob <- mean(post[, "mu_R2"] > 0.1)
  expect_equal(as.numeric(region$BF), (post_prob / (1 - post_prob)) / (prior_prob / (1 - prior_prob)), tolerance = 1e-8)

  # gated coefficients: an atom at zero in the catalog law, conditional
  # posteriors on the open draws, and region tests with the atom
  coefficient <- parameter_catalog_resolve(catalog, alias = "x")
  expect_equal(parameter_prior_density(fit, coefficient)$points$p, 1 - 0.5 * 0.5 / 0.75)
  conditional_x <- parameter_mixed_posterior(fit, coefficient, conditional = TRUE)
  expect_equal(length(conditional_x), sum(open[, 1] == 1))
  expect_equal(length(as_mixed_posteriors(fit, parameters = "mu_x", conditional = "mu_x")$mu_x), sum(open[, 1] == 1))
  region_x <- hypothesis_BF_model(fit, "x > 0.1")
  prior_prob <- mccdf(attr(fit, "prior_list")$mu_x, 0.1)
  post_prob <- mean(post[, "mu_x"] > 0.1)
  expect_equal(as.numeric(region_x$BF), (post_prob / (1 - post_prob)) / (prior_prob / (1 - prior_prob)), tolerance = 1e-4)
  expect_error(parameter_prior_density(fit, parameter_catalog_resolve(catalog, alias = "f[a]")),
    class = "BayesTools_R2_density_unavailable")

  expect_error(do.call(JAGS_bridgesampling, c(list(fit = fit, log_posterior = function(parameters, data) 0), args[-1])),
    class = "BayesTools_R2_gated_unavailable")
})

test_that("gated budgets are exact at the edges of their laws", {

  set.seed(3)
  # fixed_total with one retained unit: the linked total is Beta(total, b)
  data <- .R2_gate_data()
  budget <- prior_R2(terms = c("x", "z"), total = .R2_spike(0.6, prior("beta", list(2, 2))), law = "exact",
    convention = "linked", concentration = "fixed_total", sigma = 1, inclusion = prior("spike", list(0.5)))
  core <- JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)), R2 = budget))$prior_list$mu_R2
  draws <- .bt_R2_core_draws(core, 40000)
  on <- draws[, "mu_R2"] > 0
  expect_gt(stats::ks.test(draws[on, "mu_R2"], stats::pbeta, 2, 2)$p.value, 1e-3)
  expect_match(.bt_R2_core_syntax(core), "prior_par_eta_mu_R2_allocation[1] ~ dgamma(ifelse(mu_R2_indicator * mu_R2_gate[1] > 0, 2 / (", fixed = TRUE)

  # tiny inclusion probabilities keep the truncated law (no cancellation)
  tiny <- .bt_R2_gate_truncated_probabilities(list(type = "fixed", pi = c(1e-20, 1e-20)), 2L)
  expect_equal(tiny, c(0.5, 1))

  # supports need no enumeration of the configurations (17 gated terms)
  wide <- as.data.frame(matrix(stats::rnorm(40 * 17), 40, 17, dimnames = list(NULL, paste0("v", 1:17))))
  many <- prior_R2(terms = names(wide), total = .R2_spike(0.5, prior("beta", list(1, 1))), law = "exact",
    convention = "common", concentration = "inherited", sigma = 1, inclusion = prior("spike", list(0.5)))
  priors <- JAGS_formula(stats::reformulate(names(wide)), "mu", wide, list(intercept = prior("normal", list(0, 1)), R2 = many))$prior_list
  support <- .bt_R2_quantity_support(structure(list(), prior_list = priors), list(type = "coordinate", dependencies = "mu_R2"))
  expect_identical(support$points, 0)
  expect_error(.bt_R2_core_total_law(priors$mu_R2), class = "BayesTools_R2_density_unavailable")

  # the total of a gated budget has prior methods with its atom at zero
  spiked <- .R2_gate_bind(.R2_gate_budget(total = .R2_spike(0.6)))
  total <- spiked$mu_R2
  expect_equal(mcdf(total, 0), 0.4)
  draws <- rng(total, 100000)
  expect_equal(mean(draws == 0), 0.4, tolerance = 0.01)
  expect_equal(mean(total), mean(draws), tolerance = 0.01)
  expect_equal(var(total), stats::var(draws), tolerance = 0.02)
  expect_equal(mquant(total, 0.3), 0)

  # log densities of gated terms do not underflow; the mean needs no marginal
  expect_true(is.finite(mlpdf(spiked$mu_x, 1e100)))
  expect_identical(mean(spiked$mu_f), 0)
  set.seed(4)
  unequal <- data.frame(f = factor(rep(c("a", "b", "c"), times = c(4, 8, 12))), x = stats::rnorm(24))
  plain <- prior_R2(terms = c("x", "f"), total = prior("beta", list(1, 2)), law = "exact", convention = "common",
    concentration = "inherited", sigma = 1, contrast = "orthonormal")
  expect_identical(mean(JAGS_formula(~ x + f, "mu", unequal, list(intercept = prior("normal", list(0, 1)), R2 = plain))$prior_list$mu_f), 0)

  # a sample-based plot draws the continuous part given inclusion
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(spiked$mu_x, force_samples = TRUE))
})

test_that("conditioning on gated quantities and their scales is explicit", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  set.seed(8)
  data <- .R2_gate_data()
  data$x <- data$x * 2
  data$y <- 0.8 * (data$f == "b") + 0.25 * data$x + stats::rnorm(30)
  fit <- JAGS_fit(
    model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
    data = list(y = data$y, N = nrow(data)),
    prior_list = list(sigma = prior("invgamma", list(1, 1)),
                      theta = prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("spike", list(0.5)))),
    formula_list = list(mu = ~ x + z + f), formula_data_list = list(mu = data),
    formula_scale_list = list(mu = list(x = TRUE, z = TRUE)),
    formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 5)),
      R2 = prior_R2(terms = c("x", "z", "f"), total = .R2_spike(0.5, prior("beta", list(1, 2))), law = "exact",
        convention = "common", concentration = "inherited", sigma = "sigma", contrast = "orthonormal",
        inclusion = prior("beta", list(1, 1))))),
    chains = 1, adapt = 200, burnin = 200, sample = 1000, seed = 8, silent = TRUE)

  # a gated coefficient is tested on the fitted scale of its predictor
  expect_error(hypothesis_BF_model(fit, "x > 0.1", transform_scaled = TRUE), "Test it with 'transform_scaled = FALSE'", fixed = TRUE)
  expect_s3_class(hypothesis_BF_model(fit, "x > 0.1"), "data.frame")

  # the prior side of an event across independent priors; several terms of
  # one budget are dependent and refused
  priors <- attr(fit, "prior_list")
  open <- .bt_R2_gate_marginals(priors$mu_R2$gates, priors$mu_R2$terms)$open[["x"]]
  options <- .condition_event_model_options(priors, .condition_event(priors, c("mu_x", "theta"), "OR"))
  expect_equal(options$event_probability, 1 - (1 - open) * 0.5)
  expect_error(as_mixed_posteriors(fit, parameters = c("mu_x", "mu_z"), conditional = c("mu_x", "mu_z")),
    class = "BayesTools_R2_gated_unavailable")

  # mixed posteriors carry the prior density with the atom (or the continuous
  # part when the quantity's own inclusion is conditioned)
  total <- as_mixed_posteriors(fit, parameters = "mu_R2")$mu_R2
  expect_equal(.bt_meta_get(total, "prior_density")$points$p, 0.5)
  conditioned <- as_mixed_posteriors(fit, parameters = "mu_R2", conditional = "mu_R2")$mu_R2
  expect_true(all(conditioned > 0))
  expect_equal(nrow(.bt_meta_get(conditioned, "prior_density")$points), 0L)
  marginal <- marginal_posterior(as_mixed_posteriors(fit, parameters = "mu_R2"), parameter = "mu_R2",
    use_formula = FALSE, prior_samples = TRUE)
  expect_equal(.bt_meta_get(marginal, "prior_density")$points$p, 0.5)
  # under 'OR' with an independent prior the budget keeps zero mass 1/3
  either <- as_mixed_posteriors(fit, parameters = c("mu_R2", "theta"), conditional = c("mu_R2", "theta"), conditional_rule = "OR")$mu_R2
  expect_equal(.bt_meta_get(either, "prior_density")$points$p, 1 / 3, tolerance = 1e-10)
  # a linked slab is a mixture of Beta totals: flattened beside the atom
  linked <- .R2_gate_bind(.R2_gate_budget(convention = "linked", total = .R2_spike(0.5, prior("beta", list(3.5, 2)))))$mu_R2
  law <- .bt_R2_total_law_given(linked, 2 / 3)
  expect_true(is.prior.mixture(law))
  weights <- vapply(law, `[[`, numeric(1), "prior_weights")
  atom <- vapply(law, is.prior.point, logical(1))
  expect_true(sum(!atom) > 1)
  expect_equal(sum(weights[atom]) / sum(weights), 1 / 3, tolerance = 1e-10)
})

test_that("draws of gated budgets are extracted beyond the analytic limits and for two-level factors", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  set.seed(9)
  # a two-level orthonormal factor has one coded coefficient
  data <- data.frame(g = factor(rep(c("a", "b"), each = 15)), x = stats::rnorm(30))
  data$y <- 0.5 * (data$g == "b") + stats::rnorm(30)
  fit <- JAGS_fit(
    model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
    data = list(y = data$y, N = nrow(data)), prior_list = list(sigma = prior("invgamma", list(1, 1))),
    formula_list = list(mu = ~ x + g), formula_data_list = list(mu = data),
    formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 5)),
      R2 = prior_R2(terms = c("x", "g"), total = prior("beta", list(1, 2)), law = "exact", convention = "common",
        concentration = "inherited", sigma = "sigma", contrast = "orthonormal", inclusion = prior("spike", list(0.5))))),
    chains = 1, adapt = 100, burnin = 100, sample = 300, seed = 9, silent = TRUE)
  expect_s3_class(as_mixed_posteriors(fit, parameters = "mu_g", conditional = "mu_g")$mu_g, "mixed_posteriors")

  # 17 gated terms: the fit, its catalog and the extraction of draws need no
  # enumeration of the configurations
  wide <- as.data.frame(matrix(stats::rnorm(40 * 17), 40, 17, dimnames = list(NULL, paste0("v", 1:17))))
  wide$y <- stats::rnorm(40)
  many <- prior_R2(terms = paste0("v", 1:17), total = prior("beta", list(1, 1)), law = "exact",
    convention = "common", concentration = "inherited", sigma = "sigma", inclusion = prior("spike", list(0.5)))
  fit <- JAGS_fit(
    model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
    data = list(y = wide$y, N = nrow(wide)), prior_list = list(sigma = prior("invgamma", list(1, 1))),
    formula_list = list(mu = stats::reformulate(paste0("v", 1:17))), formula_data_list = list(mu = wide),
    formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 5)), R2 = many)),
    chains = 1, adapt = 100, burnin = 100, sample = 200, seed = 9, silent = TRUE)
  mixed <- as_mixed_posteriors(fit, parameters = c("mu_v1", "mu_R2"))
  expect_length(mixed$mu_v1, 200)
  expect_null(.bt_meta_get(mixed$mu_v1, "prior_density"))
  expect_true(inherits(runjags_inference_table(fit), "BayesTools_table"))
})
