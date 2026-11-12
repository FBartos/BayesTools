skip_if_not_test_profile("unit")

# Inclusion gates and spike-and-slab totals of R-squared budgets.

.R2_gate_data <- function(){
  set.seed(7)
  data.frame(f = factor(rep(c("a", "b", "c"), each = 10)), x = stats::rnorm(30), z = stats::rnorm(30))
}
.R2_gate_budget <- function(inclusion = prior("beta", list(1, 1)), total = prior("beta", list(1.5, 2)),
                            convention = "common", law = "exact", concentration = "inherited"){
  .R2_prior(terms = c("x", "z", "f"), total = total, allocation = prior("dirichlet", list(alpha = c(1, 0.5, 2))),
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
  expect_error(.R2_gate_budget(inclusion = list(w = prior("spike", list(0.5)))), "not terms of the budget: 'w'", fixed = TRUE)
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
  states <- .bt_R2_core_states(core)
  expect_equal(mean(rowSums(open) == 0), states$off + states$empty, tolerance = 0.02)
  expect_equal(mean(draws[, "mu_R2"] == 0), states$off + states$empty, tolerance = 0.02)
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

  # the atom of the total is the budget null: its evidence is the inclusion
  # Bayes factor (the region excluding zero gives it as well)
  expect_error(hypothesis_BF_model(fit, "R2{total} = 0"), "inclusion Bayes factor of runjags_inference_table()",
    fixed = TRUE, class = "BayesTools_point_mass_at_null")
  # a bound of a share that is no atom has no exact posterior ordinate
  expect_error(hypothesis_BF_model(fit, "R2[share: x] = 0"), "is a bound of its support", fixed = TRUE)
  # the atom of a share at one is not an inclusion event
  expect_error(hypothesis_BF_model(fit, "R2[share: x] = 1"), "the region hypothesis '(mu) R2[share: x] >= 1'",
    fixed = TRUE, class = "BayesTools_point_mass_at_null")
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
  # factor levels: the atom of the closed term and the law of each level
  level <- parameter_prior_density(fit, parameter_catalog_resolve(catalog, alias = "f[a]"))
  expect_equal(level$points$p, 1 - 0.5 * 0.5 / 0.75)
  expect_s3_class(hypothesis_BF_model(fit, "f[a] > 0.2"), "data.frame")

  # the shared inclusion rate: Beta(1 + s, 1 + 3 - s) given s open gates, its
  # prior the Beta(1, 1) rate given a nonempty configuration
  rate <- parameter_catalog_resolve(catalog, alias = "R2{inclusion rate}")
  expect_equal(table["(mu) R2{inclusion rate}", "Mean"], mean(post[, "mu_R2_rate"]))
  s <- rowSums(post[, c("mu_R2_gate[1]", "mu_R2_gate[2]", "mu_R2_gate[3]")])
  expect_equal(mean(post[, "mu_R2_rate"]), mean((1 + s) / 5), tolerance = 0.02)
  rate_density <- parameter_prior_density(fit, rate)
  at <- rate_density$density$x[c(200, 600)]
  expect_equal(rate_density$density$y[c(200, 600)],
    ((1 - (1 - at)^3)) / (1 - 1 / 4), tolerance = 1e-3)

  # predictions use the realized coefficients, zero for closed terms
  draws <- .extract_posterior_samples(fit, as_list = FALSE)
  f_closed <- draws[, "mu_R2_indicator"] * draws[, "mu_R2_gate[3]"] == 0
  expect_true(any(f_closed & draws[, "mu_x"] != 0) && any(draws[, "mu_R2_indicator"] == 0))
  newdata <- data[c(1, 11, 21), ]
  prediction <- JAGS_predict_formula(fit, parameter = "mu", data = newdata)$mean
  expected <- t(draws[, "mu_intercept"] + draws[, "mu_x"] %o% newdata$x + draws[, "mu_z"] %o% newdata$z)
  expect_equal(unname(prediction[, f_closed]), unname(expected[, f_closed]), tolerance = 1e-10)

  expect_error(do.call(JAGS_bridgesampling, c(list(fit = fit, log_posterior = function(parameters, data) 0), args[-1])),
    class = "BayesTools_R2_gated_unavailable")
})

test_that("gated budgets are exact at the edges of their laws", {

  set.seed(3)
  # fixed_total with one retained unit: the linked total is Beta(total, b)
  data <- .R2_gate_data()
  budget <- .R2_prior(terms = c("x", "z"), total = .R2_spike(0.6, prior("beta", list(2, 2))), law = "exact",
    convention = "linked", concentration = "fixed_total", sigma = 1, inclusion = prior("spike", list(0.5)))
  core <- JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)), R2 = budget))$prior_list$mu_R2
  draws <- .bt_R2_core_draws(core, 40000)
  on <- draws[, "mu_R2"] > 0
  expect_gt(stats::ks.test(draws[on, "mu_R2"], stats::pbeta, 2, 2)$p.value, 1e-3)
  # (one gate and one unit are scalar nodes)
  expect_match(.bt_R2_core_syntax(core), "prior_par_eta_mu_R2_allocation ~ dgamma(ifelse(mu_R2_indicator * mu_R2_gate > 0, 2 / (", fixed = TRUE)

  # tiny inclusion probabilities keep the truncated law (no cancellation)
  tiny <- .bt_R2_gate_truncated_probabilities(list(type = "fixed", pi = c(1e-20, 1e-20)), 2L)
  expect_equal(tiny, c(0.5, 1))

  # supports need no configurations; the laws of 17 gated terms with equal
  # concentrations have 17 open states
  wide <- as.data.frame(matrix(stats::rnorm(40 * 17), 40, 17, dimnames = list(NULL, paste0("v", 1:17))))
  many <- .R2_prior(terms = names(wide), total = .R2_spike(0.5, prior("beta", list(1, 1))), law = "exact",
    convention = "common", concentration = "inherited", sigma = 1, inclusion = prior("spike", list(0.5)))
  priors <- JAGS_formula(stats::reformulate(names(wide)), "mu", wide, list(intercept = prior("normal", list(0, 1)), R2 = many))$prior_list
  support <- .bt_R2_quantity_support(structure(list(), prior_list = priors), list(type = "coordinate", dependencies = "mu_R2"))
  expect_identical(support$points, 0)
  expect_length(.bt_R2_core_states(priors$mu_R2)$prob, 17L)
  expect_equal(mcdf(priors$mu_R2, 0), 0.5)

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
  plain <- .R2_prior(terms = c("x", "f"), total = prior("beta", list(1, 2)), law = "exact", convention = "common",
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
      R2 = .R2_prior(terms = c("x", "z", "f"), total = .R2_spike(0.5, prior("beta", list(1, 2))), law = "exact",
        convention = "common", concentration = "inherited", sigma = "sigma", contrast = "orthonormal",
        inclusion = prior("beta", list(1, 1))))),
    chains = 1, adapt = 200, burnin = 200, sample = 1000, seed = 8, silent = TRUE)

  # a gated coefficient on the original scale of its predictor: its fitted
  # law and draws rescaled, the atom of its gate at zero
  scale <- 1 / stats::sd(data$x)
  post <- as.matrix(.fit_to_posterior(fit))
  original <- hypothesis_BF_model(fit, "x > 0.1", transform_scaled = TRUE)
  prior_prob <- mccdf(attr(fit, "prior_list")$mu_x, 0.1 / scale)
  post_prob <- mean(post[, "mu_x"] * scale > 0.1)
  expect_equal(as.numeric(original$BF), (post_prob / (1 - post_prob)) / (prior_prob / (1 - prior_prob)), tolerance = 1e-4)
  expect_s3_class(hypothesis_BF_model(fit, "x > 0.1"), "data.frame")

  # the prior side of an event across independent priors; the terms of one
  # budget form one joint family with their exact joint law
  priors <- attr(fit, "prior_list")
  open <- .bt_R2_gate_marginals(priors$mu_R2$gates, priors$mu_R2$terms)$open[["x"]]
  options <- .condition_event_model_options(priors, .condition_event(priors, c("mu_x", "theta"), "OR"))
  expect_equal(options$event_probability, 1 - (1 - open) * 0.5)
  # shared Beta(1, 1) rate, three terms, truncated under the spike 0.5:
  # P(x and z open) = 0.5 * B(3, 1) / B(1, 1) / (1 - B(1, 4) / B(1, 1))
  both <- .condition_event_model_options(priors, .condition_event(priors, c("mu_x", "mu_z"), "AND"))
  expect_equal(both$event_probability, 0.5 * (1 / 3) / (3 / 4), tolerance = 1e-12)
  joint <- as_mixed_posteriors(fit, parameters = c("mu_x", "mu_z"), conditional = c("mu_x", "mu_z"))
  expect_equal(length(joint$mu_x), sum(post[, "mu_R2_indicator"] * post[, "mu_R2_gate[1]"] * post[, "mu_R2_gate[2]"]))
  # an unlabelled term of the budget follows the event: given x open, z is
  # open with probability B(3, 1) / B(2, 1) = 2/3 (the truncation cancels)
  given_x <- as_mixed_posteriors(fit, parameters = c("mu_x", "mu_z"), conditional = "mu_x")$mu_z
  expect_equal(.bt_meta_get(given_x, "prior_density")$points$p, 1 / 3, tolerance = 1e-10)

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
  linked <- .R2_gate_bind(.R2_gate_budget(convention = "linked", total = .R2_spike(0.5, prior("beta", list(3.5, 2)))))
  linked$theta <- prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("spike", list(0.5)))
  law <- .bt_R2_conditioned_law(linked, .condition_event(linked, c("mu_R2", "theta"), "OR"), "mu_R2")
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
      R2 = .R2_prior(terms = c("x", "g"), total = prior("beta", list(1, 2)), law = "exact", convention = "common",
        concentration = "inherited", sigma = "sigma", contrast = "orthonormal", inclusion = prior("spike", list(0.5))))),
    chains = 1, adapt = 100, burnin = 100, sample = 300, seed = 9, silent = TRUE)
  expect_s3_class(as_mixed_posteriors(fit, parameters = "mu_g", conditional = "mu_g")$mu_g, "mixed_posteriors")

  # 17 gated terms: the fit, its catalog and the extraction of draws need no
  # enumeration of the configurations
  wide <- as.data.frame(matrix(stats::rnorm(40 * 17), 40, 17, dimnames = list(NULL, paste0("v", 1:17))))
  wide$y <- stats::rnorm(40)
  many <- .R2_prior(terms = paste0("v", 1:17), total = prior("beta", list(1, 1)), law = "exact",
    convention = "common", concentration = "inherited", sigma = "sigma", inclusion = prior("spike", list(0.5)))
  fit <- JAGS_fit(
    model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
    data = list(y = wide$y, N = nrow(wide)), prior_list = list(sigma = prior("invgamma", list(1, 1))),
    formula_list = list(mu = stats::reformulate(paste0("v", 1:17))), formula_data_list = list(mu = wide),
    formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 5)), R2 = many)),
    chains = 1, adapt = 100, burnin = 100, sample = 200, seed = 9, silent = TRUE)
  mixed <- as_mixed_posteriors(fit, parameters = c("mu_v1", "mu_R2"))
  expect_length(mixed$mu_v1, 200)
  expect_equal(.bt_meta_get(mixed$mu_v1, "prior_density")$points$p, 0.5)
  expect_true(inherits(runjags_inference_table(fit), "BayesTools_table"))
})

test_that("the prior laws of gated budgets aggregate the configurations exactly", {

  # reference: every configuration of the gates with its Beta total and the
  # functional of the term in that model
  reference <- function(core, term, q_total, q_term, prior){
    K <- core$gates$K
    grid <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), K)))
    p_gate <- .bt_R2_gate_configuration_probability(core$gates, grid)
    on <- if(is.null(core$gates$spike)) 1 else core$gates$spike$p
    ungated <- core$terms[core$term_gate == 0L]
    total <- (1 - on) * (q_total >= 0)
    coefficient <- (1 - on) * (q_term >= 0)
    for(i in seq_len(nrow(grid))){
      open <- c(ungated, core$gates$gated[grid[i, ]])
      prob <- on * p_gate[[i]]
      if(length(open) == 0L){
        total <- total + prob * (q_total >= 0)
        coefficient <- coefficient + prob * (q_term >= 0)
        next
      }
      units <- core$unit_term %in% open
      alpha <- core$unit_alpha
      if(identical(core$concentration, "fixed_total")) alpha[units] <- alpha[units] * core$fixed_total / sum(alpha[units])
      A <- sum(alpha[units])
      A_ref <- if(identical(core$concentration, "fixed_total")) core$fixed_total else core$A_reference
      shapes <- switch(core$convention, common = c(core$a, core$b), linked = c(A, core$b),
                       conditioning = c(core$a - (A_ref - A), core$b + (A_ref - A)))
      total <- total + prob * stats::pbeta(q_total, shapes[1], shapes[2])
      if(term %in% open){
        model <- prior
        model$core <- NULL
        model$shapes <- list(a = shapes[1], b = shapes[2], alpha = alpha, A = A)
        model$n_units <- sum(units)
        coefficient <- coefficient + prob * .bt_R2_functional_cdf(.bt_R2_term_functional(model, 1), q_term)
      }else{
        coefficient <- coefficient + prob * (q_term >= 0)
      }
    }
    list(total = total, coefficient = coefficient)
  }
  q_total <- c(0, 0.2, 0.7)
  q_term <- c(-0.5, 0, 0.8)
  designs <- list(
    list(inclusion = prior("beta", list(2, 3)), total = .R2_spike(0.6), convention = "common", concentration = "inherited"),
    list(inclusion = list(x = prior("spike", list(0.3)), f = prior("spike", list(0.6))), total = .R2_spike(0.7),
         convention = "conditioning", concentration = "inherited"),
    list(inclusion = prior("spike", list(0.4)), total = prior("beta", list(3.5, 2)), convention = "linked", concentration = "fixed_total")
  )
  for(design in designs){
    budget <- .R2_gate_budget(inclusion = design$inclusion, convention = design$convention,
      concentration = design$concentration, total = if(identical(design$convention, "conditioning")) .R2_spike(0.7, prior("beta", list(4, 2))) else design$total)
    priors <- .R2_gate_bind(budget)
    expected <- reference(priors$mu_R2, "x", q_total, q_term, priors$mu_x)
    expect_equal(mcdf(priors$mu_R2, q_total), expected$total, tolerance = 1e-12)
    expect_equal(mcdf(priors$mu_x, q_term), expected$coefficient, tolerance = 1e-10)
  }

  # 30 gated terms with equal concentrations: 30 open states (one per number
  # of open gates), and the laws are available
  set.seed(1)
  wide <- as.data.frame(matrix(stats::rnorm(40 * 30), 40, 30, dimnames = list(NULL, paste0("v", 1:30))))
  expect_error(.R2_prior(terms = names(wide), total = .R2_spike(0.5, prior("beta", list(1, 1))), law = "exact",
    convention = "linked", concentration = "inherited", sigma = 1, inclusion = prior("beta", list(1, 1))),
    "must equal the total concentration of the full model, 30", fixed = TRUE)
  many <- .R2_prior(terms = names(wide), total = .R2_spike(0.5, prior("beta", list(30, 1))), law = "exact",
    convention = "linked", concentration = "inherited", sigma = 1, inclusion = prior("beta", list(1, 1)))
  priors <- JAGS_formula(stats::reformulate(names(wide)), "mu", wide, list(intercept = prior("normal", list(0, 1)), R2 = many))$prior_list
  states <- .bt_R2_core_states(priors$mu_R2)
  expect_length(states$prob, 30L)
  expect_equal(sum(states$prob), 0.5)
  # linked: given s open terms the total is Beta(s, 1), and s is uniform on
  # 1..30 under the truncated beta-binomial law
  expect_equal(mean(priors$mu_R2), 0.5 * mean((1:30) / (1:30 + 1)), tolerance = 1e-12)
  expect_true(is.finite(mcdf(priors$mu_v1, 0.1)))

  # many distinct concentrations: the number of states is limited
  set.seed(2)
  distinct <- .R2_prior(terms = names(wide)[1:17], total = prior("beta", list(1, 1)),
    allocation = prior("dirichlet", list(alpha = sqrt(c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59)))),
    law = "exact", convention = "common", concentration = "inherited", sigma = 1, inclusion = prior("spike", list(0.5)))
  priors <- JAGS_formula(stats::reformulate(names(wide)[1:17]), "mu", wide, list(intercept = prior("normal", list(0, 1)), R2 = distinct))$prior_list
  expect_error(mcdf(priors$mu_R2, 0.5), "distinct states", class = "BayesTools_R2_density_unavailable")
})

test_that("factor levels, named allocations and the last truncated gate", {

  # ungated factor levels: the catalog law of a level is its R2 law
  set.seed(3)
  data <- .R2_gate_data()
  priors <- JAGS_formula(~ x + f, "mu", data, list(intercept = prior("normal", list(0, 1)),
    R2 = .R2_prior(terms = c("x", "f"), total = prior("beta", list(1, 2)), law = "exact", convention = "common",
      concentration = "inherited", sigma = 1, contrast = "orthonormal")))$prior_list
  key <- list(type = "factor_level", dependencies = c("mu_f[1]", "mu_f[2]"), weights = contr.orthonormal(3)[2, ])
  law <- .bt_R2_quantity_law(structure(list(), prior_list = priors), key)
  expect_equal(mcdf(law, c(-0.3, 0.4)), mcdf(priors$mu_f, c(-0.3, 0.4)), tolerance = 1e-10)
  expect_equal(var(law), var(priors$mu_f), tolerance = 1e-10)

  # named concentrations follow the terms
  named <- .R2_prior(terms = c("x", "f"), total = prior("beta", list(1, 2)),
    allocation = prior("dirichlet", list(alpha = c(f = 5, x = 1))), law = "exact", convention = "common",
    concentration = "inherited", sigma = 1, contrast = "orthonormal")
  expect_equal(named$alpha, c(x = 1, f = 5))
  expect_error(.R2_prior(terms = c("x", "f"), total = prior("beta", list(1, 2)),
    allocation = prior("dirichlet", list(alpha = c(f = 5, z = 1))), law = "exact", convention = "common",
    concentration = "inherited", sigma = 1), "each once", fixed = TRUE)

  # given every earlier gate closed, the last gate opens with probability one
  gates <- .R2_gate_bind(.R2_gate_budget(inclusion = prior("beta", list(2, 3)), total = .R2_spike(0.5)))$mu_R2$gates
  expect_identical(gates$p_trunc[[3L]], 1)
})

test_that("gated laws keep distinct states, conditioned slabs and dependent members apart", {

  set.seed(5)
  data <- data.frame(x = stats::rnorm(30), z = stats::rnorm(30), w = stats::rnorm(30))
  bind <- function(budget, formula = ~ x + z){
    JAGS_formula(formula, "mu", data, list(intercept = prior("normal", list(0, 1)), R2 = budget))$prior_list
  }

  # close raw concentrations with different shapes are different states
  close <- .R2_prior(terms = c("x", "z"), total = prior("beta", list(1e14 + 2, 2)),
    allocation = prior("dirichlet", list(alpha = c(1e14, 1e14 + 1))), law = "exact", convention = "conditioning",
    concentration = "inherited", sigma = 1, inclusion = prior("spike", list(0.5)))
  q <- 1e-14
  expected <- 0.25 + 0.25 * (stats::pbeta(q, 1, 1e14 + 3) + stats::pbeta(q, 2, 1e14 + 2) + stats::pbeta(q, 1e14 + 2, 2))
  expect_equal(mcdf(bind(close)$mu_R2, q), expected, tolerance = 1e-10)

  # a core conditioned on a pattern: given x and z open (flat allocation,
  # common Beta(2, 4) total) x has the variance E[W] / 2 * scale^2 of the
  # model x + z, half of that given z closed, and its draws follow the pattern
  plain <- .R2_prior(terms = c("x", "z"), total = prior("beta", list(2, 4)), law = "exact", convention = "common",
    concentration = "inherited", sigma = 1, inclusion = prior("spike", list(0.5)))
  conditioned <- function(prior, flags){
    prior$core <- .bt_R2_condition_core(prior$core, flags)
    prior
  }
  both <- conditioned(bind(plain)$mu_x, c(x = TRUE, z = TRUE))
  draws <- rng(both, 40000)
  expect_true(all(draws != 0))
  expect_equal(stats::var(draws), var(both), tolerance = 0.05)
  alone <- conditioned(bind(plain)$mu_x, c(x = TRUE, z = FALSE))
  expect_equal(var(alone), 2 * var(both), tolerance = 1e-8)
  # given z closed, the conditioning total of x alone is Beta(1, 1.25): mean 0
  small <- .R2_prior(terms = c("x", "z"), total = prior("beta", list(2, 0.25)), law = "exact", convention = "conditioning",
    concentration = "inherited", sigma = 1, inclusion = prior("spike", list(0.5)))
  expect_identical(mean(bind(small)$mu_x), NaN)
  expect_identical(mean(conditioned(bind(small)$mu_x, c(z = FALSE))), 0)

  # the joint prior draws of a conditioned prior list follow the pattern
  priors <- bind(plain)
  options <- .condition_event_model_options(priors, .condition_event(priors, c("mu_x", "mu_z"), "AND"))
  expect_equal(options$event_probability, 0.25)
  samples <- .generate_prior_sample_matrix(options$prior_lists[[1L]], 2000, column_names = c("mu_x", "mu_z", "mu_R2"), seed = 1)
  expect_true(all(samples[, "mu_x"] != 0 & samples[, "mu_z"] != 0 & samples[, "mu_R2"] > 0))

  # a closed factor keeps its coded coefficients and a closed budget samples
  # no total, in the joint draws of the options of an event
  ff <- data.frame(f = factor(rep(c("a", "b", "c"), each = 10)), x = stats::rnorm(30))
  gated_f <- .R2_prior(terms = c("x", "f"), total = prior("beta", list(2, 4)), law = "exact", convention = "common",
    concentration = "inherited", sigma = 1, contrast = "orthonormal",
    inclusion = list(x = prior("spike", list(0.5)), f = prior("spike", list(0.5))))
  priors_f <- JAGS_formula(~ x + f, "mu", ff, list(intercept = prior("normal", list(0, 1)), R2 = gated_f))$prior_list
  either <- .condition_event_model_options(priors_f, .condition_event(priors_f, c("mu_x", "mu_f"), "OR"))
  for(option in either$prior_lists){
    sampled <- .generate_prior_sample_matrix(option, 100, column_names = c("mu_x", "mu_f[1]", "mu_f[2]", "mu_R2"), seed = 1)
    expect_false(anyNA(sampled))
  }
  single <- .R2_prior(terms = "x", total = prior("beta", list(2, 4)), law = "exact", convention = "common",
    concentration = "inherited", sigma = 1, inclusion = list(x = prior("spike", list(0.5))))
  priors_x <- JAGS_formula(~ x, "mu", ff, list(intercept = prior("normal", list(0, 1)), R2 = single))$prior_list
  priors_x$theta <- prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("spike", list(0.5)))
  closed <- .condition_event_model_options(priors_x, .condition_event(priors_x, c("mu_x", "theta"), "OR"))$prior_lists[[1L]]
  expect_true(all(.generate_prior_sample_matrix(closed, 100, column_names = c("mu_x", "mu_R2", "theta"), seed = 1)[, "mu_R2"] == 0))

  # the shape of the inclusion rate keeps a tiny second shape
  tiny <- .R2_prior(terms = c("x", "z"), total = prior("beta", list(2, 4)), law = "exact", convention = "common",
    concentration = "inherited", sigma = 1, inclusion = prior("beta", list(0.1, 1e-6)))
  expect_match(.bt_R2_core_syntax(bind(tiny)$mu_R2), "1e-06 + (2 - sum(mu_R2_gate[1:2])))", fixed = TRUE)

  # event probabilities need no concentrations (18 distinct concentrations
  # exceed the states of the laws, not those of the gates)
  wide <- as.data.frame(matrix(stats::rnorm(30 * 18), 30, 18, dimnames = list(NULL, paste0("v", 1:18))))
  distinct <- .R2_prior(terms = names(wide), total = prior("beta", list(1, 1)),
    allocation = prior("dirichlet", list(alpha = 2^(0:17))), law = "exact", convention = "common",
    concentration = "inherited", sigma = 1, inclusion = prior("spike", list(0.5)))
  many <- JAGS_formula(stats::reformulate(names(wide)), "mu", wide, list(intercept = prior("normal", list(0, 1)), R2 = distinct))$prior_list
  expect_equal(.condition_event_model_options(many, .condition_event(many, "mu_v1", "AND"))$event_probability, 0.5)

  # given z open, the other members of a shared-rate budget are dependent:
  # combinations of them have no product law
  shared <- .R2_prior(terms = c("x", "z", "w"), total = prior("beta", list(2, 2)), law = "exact", convention = "common",
    concentration = "inherited", sigma = 1, inclusion = prior("beta", list(1, 1)))
  priors <- bind(shared, ~ x + z + w)
  context <- .prior_density_build_context(priors[c("mu_x", "mu_z", "mu_w", "mu_R2")], c("mu_x", "mu_z", "mu_w"),
    conditional = "mu_z")
  expect_error(.prior_density_from_context(context, c(mu_x = 1, mu_z = 0, mu_w = 1)),
    class = "BayesTools_formula_measure_unavailable")
})

test_that("gate shapes keep tiny rates and conditioned supports drop impossible atoms", {

  set.seed(2)
  data <- data.frame(x = stats::rnorm(30))
  tiny <- prior_R2(terms = "x", total = prior("beta", list(2, 2)), inclusion = prior("beta", list(2^-53, 2^-52)))
  core <- JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)), R2 = tiny))$prior_list$mu_R2
  expect_equal(mcdf(core, 0), 2 / 3)
  expect_equal(mean(core), 1 / 6)
  on <- .bt_R2_condition_core(core, flags = logical(), total = TRUE)
  support <- .bt_R2_quantity_support(structure(list(), prior_list = list(mu_R2 = on)),
    list(type = "coordinate", dependencies = "mu_R2"))
  expect_length(support$points, 0L)
  support <- .bt_R2_quantity_support(structure(list(), prior_list = list(mu_R2 = core)),
    list(type = "coordinate", dependencies = "mu_R2"))
  expect_identical(support$points, 0)
})

test_that("gate probabilities stay exact for tiny probabilities and huge rate shapes", {

  set.seed(3)
  data <- data.frame(x = stats::rnorm(20), z = stats::rnorm(20))
  bind <- function(budget, formula = ~ x + z){
    JAGS_formula(formula, "mu", data, list(intercept = prior("normal", list(0, 1)), R2 = budget))$prior_list$mu_R2
  }
  # two gates with probability p under a spike-and-slab total with weight 1/2:
  # each term is open with probability p / (1 - (1 - p)^2) / 2 -> 1/4
  for(p in c(1e-16, 2^-54, 1e-20)){
    core <- bind(prior_R2(total = .R2_spike(0.5), inclusion = prior("spike", list(p))))
    marginals <- .bt_R2_gate_marginals(core$gates, core$terms)
    expect_equal(unname(marginals$open), c(1, 1) / (4 - 2 * p), tolerance = 1e-14)
    expect_equal(mcdf(core, c(0, 1)), c(0.5, 1))
  }
  # a shared rate concentrated at 1/2: one gate open with probability 1/2
  core <- bind(prior_R2(total = .R2_spike(0.5), inclusion = prior("beta", list(1e16, 1e16))), ~ x)
  expect_equal(.bt_R2_gate_configuration_probability(utils::modifyList(core$gates, list(truncated = FALSE)),
    matrix(FALSE, 1, 1)), 0.5)
  # a shared rate near zero, Beta(1, 1e20): each singleton of two gates has
  # probability about 1e-20, and 1/2 under the slab of a spiked total
  gates <- list(K = 2L, law = list(type = "shared", a = 1, b = 1e20), truncated = FALSE)
  singleton <- .bt_R2_gate_configuration_probability(gates, matrix(c(TRUE, FALSE), 1L))
  expect_equal(singleton / 1e-20, 1, tolerance = 1e-12)
  core <- bind(prior_R2(total = .R2_spike(0.5), inclusion = prior("beta", list(1, 1e20))))
  expect_equal(unname(.bt_R2_gate_marginals(core$gates, core$terms)$open), c(0.25, 0.25), tolerance = 1e-12)
  rate <- .bt_R2_core_rate_law(core)
  expect_true(all(is.finite(vapply(rate, function(component) component$prior_weights, numeric(1)))))
  # the product form equals the beta-binomial law at ordinary shapes
  gates <- list(K = 4L, law = list(type = "shared", a = 2, b = 3), truncated = FALSE)
  open <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), 4)))
  expect_equal(.bt_R2_gate_configuration_probability(gates, open),
    exp(lbeta(2 + rowSums(open), 3 + 4 - rowSums(open)) - lbeta(2, 3)), tolerance = 1e-14)
  expect_equal(.bt_R2_gate_truncated_probabilities(gates$law, 4L)[1:3],
    vapply(1:3, function(t) (2 / (5 + t - 1)) / (1 - exp(lbeta(2, 3 + 4) - lbeta(2, 3 + t - 1))), numeric(1)),
    tolerance = 1e-14)
})

test_that("gated totals and coefficients keep their quantiles near zero and their variance", {

  data <- data.frame(x = rep(c(-1, 1), 8), z = rep(c(-1, -1, 1, 1), 4))
  bind <- function(budget, formula = ~ x){
    JAGS_formula(formula, "mu", data, list(intercept = prior("normal", list(0, 1)), R2 = budget))$prior_list
  }
  # one gate with probability 1/2 and a Beta(0.01, 1) total: F(r) = 1/2 + r^0.01 / 2
  priors <- bind(prior_R2(total = prior("beta", list(0.01, 1)), inclusion = prior("spike", list(0.5)),
    model = R2_model(sigma = 1)))
  expect_equal(log(mquant(priors$mu_R2, 0.75)), -100 * log(2), tolerance = 1e-10)
  # the coefficient (unit design variance): P(X <= q) = 3/4 + q^0.02 / 4 for small q > 0
  q <- mquant(priors$mu_x, 0.8)
  expect_equal(log(q), 50 * log(0.2), tolerance = 1e-8)
  expect_equal(mcdf(priors$mu_x, q), 0.8, tolerance = 1e-10)
  expect_equal(log(-mquant(priors$mu_x, 0.2)), 50 * log(0.2), tolerance = 1e-8)
  # a gated term beside an ungated term leaves the total a Beta(1e9, 2)
  priors <- bind(prior_R2(total = prior("beta", list(1e9, 2)), inclusion = list(z = prior("spike", list(0.5)))), ~ x + z)
  expect_gt(var(priors$mu_R2), 0)
  concentrated <- bind(prior_R2(total = prior("beta", list(1e105, 1e105)),
    inclusion = list(z = prior("spike", list(0.5)))), ~ x + z)
  expect_equal(var(concentrated$mu_R2) / (0.25 / (2e105 + 1)), 1, tolerance = 1e-10)
  # three gated terms beside an ungated one: the merged states keep one Beta
  # component whatever the rounding of its weight
  three <- JAGS_formula(~ x + z1 + z2 + z3, "mu", data.frame(x = rep(c(-1, 1), 8), z1 = rep(c(-1, -1, 1, 1), 4),
      z2 = rep(c(1, -1, -1, 1), 4), z3 = rep(c(-1, -1, -1, -1, 1, 1, 1, 1), 2)),
    list(intercept = prior("normal", list(0, 1)), R2 = prior_R2(total = prior("beta", list(1e105, 1e105)),
      inclusion = list(z1 = prior("spike", list(0.1)), z2 = prior("spike", list(0.1)), z3 = prior("spike", list(0.1))))))$prior_list
  expect_equal(var(three$mu_R2) / (0.25 / (2e105 + 1)), 1, tolerance = 1e-10)
  # shapes whose sum is not a double are refused
  expect_error(prior_R2(inclusion = prior("beta", list(1e308, 1e308))), class = "BayesTools_R2_invalid")
  expect_error(prior_R2(total = prior("beta", list(1e308, 1e308))), class = "BayesTools_R2_invalid")
  expect_equal(.bt_R2_log_fraction(c(1, 1e300, 1e-300), c(1, 1e300, 1)), c(log(0.5), log(0.5), log(1e-300) - log1p(1e-300)))
  # a continuous part of 1e-20 keeps the upper end of the total
  tiny <- bind(prior_R2(inclusion = prior("spike", list(1e-20))))
  expect_identical(mquant(tiny$mu_R2, 1), 1)
  expect_equal(var(priors$mu_R2) / (1e9 * 2 / ((1e9 + 2)^2 * (1e9 + 3))), 1, tolerance = 1e-8)
  # a coefficient scale of 1e20 (design values of 1e-20): Beta(1, 1) total,
  # P(X <= q) = 3/4 + (q^2 / (1 + q^2)) / 4 for q > 0 with unit sigma
  huge <- JAGS_formula(~ x, "mu", data.frame(x = rep(c(-1e-20, 1e-20), 8)),
    list(intercept = prior("normal", list(0, 1)), R2 = prior_R2(inclusion = prior("spike", list(0.5)),
      model = R2_model(sigma = 1))))$prior_list
  expect_equal(mquant(huge$mu_x, 0.8), 1e20 * sqrt(0.2 / 0.8), tolerance = 1e-8)
  # one gate with probability 1/4 and a Beta(2, 3) total (mean 0.4, variance
  # 0.04): 0.25 * 0.04 + 0.25 * 0.75 * 0.4^2
  expect_equal(var(bind(prior_R2(total = prior("beta", list(2, 3)), inclusion = prior("spike", list(0.25))))$mu_R2),
    0.25 * 0.04 + 0.25 * 0.75 * 0.4^2, tolerance = 1e-12)
})

test_that("gated interaction terms are keyed by their prior-list names", {

  set.seed(4)
  data <- data.frame(x = stats::rnorm(40), z = stats::rnorm(40), f = factor(rep(c("a", "b"), 20)))
  data$y <- 0.6 * data$x * data$z + stats::rnorm(40)
  gated <- prior_R2(inclusion = prior("spike", list(0.5)), model = R2_model(sigma = 1))
  priors <- JAGS_formula(~ x * z, "mu", data, list(intercept = prior("normal", list(0, 1)), R2 = gated))$prior_list
  expect_identical(priors$mu_x__xXx__z$term, "x__xXx__z")
  expect_identical(priors$mu_x__xXx__z$declared, "x:z")
  expect_output(print(priors$mu_x__xXx__z), "included with probability 0.5", fixed = TRUE)
  expect_equal(mcdf(priors$mu_x__xXx__z, 0), 0.75)
  factor <- JAGS_formula(~ f * x, "mu", data, list(intercept = prior("normal", list(0, 1)), R2 = gated))$prior_list
  expect_equal(mcdf(factor$mu_f__xXx__x, 0), 0.75)

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  fit <- JAGS_fit(
    model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
    data = list(y = data$y, N = nrow(data)), prior_list = list(sigma = prior("invgamma", list(1, 1))),
    formula_list = list(mu = ~ x * z), formula_data_list = list(mu = data),
    formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 5)), R2 = gated)),
    chains = 1, adapt = 200, burnin = 200, sample = 1000, seed = 4, silent = TRUE)
  post <- as.matrix(.fit_to_posterior(fit))
  table <- runjags_estimates_table(fit)
  expect_true(all(c("(mu) x:z (inclusion)", "(mu) R2[share: x:z]") %in% rownames(table)))
  expect_equal(table["(mu) x:z (inclusion)", "Mean"], mean(post[, "mu_x__xXx__z"] != 0))
  inference <- runjags_inference_table(fit)
  expect_equal(inference["(mu) x:z", "prior_prob"], 0.5)
  expect_equal(inference["(mu) x:z", "post_prob"], mean(post[, "mu_x__xXx__z"] != 0))
  expect_s3_class(hypothesis_BF_model(fit, "R2[share: x:z] > 0.5"), "data.frame")
  conditioned <- as_mixed_posteriors(fit, parameters = "mu_x__xXx__z", conditional = "mu_x__xXx__z")
  expect_true(all(conditioned$mu_x__xXx__z != 0))
})

test_that("mixed posteriors keep the atoms of gated coefficients, budget inclusion and model averages", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  set.seed(11)
  data <- data.frame(x = 3 + 2 * stats::rnorm(40), z = -1 + stats::rnorm(40))
  data$y <- 0.2 * data$x + stats::rnorm(40)
  fit_args <- list(
    model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
    data = list(y = data$y, N = nrow(data)), prior_list = list(sigma = prior("invgamma", list(1, 1))),
    formula_list = list(mu = ~ x + z), formula_data_list = list(mu = data),
    formula_scale_list = list(mu = list(x = TRUE, z = TRUE)),
    chains = 1, adapt = 200, burnin = 200, sample = 2000, silent = TRUE)
  gated <- prior_R2(inclusion = list(x = prior("spike", list(0.5)), z = prior("spike", list(0.5))))
  fit <- do.call(JAGS_fit, c(fit_args, list(seed = 1,
    formula_prior_list = list(mu = list(intercept = prior("spike", list(0)), R2 = gated)))))
  post <- as.matrix(.fit_to_posterior(fit))
  closed_x <- mean(post[, "mu_x"] == 0)
  both_closed <- mean(post[, "mu_x"] == 0 & post[, "mu_z"] == 0)

  # a gated coefficient keeps the atom of its gate on the original scale
  original <- as_mixed_posteriors(fit, "mu_x", transform_scaled = TRUE)$mu_x
  atoms <- .posterior_atoms_get(original)
  expect_equal(unname(atoms$locations[, 1]), 0)
  expect_equal(atoms$mass, closed_x)
  expect_identical(.bt_meta_get(original, "support")$points, 0)

  # a prediction weighting two gated terms is zero when both are closed
  samples <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x", "mu_z"))
  row <- marginal_posterior(samples, "mu_intercept", formula = ~ 1 + x + z, at = list(x = 1, z = 1),
    prior_samples = TRUE)$intercept
  atoms <- .posterior_atoms_get(row)
  expect_equal(unname(atoms$locations[, 1]), 0)
  expect_equal(atoms$mass, both_closed)

  # the inclusion of the budget is an atomic quantity with its own plot
  selection <- parameter_catalog_resolve(parameter_catalog(fit), "R2{total} (inclusion)")
  inclusion <- parameter_mixed_posterior(fit, selection)
  atoms <- .posterior_atoms_get(inclusion)
  expect_equal(unname(atoms$locations[, 1]), c(0, 1))
  expect_equal(atoms$mass, c(both_closed, 1 - both_closed))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot_posterior(list(inclusion = inclusion), "inclusion", prior = TRUE))

  # a model average keeps the within-model exclusions of the gates
  null_fit <- do.call(JAGS_fit, c(fit_args, list(seed = 2,
    formula_prior_list = list(mu = list(intercept = prior("spike", list(0)), x = prior("spike", list(0)),
                                        z = prior("spike", list(0)))))))
  mixed <- mix_posteriors(model_list = list(list(fit = fit, marglik = bridgesampling_object(-50), prior_weights = 1),
                                            list(fit = null_fit, marglik = bridgesampling_object(-50), prior_weights = 1)),
                          parameters = "mu_x", is_null_list = list(mu_x = c(FALSE, TRUE)), seed = 1, n_samples = 4000)
  # (equal marginal likelihoods: each model has posterior probability 1/2;
  # the gated model contributes the closed frequency of its original draws,
  # also when the mixture allocates it no draw)
  atoms <- .posterior_atoms_get(mixed$mu_x)
  expect_true(all(atoms$locations == 0))
  expect_equal(sum(atoms$mass), 0.5 + 0.5 * closed_x)
  single <- mix_posteriors(model_list = list(list(fit = fit, marglik = bridgesampling_object(-50), prior_weights = 1),
                                             list(fit = null_fit, marglik = bridgesampling_object(-50), prior_weights = 1)),
                           parameters = "mu_x", is_null_list = list(mu_x = c(FALSE, TRUE)), seed = 1, n_samples = 1)
  expect_equal(sum(.posterior_atoms_get(single$mu_x)$mass), 0.5 + 0.5 * closed_x)
})

test_that("formula atoms resolve the scale that multiplies gated R2 coefficients", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  set.seed(12)
  data <- data.frame(x = stats::rnorm(40))
  data$y <- 5 + 0.3 * data$x + stats::rnorm(40)
  # the slope is s times the gated R2 coefficient, with s zero in one
  # component of its mixture: the prediction at x = 1 is the intercept 5
  # whenever the gate is closed or s is zero
  fit_with <- function(s_prior, scale = list()){
    JAGS_fit(
      model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1) } }",
      data = list(y = data$y, N = nrow(data)), prior_list = list(s = s_prior),
      formula_list = list(mu = ~ x), formula_data_list = list(mu = data), formula_scale_list = list(mu = scale),
      formula_prior_list = list(mu = list(intercept = prior("spike", list(5)),
        R2 = prior_R2(inclusion = prior("spike", list(0.5)), model = R2_model(sigma = "s")))),
      chains = 1, adapt = 200, burnin = 200, sample = 1000, seed = 12, silent = TRUE)
  }
  fit <- fit_with(prior_mixture(list(prior("spike", list(0)), prior("normal", list(1, 0.2), list(0, Inf))),
                                is_null = c(TRUE, FALSE)))
  post <- as.matrix(.fit_to_posterior(fit))
  samples <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x", "s"))
  row <- marginal_posterior(samples, "mu_intercept", formula = ~ 1 + x, at = list(x = 1), prior_samples = TRUE)$intercept
  atoms <- .posterior_atoms_get(row)
  expect_equal(unname(atoms$locations[, 1]), 5)
  expect_equal(atoms$mass, mean(post[, "mu_x"] == 0 | post[, "s"] == 0))

  # the fitted coefficient itself is not scaled by s: on the original scale of
  # a scaled predictor its atom is the closed gate only
  scaled <- fit_with(prior_mixture(list(prior("spike", list(0)), prior("normal", list(1, 0.2), list(0, Inf))),
                                   is_null = c(TRUE, FALSE)), scale = list(x = TRUE))
  post <- as.matrix(.fit_to_posterior(scaled))
  slope <- as_mixed_posteriors(scaled, "mu_x", transform_scaled = TRUE)$mu_x
  expect_equal(.posterior_atoms_get(slope)$mass, mean(post[, "mu_x"] == 0))

  # a discrete multiplier (its zero state is not resolved by the plan) leaves
  # the prediction without an atom certificate instead of declaring it
  # atom-free
  bernoulli <- fit_with(prior("bernoulli", list(0.5)))
  samples <- as_mixed_posteriors(bernoulli, c("mu_intercept", "mu_x", "s"))
  row <- marginal_posterior(samples, "mu_intercept", formula = ~ 1 + x, at = list(x = 1), prior_samples = TRUE)$intercept
  expect_true("atoms" %in% .bt_meta_get(row, "measure_unavailable")$measure)
})

test_that("a budget with one gated term names its gate and allocation coordinates as scalars", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  set.seed(13)
  data <- data.frame(x = stats::rnorm(40), z = stats::rnorm(40))
  data$y <- 0.4 * data$x + stats::rnorm(40)
  # (JAGS names a vector of length one by its scalar name)
  budget <- prior_R2(inclusion = list(z = prior("spike", list(0.5))), model = R2_model(sigma = 1))
  priors <- JAGS_formula(~ x + z, "mu", data, list(intercept = prior("normal", list(0, 5)), R2 = budget))$prior_list
  expect_identical(priors$mu_z$gating$gate_names, "mu_R2_gate")
  expect_match(.bt_R2_core_syntax(priors$mu_R2), "mu_R2_gate ~ dbern(0.5)", fixed = TRUE)
  fit <- JAGS_fit(
    model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
    data = list(y = data$y, N = nrow(data)), prior_list = list(sigma = prior("invgamma", list(1, 1))),
    formula_list = list(mu = ~ x + z), formula_data_list = list(mu = data),
    formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 5)), R2 = budget)),
    chains = 1, adapt = 200, burnin = 200, sample = 1000, seed = 13, silent = TRUE)
  post <- as.matrix(.fit_to_posterior(fit))
  expect_true(all(c("mu_R2_gate", "mu_R2_open") %in% colnames(post)))
  table <- runjags_estimates_table(fit)
  expect_equal(table["(mu) z (inclusion)", "Mean"], mean(post[, "mu_R2_gate"]))
  z <- as_mixed_posteriors(fit, "mu_z")$mu_z
  expect_equal(.posterior_atoms_get(z)$mass, mean(post[, "mu_z"] == 0))
  expect_s3_class(hypothesis_BF_model(fit, "R2[share: x] > 0.5"), "data.frame")
})
