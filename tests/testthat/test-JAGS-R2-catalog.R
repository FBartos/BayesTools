skip_if_not_test_profile("unit")

# The totals and shares of R-squared budgets in the parameter catalog,
# estimates tables, hypothesis tests and model averaging.

.R2_catalog_fit <- function(formula, law = "exact", terms = c("x", "z", "f"), allocation = NULL, sample = 400){
  set.seed(41)
  data <- data.frame(f = factor(rep(c("a", "b", "c"), each = 10)), x = stats::rnorm(30), z = stats::rnorm(30))
  data$y <- 0.5 * (data$f == "b") + 0.3 * data$x + stats::rnorm(30)
  budget <- .R2_prior(terms = terms, total = prior("beta", list(1.5, 3)), allocation = allocation, law = law,
    convention = "common", concentration = "inherited", sigma = "sigma", contrast = "orthonormal")
  JAGS_fit(
    model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i], 1 / pow(sigma, 2)) } }",
    data = list(y = data$y, N = nrow(data)), prior_list = list(sigma = prior("invgamma", list(1, 1))),
    formula_list = list(mu = formula), formula_data_list = list(mu = data),
    formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 5)), R2 = budget)),
    chains = 1, adapt = 100, burnin = 100, sample = sample, seed = 41, silent = TRUE)
}

test_that("budget totals and term shares are catalog quantities with their Beta laws", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  fit <- .R2_catalog_fit(~ x + z + f, allocation = prior("dirichlet", list(alpha = c(1, 0.5, 2))))
  catalog <- parameter_catalog(fit)
  quantities <- catalog$quantities

  total <- quantities[quantities$canonical_name == "mu_R2", , drop = FALSE]
  expect_identical(total$display_label, "(mu) R2{total}")
  shares <- quantities[startsWith(quantities$canonical_name, "mu_R2_allocation"), , drop = FALSE]
  expect_identical(shares$canonical_name, c("mu_R2_allocation[x]", "mu_R2_allocation[z]", "mu_R2_allocation[f]"))
  expect_identical(shares$display_label, c("(mu) R2[share: x]", "(mu) R2[share: z]", "(mu) R2[share: f]"))
  expect_true(all(shares$status == "derived"))
  for(support in shares$support){
    expect_identical(support$bounds, c(0, 1))
  }
  # the allocation coordinates are backend coordinates only
  expect_false(any(grepl("^mu_R2_allocation\\[[0-9]+\\]$", quantities$canonical_name)))
  expect_false(any(grepl("R2_allocation\\[1\\]", catalog$aliases$alias)))

  # draws: the normalized gamma coordinates of the allocation
  posterior <- as.matrix(.fit_to_posterior(fit))
  eta <- posterior[, paste0("prior_par_eta_mu_R2_allocation[", 1:3, "]")]
  phi <- eta / rowSums(eta)
  selection <- parameter_catalog_resolve(catalog, alias = "R2[share: z]")
  expect_equal(as.numeric(as.matrix(parameter_draws(fit, selection))), unname(phi[, 2]), tolerance = 1e-12)

  # prior densities: the Beta total and the Beta marginal of a share
  share_density <- parameter_prior_density(fit, selection)
  expect_equal(share_density$density$y[c(100, 1000)],
    stats::dbeta(share_density$density$x[c(100, 1000)], 0.5, 3), tolerance = 1e-2)
  total_density <- parameter_prior_density(fit, parameter_catalog_resolve(catalog, alias = "R2{total}"))
  expect_equal(total_density$density$y[c(100, 1000)],
    stats::dbeta(total_density$density$x[c(100, 1000)], 1.5, 3), tolerance = 1e-2)
})

test_that("estimates tables show the budget total and the shares of its terms", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  fit <- .R2_catalog_fit(~ x + f, law = "normal")
  posterior <- as.matrix(.fit_to_posterior(fit))
  for(transform_factors in c(FALSE, TRUE)){
    table <- runjags_estimates_table(fit, transform_factors = transform_factors)
    rows <- rownames(table)
    expect_false(any(grepl("allocation", rows)))
    position <- match("(mu) R2{total}", rows)
    expect_identical(rows[position + 0:2], c("(mu) R2{total}", "(mu) R2[share: x]", "(mu) R2[share: f]"))
    # a kernel share is the sum of the allocations of the term's coded coefficients
    expect_equal(table["(mu) R2[share: f]", "Mean"],
      mean(posterior[, "mu_R2_allocation[2]"] + posterior[, "mu_R2_allocation[3]"]), tolerance = 1e-12)
  }
  raw <- runjags_estimates_table(fit, random_effects_summary = "raw")
  expect_true(all(c("(mu) R2_allocation[1]", "(mu) R2_allocation[3]") %in% rownames(raw)) ||
              any(grepl("R2_allocation\\[1\\]", rownames(raw))))

  # the display filters act on the quantities: the shares follow their total
  without_x <- rownames(runjags_estimates_table(fit, remove_parameters = "mu_x"))
  expect_true(all(c("(mu) R2{total}", "(mu) R2[share: x]", "(mu) R2[share: f]") %in% without_x))
  without_total <- rownames(runjags_estimates_table(fit, remove_parameters = "mu_R2"))
  expect_false(any(grepl("R2", without_total)))

  # a budget that retains one term allocates it the whole budget: no shares
  single <- .R2_catalog_fit(~ x, sample = 200)
  expect_false(any(grepl("share", rownames(runjags_estimates_table(single)))))
  expect_false(any(grepl("share", parameter_catalog(single)$quantities$display_label)))
})

test_that("hypothesis tests of budget totals and shares use their Beta priors", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  fit <- .R2_catalog_fit(~ x + z + f, sample = 1000)
  posterior <- as.matrix(.fit_to_posterior(fit))
  odds <- function(p) p / (1 - p)

  total <- hypothesis_BF_model(fit, "R2{total} > 0.1")
  expected <- odds(mean(posterior[, "mu_R2"] > 0.1)) /
    odds(stats::pbeta(0.1, 1.5, 3, lower.tail = FALSE))
  expect_equal(as.numeric(total$BF), expected, tolerance = 1e-8)
  expect_equal(as.numeric(hypothesis_BF_model(fit, "mu_R2 > 0.1")$BF), as.numeric(total$BF))

  eta <- posterior[, paste0("prior_par_eta_mu_R2_allocation[", 1:3, "]")]
  share_x <- eta[, 1] / rowSums(eta)
  share <- hypothesis_BF_model(fit, "(mu) R2[share: x] > 0.5")
  expected <- odds(mean(share_x > 0.5)) / odds(stats::pbeta(0.5, 1, 2, lower.tail = FALSE))
  expect_equal(as.numeric(share$BF), expected, tolerance = 1e-8)

  quantities <- hypothesis_quantities(fit)
  expect_true(all(c("(mu) R2{total}", "(mu) R2[share: x]", "(mu) R2[share: f]") %in% quantities$alias))
})

test_that("model-averaged budget totals keep the label of their quantity", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  fit_x <- .R2_catalog_fit(~ x, sample = 200)
  fit_xf <- .R2_catalog_fit(~ x + f, sample = 200)
  models <- models_inference(list(
    list(fit = fit_x, marglik = bridgesampling_object(-40), prior_weights = 1),
    list(fit = fit_xf, marglik = bridgesampling_object(-39), prior_weights = 1)))
  mixed <- mix_posteriors(models, parameters = "mu_R2", is_null_list = list(mu_R2 = c(FALSE, FALSE)), seed = 1)
  expect_identical(rownames(ensemble_estimates_table(mixed, parameters = "mu_R2")), "(mu) R2{total}")
})

test_that("share densities and budget names are safe at the edges", {

  set.seed(41)
  data <- data.frame(x = stats::rnorm(20), z = stats::rnorm(20))
  budget <- .R2_prior(terms = c("x", "z"), total = prior("beta", list(1, 1)),
    allocation = prior("dirichlet", list(alpha = c(1e16, 1))), law = "exact",
    convention = "common", concentration = "inherited", sigma = "sigma")
  info <- JAGS_formula(~ x + z, "mu", data, list(intercept = prior("normal", list(0, 1)), R2 = budget))
  object <- structure(list(), prior_list = info$prior_list)
  key <- list(type = "R2_share", dependencies = paste0("prior_par_eta_mu_R2_allocation[", 1:2, "]"),
    total = "mu_R2", term = "mu_x", units = 1L)
  # the complementary concentration is summed directly: Beta(1e16, 1), not Beta(1e16, 0)
  # (qbeta warns about its accuracy at such concentrations)
  density <- suppressWarnings(
    .bt_parameter_prior_density_R2_share(object, list(quantity_id = "q"), key, n_grid = 256, tail_prob = 1e-6))
  expect_false(is.null(density))

  # budget names cannot end in the suffixes of generated nodes
  for(name in c("f_ordered", "x_indicator", "w_inclusion", "w_variable", "w_allocation")){
    priors <- list(intercept = prior("normal", list(0, 1)))
    priors[[name]] <- budget
    expect_error(JAGS_formula(~ x + z, "mu", data, priors), "cannot end in", fixed = TRUE,
      class = "BayesTools_R2_invalid")
  }
})
