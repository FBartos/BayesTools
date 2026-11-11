skip_if_not_test_profile("unit")

.R2_bridge_priors <- function(law){
  set.seed(5)
  data <- data.frame(f = factor(rep(c("a", "b", "c"), each = 6)), x = stats::rnorm(18))
  budget <- prior_R2(alpha = c(x = 1, f = 0.5), a = 1.5, b = 2, law = law, convention = "common",
    concentration = "inherited", sigma = "sigma", contrast = "meandif")
  info <- JAGS_formula(~ x + f, "mu", data, list(intercept = prior("normal", list(0, 1)), "__R2" = budget))
  info$prior_list
}

# One bridge draw: the budget total, the gamma coordinates of the allocation
# (not the normalized weights), and the direction coordinates.
.R2_bridge_row <- function(priors){
  set.seed(9)
  allocation <- priors$mu_R2_allocation
  eta_names <- paste0(.JAGS_prior_dirichlet_eta_name("mu_R2_allocation"), "[", seq_len(allocation$parameters$K), "]")
  directions <- c(.bt_R2_direction_coordinates(priors$mu_x, "mu_x"), .bt_R2_direction_coordinates(priors$mu_f, "mu_f"))
  c(stats::setNames(0.3, "mu_intercept"), stats::setNames(0.2, "mu_R2"),
    stats::setNames(stats::rgamma(length(eta_names), 2), eta_names),
    stats::setNames(stats::rnorm(length(directions)), directions))
}

test_that("bridge coordinates of R2 terms are their direction coordinates", {

  for(law in c("exact", "normal", "laplace")){
    priors <- .R2_bridge_priors(law)
    info <- .JAGS_bridgesampling_posterior_info(priors)
    directions <- c(.bt_R2_direction_coordinates(priors$mu_x, "mu_x"), .bt_R2_direction_coordinates(priors$mu_f, "mu_f"))
    expect_true(all(directions %in% info))
    expect_false(any(c("mu_x", "mu_f[1]", "mu_f[2]") %in% info))
    expect_true(all(is.infinite(attr(info, "lb")[directions]) & attr(info, "lb")[directions] < 0))
    expect_equal(unname(attr(info, "lb")["mu_R2"]), 0)
    expect_equal(unname(attr(info, "ub")["mu_R2"]), 1)
    owned <- .bt_JAGS_bridge_owned_parameter_names(priors)
    expect_true(all(c("mu_x", "mu_f[1]", directions) %in% owned))
  }
})

test_that("the bridge log prior and coefficients of R2 terms follow the declared hierarchy", {

  for(law in c("exact", "normal", "laplace")){
    priors <- .R2_bridge_priors(law)
    row <- .R2_bridge_row(priors)
    evaluator <- .bt_JAGS_bridge_compile_prior_list_evaluator(priors)

    eta <- row[grepl("prior_par_eta", names(row))]
    z_x <- row[.bt_R2_direction_coordinates(priors$mu_x, "mu_x")]
    z_f <- row[.bt_R2_direction_coordinates(priors$mu_f, "mu_f")]
    log_z <- function(z) if(law == "laplace") sum(log(sqrt(2) / 2) - sqrt(2) * abs(z)) else sum(stats::dnorm(z, log = TRUE))
    expected_log_prior <- stats::dnorm(0.3, 0, 1, log = TRUE) + stats::dbeta(0.2, 1.5, 2, log = TRUE) +
      sum(stats::dgamma(eta, priors$mu_R2_allocation$parameters$alpha, 1, log = TRUE)) + log_z(z_x) + log_z(z_f)
    expect_equal(evaluator$log_prior(row), expected_log_prior, tolerance = 1e-12)
    expect_equal(.bt_JAGS_marglik_priors(as.list(row), priors), expected_log_prior, tolerance = 1e-12)
    rows <- rbind(row, row)
    expect_equal(JAGS_marglik_priors_rows(rows, priors), rep(expected_log_prior, 2), tolerance = 1e-12)

    W <- 0.2 / 0.8
    phi <- eta / sum(eta)
    parameters <- evaluator$parameters(row)
    if(law == "exact"){
      expected_x <- sqrt(W * phi[1]) * sign(z_x) * priors$mu_x$map[1, 1]
      expected_f <- sqrt(W * phi[2]) * as.vector(priors$mu_f$map %*% (z_f / sqrt(sum(z_f^2))))
    }else{
      expected_x <- sqrt(W * phi[1]) * z_x * priors$mu_x$map[1, 1]
      expected_f <- sqrt(W * phi[2:3]) * z_f * diag(priors$mu_f$map)
    }
    expect_equal(unname(parameters$mu_x), unname(expected_x), tolerance = 1e-13)
    expect_equal(unname(parameters$mu_f), unname(expected_f), tolerance = 1e-13)
    expect_equal(JAGS_marglik_parameters(row, priors)$mu_f, parameters$mu_f)
  }
})

test_that("the fixed linear predictor of a bridge draw ignores the budget priors and applies sigma", {

  set.seed(5)
  data <- data.frame(f = factor(rep(c("a", "b", "c"), each = 6)), x = stats::rnorm(18))
  budget <- prior_R2(alpha = c(x = 1, f = 0.5), a = 1.5, b = 2, law = "exact", convention = "common",
    concentration = "inherited", sigma = "sigma", contrast = "meandif")
  info <- JAGS_formula(~ x + f, "mu", data, list(intercept = prior("normal", list(0, 1)), "__R2" = budget))
  plan <- .bt_dnode_linear_predictor_fixed_plan(info$formula_design, info$prior_list)
  row <- .R2_bridge_row(info$prior_list)
  coefficients <- .bt_JAGS_bridge_compile_prior_list_evaluator(info$prior_list)$parameters(row)
  sigma <- 1.7
  value <- plan$value(row, c(list(sigma = sigma), coefficients))
  expected <- 0.3 + sigma * (info$data$mu_data_x * coefficients$mu_x + info$data$mu_data_f %*% coefficients$mu_f)
  expect_equal(as.vector(value), as.vector(expected), tolerance = 1e-12)
})

test_that("coefficient-only draws and out-of-support budget states are handled", {

  priors <- .R2_bridge_priors("exact")
  values <- .bt_JAGS_marglik_compile_R2_parameter_values(priors$mu_x, "mu_x")

  # a monitored coefficient with an extra direction coordinate but no budget
  for(row in list(c(mu_x = 0.4, prior_par_z_mu_x = -1), list(mu_x = 0.4, prior_par_z_mu_x = -1))){
    expect_equal(values(row), 0.4)
    expect_true(.bt_R2_values_available(priors$mu_x, "mu_x", names(row)))
  }
  expect_false(.bt_R2_values_available(priors$mu_x, "mu_x", "prior_par_z_mu_x"))
  expect_error(values(c(prior_par_z_mu_x = -1)), class = "BayesTools_missing_monitored_columns")

  row <- .R2_bridge_row(priors)
  row[["mu_R2"]] <- 1
  expect_error(values(row), class = "BayesTools_marglik_out_of_support")
  row <- .R2_bridge_row(priors)
  row[grepl("prior_par_eta", names(row))][1] <- -1
  expect_error(values(row), class = "BayesTools_marglik_out_of_support")
})
