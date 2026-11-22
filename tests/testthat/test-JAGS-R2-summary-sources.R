skip_if_not_test_profile("unit")

# One draw of each joint gate state, including a raw gate which is on
# under the spike. The semantic table needs only the declared open flags.
.R2_summary_source_fit <- function(parameters = "mu"){

  data <- data.frame(x = c(-1, -1, 1, 1), z = c(-1, 1, -1, 1))
  inclusion <- prior("spike", list(0.5))
  priors <- list(intercept = prior("normal", list(0, 1)),
                 signal = prior_R2(
                   total = prior_spike_and_slab(prior("beta", list(2, 2)), prior_inclusion = inclusion),
                   inclusion = inclusion, model = R2_model(sigma = 1)
                 ))
  compiled <- lapply(parameters, function(parameter) JAGS_formula(~ x + z, parameter, data, priors))
  rows <- cbind(
    mu_intercept = rep(0.1, 4), mu_signal = c(0, 0.2, 0.4, 0.6),
    mu_x = c(0, 2, 0, 4), mu_z = c(0, 0, 1.5, 2.5),
    mu_signal_indicator = c(0, 1, 1, 1),
    "mu_signal_gate[1]" = c(1, 1, 0, 1), "mu_signal_gate[2]" = c(0, 0, 1, 1),
    "mu_signal_open[1]" = c(0, 1, 0, 1), "mu_signal_open[2]" = c(0, 0, 1, 1),
    "mu_signal_allocation[1]" = c(0, 1, 0, 0.25), "mu_signal_allocation[2]" = c(0, 0, 1, 0.75),
    "prior_par_eta_mu_signal_allocation[1]" = rep(1, 4),
    "prior_par_eta_mu_signal_allocation[2]" = rep(3, 4)
  )
  sources <- lapply(parameters, function(parameter){
    out <- rows
    colnames(out) <- sub("mu_", paste0(parameter, "_"), colnames(out), fixed = TRUE)
    out
  })
  chains <- coda::mcmc.list(coda::mcmc(do.call(cbind, sources), start = 2, thin = 2))
  .parameter_catalog_test_fit(
    chains, do.call(c, lapply(compiled, `[[`, "prior_list")),
    formula_design = stats::setNames(lapply(compiled, `[[`, "formula_design"), parameters)
  )
}

test_that("R2 values, inclusion events and tables use the declared open states in the supplied view", {

  fit <- .R2_summary_source_fit()
  chains <- JAGS_materialize_draws(fit, include_internal = TRUE)
  cols <- colnames(chains[[1]])
  # Neither raw gates nor monitored normalized weights are public value inputs.
  keep <- !startsWith(cols, "mu_signal_gate") & !startsWith(cols, "mu_signal_allocation")
  view <- JAGS_with_draws(fit, coda::mcmc.list(lapply(chains, function(ch) ch[, keep, drop = FALSE])))
  share <- parameter_catalog_resolve(parameter_catalog(view), "signal[share: x]", namespace = "mu")
  source_draws <- parameter_draws(view, share)
  expect_identical(as.numeric(as.matrix(source_draws)), c(NA_real_, 1, NA_real_, 0.25))
  expect_identical(attr(source_draws[[1]], "mcpar"), attr(chains[[1]], "mcpar"))
  states <- parameter_gate_states(view, share)
  expect_identical(states$defined, c(FALSE, TRUE, FALSE, TRUE))
  expect_equal(states$atom, c(NA_real_, 1, NA_real_, NA_real_))

  estimates <- JAGS_estimates_table(view, conditional = TRUE, remove_diagnostics = TRUE)
  expect_equal(estimates["(mu) x", "Mean"], 3)
  expect_equal(estimates["(mu) signal{total}", "Mean"], 0.4)
  expect_equal(estimates["(mu) signal[share: x]", "Mean"], 0.625)
  raw_estimates <- JAGS_estimates_table(view, conditional = TRUE, random_effects_summary = "raw",
                                        remove_diagnostics = TRUE)
  expect_equal(raw_estimates["(mu) x", "Mean"], 3)
  expect_equal(raw_estimates["(mu) signal{total}", "Mean"], 0.4)
  note <- "(mu) signal[share: x]: summarized over 2 of 4 draws where the term is included in its R-squared budget."
  expect_identical(unname(attr(estimates, "footnotes")[["(mu) signal[share: x]"]]), note)
  samples <- JAGS_estimates_table(view, return_samples = TRUE)
  expect_identical(posterior_metadata(samples, "undefined_draws")[["(mu) signal[share: x]"]], "term_included")
  expect_identical(unname(samples[, "(mu) signal[share: x]"]), c(NA_real_, 1, NA_real_, 0.25))
  inclusion <- JAGS_inference_table(view)
  expect_equal(inclusion["(mu) x", "post_prob"], 0.5)
  expect_equal(inclusion["(mu) signal{total}", "post_prob"], 0.75)
  expect_identical(attr(inclusion, "parameter_roles"), c("r2_budget_inclusion", "r2_term_inclusion", "r2_term_inclusion"))
  events <- c("signal{total} (inclusion)", "x (inclusion)", "z (inclusion)")
  ids <- vapply(events, function(alias){
    parameter_catalog_resolve(parameter_catalog(view), alias, namespace = "mu")$quantity_id
  }, character(1), USE.NAMES = FALSE)
  expect_identical(attr(inclusion, "event_quantity_ids"), ids)
  expect_identical(attr(inclusion[3:2, , drop = FALSE], "event_quantity_ids"), ids[3:2])
  expect_identical(attr(inclusion[c(3,1), , drop = FALSE], "parameter_roles"),
                   c("r2_term_inclusion", "r2_budget_inclusion"))
  estimate_metadata <- attr(estimates, "quantities")
  expect_identical(estimate_metadata$quantity_id[match(paste("(mu)", events), estimate_metadata$row)], ids)

  mixed <- parameter_mixed_posterior(view, share)
  expect_equal(as.numeric(mixed), c(1, 0.25))
  name <- share$quantities$canonical_name
  expect_equal(posterior_metadata(mixed, "omitted_undefined_draws"), 2)
  expect_null(posterior_metadata(.bt_draws_subset_rows(mixed, 1L), "omitted_undefined_draws"))
  summary <- ensemble_estimates_table(stats::setNames(list(mixed), name), name)
  expect_identical(unname(attr(summary, "footnotes")[["(mu) signal[share: x]"]]), note)
})

test_that("a missing R2 event refuses dependent operations and permits unrelated terms", {

  fit <- .R2_summary_source_fit()
  chains <- JAGS_materialize_draws(fit, include_internal = TRUE)
  keep <- colnames(chains[[1]]) != "mu_signal_open[2]"
  view <- JAGS_with_draws(fit, coda::mcmc.list(lapply(chains, function(ch) ch[, keep, drop = FALSE])))
  share <- parameter_catalog_resolve(parameter_catalog(view), "signal[share: x]", namespace = "mu")
  expect_error(parameter_draws(view, share), class = "BayesTools_draws_view_coordinate_unavailable")
  expect_error(parameter_gate_states(view, share), class = "BayesTools_R2_coordinates_unavailable")
  expect_error(JAGS_estimates_table(view, conditional = TRUE, remove_diagnostics = TRUE),
               class = "BayesTools_draws_view_coordinate_unavailable")
  expect_error(JAGS_inference_table(view), class = "BayesTools_R2_coordinates_unavailable")
  selected <- JAGS_estimates_table(view, conditional = TRUE, keep_parameters = "mu_x", remove_diagnostics = TRUE)
  expect_equal(selected["(mu) x", "Mean"], 3)
  expect_false(any(grepl("share", rownames(selected))))
})

test_that("a wholly undefined R2 share retains its denominator and a local refusal", {

  fit <- .R2_summary_source_fit()
  raw <- as.matrix(JAGS_materialize_draws(fit, include_internal = TRUE)[[1]])
  closed <- JAGS_with_draws(fit, coda::mcmc.list(coda::mcmc(raw[c(1, 1), , drop = FALSE])))
  share <- parameter_catalog_resolve(parameter_catalog(closed), "signal[share: x]", namespace = "mu")
  table <- JAGS_estimates_table(closed, remove_diagnostics = TRUE)
  label <- "(mu) signal[share: x]"
  expect_identical(table[label, "Mean"], NA_real_)
  expect_match(attr(table, "footnotes")[[label]], "0 of 2 draws", fixed = TRUE)
  expect_error(parameter_mixed_posterior(closed, share), class = "BayesTools_no_defined_posterior_draws")
  bad <- raw
  bad[1, "mu_signal_open[1]"] <- 0.5
  expect_error(parameter_gate_states(fit, share, draws = bad), class = "BayesTools_R2_unavailable")
})

test_that("inference event IDs refer to Bernoulli quantities", {

  priors <- list(
    a = prior_spike_and_slab(prior("normal", list(0, 1))),
    b = prior_mixture(list(prior("normal", list(0, 1)), prior("normal", list(1, 1))),
                      components = c("one", "two"))
  )
  rows <- cbind(a = c(0, 1, 0, 2), a_indicator = c(0, 1, 0, 1),
                b = c(0.4, 1.3, 0.7, 0.9), b_indicator = c(1, 2, 1, 2))
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(rows)), priors)
  table <- JAGS_inference_table(fit)
  expected <- parameter_catalog_resolve(parameter_catalog(fit), "a_indicator")$quantity_id
  expect_identical(attr(table, "event_quantity_ids"), c(expected, "", ""))
})

test_that("budgets with the same name in two formulas keep distinct table event IDs", {

  fit <- .R2_summary_source_fit(c("mu", "eta"))
  table <- JAGS_inference_table(fit, formula_prefix = FALSE)
  expect_identical(rownames(table), c("signal{total}", "x", "z", "signal{total}.1", "x.1", "z.1"))
  expected <- vapply(c("mu", "eta"), function(parameter){
    parameter_catalog_resolve(parameter_catalog(fit), "x (inclusion)", namespace = parameter)$quantity_id
  }, character(1), USE.NAMES = FALSE)
  expect_identical(attr(table, "event_quantity_ids")[c(2, 5)], expected)
  expect_length(unique(expected), 2L)
  summary <- JAGS_estimates_table(fit, formula_prefix = FALSE, remove_diagnostics = TRUE)
  expect_true(all(c("signal[share: x]", "signal[share: x].1") %in% names(attr(summary, "footnotes"))))
})
