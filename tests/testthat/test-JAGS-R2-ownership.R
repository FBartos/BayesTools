skip_if_not_test_profile("unit")

test_that("R2 ownership leaves the coefficient selector and factor cells intact", {

  data <- data.frame(x = seq_len(12), f = factor(rep(letters[1:3], 4)))
  half <- prior("spike", list(0.5))
  budget <- prior_R2(
    total = prior_spike_and_slab(prior("beta", list(1, 1)), prior_inclusion = half),
    inclusion = prior("beta", list(2, 2)),
    model = R2_model(sigma = 1, contrast = "treatment")
  )
  priors <- list(intercept = prior("normal", list(0, 1)), x_signal = budget)
  formula <- JAGS_formula(~ x + f, "mu", data, priors)
  name_map <- formula$formula_design$name_map
  total_node <- name_map[name_map$jags_name == "mu_x_signal", ]
  expect_identical(total_node$kind, "generated")
  expect_identical(total_node$term, "x_signal")
  expect_identical(total_node$role, "r2_total")

  columns <- c("mu_intercept", "mu_x", "mu_f[1]", "mu_f[2]", "mu_x_signal",
               "mu_x_signal_rate", "mu_x_signal_rate_u", "mu_x_signal_indicator", "mu_x_signal_open[1]",
               "mu_x_signal_open[2]", "mu_x_signal_gate[1]", "mu_x_signal_gate[2]",
               paste0("prior_par_eta_mu_x_signal_allocation[", 1:2, "]"))
  values <- matrix(0.5, 3, length(columns), dimnames = list(NULL, columns))
  values[, "mu_x_signal_indicator"] <- c(0, 1, 1)
  values[, "mu_x_signal_gate[1]"] <- c(1, 1, 0)
  values[, "mu_x_signal_open[1]"] <- c(0, 1, 0)
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(values)), formula$prior_list,
                                     list(mu = formula$formula_design))
  catalog <- parameter_catalog(fit)
  row <- function(label) parameter_catalog_resolve(catalog, label, namespace = "mu")$quantities
  expect_identical(row("x")$canonical_name, "mu_x")
  expect_identical(row("x")$role, "fixed_coefficient")
  expect_identical(row("x")$quantity_id, .bt_parameter_catalog_quantity_id("mu_x", "mu", "fixed_coefficient"))
  expect_identical(row("x_signal{total}")$quantity, "total")
  expect_identical(row("x_signal{total}")$role, "parameter")
  expect_identical(row("x_signal{total} (inclusion)")$quantity, "budget_inclusion")
  expect_identical(row("x_signal{inclusion rate}")$quantity, "inclusion_rate")
  expect_identical(row("x_signal[share: x]")$arguments[[1L]], "x")
  expect_identical(row("x_signal[share: x]")$definedness, "term_included")
  expect_identical(row("x_signal[share: x]")$quantity, "share")
  for(level in levels(data$f)){
    cell <- row(paste0("f[", level, "]"))
    expect_identical(cell$role, "fixed_coefficient")
    expect_identical(cell$quantity, "coefficient")
    expect_identical(cell$term, "f")
  }
  inclusion <- parameter_catalog_resolve(catalog, "x (inclusion)", namespace = "mu")
  expect_identical(inclusion$quantities$canonical_name, "mu_x_signal_open[x]")
  expect_identical(inclusion$quantities$quantity, "term_inclusion")
  expect_identical(inclusion$quantities$source_type, "identity")
  expect_error(parameter_catalog_resolve(catalog, "mu_x_signal_rate_u"), class = "BayesTools_parameter_not_found")
  expect_identical(row("mu_x_signal_gate[x]")$quantity_id, inclusion$quantity_id)
  expect_equal(as.numeric(as.matrix(parameter_draws(fit, inclusion))), c(0, 1, 0))
  expect_error(parameter_draws(fit, inclusion, model_samples = values[, "mu_x", drop = FALSE]),
               class = "BayesTools_R2_coordinates_unavailable")

  owned <- catalog$quantities[catalog$quantities$owner_type == "r2_budget", ]
  expect_true(all(owned$owner_name == "x_signal" & owned$formula_parameter == "mu"))
  expect_true(all(owned$parent_quantity_id[owned$quantity != "total"] == row("x_signal{total}")$quantity_id))
  expect_identical(row("intercept")$owner_type, "")
  expect_false(any(grepl("_gate\\[[12]\\]|prior_par_", owned$canonical_name)))
})

test_that("a compulsory share is defined and missing R2 sources retain their identity", {

  data <- data.frame(x = seq_len(10), z = rep(c(0, 1), 5))
  formula <- JAGS_formula(~ x + z, "mu", data, list(
    intercept = prior("normal", list(0, 1)),
    group = prior_R2(inclusion = list(z = prior("spike", list(0.3))), model = R2_model(sigma = 1))
  ))
  fit <- .prior_monitor_test_fit(formula, c("mu_intercept", "mu_x", "mu_z", "mu_group"))
  catalog <- parameter_catalog(fit)
  compulsory <- parameter_catalog_resolve(catalog, "group[share: x]")
  gated <- parameter_catalog_resolve(catalog, "group[share: z]")
  event <- parameter_catalog_resolve(catalog, "z (inclusion)")
  expect_identical(compulsory$quantities$definedness, "always")
  expect_identical(gated$quantities$definedness, "term_included")
  expect_identical(event$quantities$status, "unavailable")
  expect_identical(gated$quantities$status, "unavailable")
  coordinates <- parameter_coordinates(fit)
  sources <- event$quantities$extraction_key[[1L]]$dependencies
  expect_identical(coordinates$monitor_status[match(sources, coordinates$coordinate_name)], "unavailable")
  expect_error(parameter_draws(fit, gated), class = "BayesTools_refit_monitoring")
  expect_error(parameter_draws(fit, event), class = "BayesTools_refit_monitoring")

  other <- JAGS_formula(~ x + z, "eta", data, list(
    intercept = prior("normal", list(0, 1)), group = prior_R2(model = R2_model(sigma = 1))
  ))
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(cbind(mu_x = 1:3, eta_x = 3:1))),
    c(formula$prior_list, other$prior_list), list(mu = formula$formula_design, eta = other$formula_design))
  catalog <- parameter_catalog(fit)
  expect_error(parameter_catalog_resolve(catalog, "group{total}"), class = "BayesTools_parameter_ambiguous")
  for(parameter in c("mu", "eta")){
    total <- parameter_catalog_resolve(catalog, "group{total}", namespace = parameter)
    share <- parameter_catalog_resolve(catalog, "group[share: x]", namespace = parameter)
    expect_identical(total$quantities$status, "unavailable")
    expect_identical(share$quantities$parent_quantity_id, total$quantity_id)
  }
})

test_that("a compiled R2 term cannot enter a new formula by itself or in a mixture", {

  data <- data.frame(x = seq_len(5))
  intercept <- prior("normal", list(0, 1))
  info <- JAGS_formula(~ x, "mu", data, list(intercept = intercept,
    budget = prior_R2(model = R2_model(sigma = 1))))
  bound <- info$prior_list$mu_x
  candidates <- list(bound, prior_mixture(list(prior("spike", list(0)), bound)),
                     prior_spike_and_slab(bound))
  for(candidate in candidates){
    expect_error(JAGS_formula(~ x, "other", data, list(intercept = intercept, x = candidate)),
      "original named 'prior_R2()' budget", fixed = TRUE, class = "BayesTools_R2_invalid")
  }
})

test_that("R2 relationship edges retain external scale ownership, member units and realized events", {

  data <- data.frame(x = seq_len(12), f = factor(rep(letters[1:3], 4)))
  budget <- prior_R2(
    inclusion = prior("beta", list(2, 3)),
    total = prior_spike_and_slab(prior("beta", list(1, 1))),
    model = R2_model(law = "normal", sigma = "tau_total", contrast = "treatment")
  )
  info <- JAGS_formula(~ x + f, "mu", data, list(intercept = prior("normal", list(0, 1)), signal = budget))
  chains <- coda::mcmc.list(coda::mcmc(cbind(tau_total = c(1, 2), mu_x = c(-0.1, 0.4))))
  fit <- .parameter_catalog_test_fit(chains, c(info$prior_list, list(tau_total = prior("gamma", list(2, 2)))),
                                     list(mu = info$formula_design))
  edges <- parameter_relationships(fit)
  edges <- edges[edges$allocation == "signal", , drop = FALSE]
  catalog <- parameter_catalog(fit)
  qid <- function(alias) parameter_catalog_resolve(catalog, alias)$quantity_id
  vertex <- edges$consumer_id[edges$relation == "budget_total" & edges$source_id == qid("signal{total}")]
  expect_length(vertex, 1L)
  expect_identical(edges$source_type[edges$relation == "residual_sd"], "quantity")
  expect_identical(edges$source_id[edges$relation == "residual_sd"], qid("tau_total"))
  expect_identical(edges$consumer_id[edges$relation == "residual_sd"], vertex)
  expect_true(all(edges$normalization == "normal" & edges$formula_parameter == "mu"))
  branches <- edges[edges$relation == "allocates" & edges$source_id == vertex, ]
  expect_setequal(branches$component, c("x", "f"))
  factor_id <- branches$consumer_id[branches$component == "f"]
  factor <- edges[edges$consumer_id == factor_id, ]
  expect_setequal(factor$source_id[factor$relation == "allocation_weight"], paste0("mu_signal_allocation[", 2:3, "]"))
  expect_setequal(factor$source_id[factor$relation == "coefficient_direction"], paste0("prior_par_z_mu_f[", 1:2, "]"))
  expect_identical(factor$source_id[factor$relation == "inclusion_gate"], "mu_signal_open[2]")
  reports <- edges[edges$source_id == factor_id, ]
  expect_identical(reports$consumer_id[reports$relation == "share"], qid("signal[share: f]"))
  expect_identical(reports$consumer_id[reports$relation == "term_inclusion"], qid("f (inclusion)"))
  expect_setequal(reports$consumer_id[reports$relation == "budget_coefficient"],
                  vapply(paste0("f[", levels(data$f), "]"), qid, character(1), USE.NAMES = FALSE))
  event <- edges[edges$consumer_id == qid("f (inclusion)"), ]
  expect_identical(event$source_id[event$relation == "raw_gate"], "mu_signal_gate[2]")
  expect_identical(event$source_id[event$relation == "budget_gate"], "mu_signal_indicator")
  rate <- edges[edges$consumer_id == qid("signal{inclusion rate}"), ]
  expect_identical(rate$source_id[rate$relation == "rate_draw"], "mu_signal_rate_u")
  expect_setequal(rate$source_id[rate$relation == "conditions_rate"], paste0("mu_signal_gate[", 1:2, "]"))
  expect_identical(parameter_relationships(unserialize(serialize(fit, NULL))), parameter_relationships(fit))

  one <- JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)),
    signal = prior_R2(model = R2_model(sigma = 1.25))))
  single <- parameter_relationships(.prior_monitor_test_fit(one, c("mu_intercept", "mu_x", "mu_signal")))
  expect_identical(single$source_type[single$relation == "residual_sd"], "constant")
  expect_identical(single$source_id[single$relation == "residual_sd"], "1.25")
  expect_equal(sum(single$relation == "allocates"), 1)
  expect_false("share" %in% single$relation)

  # A residual scale may itself be a realized model quantity. The shared rate
  # is drawn conditional on the raw gates, independently of this scale use.
  own_rate <- JAGS_formula(~ x + f, "mu", data, list(intercept = prior("normal", list(0, 1)),
    signal = prior_R2(inclusion = prior("beta", list(2, 2)), model = R2_model(sigma = "mu_signal_rate"))))
  rate_fit <- .prior_monitor_test_fit(own_rate, c("mu_x", "mu_signal", "mu_signal_rate"))
  rate_id <- parameter_catalog_resolve(parameter_catalog(rate_fit), "signal{inclusion rate}")$quantity_id
  reused <- parameter_relationships(rate_fit)
  expect_identical(reused$source_id[reused$relation == "residual_sd"], rate_id)
})
