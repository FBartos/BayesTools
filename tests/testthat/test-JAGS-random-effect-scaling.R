skip_if_not_test_profile("unit")

make_random_scale_table_fit <- function(formula_result, posterior){

  fit <- structure(
    list(
      mcmc = coda::mcmc.list(coda::mcmc(posterior)),
      sample = nrow(posterior),
      summary.pars = list(mutate = NULL),
      monitor = colnames(posterior)
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  attach_test_parameter_map(fit)
}

test_that("JAGS_estimates_table suppresses fixed warnings for random-only scaled slopes", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  posterior <- matrix(
    rep(c(0, 2), 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", "mu__xREx__id_x"))
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  samples <- expect_silent(JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  ))

  expect_equal(
    unname(samples[, "(mu) sd(x)"]),
    rep(2 / formula_result$formula_scale$mu_x$sd, nrow(samples)),
    tolerance = 1e-12
  )
})

test_that("JAGS_estimates_table suppresses fixed warnings for homogeneous random-only slopes", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + id(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  posterior <- matrix(
    rep(c(0, 2), 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", random_term$sd_parameter_names))
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  samples <- expect_silent(JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  ))

  expect_equal(
    unname(samples[, "(mu) sd"]),
    rep(2 / formula_result$formula_scale$mu_x$sd, nrow(samples)),
    tolerance = 1e-12
  )
})

test_that("JAGS_estimates_table treats sd as a valid random-only predictor name", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    sd = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + sd | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(sd = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  posterior <- matrix(
    rep(c(0, 2), 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", random_term$sd_parameter_names))
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  samples <- expect_silent(JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  ))

  expect_equal(
    unname(samples[, "(mu) sd(sd)"]),
    rep(2 / formula_result$formula_scale$mu_sd$sd, nrow(samples)),
    tolerance = 1e-12
  )
})

test_that("JAGS_estimates_table still warns about genuinely unused scale entries", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  formula_result$formula_scale$mu_z <- list(mean = 0, sd = 1)
  random_term <- formula_result$formula_design$random_effects[[1]]
  posterior <- matrix(
    rep(c(0, 2), 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", random_term$sd_parameter_names))
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  warnings <- character()
  samples <- withCallingHandlers(
    JAGS_estimates_table(
      fit,
      transform_scaled = TRUE,
      random_effects_summary = "standard",
      remove_diagnostics = TRUE,
      return_samples = TRUE
    ),
    warning = function(w){
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_length(warnings, 1L)
  expect_match(warnings, "mu_z", fixed = TRUE)
  expect_false(grepl("mu_x", warnings, fixed = TRUE))
  expect_equal(
    unname(samples[, "(mu) sd(x)"]),
    rep(2 / formula_result$formula_scale$mu_x$sd, nrow(samples)),
    tolerance = 1e-12
  )
})

test_that("JAGS_estimates_table unscales diagonal random intercept-slope blocks without correlations", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  posterior <- matrix(
    rep(c(0, 1, 2), 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c("mu_intercept", random_term$sd_parameter_names)
    )
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  samples <- expect_silent(JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  ))

  scale_info <- formula_result$formula_scale$mu_x
  expected_intercept_sd <- sqrt(1^2 + (scale_info$mean / scale_info$sd)^2 * 2^2)
  expected_slope_sd <- 2 / scale_info$sd

  expect_equal(
    unname(samples[, "(mu) sd(intercept)"]),
    rep(expected_intercept_sd, nrow(samples)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(samples[, "(mu) sd(x)"]),
    rep(expected_slope_sd, nrow(samples)),
    tolerance = 1e-12
  )
})

test_that("JAGS_estimates_table keeps fixed and random scaled slope transforms together", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + diag(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  posterior <- matrix(
    rep(c(10, 4, 2), 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c("mu_intercept", "mu_x", random_term$sd_parameter_names)
    )
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  samples <- expect_silent(JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  ))

  scale_info <- formula_result$formula_scale$mu_x
  expect_equal(
    unname(samples[, "(mu) intercept"]),
    rep(10 - scale_info$mean / scale_info$sd * 4, nrow(samples)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(samples[, "(mu) x"]),
    rep(4 / scale_info$sd, nrow(samples)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(samples[, "(mu) sd(x)"]),
    rep(2 / scale_info$sd, nrow(samples)),
    tolerance = 1e-12
  )
})

test_that("original-scale accessors exclude internal random coordinates", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(
          latent = TRUE,
          coefficients = TRUE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  internal_names <- c(
    "mu__xREx__id_xRE_Zx[1,1]",
    "mu__xREx__id_xRE_COEFx[1,1]"
  )
  posterior_names <- c(
    "mu_intercept",
    "mu_x",
    random_term$sd_parameter_names,
    internal_names
  )
  posterior <- matrix(
    rep(seq_along(posterior_names), each = 3L),
    nrow = 3L,
    dimnames = list(NULL, posterior_names)
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  transformed_fit <- transform_scale_samples(fit)
  expect_false(any(internal_names %in% colnames(transformed_fit)))
  expect_true(all(c("mu_intercept", "mu_x") %in% colnames(transformed_fit)))

  transformed_matrix <- transform_scale_samples(
    posterior,
    list(mu = formula_result$formula_scale)
  )
  expect_false(any(internal_names %in% colnames(transformed_matrix)))

  raw_fitted <- JAGS_estimates_table(
    fit,
    transform_scaled = FALSE,
    random_effects_summary = "raw",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  )
  expect_true(any(grepl("z\\(", colnames(raw_fitted))))
  expect_true(any(grepl("coef\\(", colnames(raw_fitted))))

  raw_original <- JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "raw",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  )
  expect_false(any(grepl("z\\(", colnames(raw_original))))
  expect_false(any(grepl("coef\\(", colnames(raw_original))))
})

.undefined_correlation_table_fit <- function(source_sd, source_rho){

  df <- data.frame(
    x = c(4, 5, 6, 8),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "id", covariance = "us"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("normal", list(0, 1), truncation = list(0, Inf)),
        cor = prior_lkj(eta = 1, include_primitives = TRUE),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = TRUE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  R_names <- as.vector(outer(
    seq_len(2),
    seq_len(2),
    Vectorize(function(row, column){
      paste0(random_term$parameter_stem, "_xRE_CORx_R[", row, ",", column, "]")
    })
  ))
  L_names <- as.vector(
    BayesTools:::.bt_random_effect_cholesky_names(random_term, 2L)
  )
  rows <- lapply(seq_along(source_rho), function(draw){
    source_cor <- matrix(c(1, source_rho[draw], source_rho[draw], 1), 2, 2)
    c(0, 0, source_sd[draw, ], as.vector(source_cor),
      as.vector(t(chol(source_cor))), (source_rho[draw] + 1) / 2,
      source_rho[draw])
  })
  posterior <- do.call(rbind, rows)
  colnames(posterior) <- c(
    "mu_intercept", "mu_x", random_term$sd_parameter_names, R_names, L_names,
    random_term$correlation$primitive_names, random_term$correlation$cpc_names
  )

  make_random_scale_table_fit(formula_result, posterior)
}

test_that("estimates tables report the share of draws with a defined original-scale correlation", {

  skip_if_not_installed("runjags")

  # A zero slope SD (e.g., an excluded spike-and-slab SD) leaves the
  # original-scale intercept-slope correlation undefined (0/0): draws 3 and 4
  # are missing and the row summarises the 3 defined draws of 5.
  fit <- .undefined_correlation_table_fit(
    source_sd = rbind(c(1, 2), c(0.5, 1.5), c(0.7, 0), c(0.3, 0), c(1, 1)),
    source_rho = c(0.8, -0.3, 0.4, 0.1, 0.2)
  )
  expected_footnote <- c(
    "(mu) cor(intercept,x)" = paste0(
      "(mu) cor(intercept,x): summarised over 3 of 5 draws where the ",
      "correlation is defined."
    )
  )
  for(transform_scaled in c(FALSE, TRUE)){
    samples <- JAGS_estimates_table(
      fit,
      transform_scaled = transform_scaled,
      return_samples = TRUE
    )
    expect_equal(
      which(is.na(samples[, "(mu) cor(intercept,x)"])),
      c(3L, 4L)
    )
    for(mode in c("standard", "full")){
      table <- JAGS_estimates_table(
        fit,
        transform_scaled = transform_scaled,
        random_effects_summary = mode
      )
      expect_identical(
        attr(table, "footnotes"),
        expected_footnote,
        info = paste(mode, transform_scaled)
      )
    }
  }

  table <- JAGS_estimates_table(fit, footnotes = "User note.")
  expect_identical(
    attr(table, "footnotes"),
    c("User note.", expected_footnote)
  )
  expect_output(
    print(table),
    "summarised over 3 of 5 draws where the correlation is defined.",
    fixed = TRUE
  )
  # The row footnote follows its row when the table is subset.
  expect_identical(
    attr(
      update(table, remove_parameters = "(mu) cor(intercept,x)"),
      "footnotes"
    ),
    "User note."
  )
  expect_identical(
    attr(table[c("(mu) sd(x)", "(mu) cor(intercept,x)"), ], "footnotes"),
    c("User note.", expected_footnote)
  )
  expect_null(attr(
    JAGS_estimates_table(fit)[c("(mu) sd(intercept)", "(mu) sd(x)"), ],
    "footnotes"
  ))
  simplified <- JAGS_estimates_table(
    fit,
    formula_prefix = FALSE,
    simplify_names = TRUE
  )
  expect_identical(
    unname(attr(simplified, "footnotes")),
    paste0(
      "cor(intercept,x): summarised over 3 of 5 draws where the correlation ",
      "is defined."
    )
  )

  # Raw original-scale correlation coordinates are reported per row as well;
  # on the fitted scale every draw is defined.
  raw <- JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "raw"
  )
  expect_true(
    "(mu) cor(intercept,x | id)" %in% names(attr(raw, "footnotes"))
  )
  expect_true(all(grepl(
    ": summarised over 3 of 5 draws where the correlation is defined.$",
    attr(raw, "footnotes")
  )))
  expect_null(attr(
    JAGS_estimates_table(fit, random_effects_summary = "raw"),
    "footnotes"
  ))

  # Every draw defined: no footnote.
  defined <- .undefined_correlation_table_fit(
    source_sd = rbind(c(1, 2), c(0.5, 1.5), c(0.7, 0.2)),
    source_rho = c(0.8, -0.3, 0.4)
  )
  for(transform_scaled in c(FALSE, TRUE)){
    expect_null(attr(
      JAGS_estimates_table(defined, transform_scaled = transform_scaled),
      "footnotes"
    ))
  }
})

test_that("correlation footnotes count only draws retained by conditioning", {

  correlation_prior <- prior_none()
  attr(correlation_prior, "random_summary") <- "cor"
  model_samples <- cbind(
    "(mu) cor(a,b)" = c(0.1, NA, NA, 0.3, NA),
    "(mu) sd(a)"    = c(NA, 1, 1, 1, 1)
  )
  coordinates <- data.frame(
    coordinate_name = character(),
    role = character(),
    stringsAsFactors = FALSE
  )
  footnotes <- BayesTools:::.bt_random_effect_summary_correlation_footnotes(
    model_samples = model_samples,
    parameter_names = colnames(model_samples),
    prior_list = list("(mu) cor(a,b)" = correlation_prior),
    coordinates = coordinates,
    included = list("(mu) cor(a,b)" = c(TRUE, TRUE, FALSE, TRUE, FALSE))
  )
  expect_identical(
    footnotes,
    c("(mu) cor(a,b)" = paste0(
      "(mu) cor(a,b): summarised over 2 of 3 draws where the correlation is ",
      "defined."
    ))
  )
  expect_null(BayesTools:::.bt_random_effect_summary_correlation_footnotes(
    model_samples = model_samples,
    parameter_names = colnames(model_samples),
    prior_list = list("(mu) cor(a,b)" = correlation_prior),
    coordinates = coordinates,
    included = list("(mu) cor(a,b)" = c(TRUE, FALSE, FALSE, TRUE, FALSE))
  ))
})

test_that("scaled predictors are centred in terms without a free intercept", {

  skip_if_not_installed("runjags")

  # Maintainer decision M52 (c): `~ 0 + x` with scaled x fits
  # mu = b (x - m) / s, so the original-scale intercept is -b m / s.
  df <- data.frame(
    x = c(3, 7, 3, 7),
    id = factor(c("a", "a", "b", "b"))
  )
  m <- mean(df$x)
  s <- stats::sd(df$x)
  formula_result <- JAGS_formula(
    formula = ~ 0 + x,
    parameter = "mu",
    data = df,
    prior_list = list(x = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE)
  )
  b <- c(0.75, 0.8, 0.7)
  fit <- make_random_scale_table_fit(formula_result, cbind(mu_x = b))

  expect_equal(
    unname(drop(JAGS_evaluate_formula(fit, parameter = "mu", data = data.frame(x = 0)))),
    -b * m / s,
    tolerance = 1e-12
  )
  samples <- JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    return_samples = TRUE
  )
  expect_equal(unname(samples[, "(mu) intercept"]), -b * m / s, tolerance = 1e-12)
  expect_equal(unname(samples[, "(mu) x"]), b / s, tolerance = 1e-12)
  table <- JAGS_estimates_table(fit, transform_scaled = TRUE)
  expect_true("(mu) intercept" %in% rownames(table))
  expect_false("(mu) intercept" %in% rownames(JAGS_estimates_table(fit)))

  # The random analogue `(0 + x | id)` implies the group intercept -u m / s.
  random_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = TRUE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- random_result$formula_design$random_effects[[1]]
  coefficient_names <- BayesTools:::.bt_random_effect_coefficient_names(
    random_term = random_term,
    n_groups = 2L,
    n_columns = 1L
  )
  u <- c(1, -1)
  posterior <- matrix(
    c(0, 2, u),
    nrow = 1L,
    dimnames = list(NULL, c(
      "mu_intercept",
      random_term$sd_parameter_names,
      as.vector(coefficient_names)
    ))
  )
  random_fit <- make_random_scale_table_fit(random_result, posterior)
  offsets <- JAGS_evaluate_formula(
    random_fit,
    parameter = "mu",
    data = data.frame(x = c(0, 0), id = factor(c("a", "b"))),
    formula_target = "conditional"
  )
  expect_equal(unname(drop(offsets)), -u * m / s, tolerance = 1e-12)
})

test_that("independent scaled random slopes imply the documented original-scale correlation", {

  # `(1 + x || id)` is independent on the centred scale. The original-scale
  # intercept u0 - u1 m / s and slope u1 / s have correlation
  # -(t1 m / s) / sqrt(t0^2 + (t1 m / s)^2). Reference: the package's marginal
  # covariance Z G Z' for one group at x = 0 and x = 1, which gives
  # Var(a) = c00, Cov(a, b) = c01 - c00, and Var(b) = c11 - 2 c01 + c00.
  df <- data.frame(
    x = c(4, 6, 5, 7, 3, 5),
    id = factor(c("a", "a", "b", "b", "c", "c"))
  )
  m <- mean(df$x)
  s <- stats::sd(df$x)
  formula_result <- JAGS_formula(
    formula = ~ 1 + (1 + x || id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  sd_draws <- rbind(c(1, 1), c(0.5, 2), c(2, 0.3))
  posterior <- cbind(0, sd_draws)
  colnames(posterior) <- c("mu_intercept", random_term$sd_parameter_names)
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  vcov <- random_effects_marginal_vcov(
    fit = fit,
    parameter = "mu",
    data = data.frame(x = c(0, 1), id = factor(c("a", "a"), levels = c("a", "b", "c"))),
    posterior_samples = posterior,
    prior_list = formula_result$prior_list
  )$samples
  c00 <- vcov[, 1L, 1L]
  c01 <- vcov[, 1L, 2L]
  c11 <- vcov[, 2L, 2L]
  implied <- (c01 - c00) / sqrt(c00 * (c11 - 2 * c01 + c00))
  documented <- -(sd_draws[, 2L] * m / s) /
    sqrt(sd_draws[, 1L]^2 + (sd_draws[, 2L] * m / s)^2)
  expect_equal(unname(implied), documented, tolerance = 1e-12)
  # The worked magnitude in the documentation: equal SDs and m / s = 2.5.
  expect_equal(-2.5 / sqrt(1 + 2.5^2), -0.9284767, tolerance = 1e-7)
})
