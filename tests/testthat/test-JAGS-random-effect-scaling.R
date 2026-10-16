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

test_that("D5 unscaled samples omit only internal latent and group coefficients", {

  skip_if_not_installed("runjags")

  formula_result <- JAGS_formula(
    ~ 1 + x + us(1 + x | id), "mu",
    data.frame(x = c(1, 2, 3, 4), id = factor(c("a", "a", "b", "b"))),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior_spike_and_slab(
        prior("normal", list(0, 1)),
        prior_inclusion = prior("beta", list(2, 2))
      )
    ),
    prior_random = prior_random(id = random_block(
      sd = prior("gamma", list(2, 2)),
      cor = prior_lkj(eta = 1, include_primitives = TRUE),
      monitor = random_monitor(latent = TRUE, coefficients = TRUE,
                               correlation = TRUE)
    ))
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  internal_names <- c("mu__xREx__id_xRE_Zx[1,1]",
                      "mu__xREx__id_xRE_COEFx[1,1]")
  retained_names <- c("mu_x_inclusion", "mu_intercept",
                      random_term$sd_parameter_names,
                      random_term$correlation$primitive_names,
                      "mu_x", "mu_x_indicator")
  posterior_names <- c(internal_names[1L], retained_names,
                       internal_names[2L])
  posterior <- matrix(
    c(.1, .2, .3, .4, .5, .6, .7, .8, .9,
      1, 1.1, 1.2, 1.3, 1.4, 1.5, .2, .4, .6,
      0, 2, 0, 0, 1, 0, .1, .22, .36),
    nrow = 3L, dimnames = list(NULL, posterior_names)
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)
  attr(fit, "formula_scale") <- NULL
  original_fit <- serialize(fit, NULL)
  expected <- matrix(
    c(.4, .5, .6, .7, .8, .9, 1, 1.1, 1.2,
      1.3, 1.4, 1.5, .2, .4, .6, 0, 2, 0, 0, 1, 0),
    nrow = 3L, dimnames = list(NULL, retained_names)
  )

  coordinates <- parameter_coordinates(fit)
  retained_rows <- match(retained_names, coordinates$coordinate_name)
  expect_true(any(coordinates$internal[retained_rows]))
  for(formula_scale in list(NULL, list())){
    samples <- transform_scale_samples(fit, formula_scale)
    expect_true(is.matrix(samples))
    expect_true(is.numeric(samples))
    expect_identical(dim(samples), c(3L, 7L))
    expect_identical(colnames(samples), retained_names)
    expect_identical(samples, expected)
  }
  expect_identical(serialize(fit, NULL), original_fit)
})

test_that("D5 removing every random coordinate retains a zero-column matrix", {

  skip_if_not_installed("runjags")

  formula_result <- JAGS_formula(
    ~ 1 + diag(1 | id), "mu", data.frame(id = factor(c("a", "a", "b", "b"))),
    prior_list = list(intercept = prior("point", list(0))),
    prior_random = prior_random(id = random_block(
      sd = prior("point", list(1)),
      monitor = random_monitor(latent = TRUE, coefficients = TRUE,
                               correlation = FALSE)
    ))
  )
  posterior <- matrix(c(1, 2, 3, 1, 2, 3), nrow = 3L,
                      dimnames = list(NULL, c("mu__xREx__id_xRE_Zx[1,1]",
                                             "mu__xREx__id_xRE_COEFx[1,1]")))
  fit <- make_random_scale_table_fit(formula_result, posterior)
  original_fit <- serialize(fit, NULL)
  samples <- transform_scale_samples(fit, formula_scale = list())

  expect_true(is.matrix(samples))
  expect_true(is.numeric(samples))
  expect_identical(dim(samples), c(3L, 0L))
  expect_identical(samples, matrix(numeric(), nrow = 3L, ncol = 0L,
                                  dimnames = list(NULL, character())))
  expect_identical(serialize(fit, NULL), original_fit)
})

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

test_that("contrast-coefficient random slopes `{j}` unscale with their scaled interactions", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 2, 3, 4, 5, 6),
    g = factor(c("u", "v", "w", "u", "v", "w")),
    id = factor(c("a", "a", "a", "b", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 + g * x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        contrasts = c(g = "meandif"),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  # The factor comes first in the interaction term, so its coefficient index
  # `{j}` trails the scaled covariate: `g__xXx__x{j}`.
  expect_identical(
    unname(random_term$sd_leaves$leaf_terms),
    c("intercept", "g{1}", "g{2}", "x", "g__xXx__x{1}", "g__xXx__x{2}")
  )
  fitted_sd <- c(1, 2, 3, 4, 5, 6)
  posterior <- matrix(
    rep(c(0, fitted_sd), 3),
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

  # Diagonal block: original-scale SDs combine each coefficient with its
  # scaled interaction through the unscaling coefficients -m/s and 1/s.
  scale_info <- formula_result$formula_scale$mu_x
  ratio <- scale_info$mean / scale_info$sd
  expected <- c(
    "(mu) sd(intercept)" = sqrt(1^2 + ratio^2 * 4^2),
    "(mu) sd(g{1})"      = sqrt(2^2 + ratio^2 * 5^2),
    "(mu) sd(g{2})"      = sqrt(3^2 + ratio^2 * 6^2),
    "(mu) sd(x)"         = 4 / scale_info$sd,
    "(mu) sd(g:x{1})"    = 5 / scale_info$sd,
    "(mu) sd(g:x{2})"    = 6 / scale_info$sd
  )
  for(row in names(expected)){
    expect_equal(
      unname(samples[, row]),
      rep(expected[[row]], nrow(samples)),
      tolerance = 1e-12,
      info = row
    )
  }
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

  transformed_matrix <- BayesTools:::.bt_transform_scale_posterior(
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
  # are missing and the row summarizes the 3 defined draws of 5.
  fit <- .undefined_correlation_table_fit(
    source_sd = rbind(c(1, 2), c(0.5, 1.5), c(0.7, 0), c(0.3, 0), c(1, 1)),
    source_rho = c(0.8, -0.3, 0.4, 0.1, 0.2)
  )
  expected_footnote <- c(
    "(mu) cor(intercept,x)" = paste0(
      "(mu) cor(intercept,x): summarized over 3 of 5 draws where the ",
      "correlation is defined, i.e. both SDs are positive."
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
    "summarized over 3 of 5 draws where the correlation is defined, i.e. both SDs are positive.",
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
      "cor(intercept,x): summarized over 3 of 5 draws where the correlation ",
      "is defined, i.e. both SDs are positive."
    )
  )

  # Raw original-scale correlation coordinates are reported per row as well;
  # the Cholesky and LKJ coordinates exist only for positive-definite
  # correlation matrices. On the fitted scale every draw is defined.
  raw <- JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "raw"
  )
  raw_footnotes <- attr(raw, "footnotes")
  correlation_rows <- c(
    "(mu) cor(x,intercept | id)", "(mu) cor(intercept,x | id)",
    "(mu) cor(x,x | id)"
  )
  expect_identical(
    unname(raw_footnotes[correlation_rows]),
    paste0(
      correlation_rows, ": summarized over 3 of 5 draws where the ",
      "correlation is defined, i.e. both SDs are positive."
    )
  )
  # cor(intercept,intercept) needs only the intercept SD, positive in all draws.
  expect_false("(mu) cor(intercept,intercept | id)" %in% names(raw_footnotes))
  coordinate_rows <- setdiff(names(raw_footnotes), correlation_rows)
  expect_setequal(
    coordinate_rows,
    c(
      "(mu) cor_chol(intercept,intercept | id)",
      "(mu) cor_chol(x,intercept | id)",
      "(mu) cor_chol(intercept,x | id)",
      "(mu) cor_chol(x,x | id)",
      "(mu) lkj_u(intercept,x | id)",
      "(mu) lkj_cpc(intercept,x | id)"
    )
  )
  expect_true(all(grepl(
    ": summarized over 3 of 5 draws where the correlation matrix is positive definite.$",
    raw_footnotes[coordinate_rows]
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

test_that("original-scale correlations of singular draws with positive SDs are +-1", {

  skip_if_not_installed("runjags")

  # With a zero scaled intercept SD and a positive scaled slope SD, the
  # original-scale intercept -u1 m / s and slope u1 / s are perfectly
  # correlated: cor = -sign(m) = -1 (m = mean(x) > 0), although the covariance
  # is singular. Only a zero original-scale SD leaves the correlation missing.
  source_sd <- rbind(c(1, 2), c(0, 1), c(0, 2), c(0.7, 0), c(0, 0))
  source_rho <- c(0.8, 0.5, -0.7, 0.4, 0.3)
  fit <- .undefined_correlation_table_fit(source_sd, source_rho)
  scale_info <- attr(fit, "formula_scale")$mu$mu_x
  expect_gt(scale_info$mean, 0)
  M <- matrix(
    c(1, -scale_info$mean / scale_info$sd, 0, 1 / scale_info$sd),
    nrow = 2,
    byrow = TRUE
  )
  source_cov <- diag(source_sd[1, ]) %*%
    matrix(c(1, source_rho[1], source_rho[1], 1), 2, 2) %*%
    diag(source_sd[1, ])
  expected_cov <- M %*% source_cov %*% t(M)
  expected <- c(
    expected_cov[1, 2] / sqrt(prod(diag(expected_cov))),
    -1, -1, NA, NA
  )

  draws <- parameter_draws(
    fit,
    parameter_catalog_resolve(parameter_catalog(fit), "(mu) cor(intercept,x)")
  )
  values <- as.numeric(as.matrix(draws[[1]]))
  expect_equal(values, expected, tolerance = 1e-12)
  expect_true(all(abs(values) <= 1, na.rm = TRUE))

  for(transform_scaled in c(FALSE, TRUE)){
    table <- JAGS_estimates_table(fit, transform_scaled = transform_scaled)
    expect_identical(
      unname(attr(table, "footnotes")),
      paste0(
        "(mu) cor(intercept,x): summarized over 3 of 5 draws where the ",
        "correlation is defined, i.e. both SDs are positive."
      )
    )
  }

  # The monitored correlation matrix follows the same rule entrywise; the
  # Cholesky factor and LKJ primitives exist only for positive-definite
  # correlation matrices.
  transformed <- transform_scale_samples(fit)
  R_21 <- "mu__xREx__id_xRE_CORx_R[2,1]"
  R_11 <- "mu__xREx__id_xRE_CORx_R[1,1]"
  L_21 <- "mu__xREx__id_xRE_CORx_L[2,1]"
  u_1  <- "mu__xREx__id_xRE_CORx_lkj_u[1]"
  expect_equal(unname(transformed[, R_21]), expected, tolerance = 1e-12)
  expect_equal(unname(transformed[, R_11]), c(1, 1, 1, 1, NA))
  expect_equal(which(is.na(transformed[, L_21])), 2:5)
  expect_equal(which(is.na(transformed[, u_1])), 2:5)

  # Depending on the slope SD, the perfect correlation rounds to exactly -1 or
  # to just inside it, where chol() succeeds. Singularity is decided from the
  # zero scaled SD, so the Cholesky and LKJ coordinates are missing in every
  # singular draw.
  slope_sd <- exp(seq(log(0.01), log(100), length.out = 100))
  singular <- transform_scale_samples(.undefined_correlation_table_fit(
    source_sd = cbind(0, slope_sd),
    source_rho = rep(0.37, length(slope_sd))
  ))
  expect_lte(max(abs(singular[, R_21] + 1)), 1e-12)
  expect_true(all(is.na(singular[, L_21])))
  expect_true(all(is.na(singular[, u_1])))
})

test_that("original-scale random-effect SDs beyond the double range are refused", {

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
  sd_names <- formula_result$formula_design$random_effects[[1L]]$sd_parameter_names
  transformed <- function(sd_values){
    posterior <- matrix(
      c(0, sd_values), nrow = 1L,
      dimnames = list(NULL, c("mu_intercept", sd_names))
    )
    transform_scale_samples(make_random_scale_table_fit(formula_result, posterior))
  }

  # The SDs are squared on the way: 1e100 stays representable, 1e200 and 1e-200
  # do not, and a zero SD stays a zero SD.
  large <- transformed(c(1e100, 1e100))
  expect_true(all(is.finite(large[, sd_names])))
  expect_true(all(large[, sd_names] > 0))
  expect_equal(unname(transformed(c(0, 0))[, sd_names]), c(0, 0))
  for(sd_value in c(1e200, 1e-200)){
    refused <- tryCatch(transformed(c(sd_value, sd_value)), error = identity)
    expect_s3_class(refused, "BayesTools_formula_transform_unavailable")
    expect_identical(refused$reason, "nonfinite_transform")
    expect_match(
      conditionMessage(refused),
      paste0("Original-scale random-effect SDs of block 'id' are unavailable: ",
             "the transformed SD of 1 draw is not representable in double ",
             "precision although the fitted-scale SDs are finite."),
      fixed = TRUE
    )
  }

  # With a correlation the same overflow is refused; a representable draw next
  # to it is not.
  expect_s3_class(
    tryCatch(
      transform_scale_samples(.undefined_correlation_table_fit(
        source_sd = rbind(c(1, 1), c(1e200, 1e200)), source_rho = c(0.3, 0.3)
      )),
      error = identity
    ),
    "BayesTools_formula_transform_unavailable"
  )
  expect_no_error(transform_scale_samples(.undefined_correlation_table_fit(
    source_sd = rbind(c(1, 1), c(1e100, 1e100)), source_rho = c(0.3, 0.3)
  )))
})

test_that("larger blocks define correlations entrywise in singular draws", {

  skip_if_not_installed("runjags")

  # (1 + x + z | id) with a zero scaled intercept SD. Draw 1: the
  # original-scale intercept is a linear combination of both slopes, so the
  # covariance has rank 2, yet every SD is positive and every correlation is
  # defined. Draw 2 also has a zero z-slope SD: correlations involving z are
  # missing, while the intercept is -u1 m_x / s_x, perfectly correlated with
  # the x slope (cor = -1).
  df <- data.frame(
    x = c(4, 5, 6, 8, 3, 7),
    z = c(10, 14, 11, 16, 12, 9),
    id = factor(c("a", "a", "b", "b", "c", "c"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + z + random(1 + x + z | id, name = "id", covariance = "us"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      z = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE, z = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("normal", list(0, 1), truncation = list(0, Inf)),
        cor = prior_lkj(eta = 1),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = TRUE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  R_names <- outer(
    seq_len(3),
    seq_len(3),
    Vectorize(function(row, column){
      paste0(random_term$parameter_stem, "_xRE_CORx_R[", row, ",", column, "]")
    })
  )
  L_names <- BayesTools:::.bt_random_effect_cholesky_names(random_term, 3L)
  source_sd <- rbind(c(0, 1, 2), c(0, 1, 0))
  source_cor <- matrix(
    c(1, 0.3, -0.2, 0.3, 1, 0.4, -0.2, 0.4, 1),
    nrow = 3
  )
  posterior <- t(vapply(seq_len(2), function(draw){
    c(0, 0, 0, source_sd[draw, ], as.vector(source_cor),
      as.vector(t(chol(source_cor))))
  }, numeric(24)))
  colnames(posterior) <- c(
    "mu_intercept", "mu_x", "mu_z", random_term$sd_parameter_names,
    as.vector(R_names), as.vector(L_names)
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  x_scale <- formula_result$formula_scale$mu_x
  z_scale <- formula_result$formula_scale$mu_z
  M <- rbind(
    c(1, -x_scale$mean / x_scale$sd, -z_scale$mean / z_scale$sd),
    c(0, 1 / x_scale$sd, 0),
    c(0, 0, 1 / z_scale$sd)
  )
  expected_cov <- M %*% diag(source_sd[1, ]) %*% source_cor %*%
    diag(source_sd[1, ]) %*% t(M)
  expected_cor <- stats::cov2cor(expected_cov)
  expect_lt(abs(det(expected_cor)), 1e-12)
  expected <- list(
    "intercept,x" = c(expected_cor[1, 2], -1),
    "intercept,z" = c(expected_cor[1, 3], NA),
    "x,z"         = c(expected_cor[2, 3], NA)
  )

  for(pair in names(expected)){
    name <- paste0("(mu) cor(", pair, ")")
    value <- as.numeric(as.matrix(parameter_draws(
      fit,
      parameter_catalog_resolve(parameter_catalog(fit), name)
    )[[1]]))
    expect_equal(value, expected[[pair]], tolerance = 1e-12, info = name)
  }
  footnotes <- attr(
    JAGS_estimates_table(fit, remove_diagnostics = TRUE),
    "footnotes"
  )
  expect_identical(
    unname(footnotes),
    paste0(
      c("(mu) cor(intercept,z)", "(mu) cor(x,z)"),
      ": summarized over 1 of 2 draws where the correlation is defined, ",
      "i.e. both SDs are positive."
    )
  )
})

test_that("original-scale correlations are defined entrywise and clamp rounding only", {

  correlation <- BayesTools:::.random_sd_transformed_correlation(
    covariance = matrix(c(4, 2, 0, 2, 1, 0, 0, 0, 0), nrow = 3),
    sd = c(2, 1, 0),
    group_key = "id"
  )
  expect_identical(correlation[1:2, 1:2], matrix(1, 2, 2))
  expect_true(all(is.na(correlation[3, ])))
  expect_true(all(is.na(correlation[, 3])))

  # A rounding excess of at most 1e-8 is set to +-1.
  rounded <- BayesTools:::.random_sd_transformed_correlation(
    covariance = matrix(c(1, -(1 + 1e-12), -(1 + 1e-12), 1), nrow = 2),
    sd = c(1, 1),
    group_key = "id"
  )
  expect_identical(rounded[1, 2], -1)
  expect_error(
    BayesTools:::.random_sd_transformed_correlation(
      covariance = matrix(c(1, 1.5, 1.5, 1), nrow = 2),
      sd = c(1, 1),
      group_key = "id"
    ),
    paste0(
      "Internal error: an original-scale random-effect correlation of block ",
      "'id' exceeds 1 in absolute value by 0.5."
    ),
    fixed = TRUE
  )
})

test_that("parameter_draws declares undefined correlation draws for ensemble tables", {

  skip_if_not_installed("runjags")

  # Zero slope SDs in draws 3 and 4 leave the original-scale correlation
  # undefined (NA).
  fit <- .undefined_correlation_table_fit(
    source_sd = rbind(c(1, 2), c(0.5, 1.5), c(0.7, 0), c(0.3, 0), c(1, 1)),
    source_rho = c(0.8, -0.3, 0.4, 0.1, 0.2)
  )

  # Built as RoBMA's multivariate heterogeneity summary does
  # (RoBMA 04d283e9, R/heterogeneity-mv.R:491-528,
  # .brma_mv_correlation_sample_lists(): catalog rows with quantity "cor"
  # owned by a random block, resolved and drawn with
  # parameter_draws(model_samples = posterior_samples) at :515, one numeric
  # vector per correlation at :525; .summary_heterogeneity_brma_mv_one() at
  # :459-487 appends them to tau and tau2 and calls ensemble_estimates_table()
  # at :473).
  catalog <- parameter_catalog(fit)
  quantities <- catalog$quantities
  selected <- quantities[
    quantities$namespace == "mu" &
      quantities$owner_type == "random_block" &
      quantities$quantity == "cor" &
      !quantities$internal &
      quantities$status != "unavailable",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(selected), 1L)
  posterior_samples <- as.matrix(fit$mcmc[[1L]])
  selection <- parameter_catalog_resolve(
    catalog = catalog,
    alias = selected$canonical_name[1L],
    namespace = selected$namespace[1L]
  )
  draws <- parameter_draws(fit, selection, model_samples = posterior_samples)
  expect_identical(
    .bt_meta_get(draws, "undefined_draws"),
    c("(mu) cor(intercept,x)" = "correlation")
  )
  expect_identical(
    .bt_meta_get(parameter_draws(fit, selection), "undefined_draws"),
    .bt_meta_get(draws, "undefined_draws")
  )
  sd_samples <- JAGS_estimates_table(fit, return_samples = TRUE)[
    , c("(mu) sd(intercept)", "(mu) sd(x)")
  ]
  var_samples <- rowMeans(sd_samples^2)

  # RoBMA's current extraction, as.numeric(draws[[1L]][, 1L])
  # (heterogeneity-mv.R:525 and :749), drops the declaration: the table stops clearly.
  undeclared <- as.numeric(draws[[1L]][, 1L])
  expect_error(
    ensemble_estimates_table(
      samples = list(tau = sqrt(var_samples), tau2 = var_samples,
                     "cor(intercept,x)" = undeclared),
      parameters = c("tau", "tau2", "cor(intercept,x)")
    ),
    "The posterior draws of 'cor(intercept,x)' contain missing values.",
    fixed = TRUE
  )

  # Keeping the declaration on the extracted vector summarizes the defined
  # draws and footnotes their share.
  declared <- as.numeric(draws[[1L]][, 1L])
  declared <- .bt_meta_set(declared, "undefined_draws", .bt_meta_get(draws, "undefined_draws")[[1L]])
  samples_list <- list(tau = sqrt(var_samples), tau2 = var_samples,
                       "cor(intercept,x)" = declared)
  estimates <- ensemble_estimates_table(
    samples = samples_list,
    parameters = names(samples_list),
    probs = c(.025, .975),
    title = "Heterogeneity Estimates (id):"
  )
  defined <- declared[!is.na(declared)]
  expect_equal(length(defined), 3L)
  expect_equal(estimates["cor(intercept,x)", "Mean"], mean(defined))
  expect_equal(estimates["cor(intercept,x)", "Median"], stats::median(defined))
  expect_identical(
    unname(attr(estimates, "footnotes")),
    paste0(
      "cor(intercept,x): summarized over 3 of 5 draws where the correlation ",
      "is defined, i.e. both SDs are positive."
    )
  )
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
      "(mu) cor(a,b): summarized over 2 of 3 draws where the correlation is ",
      "defined, i.e. both SDs are positive."
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

test_that("scaled predictors are centered in terms without a free intercept", {

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

  # `(1 + x || id)` is independent on the centered scale. The original-scale
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

test_that("a hand-built formula_scale unscales random-effect SDs like the automatic fit", {

  set.seed(3)
  df <- data.frame(id = factor(rep(letters[1:8], each = 6)), x = rnorm(48, mean = 10, sd = 3))
  m <- mean(df$x)
  s <- stats::sd(df$x)
  df_manual <- df
  df_manual$x <- (df$x - m) / s
  manual_scale <- list(mu = list(mu_x = list(mean = m, sd = s)))
  priors <- list(intercept = prior("normal", list(0, 5)), x = prior("normal", list(0, 2)))
  fit_of <- function(carrier){
    fit <- coda::mcmc.list(coda::mcmc(as.matrix(carrier)))
    class(fit) <- c("BayesTools_fit", class(fit))
    for(name in c("prior_list", "formula_design", "formula_scale")){
      attr(fit, name) <- attr(carrier, name)
    }
    fit <- .bt_attach_parameter_map(fit)
    fit <- .bt_attach_draw_geometry(fit)
    .bt_attach_fit_contract(fit)
  }

  # independent block: the standardized SDs map to the original-scale SDs
  formula <- ~ 1 + x + diag(1 + x | id)
  random_priors <- prior_random(id = random_block(
    sd = prior("normal", list(0, 1), list(0, Inf)),
    monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = FALSE)
  ))
  compiled <- JAGS_formula(formula, "mu", df_manual, priors, prior_random = random_priors)
  draws <- .generate_prior_sample_matrix(compiled$prior_list, 4L, seed = 7)
  manual_fit <- fit_of(JAGS_formula_draws(draws, formula, "mu", df_manual, priors,
    prior_random = random_priors))
  auto_fit <- fit_of(JAGS_formula_draws(draws, formula, "mu", df, priors,
    prior_random = random_priors, formula_scale = list(x = TRUE)))
  transformed <- transform_scale_samples(manual_fit, formula_scale = manual_scale)
  sd_columns <- grep("__xREx__", colnames(draws), value = TRUE)
  expect_length(sd_columns, 2L)
  expect_false(isTRUE(all.equal(unname(transformed[, sd_columns]), unname(draws[, sd_columns]))))
  expect_equal(transformed, transform_scale_samples(auto_fit), tolerance = 1e-12)

  # correlated block with a point intercept: the completion carries the same
  # structure that the automatic scaling attaches
  formula <- ~ 1 + x + us(1 + x | id)
  priors$intercept <- prior("point", list(1.5))
  random_priors <- prior_random(id = random_block(
    sd = prior("gamma", list(2, 2)), cor = prior_lkj(eta = 1.5)
  ))
  automatic <- JAGS_formula(formula, "mu", df, priors, prior_random = random_priors,
    formula_scale = list(x = TRUE))
  manual <- JAGS_formula(formula, "mu", df_manual, priors, prior_random = random_priors)
  scale <- .bt_formula_scale_finalize(manual$formula_scale, manual$formula_design,
    prior_list = manual$prior_list, owner_scope = "fit")
  design <- manual$formula_design
  design$formula_scale <- scale
  fit <- structure(list(), class = "BayesTools_fit")
  attr(fit, "formula_scale") <- list(mu = scale)
  attr(fit, "formula_design") <- list(mu = design)
  completed <- .bt_formula_scale_list_complete(manual_scale, fit)$mu
  for(key in c("parameter", "log_intercept", "point_terms",
               "random_effect_sd_leaves", "random_effect_correlation_required")){
    expect_false(is.null(attr(automatic$formula_scale, key, exact = TRUE)))
    expect_identical(attr(completed, key, exact = TRUE), attr(automatic$formula_scale, key, exact = TRUE))
  }
  expect_identical(attr(completed, "random_effect_correlation_required"), "__xREx__id")
})

# Independent reference for the original-scale map of a random-effect block:
# the block's design on a grid of original-scale predictor values and on the
# standardized values satisfies Z_s = Z_o A, so the block coefficients transform
# as u_o = A u_s and their covariance as A Sigma A'.
random_block_unscale_reference <- function(rhs, grid, m, s, contrasts = list()){

  grid_scaled <- grid
  grid_scaled$x <- (grid$x - m) / s
  for(factor_name in names(contrasts)){
    stats::contrasts(grid[[factor_name]]) <- contrasts[[factor_name]]
    stats::contrasts(grid_scaled[[factor_name]]) <- contrasts[[factor_name]]
  }
  Z_original <- stats::model.matrix(rhs, grid)
  Z_scaled <- stats::model.matrix(rhs, grid_scaled)
  A <- qr.solve(Z_original, Z_scaled)
  stopifnot(max(abs(Z_original %*% A - Z_scaled)) < 1e-10)

  A
}

random_block_scale_fit <- function(formula, data, prior_random, extra = NULL, tau = NULL){

  formula_result <- JAGS_formula(
    formula, "mu", data,
    prior_list = list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random
  )
  sd_names <- unlist(lapply(formula_result$formula_design$random_effects, function(term){
    unique(term$sd_parameter_names)
  }))
  if(is.null(tau)) tau <- seq_along(sd_names) / 2
  posterior <- matrix(
    rep(c(0, 0, tau, extra$values), 3L), nrow = 3L, byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", "mu_x", sd_names, extra$names))
  )

  list(
    formula_result = formula_result,
    fit = make_random_scale_table_fit(formula_result, posterior),
    sd_names = sd_names,
    tau = tau
  )
}

random_block_scale_data <- function(){

  set.seed(5)
  data.frame(
    x  = rnorm(36, mean = 10, sd = 2),
    f  = factor(rep(c("a", "b", "c"), length.out = 36)),
    id = factor(rep(paste0("g", 1:6), each = 6))
  )
}

random_block_scale_monitor <- function(correlation = FALSE){

  random_monitor(latent = FALSE, coefficients = FALSE, correlation = correlation)
}

test_that("tiny predictor units retain the original-scale SD map and unit conversion", {

  skip_if_not_installed("runjags")

  x <- as.numeric(scale(rep(c(-3, -2, -1, 1, 2, 3), 6L))) * 1e-12 + 5e-11
  results <- lapply(c(1, 1e12), function(units){
    data <- data.frame(x = x * units, id = factor(rep(1:6, each = 6L)))
    m <- mean(data$x)
    s <- stats::sd(data$x)
    setup <- random_block_scale_fit(
      ~ 1 + x + diag(1 + x | id), data,
      prior_random(id = random_block(
        sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
      )), tau = c(1, 1)
    )
    term <- setup$formula_result$formula_design$random_effects[[1L]]
    A <- BayesTools:::.random_sd_unscale_matrix(
      term$sd_parameter_names, setup$formula_result$formula_scale, "mu",
      random_term = term, source_data = data
    )
    reference <- matrix(c(1, 0, -m / s, 1 / s), 2L)
    expect_equal(A, reference, tolerance = 1e-14)
    expect_equal(cbind(1, data$x) %*% A, cbind(1, (data$x - m) / s),
                 tolerance = 1e-13, ignore_attr = TRUE)
    expected <- c(sqrt(1 + (m / s)^2), 1 / s)
    transformed <- transform_scale_samples(setup$fit)
    mixed <- as_mixed_posteriors(setup$fit, setup$sd_names, transform_scaled = TRUE)
    expect_equal(unname(transformed[1L, setup$sd_names]), expected, tolerance = 1e-13)
    expect_equal(unname(vapply(mixed, `[`, numeric(1), 1L)), expected, tolerance = 1e-13)
    unname(transformed[1L, setup$sd_names])
  })
  expect_equal(results[[1L]][1L], sqrt(2501), tolerance = 1e-13)
  expect_equal(results[[1L]][1L], results[[2L]][1L], tolerance = 1e-13)
  expect_equal(results[[1L]][2L] / 1e12, results[[2L]][2L], tolerance = 1e-13)
})

test_that("tiny-unit correlated blocks retain the covariance cross term", {

  skip_if_not_installed("runjags")

  data <- data.frame(
    x = as.numeric(scale(rep(c(-3, -2, -1, 1, 2, 3), 6L))) * 1e-12 + 5e-11,
    id = factor(rep(1:6, each = 6L))
  )
  rho <- 0.4
  tau <- c(0.8, 1.3)
  correlation_names <- outer(1:2, 1:2, Vectorize(function(row, column){
    paste0("mu__xREx__id_xRE_CORx_R[", row, ",", column, "]")
  }))
  setup <- random_block_scale_fit(
    ~ 1 + x + us(1 + x | id), data,
    prior_random(id = random_block(
      sd = prior("gamma", list(2, 2)), cor = prior_lkj(eta = 1),
      monitor = random_block_scale_monitor(correlation = TRUE)
    )),
    extra = list(values = c(1, rho, rho, 1), names = as.vector(correlation_names)),
    tau = tau
  )
  s <- stats::sd(data$x)
  ratio <- mean(data$x) / s
  # Direct entries of A Sigma A', independent of the implementation's map.
  v00 <- tau[1L]^2 - 2 * ratio * rho * tau[1L] * tau[2L] + ratio^2 * tau[2L]^2
  v11 <- tau[2L]^2 / s^2
  v01 <- (rho * tau[1L] * tau[2L] - ratio * tau[2L]^2) / s
  expected_sd <- sqrt(c(v00, v11))
  expected_rho <- v01 / prod(expected_sd)
  transformed <- transform_scale_samples(setup$fit)
  expect_equal(unname(transformed[1L, setup$sd_names]), expected_sd, tolerance = 1e-13)
  expect_equal(unname(transformed[1L, correlation_names[1L, 2L]]), expected_rho,
               tolerance = 1e-13)
  mixed <- as_mixed_posteriors(setup$fit, setup$sd_names, transform_scaled = TRUE)
  expect_equal(unname(vapply(mixed, `[`, numeric(1), 1L)), expected_sd, tolerance = 1e-13)
})

test_that("random-effect SDs of a block with a factor-by-scaled-slope term follow the block design", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  m <- mean(df$x)
  s <- stats::sd(df$x)
  grid <- expand.grid(f = factor(c("a", "b", "c")), x = c(-1, 0, 1, 2))
  block <- function(...){
    prior_random(id = random_block(
      sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor(), ...
    ))
  }

  # `f` enters the block through its main effect and through its interaction
  # with the scaled slope: the induced shift of the intercept and of the factor
  # levels (u0 - m / s times the slopes by level) is not carried by coefficient
  # names.
  cases <- list(
    list(~ 1 + x + diag(1 + f + f:x | id), ~ f + f:x, block(), list()),
    list(~ 1 + x + diag(1 + f + f:x | id), ~ f + f:x, block(contrasts = c(f = "meandif")),
         list(f = contr.meandif(3))),
    list(~ 1 + x + diag(1 + f * x | id), ~ f * x, block(), list()),
    list(~ 1 + x + diag(1 + x | id), ~ x, block(), list())
  )
  for(case in cases){
    setup <- random_block_scale_fit(case[[1L]], df, case[[3L]])
    A <- random_block_unscale_reference(case[[2L]], grid, m, s, case[[4L]])
    expected <- unname(sqrt(diag(A %*% diag(setup$tau^2, nrow = length(setup$tau)) %*% t(A))))
    transformed <- transform_scale_samples(setup$fit)
    expect_equal(
      unname(transformed[1L, setup$sd_names]), expected, tolerance = 1e-10,
      info = paste(deparse(case[[1L]]), paste(names(case[[4L]]), collapse = ","))
    )
  }

  # the analytic intercept SD of centered sum-to-zero level contrasts is
  # sqrt(tau0^2 + (m / s)^2 * sum(tau_slopes_by_level^2) / 9)
  setup <- random_block_scale_fit(~ 1 + x + diag(1 + f + f:x | id), df, block(contrasts = c(f = "meandif")))
  expect_equal(
    unname(transform_scale_samples(setup$fit)[1L, setup$sd_names[1L]]),
    sqrt(setup$tau[1L]^2 + (m / s)^2 * sum(setup$tau[4:6]^2) / 9),
    tolerance = 1e-12
  )

  # hand-checked case: six rows, SDs 1:6, intercept SD 5.5627730894
  df_small <- data.frame(
    x = 1:6, f = factor(c("a", "b", "c", "a", "b", "c")),
    id = factor(c("u", "u", "u", "v", "v", "v"))
  )
  small <- random_block_scale_fit(
    ~ 1 + x + diag(1 + f + f:x | id), df_small, block(contrasts = c(f = "meandif")), tau = 1:6
  )
  expect_equal(
    unname(transform_scale_samples(small$fit)[1L, small$sd_names[1L]]),
    5.562773089426212, tolerance = 1e-12
  )
})

test_that("every consumer of the original-scale random-effect SDs uses the block design", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  grid <- expand.grid(f = factor(c("a", "b", "c")), x = c(-1, 0, 1, 2))
  setup <- random_block_scale_fit(
    ~ 1 + x + diag(1 + f + f:x | id), df,
    prior_random(id = random_block(
      sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
    ))
  )
  A <- random_block_unscale_reference(~ f + f:x, grid, mean(df$x), stats::sd(df$x))
  expected <- unname(sqrt(diag(A %*% diag(setup$tau^2) %*% t(A))))
  expect_gt(abs(expected[1L] - setup$tau[1L]), 5)

  # JAGS_estimates_table(transform_scaled = TRUE)
  samples <- JAGS_estimates_table(
    setup$fit, transform_scaled = TRUE, remove_diagnostics = TRUE, return_samples = TRUE
  )
  expect_equal(
    unname(samples[1L, grep("sd(", colnames(samples), fixed = TRUE)]),
    expected, tolerance = 1e-10
  )

  # random-effect summaries through the parameter catalog
  catalog <- parameter_catalog(setup$fit)
  draws <- parameter_draws(setup$fit, parameter_catalog_resolve(catalog, "(mu) sd(f[b])"))
  expect_equal(unname(drop(as.matrix(draws))[1L]), expected[2L], tolerance = 1e-10)

  # mixed posteriors of the SD parameters
  mixed <- as_mixed_posteriors(setup$fit, "mu__xREx__id_f", transform_scaled = TRUE)
  expect_equal(unname(mixed$mu__xREx__id_f[1L, ]), expected[2:3], tolerance = 1e-10)

  # prior draws on the original scale
  prior_draws <- transform_prior_samples(setup$fit, n_samples = 6L, seed = 3)
  raw <- .generate_prior_sample_matrix(setup$formula_result$prior_list, 6L, column_names = NULL, seed = 3)
  expect_equal(
    unname(prior_draws[, setup$sd_names]),
    unname(sqrt(raw[, setup$sd_names]^2 %*% t(A^2))),
    tolerance = 1e-10
  )
})

test_that("a correlated block with a factor-by-scaled-slope term transforms its covariance by the block map", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  grid <- expand.grid(f = factor(c("a", "b", "c")), x = c(-1, 0, 1, 2))
  K <- 6L
  correlation <- matrix(0.3, K, K)
  diag(correlation) <- 1
  correlation_names <- outer(seq_len(K), seq_len(K), Vectorize(function(row, column){
    paste0("mu__xREx__id_xRE_CORx_R[", row, ",", column, "]")
  }))
  setup <- random_block_scale_fit(
    ~ 1 + x + us(1 + f + f:x | id), df,
    prior_random(id = random_block(
      sd = prior("gamma", list(2, 2)), cor = prior_lkj(eta = 1),
      monitor = random_block_scale_monitor(correlation = TRUE)
    )),
    extra = list(values = as.vector(correlation), names = as.vector(correlation_names))
  )
  A <- random_block_unscale_reference(~ f + f:x, grid, mean(df$x), stats::sd(df$x))
  covariance <- A %*% diag(setup$tau) %*% correlation %*% diag(setup$tau) %*% t(A)
  transformed <- transform_scale_samples(setup$fit)

  expect_equal(unname(transformed[1L, setup$sd_names]), unname(sqrt(diag(covariance))), tolerance = 1e-10)
  expect_equal(
    unname(matrix(transformed[1L, as.vector(correlation_names)], K, K)),
    unname(stats::cov2cor(covariance)),
    tolerance = 1e-10
  )
})

test_that("blocks that share a grouping variable use their own design", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  grid <- expand.grid(f = factor(c("a", "b", "c")), x = c(-1, 0, 1, 2))
  setup <- random_block_scale_fit(
    ~ 1 + x + random(1 + x | id, name = "slope", covariance = "diag") +
      random(1 + f + f:x | id, name = "level", covariance = "diag"),
    df,
    prior_random(
      slope = random_block(sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()),
      level = random_block(sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor())
    )
  )
  m <- mean(df$x)
  s <- stats::sd(df$x)
  A_slope <- random_block_unscale_reference(~ x, grid, m, s)
  A_level <- random_block_unscale_reference(~ f + f:x, grid, m, s)
  transformed <- transform_scale_samples(setup$fit)

  expect_equal(
    unname(transformed[1L, setup$sd_names]),
    unname(c(sqrt(diag(A_slope %*% diag(setup$tau[1:2]^2) %*% t(A_slope))),
      sqrt(diag(A_level %*% diag(setup$tau[3:8]^2) %*% t(A_level))))),
    tolerance = 1e-10
  )
})

test_that("a block with no original-scale representation keeps the documented drop policy", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  s <- stats::sd(df$x)
  block <- prior_random(id = random_block(
    sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
  ))

  # `f:x` needs the per-level indicators, which no column of the block carries:
  # the SD of the realized intercept stays the fitted one (documented), and
  # the slopes by level are divided by the SD of the predictor.
  setup <- random_block_scale_fit(~ 1 + x + diag(1 + f:x | id), df, block)
  expect_equal(
    unname(transform_scale_samples(setup$fit)[1L, setup$sd_names]),
    c(setup$tau[1L], setup$tau[-1L] / s), tolerance = 1e-12
  )
  # `(0 + x | g)` implies an intercept that no column carries
  setup <- random_block_scale_fit(~ 1 + x + diag(0 + x | id), df, block)
  expect_equal(
    unname(transform_scale_samples(setup$fit)[1L, setup$sd_names]),
    setup$tau / s, tolerance = 1e-12
  )

  # any other failure of the block design propagates
  setup <- random_block_scale_fit(~ 1 + x + diag(1 + f + f:x | id), df, block)
  design <- attr(setup$fit, "formula_design")
  term <- design$mu$random_effects[[1L]]
  term$raw_column_names <- rev(term$raw_column_names)
  design$mu$random_effects[[1L]] <- term
  attr(setup$fit, "formula_design") <- design
  expect_error(
    transform_scale_samples(setup$fit),
    "do not reproduce the fitted design", fixed = TRUE
  )
})

test_that("a logical predictor in a block with a scaled slope gets design-derived SDs", {

  skip_if_not_installed("runjags")

  # `flag` enters the block through its main effect and through `flag:x`; the
  # interaction is indicator coded (flagFALSE:x, flagTRUE:x) and shares one SD.
  # The least-squares map gives the two columns equal weights only up to
  # rounding.
  df <- random_block_scale_data()
  df$flag <- rep(c(TRUE, FALSE), 18)
  grid <- expand.grid(flag = c(FALSE, TRUE), x = c(-1, 0, 1, 2))
  setup <- random_block_scale_fit(
    ~ 1 + x + diag(1 + flag + flag:x | id), df,
    prior_random(id = random_block(
      sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
    ))
  )
  A <- random_block_unscale_reference(~ flag + flag:x, grid, mean(df$x), stats::sd(df$x))
  expect_identical(nrow(A), 4L)
  expect_identical(length(setup$sd_names), 3L)
  tau_columns <- setup$tau[c(1L, 2L, 3L, 3L)]
  expected <- unname(sqrt(diag(A %*% diag(tau_columns^2) %*% t(A))))
  # the two indicator columns of the shared SD transform to the same SD
  expect_equal(expected[3L], expected[4L], tolerance = 1e-12)
  expect_gt(abs(expected[1L] - setup$tau[1L]), 5)

  transformed <- transform_scale_samples(setup$fit)
  expect_equal(unname(transformed[1L, setup$sd_names]), expected[1:3], tolerance = 1e-10)

  # the default estimates table returns the design-derived SDs
  samples <- JAGS_estimates_table(setup$fit, remove_diagnostics = TRUE, return_samples = TRUE)
  expect_equal(
    unname(samples[1L, grep("sd(", colnames(samples), fixed = TRUE)]),
    expected[1:3], tolerance = 1e-10
  )
  table <- as.data.frame(JAGS_estimates_table(setup$fit))
  expect_equal(
    unname(table[grep("sd(", rownames(table), fixed = TRUE), "Mean"]),
    expected[1:3], tolerance = 1e-10
  )
})

test_that("the weights of a shared random-effect SD leaf are compared relative to their magnitude", {

  skip_if_not_installed("runjags")

  # two columns that share leaf "a" (target rows 1 and 2) next to leaf "b" (column 3);
  # `first` and `second` are the weights of the two target rows on the two leaves
  invariant <- function(first, second){
    M <- rbind(
      c(sqrt(first[1L]), 0, sqrt(first[2L])),
      c(sqrt(second[1L]), 0, sqrt(second[2L])),
      0
    )
    .random_sd_shared_leaf_is_invariant(
      M, c("a", "a", "b"), target_columns = 1:2, correlated = FALSE
    )
  }

  # an absolute tolerance would merge weights of order 1e-14
  expect_false(invariant(c(1e-14, 1e-14), c(2e-14, 1e-14)))
  expect_false(invariant(c(1e-14, 2e-14), c(2e-14, 4e-14)))
  expect_true(invariant(c(1, 2) * 1e-14, c(1, 2) * 1e-14 * (1 + 1e-15)))
  # equal weights up to rounding of a few ulps, whatever their magnitude
  expect_true(invariant(c(1, 2), c(1, 2) * (1 + 1e-15)))
  expect_true(invariant(c(1, 2) * 1e6, c(1, 2) * 1e6 * (1 - 1e-15)))
  expect_false(invariant(c(1, 2), c(1, 2) * (1 + 1e-9)))
  # all-zero weights are equal
  expect_true(invariant(c(0, 0), c(0, 0)))

  # a shared leaf whose columns genuinely transform to different SDs is refused
  # also when the weights are of order 1 / sd(x)^2 for a predictor with a very
  # large SD (`flagFALSE:x` and `flagTRUE:x` share one SD, and only the TRUE
  # level picks up the induced term of `flag:x:z`)
  shared_leaf_samples <- function(x_sd){
    set.seed(11)
    df <- data.frame(
      x  = rnorm(48, 3 * x_sd, x_sd),
      z  = rnorm(48, 1, 1),
      id = factor(rep(paste0("g", 1:8), each = 6)),
      flag = rep(c(TRUE, FALSE), 24)
    )
    formula_result <- JAGS_formula(
      ~ 1 + x + z + diag(1 + flag + z + flag:z + flag:x + x:z + flag:x:z | id), "mu", df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1)),
        z = prior("normal", list(0, 1))
      ),
      formula_scale = list(x = TRUE, z = TRUE),
      prior_random = prior_random(id = random_block(
        sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
      ))
    )
    sd_names <- unique(unlist(lapply(formula_result$formula_design$random_effects, function(term){
      term$sd_parameter_names
    })))
    posterior <- matrix(
      rep(c(0, 0, 0, seq_along(sd_names) / 2), 3L), nrow = 3L, byrow = TRUE,
      dimnames = list(NULL, c("mu_intercept", "mu_x", "mu_z", sd_names))
    )
    make_random_scale_table_fit(formula_result, posterior)
  }
  shared <- "mu__xREx__id_flag__xXx__x"
  for(x_sd in c(2, 1e7, 1e8)){
    expect_error(
      transform_scale_samples(shared_leaf_samples(x_sd)),
      paste0("Random-effect SD column '", shared, "' cannot be unscaled because its shared design columns transform to different SDs."),
      fixed = TRUE, info = paste("sd(x) =", x_sd)
    )
  }
})

# Mixed posteriors of the fixed and SD coordinates of a fitted block, on the
# fitted scale, and the same without the fitted block design in their formula
# state or without any formula state.
random_block_table_samples <- function(fit, parameters, design = c("block", "none", "state")){

  design <- match.arg(design)
  samples <- as_mixed_posteriors(fit, parameters, transform_scaled = FALSE)
  state <- .bt_formula_state_get(samples)
  if(identical(design, "none")){
    state$models[[1L]]$formula_design <- NULL
    samples <- .bt_formula_state_attach(samples, state)
  }else if(identical(design, "state")){
    samples <- .bt_meta_set(samples, "formula_state", NULL)
    for(parameter in names(samples)){
      samples[[parameter]] <- .bt_meta_set(samples[[parameter]], "formula_state", NULL)
    }
  }

  samples
}

random_block_table_stub <- function(parameters){

  stats::setNames(lapply(parameters, function(parameter){
    structure(
      list(BF = 2, prior_probs = c(0.5, 0.5), post_probs = c(0.4, 0.6)),
      class = c("list", "marginal_inference")
    )
  }), parameters)
}

test_that("ensemble and marginal tables unscale random-effect SDs through the block design", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  grid <- expand.grid(f = factor(c("a", "b", "c")), x = c(-1, 0, 1, 2))
  setup <- random_block_scale_fit(
    ~ 1 + x + diag(1 + f + f:x | id), df,
    prior_random(id = random_block(
      sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
    ))
  )
  A <- random_block_unscale_reference(~ f + f:x, grid, mean(df$x), stats::sd(df$x))
  expected <- unname(sqrt(diag(A %*% diag(setup$tau^2) %*% t(A))))
  expect_gt(abs(expected[1L] - setup$tau[1L]), 5)
  parameters <- c("mu__xREx__id_intercept", "mu__xREx__id_f", "mu__xREx__id_f__xXx__x")
  formula_scale <- attr(setup$fit, "formula_scale")

  # the mixed posteriors that transform the SDs at their creation
  original <- as_mixed_posteriors(setup$fit, parameters, transform_scaled = TRUE)
  expect_equal(
    unname(c(mean(original[[1L]]), colMeans(original[[2L]]), colMeans(original[[3L]]))),
    expected, tolerance = 1e-10
  )

  # the tables transform the fitted-scale samples with the design their state carries
  samples <- random_block_table_samples(setup$fit, parameters)
  ensemble <- ensemble_estimates_table(
    samples, parameters, transform_scaled = TRUE, formula_scale = formula_scale
  )
  expect_equal(unname(ensemble[, "Mean"]), expected, tolerance = 1e-10)
  marginal <- marginal_estimates_table(
    samples, random_block_table_stub(parameters), parameters[1L],
    transform_scaled = TRUE, formula_scale = formula_scale
  )
  expect_equal(unname(marginal[1L, "Mean"]), expected[1L], tolerance = 1e-10)
})

test_that("a mixture unscales the random-effect SDs of each model with the design of its own block", {

  skip_if_not_installed("runjags")

  # the models share the formula and the SD parameters but not the design of
  # the block: `flag` is logical (two indicator columns of `flag:x` that share
  # one SD) in the first model and numeric in the second
  df <- random_block_scale_data()
  df$flag <- rep(c(TRUE, FALSE), 18)
  df_numeric <- df
  df_numeric$flag <- as.numeric(df$flag)
  m <- mean(df$x)
  s <- stats::sd(df$x)
  block <- prior_random(id = random_block(
    sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
  ))
  first <- random_block_scale_fit(~ 1 + x + diag(1 + flag + flag:x | id), df, block)
  second <- random_block_scale_fit(
    ~ 1 + x + diag(1 + flag + flag:x | id), df_numeric, block, tau = first$tau * 1.5
  )
  A_first <- random_block_unscale_reference(
    ~ flag + flag:x, expand.grid(flag = c(FALSE, TRUE), x = c(-1, 0, 1, 2)), m, s
  )
  A_second <- random_block_unscale_reference(
    ~ flag + flag:x, expand.grid(flag = c(0, 1), x = c(-1, 0, 1, 2)), m, s
  )
  expected_first <- unname(sqrt(diag(
    A_first %*% diag(first$tau[c(1L, 2L, 3L, 3L)]^2) %*% t(A_first)
  )))[1:3]
  expected_second <- unname(sqrt(diag(A_second %*% diag(second$tau^2) %*% t(A_second))))
  # the first model's design would give the second model's draws other SDs
  expect_gt(max(abs(expected_second - expected_first * 1.5)), 1)
  expect_identical(first$sd_names, second$sd_names)
  parameters <- first$sd_names
  formula_scale <- attr(first$fit, "formula_scale")

  models <- list(
    list(fit = first$fit, marglik = bridgesampling_object(0), prior_weights = 0.5),
    list(fit = second$fit, marglik = bridgesampling_object(0), prior_weights = 0.5)
  )
  # with and without fixed coefficients among the mixed parameters
  for(mixture_parameters in list(c("mu_intercept", "mu_x", parameters), parameters)){
    mixed <- mix_posteriors(
      models, mixture_parameters,
      stats::setNames(rep(list(c(FALSE, FALSE)), length(mixture_parameters)), mixture_parameters),
      seed = 1, n_samples = 40
    )
    state <- .bt_formula_state_get(mixed)
    expect_false(is.null(state))
    model <- state$model
    expect_setequal(model, 1:2)
    expect_false(is.null(state$models[[1L]]$formula_design$mu))
    expect_false(is.null(state$models[[2L]]$formula_design$mu))
    transformed <- .transform_scale_samples_list(mixed, formula_scale)
    got <- vapply(parameters, function(parameter) as.numeric(transformed[[parameter]]), numeric(length(model)))
    expect_equal(
      unname(got), unname(rbind(expected_first, expected_second)[model, ]), tolerance = 1e-10,
      info = paste(mixture_parameters, collapse = ",")
    )
  }

  # the model's own design is also what the tables use
  table <- ensemble_estimates_table(
    mixed, parameters, transform_scaled = TRUE, formula_scale = formula_scale
  )
  weights <- tabulate(model, 2L) / length(model)
  expect_equal(
    unname(table[, "Mean"]),
    as.numeric(weights %*% rbind(expected_first, expected_second)), tolerance = 1e-10
  )
})

test_that("tables without the fitted block design refuse random-effect SDs that pairing names cannot map", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  grid <- expand.grid(f = factor(c("a", "b", "c")), x = c(-1, 0, 1, 2))
  m <- mean(df$x)
  s <- stats::sd(df$x)
  block <- prior_random(id = random_block(
    sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
  ))
  parameters <- list(
    c("mu__xREx__id_intercept", "mu__xREx__id_f", "mu__xREx__id_f__xXx__x"),
    c("mu__xREx__id_intercept", "mu__xREx__id_f", "mu__xREx__id_x", "mu__xREx__id_f__xXx__x"),
    c("mu__xREx__id_intercept", "mu__xREx__id_x")
  )

  # `f` and `f:x` without `x`: the induced shifts of the intercept and of the
  # factor levels are not carried by coefficient names
  setup <- random_block_scale_fit(~ 1 + x + diag(1 + f + f:x | id), df, block)
  formula_scale <- attr(setup$fit, "formula_scale")
  for(design in c("none", "state")){
    samples <- random_block_table_samples(setup$fit, parameters[[1L]], design)
    refusal <- expect_error(
      ensemble_estimates_table(
        samples, parameters[[1L]], transform_scaled = TRUE, formula_scale = formula_scale
      ),
      class = "BayesTools_formula_transform_unavailable"
    )
    expect_identical(refusal$reason, "random_effect_design_unavailable")
    expect_identical(
      conditionMessage(refusal),
      paste0(
        "Cannot transform 'mu__xREx__id_intercept', 'mu__xREx__id_f[1]', 'mu__xREx__id_f[2]', ",
        "'mu__xREx__id_f__xXx__x[1]', 'mu__xREx__id_f__xXx__x[2]', 'mu__xREx__id_f__xXx__x[3]' ",
        "to the original predictor scale: the samples carry no fitted design of the ",
        "random-effect block 'id' of formula parameter 'mu', which is needed to map its ",
        "SDs to the original predictor scale. Create the samples with as_mixed_posteriors() ",
        "or mix_posteriors() from fits that carry their formula design; mixed posteriors of a ",
        "model without a fitted design cannot be mapped to the original predictor scale."
      )
    )
    expect_error(
      marginal_estimates_table(
        samples, random_block_table_stub(parameters[[1L]]), parameters[[1L]][1L],
        transform_scaled = TRUE, formula_scale = formula_scale
      ),
      class = "BayesTools_formula_transform_unavailable"
    )
  }
  # fixed coefficients do not depend on the design of the random blocks
  expect_no_error(ensemble_estimates_table(
    random_block_table_samples(setup$fit, "mu_x", "none"), "mu_x",
    transform_scaled = TRUE, formula_scale = formula_scale
  ))

  # blocks that pairing coefficient names maps exactly keep their values
  cases <- list(
    list(~ 1 + x + diag(1 + f * x | id), ~ f * x, parameters[[2L]]),
    list(~ 1 + x + diag(1 + x | id), ~ x, parameters[[3L]])
  )
  for(case in cases){
    setup <- random_block_scale_fit(case[[1L]], df, block)
    A <- random_block_unscale_reference(case[[2L]], grid, m, s)
    expected <- unname(sqrt(diag(A %*% diag(setup$tau^2, nrow = length(setup$tau)) %*% t(A))))
    for(design in c("none", "state")){
      samples <- random_block_table_samples(setup$fit, case[[3L]], design)
      table <- ensemble_estimates_table(
        samples, case[[3L]], transform_scaled = TRUE,
        formula_scale = attr(setup$fit, "formula_scale")
      )
      expect_equal(
        unname(table[, "Mean"]), expected, tolerance = 1e-10,
        info = paste(deparse(case[[1L]]), design)
      )
    }
  }
})

test_that("a logical predictor in a block with a scaled slope is refused without the block design and mapped with it", {

  skip_if_not_installed("runjags")

  # `flagFALSE:x` and `flagTRUE:x` share one term name: pairing names finds the
  # partner `flag` of both, but their induced terms are different
  df <- random_block_scale_data()
  df$flag <- rep(c(TRUE, FALSE), 18)
  grid <- expand.grid(flag = c(FALSE, TRUE), x = c(-1, 0, 1, 2))
  block <- prior_random(id = random_block(
    sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
  ))
  setup <- random_block_scale_fit(~ 1 + x + diag(1 + flag + flag:x | id), df, block)
  A <- random_block_unscale_reference(~ flag + flag:x, grid, mean(df$x), stats::sd(df$x))
  expected <- unname(sqrt(diag(A %*% diag(setup$tau[c(1L, 2L, 3L, 3L)]^2) %*% t(A))))[1:3]
  expect_gt(abs(expected[1L] - setup$tau[1L]), 5)
  parameters <- setup$sd_names
  expect_identical(length(parameters), 3L)
  formula_scale <- attr(setup$fit, "formula_scale")
  tables <- function(samples){
    ensemble_estimates_table(
      samples, parameters, transform_scaled = TRUE, formula_scale = formula_scale
    )
  }

  # without the design the pairing of names would report 0.5 / 10.7 / 1.02
  for(design in c("none", "state")){
    samples <- random_block_table_samples(setup$fit, parameters, design)
    refusal <- expect_error(tables(samples), class = "BayesTools_formula_transform_unavailable")
    expect_identical(refusal$reason, "random_effect_design_unavailable")
    expect_identical(
      conditionMessage(refusal),
      paste0(
        "Cannot transform 'mu__xREx__id_intercept', 'mu__xREx__id_flag', ",
        "'mu__xREx__id_flag__xXx__x' to the original predictor scale: the samples carry ",
        "no fitted design of the random-effect block 'id' of formula parameter 'mu', ",
        "which is needed to map its SDs to the original predictor scale. Create the ",
        "samples with as_mixed_posteriors() or mix_posteriors() from fits that carry their ",
        "formula design; mixed posteriors of a model without a fitted design cannot be mapped ",
        "to the original predictor scale."
      )
    )
    expect_error(
      marginal_estimates_table(
        samples, random_block_table_stub(parameters), parameters[1L],
        transform_scaled = TRUE, formula_scale = formula_scale
      ),
      class = "BayesTools_formula_transform_unavailable"
    )
  }

  # with the design: mixed posteriors, and a mixture of the SDs alone
  expect_equal(
    unname(tables(random_block_table_samples(setup$fit, parameters, "block"))[, "Mean"]),
    expected, tolerance = 1e-10
  )
  original <- as_mixed_posteriors(setup$fit, parameters, transform_scaled = TRUE)
  expect_equal(
    unname(vapply(original, function(x) mean(as.matrix(x)[, 1L]), numeric(1))),
    expected, tolerance = 1e-10
  )
  mixed <- mix_posteriors(
    list(list(fit = setup$fit, marglik = bridgesampling_object(0), prior_weights = 1)),
    parameters, stats::setNames(rep(list(FALSE), length(parameters)), parameters),
    seed = 1, n_samples = 3
  )
  expect_false(is.null(.bt_formula_state_get(mixed)))
  expect_equal(unname(tables(mixed)[, "Mean"]), expected, tolerance = 1e-10)

  # `(0 + x | g)` documents its drop of the implied intercept: the design path
  # and a mixture of its SDs alone report tau / sd(x); without the design the
  # block is refused, since pairing names is not its documented map either
  s <- stats::sd(df$x)
  setup <- random_block_scale_fit(~ 1 + x + diag(0 + x | id), df, block)
  expect_equal(
    unname(transform_scale_samples(setup$fit)[1L, setup$sd_names]),
    setup$tau / s, tolerance = 1e-12
  )
  mixed <- mix_posteriors(
    list(list(fit = setup$fit, marglik = bridgesampling_object(0), prior_weights = 1)),
    setup$sd_names, stats::setNames(rep(list(FALSE), 1L), setup$sd_names),
    seed = 1, n_samples = 3
  )
  formula_scale <- attr(setup$fit, "formula_scale")
  expect_equal(
    unname(ensemble_estimates_table(
      mixed, setup$sd_names, transform_scaled = TRUE, formula_scale = formula_scale
    )[, "Mean"]),
    setup$tau / s, tolerance = 1e-12
  )
  refusal <- expect_error(
    ensemble_estimates_table(
      random_block_table_samples(setup$fit, setup$sd_names, "none"), setup$sd_names,
      transform_scaled = TRUE, formula_scale = formula_scale
    ),
    class = "BayesTools_formula_transform_unavailable"
  )
  expect_identical(refusal$reason, "random_effect_design_unavailable")
  expect_match(conditionMessage(refusal), "which is needed to map its SDs to the original predictor scale", fixed = TRUE)
})

test_that("a mixture of random-effect SDs with a model that has no formula design attaches no state", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  setup <- random_block_scale_fit(
    ~ 1 + x + diag(1 + x | id), df,
    prior_random(id = random_block(
      sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
    ))
  )
  plain <- structure(
    list(
      mcmc = coda::mcmc.list(coda::mcmc(cbind(theta = seq(-1, 1, length.out = 20)))),
      sample = 20L, summary.pars = list(mutate = NULL), monitor = "theta"
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(plain, "prior_list") <- list(theta = prior("normal", list(0, 1)))
  plain <- attach_test_parameter_map(plain)

  mixed <- mix_posteriors(
    list(
      list(fit = setup$fit, marglik = bridgesampling_object(0), prior_weights = 0.5),
      list(fit = plain, marglik = bridgesampling_object(0), prior_weights = 0.5)
    ),
    setup$sd_names, stats::setNames(rep(list(c(FALSE, FALSE)), 2L), setup$sd_names),
    seed = 1, n_samples = 20
  )
  expect_s3_class(mixed, "mixed_posteriors")
  expect_null(.bt_formula_state_get(mixed))
  # the SDs of an exactly paired block are still unscaled by pairing names
  table <- ensemble_estimates_table(
    mixed, setup$sd_names, transform_scaled = TRUE,
    formula_scale = attr(setup$fit, "formula_scale")
  )
  expect_identical(nrow(table), 2L)
})

test_that("the no-design gate admits a block only when pairing names gives its design-derived map", {

  skip_if_not_installed("runjags")

  set.seed(7)
  df <- data.frame(
    x  = rnorm(60, 10, 2),
    z  = rnorm(60, -3, 1.5),
    f  = factor(rep(c("a", "b", "c"), 20)),
    id = factor(rep(paste0("g", 1:10), each = 6)),
    flag = rep(c(TRUE, FALSE), 30)
  )
  df$xc <- df$x - mean(df$x)
  df$xs <- rep(c(-3, -2, -1, 0, 1, 2, 3, -3.5, 3.5, 0), 6)

  # shapes: random-effect block, the fixed predictors, the scaled predictors,
  # and the contrasts of the block
  shape <- function(block, fixed, scaled = fixed, contrasts = NULL){
    list(block = block, fixed = fixed, scaled = scaled, contrasts = contrasts)
  }
  shapes <- list(
    shape("diag(1 + x | id)", "x"),
    shape("diag(1 + f * x | id)", "x"),
    shape("diag(1 + f + f:x | id)", "x"),
    shape("diag(1 + f:x | id)", "x"),
    shape("diag(0 + x | id)", "x"),
    shape("diag(0 + f + f:x | id)", "x"),
    shape("diag(0 + f + x + f:x | id)", "x"),
    shape("diag(1 + flag * x | id)", "x"),
    shape("diag(1 + flag + flag:x | id)", "x"),
    shape("diag(0 + flag + flag:x | id)", "x"),
    shape("diag(0 + flag:x | id)", "x"),
    shape("diag(1 + x * z | id)", c("x", "z")),
    shape("diag(1 + x:z | id)", c("x", "z")),
    shape("diag(1 + x + x:z | id)", c("x", "z")),
    shape("diag(1 + x + z + x:z | id)", c("x", "z"), "x"),
    shape("diag(1 + x:z | id)", c("x", "z"), "x"),
    shape("diag(1 + f * x * z | id)", c("x", "z")),
    shape("diag(1 + f + x + z + f:z + x:z + f:x:z | id)", c("x", "z")),
    shape("diag(1 + flag + z + flag:z + flag:x + x:z + flag:x:z | id)", c("x", "z")),
    shape("(1 + f + f:x | id)", "x"),
    shape("(1 + f * x | id)", "x"),
    shape("diag(1 + xc | id)", "xc"),
    shape("diag(1 + f + f:xc | id)", "xc"),
    shape("diag(1 + xs | id)", "xs"),
    shape("diag(1 + flag * x | id)", c("x", "flag")),
    shape("diag(1 + flag + flag:x | id)", c("x", "flag")),
    shape("diag(1 + f * x | id)", "x", contrasts = c(f = "orthonormal")),
    shape("diag(1 + f + f:x | id)", "x", contrasts = c(f = "orthonormal")),
    shape("diag(1 + f * x | id)", "x", contrasts = c(f = "meandif")),
    shape("diag(1 + f + f:x | id)", "x", contrasts = c(f = "meandif")),
    shape("diag(1 + f * x | id)", "x", contrasts = c(f = "independent")),
    shape("diag(0 + f + f:x | id)", "x", contrasts = c(f = "independent"))
  )

  problems <- character()
  admitted <- character()
  refused <- character()
  checked <- 0L
  for(case in shapes){
    label <- paste0(
      case$block, " [", paste(case$fixed, collapse = ","), "; scaled ",
      paste(case$scaled, collapse = ","), "]",
      if(!is.null(case$contrasts)) paste0(" ", unname(case$contrasts))
    )
    formula <- stats::as.formula(paste("~ 1 +", paste(case$fixed, collapse = " + "), "+", case$block))
    prior_list <- c(
      list(intercept = prior("normal", list(0, 1))),
      stats::setNames(lapply(case$fixed, function(variable) prior("normal", list(0, 1))), case$fixed)
    )
    formula_result <- JAGS_formula(
      formula, "mu", df, prior_list = prior_list,
      formula_scale = stats::setNames(rep(list(TRUE), length(case$scaled)), case$scaled),
      prior_random = prior_random(id = random_block(
        sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor(),
        contrasts = case$contrasts
      ))
    )
    design <- formula_result$formula_design
    parameter_scale <- formula_result$formula_scale
    sd_names <- unique(unlist(lapply(design$random_effects, function(term) term$sd_parameter_names)))
    groups <- .random_sd_column_unscale_groups(sd_names, parameter_scale, "mu")
    for(group in groups){
      checked <- checked + 1L
      gate <- .random_sd_name_map_is_exact(group$column_terms, parameter_scale, "mu")
      by_name <- .random_sd_unscale_matrix(group$column_terms, parameter_scale, "mu")
      by_design <- .random_sd_unscale_matrix(
        group$column_terms, parameter_scale, "mu",
        random_term = .random_sd_design_block(design, group$group_key, group$leaf_names),
        source_data = design$source_data
      )
      if(gate){
        admitted <- c(admitted, label)
        if(!isTRUE(all.equal(unname(by_name), unname(by_design), tolerance = 1e-10))){
          problems <- c(problems, label)
        }
      }else{
        refused <- c(refused, label)
      }
    }
  }

  # the gate never admits a block whose pairing of names differs from its design map
  expect_identical(problems, character())
  expect_gt(checked, 30L)
  expect_gt(length(admitted), 10L)
  expect_true("diag(1 + x | id) [x; scaled x]" %in% admitted)
  expect_true("diag(1 + f * x | id) [x; scaled x]" %in% admitted)
  # blocks whose induced terms the pairing of names misses or cannot tell apart
  expect_true("diag(1 + f + f:x | id) [x; scaled x]" %in% refused)
  expect_true("diag(1 + flag + flag:x | id) [x; scaled x]" %in% refused)
  expect_true("diag(1 + flag + z + flag:z + flag:x + x:z + flag:x:z | id) [x,z; scaled x,z]" %in% refused)
})

test_that("shared random-effect SD leaves keep one coordinate per design column in the name-paired map", {

  skip_if_not_installed("runjags")

  # the columns of one shared leaf have one term name; each is its own coordinate
  # and puts no weight on the others
  scale_list <- list(mu_x = list(mean = 10, sd = 2))
  M <- .build_unscale_matrix_by_names(
    c("mu_intercept", "mu_flag", "mu_flag__xXx__x", "mu_flag__xXx__x"),
    scale_list, "mu"
  )
  # Row mu_flag: 1 on its own column and, on each of the two interaction columns
  # of the shared leaf, -mean / sd = -10 / 2 = -5. Both columns carry the
  # constant, which is the only part the exact map drops, so the weight is the
  # same on both and neither column puts weight on the other.
  expect_equal(
    unname(M),
    rbind(
      c(1, 0, 0, 0),
      c(0, 1, -5, -5),
      c(0, 0, 0.5, 0),
      c(0, 0, 0, 0.5)
    ),
    tolerance = 1e-14
  )
  # without a duplicate the map is the same
  expect_equal(
    unname(.build_unscale_matrix_by_names(
      c("mu_intercept", "mu_flag", "mu_flag__xXx__x"), scale_list, "mu"
    )),
    rbind(c(1, 0, 0), c(0, 1, -5), c(0, 0, 0.5)),
    tolerance = 1e-14
  )

  # blocks with a shared leaf that keep the documented drop: the SD of the
  # leaf is the common per-column value tau / sd(x)
  df <- random_block_scale_data()
  df$flag <- rep(c(TRUE, FALSE), 18)
  s <- stats::sd(df$x)
  block <- prior_random(id = random_block(
    sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
  ))
  setup <- random_block_scale_fit(~ 1 + x + diag(0 + flag:x | id), df, block, tau = 0.8)
  expect_identical(length(setup$sd_names), 1L)
  expect_equal(
    unname(transform_scale_samples(setup$fit)[1L, setup$sd_names]), 0.8 / s, tolerance = 1e-12
  )

  # `0 + flag + flag:x` codes `flag` by `flagTRUE` only (a block without an
  # intercept keeps the contrast coding of its factors), so the induced term of
  # the `flagFALSE:x` column on `1 - flag` has no column: the block has no
  # original-scale representation and keeps the documented drop as well
  setup <- random_block_scale_fit(~ 1 + x + diag(0 + flag + flag:x | id), df, block, tau = c(0.5, 0.8))
  term <- setup$formula_result$formula_design$random_effects[[1L]]
  expect_identical(term$raw_column_names, c("flagTRUE", "flagFALSE:x", "flagTRUE:x"))
  refusal <- expect_error(
    .bt_formula_unscale_design_transform(
      .bt_random_effect_unscale_design_spec(term, setup$formula_result$formula_design$source_data),
      setup$formula_result$formula_scale, "mu"
    ),
    class = "BayesTools_formula_transform_unavailable"
  )
  expect_identical(refusal$reason, "original_scale_not_representable")
  expect_identical(refusal$terms, "flag:x")
  transformed <- transform_scale_samples(setup$fit)[1L, setup$sd_names]
  expect_equal(unname(transformed[2L]), 0.8 / s, tolerance = 1e-12)

  # a block without an intercept is mapped by its design when its scaled design
  # lies in the column space of its original one (full indicator coding of `f`)
  grid <- expand.grid(f = factor(c("a", "b", "c")), x = c(-1, 0, 1, 2))
  setup <- random_block_scale_fit(
    ~ 1 + x + diag(0 + f + f:x | id), df,
    prior_random(id = random_block(
      sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor(),
      contrasts = c(f = "independent")
    ))
  )
  A <- random_block_unscale_reference(~ 0 + f + f:x, grid, mean(df$x), s, list(f = diag(3)))
  expect_identical(nrow(A), 6L)
  expect_equal(
    unname(transform_scale_samples(setup$fit)[1L, setup$sd_names]),
    unname(sqrt(diag(A %*% diag(setup$tau^2) %*% t(A)))), tolerance = 1e-10
  )

  # unchanged: `(0 + x | g)` and the block of a logical predictor with a representable design
  setup <- random_block_scale_fit(~ 1 + x + diag(0 + x | id), df, block, tau = 0.5)
  expect_equal(unname(transform_scale_samples(setup$fit)[1L, setup$sd_names]), 0.5 / s, tolerance = 1e-12)
  grid <- expand.grid(flag = c(FALSE, TRUE), x = c(-1, 0, 1, 2))
  setup <- random_block_scale_fit(~ 1 + x + diag(1 + flag + flag:x | id), df, block)
  A <- random_block_unscale_reference(~ flag + flag:x, grid, mean(df$x), s)
  expect_equal(
    unname(transform_scale_samples(setup$fit)[1L, setup$sd_names]),
    unname(sqrt(diag(A %*% diag(setup$tau[c(1L, 2L, 3L, 3L)]^2) %*% t(A))))[1:3],
    tolerance = 1e-10
  )
})

# A synthetic fit whose posterior draws of the random-effect SDs are draws of their
# gamma(2, 2) priors (fitted scale), so that on the original scale they follow the
# prior law of the original-scale SDs: a hypothesis on such an SD has a Bayes factor
# near 1 exactly when the prior law is the law of the original-scale SD.
random_block_prior_draw_fit <- function(formula, data, variable = "x", n = 4000L, seed = 21,
                                        sd_prior = prior("gamma", list(2, 2)),
                                        rsd = function(n) stats::rgamma(n, 2, 2)){

  formula_result <- JAGS_formula(
    formula, "mu", data,
    prior_list = stats::setNames(
      list(prior("normal", list(0, 1)), prior("normal", list(0, 1))),
      c("intercept", variable)
    ),
    formula_scale = stats::setNames(list(TRUE), variable),
    prior_random = prior_random(id = random_block(
      sd = sd_prior, monitor = random_block_scale_monitor()
    ))
  )
  sd_names <- unlist(lapply(formula_result$formula_design$random_effects, function(term){
    unique(term$sd_parameter_names)
  }))
  set.seed(seed)
  posterior <- cbind(
    matrix(
      stats::rnorm(n * 2L, 0, 0.3), n, 2L,
      dimnames = list(NULL, c("mu_intercept", paste0("mu_", variable)))
    ),
    matrix(
      rsd(n * length(sd_names)), n, length(sd_names),
      dimnames = list(NULL, sd_names)
    )
  )

  list(
    formula_result = formula_result,
    fit = make_random_scale_table_fit(formula_result, posterior),
    sd_names = sd_names,
    n = n
  )
}

random_block_density_ordinate <- function(density, values){

  vapply(values, function(value){
    exp(prior_density_ordinate(density, value)$log_density)
  }, numeric(1))
}

# Monte Carlo bound of the log Bayes factor of an event with probability p
# under draws that follow the prior: four standard errors of the log odds.
random_block_bf_bound <- function(p, n){

  4 * sqrt(1 / (n * p * (1 - p)))
}

test_that("the prior law of a pure-rescaling random-effect SD is the fitted law rescaled and the composite SD is refused", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  s <- stats::sd(df$x)
  setup <- random_block_prior_draw_fit(~ 1 + x + diag(1 + x | id), df)
  sd_intercept <- "mu__xREx__id_intercept"
  sd_slope <- "mu__xREx__id_x"
  mixed <- as_mixed_posteriors(setup$fit, setup$sd_names, transform_scaled = TRUE)

  # the stored density of the slope SD is the fitted gamma(2, 2) rescaled: the
  # original-scale SD is tau / sd(x), a pure rescaling (exact ordinates)
  grid <- seq(0.02, 3, length.out = 60)
  densities <- posterior_metadata(mixed, "prior_densities")
  expect_equal(
    random_block_density_ordinate(densities[[sd_slope]], grid),
    stats::dgamma(grid * s, 2, 2) * s, tolerance = 1e-8
  )
  context <- posterior_metadata(mixed, "prior_context")
  expect_false(is.null(context$formula_design$mu))
  marginal <- marginal_posterior(mixed, sd_slope, use_formula = FALSE, prior_samples = TRUE)
  marginal_density <- posterior_metadata(marginal, "prior_density")
  expect_equal(
    random_block_density_ordinate(marginal_density, grid),
    stats::dgamma(grid * s, 2, 2) * s, tolerance = 1e-8
  )
  # prior mass of `sd(x) > 0.5`, and the Bayes factor of posterior = prior draws
  mass <- stats::integrate(
    function(values) random_block_density_ordinate(marginal_density, values), 0.5, Inf
  )$value
  p <- stats::pgamma(0.5 * s, 2, 2, lower.tail = FALSE)
  expect_equal(mass, p, tolerance = 1e-6)
  hypothesis <- as.data.frame(hypothesis_BF(
    marginal, hypothesis = paste0(sd_slope, " > 0.5"), parameter = sd_slope
  ))
  expect_lt(abs(log(hypothesis$BF10)), random_block_bf_bound(p, setup$n))

  # the intercept SD combines the fitted SDs of the intercept and of the slope
  # (mean(x) differs from 0): its law is not the one of a fitted SD
  expect_gt(abs(mean(df$x)), 1)
  expect_null(densities[[sd_intercept]])
  entry <- posterior_metadata(densities, "measure_unavailable")
  expect_identical(entry$column, sd_intercept)
  expect_identical(entry$cause, "random_sd_prior_law_unavailable")
  element <- mixed[[sd_intercept]]
  refusal <- expect_error(
    posterior_measure_check(element, "prior_density"),
    class = "BayesTools_formula_measure_unavailable"
  )
  expect_identical(refusal$reason, "random_sd_prior_law_unavailable")
  refusals <- list(
    marginal = function() marginal_posterior(mixed, sd_intercept, use_formula = FALSE, prior_samples = TRUE),
    plot = function() plot_posterior(mixed, sd_intercept, plot_type = "ggplot", prior = TRUE)
  )
  for(label in names(refusals)){
    refusal <- expect_error(refusals[[label]](), class = "BayesTools_formula_measure_unavailable", info = label)
    expect_identical(refusal$reason, "random_sd_prior_law_unavailable", info = label)
  }
  # what does not read the prior law of the SD keeps working
  expect_no_error(marginal_posterior(mixed, sd_intercept, use_formula = FALSE, prior_samples = FALSE))
  expect_no_error(plot_posterior(mixed, sd_intercept, plot_type = "ggplot", prior = FALSE))
  table <- ensemble_estimates_table(mixed, c(sd_intercept, sd_slope))
  expect_identical(nrow(table), 2L)
  expect_true(all(is.finite(table[, "Mean"])))

  # the fitted-scale control is unchanged: the fitted gamma(2, 2) and a Bayes factor near 1
  control <- as_mixed_posteriors(setup$fit, setup$sd_names, transform_scaled = FALSE)
  control_marginal <- marginal_posterior(control, sd_slope, use_formula = FALSE, prior_samples = TRUE)
  expect_equal(
    random_block_density_ordinate(posterior_metadata(control_marginal, "prior_density"), grid),
    stats::dgamma(grid, 2, 2), tolerance = 1e-8
  )
  p_control <- stats::pgamma(1, 2, 2, lower.tail = FALSE)
  control_bf <- as.data.frame(hypothesis_BF(
    control_marginal, hypothesis = paste0(sd_slope, " > 1"), parameter = sd_slope
  ))$BF10
  expect_lt(abs(log(control_bf)), random_block_bf_bound(p_control, setup$n))
})

test_that("exact-zero means retain pure SD rescaling and nonzero subtraction-centered means refuse it", {

  skip_if_not_installed("runjags")

  # mean(x) is exactly zero: no term is induced in the intercept SD (c = 1)
  df <- data.frame(
    x  = rep(c(-3, -2, -1, 0, 1, 2, 3, -3.5, 3.5, 0), 4),
    id = factor(rep(paste0("g", 1:8), each = 5))
  )
  expect_true(mean(df$x) == 0)
  s <- stats::sd(df$x)
  setup <- random_block_prior_draw_fit(~ 1 + x + diag(1 + x | id), df)
  mixed <- as_mixed_posteriors(setup$fit, setup$sd_names, transform_scaled = TRUE)
  densities <- posterior_metadata(mixed, "prior_densities")
  grid <- seq(0.02, 3, length.out = 40)
  expect_equal(
    random_block_density_ordinate(densities[["mu__xREx__id_intercept"]], grid),
    stats::dgamma(grid, 2, 2), tolerance = 1e-8
  )
  expect_equal(
    random_block_density_ordinate(densities[["mu__xREx__id_x"]], grid),
    stats::dgamma(grid * s, 2, 2) * s, tolerance = 1e-8
  )
  marginal <- marginal_posterior(mixed, "mu__xREx__id_intercept", use_formula = FALSE, prior_samples = TRUE)
  p <- stats::pgamma(1, 2, 2, lower.tail = FALSE)
  hypothesis <- as.data.frame(hypothesis_BF(
    marginal, hypothesis = "mu__xREx__id_intercept > 1", parameter = "mu__xREx__id_intercept"
  ))
  expect_lt(abs(log(hypothesis$BF10)), random_block_bf_bound(p, setup$n))

  # Centering by subtraction leaves a nonzero mean and a retained design entry.
  # Its size does not establish a pure rescaling.
  df_centered <- random_block_scale_data()
  df_centered$xc <- df_centered$x - mean(df_centered$x)
  expect_true(mean(df_centered$xc) != 0 && abs(mean(df_centered$xc)) < 1e-12)
  setup <- random_block_prior_draw_fit(~ 1 + xc + diag(1 + xc | id), df_centered, variable = "xc")
  mixed <- as_mixed_posteriors(setup$fit, setup$sd_names, transform_scaled = TRUE)
  term <- setup$formula_result$formula_design$random_effects[[1L]]
  A <- BayesTools:::.random_sd_unscale_matrix(
    term$sd_parameter_names, setup$formula_result$formula_scale, "mu",
    random_term = term, source_data = df_centered
  )
  expect_identical(A[1L, 2L], -mean(df_centered$xc) / stats::sd(df_centered$xc))
  expect_true(A[1L, 2L] != 0)
  densities <- posterior_metadata(mixed, "prior_densities")
  expect_null(densities[["mu__xREx__id_intercept"]])
  refusal <- expect_error(
    marginal_posterior(mixed, "mu__xREx__id_intercept", use_formula = FALSE, prior_samples = TRUE),
    class = "BayesTools_formula_measure_unavailable"
  )
  expect_identical(refusal$reason, "random_sd_prior_law_unavailable")
})

test_that("a tiny nonzero predictor mean retains composite SD draws and prior-law refusal", {

  skip_if_not_installed("runjags")

  sd_names <- c("mu__xREx__id_intercept", "mu__xREx__id_x")
  base <- rep(c(-3, -2, -1, 1, 2, 3), 6L)
  tau0 <- c(0.25, 0.5, 1, 2)
  tau1 <- c(1e13, 4e13, 1e14, 2e14)
  for(shift in c(1e-13, 0)){
    data <- data.frame(x = base + shift, id = factor(rep(1:6, each = 6L)))
    setup <- random_block_scale_fit(
      ~ 1 + x + diag(1 + x | id), data,
      prior_random(id = random_block(
        sd = prior("gamma", list(2, 2)),
        terms = list(x = prior("gamma", list(2, 2e-14))),
        monitor = random_block_scale_monitor()
      ))
    )
    posterior <- cbind(mu_intercept = 0, mu_x = 0, tau0, tau1)
    colnames(posterior)[3:4] <- sd_names
    fit <- make_random_scale_table_fit(setup$formula_result, posterior)
    mixed <- as_mixed_posteriors(fit, sd_names, transform_scaled = TRUE)
    ratio <- mean(data$x) / stats::sd(data$x)
    expected <- sqrt(tau0^2 + ratio^2 * tau1^2)
    expect_equal(as.numeric(mixed[[sd_names[1L]]]), expected, tolerance = 1e-13)
    expect_equal(unname(transform_scale_samples(fit)[, sd_names[1L]]), expected,
                 tolerance = 1e-13)
    densities <- posterior_metadata(mixed, "prior_densities")
    if(shift != 0){
      expect_true(ratio != 0)
      expect_gt(max(expected - tau0), 5)
      expect_null(densities[[sd_names[1L]]])
      # Prior construction refuses here, before any region consumer is called.
      refusal <- expect_error(
        marginal_posterior(mixed, sd_names[1L], use_formula = FALSE, prior_samples = TRUE),
        class = "BayesTools_formula_measure_unavailable"
      )
      expect_identical(refusal$reason, "random_sd_prior_law_unavailable")
      quantities <- parameter_catalog(fit)$quantities
      key <- quantities$extraction_key[[match("(mu) sd(intercept)", quantities$canonical_name)]]
      expect_identical(key$source_type, "composite")
      expect_identical(BayesTools:::.bt_parameter_catalog_own_scale_dependencies(fit, key), sd_names)
      expect_true(all(sd_names %in% key$dependencies))
    }else{
      expect_identical(ratio, 0)
      grid <- c(0.25, 0.5, 1, 2)
      expect_equal(random_block_density_ordinate(densities[[sd_names[1L]]], grid),
                   stats::dgamma(grid, 2, 2), tolerance = 1e-8)
      s <- stats::sd(data$x)
      slope_grid <- stats::qgamma(c(0.1, 0.5, 0.9), 2, 2e-14) / s
      expect_equal(
        random_block_density_ordinate(densities[[sd_names[2L]]], slope_grid) /
          (stats::dgamma(slope_grid * s, 2, 2e-14) * s),
        rep(1, length(slope_grid)), tolerance = 1e-8
      )
    }
  }
})

# Data of a block with a nearly centered predictor `u` (mean / SD = `ratio`, SD
# 1e4) and a centered predictor `v` with a small SD (1e-3), so that the row of `v`
# in the original-scale map has the large entry 1 / sd(v) on its own coordinate
# and the entry -mean(u) / (sd(u) sd(v)) on the interaction `u:v`, `ratio / sd(v)`
# in size and 1e-12 times the diagonal at `ratio = 1e-12`. `q`, a non-centered
# scaled slope, makes a categorical recoding of `f:q` necessary.
random_block_nearly_centered_data <- function(ratio){

  set.seed(11)
  n <- 60L
  df <- data.frame(
    id = factor(rep(1:10, each = 6L)),
    f  = factor(rep(c("a", "b", "c"), length.out = n)),
    q  = rnorm(n, 10, 2)
  )
  u <- rnorm(n, 0, 1e4)
  u <- u - mean(u)
  v <- rnorm(n, 0, 1e-3)
  df$u <- u + ratio * stats::sd(u)
  df$v <- v - mean(v)

  df
}

# A synthetic fit of the block `rhs` over `u`, `v` (and `q`), whose SD draws are
# draws of their priors: gamma(2, 2) for every SD but the interaction `u:v`, whose
# prior scale is 1e13 times larger, so that a tiny entry on `u:v` matters.
random_block_nearly_centered_fit <- function(rhs, scaled, data, n = 4000L, seed = 21){

  formula_result <- JAGS_formula(
    stats::as.formula(paste("~ 1 +", paste(scaled, collapse = " + "), "+ diag(", rhs, "| id)")),
    "mu", data,
    prior_list = stats::setNames(
      rep(list(prior("normal", list(0, 1))), length(scaled) + 1L), c("intercept", scaled)
    ),
    formula_scale = stats::setNames(rep(list(TRUE), length(scaled)), scaled),
    prior_random = prior_random(id = random_block(
      sd = prior("gamma", list(2, 2)),
      terms = list("u:v" = prior("gamma", list(2, 2e-13))),
      monitor = random_block_scale_monitor()
    ))
  )
  term <- formula_result$formula_design$random_effects[[1L]]
  sd_names <- unique(term$sd_parameter_names)
  leaf <- function(name) sub("[[]([0-9]+)[]]$", "", name)
  prior_draws <- function(n){
    vapply(sd_names, function(name) rng(formula_result$prior_list[[leaf(name)]], n), numeric(n))
  }
  set.seed(seed)
  posterior <- cbind(
    matrix(
      stats::rnorm(n * (length(scaled) + 1L), 0, 0.3), n, length(scaled) + 1L,
      dimnames = list(NULL, paste0("mu_", c("intercept", scaled)))
    ),
    prior_draws(n)
  )

  list(
    formula_result = formula_result,
    term = term,
    data = data,
    fit = make_random_scale_table_fit(formula_result, posterior),
    sd_names = sd_names,
    leaves = unique(leaf(sd_names)),
    posterior = posterior[, sd_names, drop = FALSE],
    prior_draws = prior_draws
  )
}

# Elementary affine expansion for the C5 block, independently of the unscaler.
random_block_c5_reference <- function(setup){

  spec <- BayesTools:::.bt_random_effect_unscale_design_spec(setup$term, setup$data)
  cn <- BayesTools:::.bt_formula_unscale_coefficient_names(spec, "mu")
  A <- diag(length(cn))
  dimnames(A) <- list(cn, cn)
  i0 <- match("mu_intercept", cn)
  iu <- match("mu_u", cn)
  iv <- match("mu_v", cn)
  iuv <- match("mu_u__xXx__v", cn)
  mu <- mean(setup$data$u); su <- stats::sd(setup$data$u)
  mv <- mean(setup$data$v); sv <- stats::sd(setup$data$v)
  A[i0, iu] <- -mu / su; A[iu, iu] <- 1 / su
  A[i0, iv] <- -mv / sv; A[iv, iv] <- 1 / sv
  A[c(i0, iu, iv, iuv), iuv] <- c(mu * mv, -mv, -mu, 1) / (su * sv)
  iq <- grep("^mu_f__xXx__q", cn)
  if(length(iq)){
    ib <- c(i0, grep("^mu_f\\[", cn))
    # f:q is indicator coded; intercept/f use the persisted factor contrast.
    B <- cbind(1, spec$contrast_matrices$f)
    A[ib, iq] <- -mean(setup$data$q) / stats::sd(setup$data$q) * solve(B)
    A[iq, iq] <- diag(1 / stats::sd(setup$data$q), length(iq))
  }
  A
}

test_that("categorical recoding preserves tiny numeric affine entries across predictor units", {

  skip_if_not_installed("runjags")
  for(units in c(1e-3, 1, 1e3)){
    data <- random_block_nearly_centered_data(4e-13)
    data$v <- data$v * units
    setup <- random_block_nearly_centered_fit(
      "1 + f + f:q + u + v + u:v", c("q", "u", "v"), data, n = 12L
    )
    spec <- BayesTools:::.bt_random_effect_unscale_design_spec(setup$term, data)
    A <- BayesTools:::.bt_formula_unscale_design_transform(spec, setup$formula_result$formula_scale, "mu")$matrix
    ref <- random_block_c5_reference(setup)
    expect_identical(A[ref == 0], ref[ref == 0])
    expect_true(all(A[ref != 0] != 0))
    expect_lt(max(abs(A[ref != 0] / ref[ref != 0] - 1)), 2e-14)
    # Full covariance, including an independent correlated control.
    tau <- as.numeric(setup$posterior[1L, ])
    for(rho in c(0, 0.2)){
      R <- matrix(rho, length(tau), length(tau)); diag(R) <- 1
      Sigma <- R * outer(tau, tau)
      got <- A %*% Sigma %*% t(A)
      expected <- ref %*% Sigma %*% t(ref)
      expect_lt(max(abs(got - expected) / sqrt(outer(diag(expected), diag(expected)))), 3e-14)
    }
    expected <- sqrt(setup$posterior^2 %*% t(ref^2))
    got <- transform_scale_samples(setup$fit)[, setup$sd_names, drop = FALSE]
    expect_lt(max(abs(unname(got) / unname(expected) - 1)), 3e-14)
    mixed <- as_mixed_posteriors(setup$fit, setup$leaves, transform_scaled = TRUE)
    mixed_draws <- do.call(cbind, lapply(mixed, function(x) matrix(as.numeric(x), nrow = nrow(expected))))
    expect_identical(dim(mixed_draws), dim(expected))
    expect_lt(max(abs(unname(mixed_draws) / unname(expected) - 1)), 3e-14)
    catalog <- parameter_catalog(setup$fit)
    selection <- parameter_catalog_resolve(catalog, "(mu) sd(v)")
    expect_equal(as.numeric(as.matrix(parameter_draws(setup$fit, selection))), expected[, match("mu__xREx__id_v", setup$sd_names)],
      tolerance = 2e-14)
  }
})

test_that("categorical hierarchy and indicator partitions keep unrelated SD laws pure", {

  skip_if_not_installed("runjags")
  data <- expand.grid(f = letters[1:4], g = LETTERS[1:3], d1 = c(FALSE, TRUE),
    d2 = c(FALSE, TRUE), x = c(2, 5), KEEP.OUT.ATTRS = FALSE)
  data$f <- factor(data$f); data$g <- factor(data$g)
  data$id <- factor(rep(1:4, length.out = nrow(data)))
  for(coding in c("meandif", "orthonormal")) for(rhs in c(
    "f*g+f:x", "1+g+f:g+f:x", "1+g+f:g+f:x+d1+d2"
  )){
    compiled <- JAGS_formula(stats::as.formula(paste("~1+x+diag(", rhs, "|id)")), "mu", data,
      list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
      formula_scale = list(x = TRUE), prior_random = prior_random(id = random_block(
        sd = prior("gamma", list(2, 2)), contrasts = c(f = coding, g = coding),
        monitor = random_block_scale_monitor())))
    term <- compiled$formula_design$random_effects[[1L]]
    spec <- BayesTools:::.bt_random_effect_unscale_design_spec(term, data)
    A <- BayesTools:::.bt_formula_unscale_design_transform(spec, compiled$formula_scale, "mu")$matrix
    untouched <- grep("^mu_g\\[|^mu_d1$|^mu_d2$", rownames(A))
    slopes <- grep("^mu_f__xXx__x", colnames(A))
    # These effects cannot be induced by a function of f alone.
    expect_true(all(A[untouched, slopes, drop = FALSE] == 0))
    expect_identical(unname(rowSums(A[untouched, , drop = FALSE] != 0)), rep(1, length(untouched)))
    laws <- BayesTools:::.random_sd_original_scale_laws(compiled$formula_scale, "mu",
      BayesTools:::.bt_formula_design_random_unscale_view(compiled$formula_design))
    leaves <- unique(term$sd_parameter_names[untouched])
    expect_true(all(leaves %in% names(laws$scale)))
    expect_equal(unname(laws$scale[leaves]), rep(1, length(leaves)), tolerance = 0)
    # Coefficients within the f basis may be tiny: no claim of idealized
    # eigenvector orthogonality is made about those persisted doubles.
  }
})

test_that("source logical functions retain the elementary Boolean map on a nonproduct grid", {

  skip_if_not_installed("runjags")
  data <- expand.grid(f = letters[1:3], d1 = c(FALSE, TRUE), d2 = c(FALSE, TRUE),
    d3 = c(FALSE, TRUE), x = c(2, 5), KEEP.OUT.ATTRS = FALSE)
  data$f <- factor(data$f); data$id <- factor(rep(1:4, length.out = nrow(data)))
  compiled <- JAGS_formula(~1+x+diag(1+d2+d1:d2+d1:x+d3+f|id), "mu", data,
    list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE), prior_random = prior_random(id = random_block(
      sd = prior("gamma", list(2, 2)), contrasts = c(f = "meandif"), monitor = random_block_scale_monitor())))
  term <- compiled$formula_design$random_effects[[1L]]
  spec <- BayesTools:::.bt_random_effect_unscale_design_spec(term, data)
  A <- BayesTools:::.bt_formula_unscale_design_transform(spec, compiled$formula_scale, "mu")$matrix
  factors <- attr(stats::terms(spec$formula), "factors")
  find_term <- function(vars) which(vapply(seq_len(ncol(factors)), function(i)
    setequal(rownames(factors)[factors[, i] > 0], vars), logical(1)))
  source <- which(spec$assign == find_term(c("d1", "x")))
  target <- which(spec$assign == find_term(c("d1", "d2")))
  expect_length(source, 2L); expect_length(target, 2L)
  ratio <- mean(data$x) / stats::sd(data$x)
  expected <- diag(nrow(A))
  expected[source, source] <- diag(1 / stats::sd(data$x), 2L)
  # FALSE/TRUE sources are 1-d1 and d1; d1 is the sum of its two d2 cells.
  expected[1L, source] <- c(-ratio, 0)
  expected[target, source] <- rbind(c(ratio, -ratio), c(ratio, -ratio))
  expect_identical(A[expected == 0], expected[expected == 0])
  expect_equal(unname(A), expected, tolerance = 2e-14)
})

test_that("random rank refusal precedes dropping and matrix support keeps its boundary", {

  skip_if_not_installed("runjags")
  data <- expand.grid(f = letters[1:3], f2 = c("u", "v"), x = c(2, 5), KEEP.OUT.ATTRS = FALSE)
  data$f <- factor(data$f); data$f2 <- factor(data$f2)
  data$id <- factor(rep(1:4, length.out = nrow(data)))
  priors <- list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1)))
  block <- prior_random(id = random_block(sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()))
  compiled <- JAGS_formula(~1+x+diag(0+x+x:f:f2|id), "mu", data, priors,
    formula_scale = list(x = TRUE), prior_random = block)
  term <- compiled$formula_design$random_effects[[1L]]
  refusal <- expect_error(BayesTools:::.random_sd_unscale_matrix(term$sd_parameter_names,
    compiled$formula_scale, "mu", term, data), class = "BayesTools_formula_transform_unavailable")
  expect_identical(refusal$reason, "original_scale_not_identified")
  data$M <- I(cbind(a = data$x, b = seq_len(nrow(data))))
  data$M1 <- I(matrix(seq_len(nrow(data)), ncol = 1L))
  for(variable in c("M", "M1")){
    formula <- stats::as.formula(paste0("~1+x+diag(1+x+", variable, "+x:", variable, "|id)"))
    compiled <- JAGS_formula(formula, "mu", data, priors,
      formula_scale = list(x = TRUE), prior_random = block)
    term <- compiled$formula_design$random_effects[[1L]]
    request <- function() BayesTools:::.random_sd_unscale_matrix(term$sd_parameter_names,
      compiled$formula_scale, "mu", term, data)
    if(variable == "M"){
      expect_error(request(), class = "BayesTools_refit_required")
    }else{
      m <- mean(data$x); s <- stats::sd(data$x)
      expected <- diag(c(1, 1/s, 1, 1/s))
      expected[1L, 2L] <- expected[3L, 4L] <- -m/s
      expect_equal(unname(request()), expected, tolerance = 2e-15)
    }
  }
})

test_that("an entry of the original-scale map that the draws keep makes the SD a composite, however small beside its row", {

  skip_if_not_installed("runjags")

  # The SD of `v` reads the SD of `u:v` through the entry -mean(u) / (sd(u) sd(v))
  # of its row, which is 4e-13 times its large diagonal entry 1 / sd(v) and still
  # retained by both design-transform paths. With SD-prior scales 1e13
  # apart the entry decides the law of the SD, so it is a composite whose prior
  # law is refused; a relative threshold on the row made it a rescaling of the
  # fitted SD, with the prior event probability 0.46 where it is above 0.9.
  designs <- list(
    categorical_basis = list(rhs = "1 + f + f:q + u + v + u:v", scaled = c("q", "u", "v")),
    name_paired   = list(rhs = "1 + u + v + u:v",           scaled = c("u", "v"))
  )
  sd_v <- "mu__xREx__id_v"
  sd_uv <- "mu__xREx__id_u__xXx__v"
  for(method in names(designs)){
    design <- designs[[method]]
    setup <- random_block_nearly_centered_fit(
      design$rhs, design$scaled, random_block_nearly_centered_data(4e-13)
    )
    formula_scale <- setup$formula_result$formula_scale
    spec <- BayesTools:::.bt_random_effect_unscale_design_spec(setup$term, setup$data)
    expect_identical(
      BayesTools:::.bt_formula_unscale_design_transform(spec, formula_scale, "mu")$method,
      method, info = method
    )

    A <- BayesTools:::.random_sd_unscale_matrix(
      column_terms = setup$term$sd_parameter_names, formula_scale = formula_scale,
      prefix = "mu", random_term = setup$term, source_data = setup$data
    )
    i_v <- match(sd_v, setup$term$sd_parameter_names)
    i_uv <- match(sd_uv, setup$term$sd_parameter_names)
    expect_gt(abs(A[i_v, i_v]), 900)
    expect_true(A[i_v, i_uv] != 0, info = method)
    expect_lt(abs(A[i_v, i_uv]) / abs(A[i_v, i_v]), 1e-12)
    expect_identical(which(A[i_v, ] != 0), c(i_v, i_uv), info = method)

    # the classification follows the matrix of the draws: a composite, refused
    laws <- BayesTools:::.random_sd_original_scale_laws(
      formula_scale, "mu",
      BayesTools:::.bt_formula_design_random_unscale_view(setup$formula_result$formula_design)
    )
    expect_identical(unname(laws$refused[sd_v]), "random_sd_prior_law_unavailable", info = method)
    expect_false(sd_v %in% names(laws$scale), info = method)
    expect_equal(laws$scale[[sd_uv]], abs(A[i_uv, i_uv]), tolerance = 1e-10)

    mixed <- as_mixed_posteriors(setup$fit, setup$leaves, transform_scaled = TRUE)
    densities <- posterior_metadata(mixed, "prior_densities")
    expect_null(densities[[sd_v]])
    expect_true(sd_v %in% posterior_metadata(densities, "measure_unavailable")$column, info = method)
    for(request in list(
      marginal = function() marginal_posterior(mixed, sd_v, use_formula = FALSE, prior_samples = TRUE),
      plot = function() plot_posterior(mixed, sd_v, plot_type = "ggplot", prior = TRUE)
    )){
      refusal <- expect_error(request(), class = "BayesTools_formula_measure_unavailable", info = method)
      expect_identical(refusal$reason, "random_sd_prior_law_unavailable", info = method)
    }

    # the draws keep the entry (the composite of the row), and it matters: the
    # prior probability that sd(v) exceeds 0.9 times its rescaling is far above
    # the one of the fitted law that the rescaling would give
    # Independent two-term expansion for this row, including the stored mean.
    composite <- function(tau) sqrt(tau[, sd_v]^2 +
      (mean(setup$data$u) / stats::sd(setup$data$u))^2 * tau[, sd_uv]^2) / stats::sd(setup$data$v)
    expect_equal(as.numeric(mixed[[sd_v]]), composite(setup$posterior), tolerance = 1e-10)
    threshold <- 0.9 * abs(A[i_v, i_v])
    expect_gt(mean(composite(setup$prior_draws(20000L)) > threshold), 0.9)
    expect_lt(stats::pgamma(0.9, 2, 2, lower.tail = FALSE), 0.5)

    # the catalog reads the same entry: the SD depends on the SD of `u:v`
    catalog <- parameter_catalog(setup$fit)$quantities
    key <- catalog$extraction_key[[match("(mu) sd(v)", catalog$canonical_name)]]
    expect_identical(key$source_type, "composite", info = method)
    expect_true(all(c(sd_v, sd_uv) %in%
      BayesTools:::.bt_parameter_catalog_own_scale_dependencies(setup$fit, key)), info = method)

    # Both paths retain the nonzero mean left by subtraction-centering.
    centered <- random_block_nearly_centered_fit(
      design$rhs, design$scaled, random_block_nearly_centered_data(0)
    )
    A0 <- BayesTools:::.random_sd_unscale_matrix(
      column_terms = centered$term$sd_parameter_names,
      formula_scale = centered$formula_result$formula_scale,
      prefix = "mu", random_term = centered$term, source_data = centered$data
    )
    mixed0 <- as_mixed_posteriors(centered$fit, centered$leaves, transform_scaled = TRUE)
    densities0 <- posterior_metadata(mixed0, "prior_densities")
    expected_entry <- -mean(centered$data$u) / (stats::sd(centered$data$u) * stats::sd(centered$data$v))
    expect_equal(A0[i_v, i_uv], expected_entry, tolerance = 1e-15)
    expect_true(expected_entry != 0)
    expect_null(densities0[[sd_v]])
    refusal <- expect_error(
      marginal_posterior(mixed0, sd_v, use_formula = FALSE, prior_samples = TRUE),
      class = "BayesTools_formula_measure_unavailable"
    )
    expect_identical(refusal$reason, "random_sd_prior_law_unavailable")
    catalog0 <- parameter_catalog(centered$fit)$quantities
    key0 <- catalog0$extraction_key[[match("(mu) sd(v)", catalog0$canonical_name)]]
    expect_identical(key0$source_type, "composite")
    expect_identical(BayesTools:::.bt_parameter_catalog_own_scale_dependencies(centered$fit, key0),
                     c(sd_v, sd_uv))
  }
})

test_that("the SDs of a block with a factor-by-slope term are rescalings or composites by their design row", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  s <- stats::sd(df$x)
  setup <- random_block_prior_draw_fit(~ 1 + x + diag(1 + f + f:x | id), df)
  parameters <- c("mu__xREx__id_intercept", "mu__xREx__id_f", "mu__xREx__id_f__xXx__x")
  mixed <- as_mixed_posteriors(setup$fit, parameters, transform_scaled = TRUE)
  densities <- posterior_metadata(mixed, "prior_densities")
  slopes <- grep("__xXx__x", setup$sd_names, value = TRUE)
  composites <- setdiff(setup$sd_names, slopes)
  expect_identical(length(slopes), 3L)
  expect_identical(length(composites), 3L)

  # the per-level slopes are exact rescalings; the intercept SD and the level SDs
  # absorb the induced shifts and are composites
  grid <- seq(0.02, 3, length.out = 30)
  for(slope in slopes){
    expect_equal(
      random_block_density_ordinate(densities[[slope]], grid),
      stats::dgamma(grid * s, 2, 2) * s, tolerance = 1e-8, info = slope
    )
  }
  for(composite in composites){
    expect_null(densities[[composite]])
  }
  unavailable <- posterior_metadata(densities, "measure_unavailable")
  expect_setequal(unavailable$column, composites)
  expect_true(all(unavailable$cause == "random_sd_prior_law_unavailable"))
  # the elements carry the refusal of their columns
  for(parameter in parameters[1:2]){
    expect_error(
      posterior_measure_check(mixed[[parameter]], "prior_density"),
      class = "BayesTools_formula_measure_unavailable"
    )
  }
  expect_no_error(posterior_measure_check(mixed[[parameters[3L]]], "prior_density"))
  # and the mixed posteriors of the composite SDs alone are still created
  expect_no_error(as_mixed_posteriors(setup$fit, parameters[2L], transform_scaled = TRUE))
})

test_that("a context without the fitted block design refuses the SDs of a scaled block and keeps the others", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  setup <- random_block_prior_draw_fit(~ 1 + x + diag(1 + x | id), df)
  mixed <- as_mixed_posteriors(setup$fit, setup$sd_names, transform_scaled = TRUE)
  context <- posterior_metadata(mixed, "prior_context")
  weights <- function(context, column){
    stats::setNames(as.numeric(context$column_names == column), context$column_names)
  }

  # the design is gone (a context stored by an older version carries none)
  context$formula_design <- NULL
  for(column in setup$sd_names){
    refusal <- expect_error(
      .prior_density_from_context(context, weights(context, column)),
      class = "BayesTools_formula_measure_unavailable", info = column
    )
    expect_identical(refusal$reason, "random_effect_design_unavailable", info = column)
  }
  # the context is not refused where no SD is requested (fixed coefficients)
  expect_no_error(.prior_density_from_context(context, weights(context, "mu_x")))

  # the exported plot of a transformed prior has no design either
  prior_list <- attr(setup$fit, "prior_list")
  formula_scale <- attr(setup$fit, "formula_scale")
  columns <- c("mu_intercept", "mu_x", setup$sd_names)
  refusal <- expect_error(
    plot_transformed_prior(prior_list, columns, formula_scale, parameter = "mu__xREx__id_x"),
    class = "BayesTools_formula_measure_unavailable"
  )
  expect_identical(refusal$reason, "random_effect_design_unavailable")

  # a block without a scaled column keeps its fitted law, which is the law on the original scale
  setup <- random_block_prior_draw_fit(~ 1 + x + (1 | id), df)
  expect_identical(setup$sd_names, "mu__xREx__id_intercept")
  mixed <- as_mixed_posteriors(setup$fit, setup$sd_names, transform_scaled = TRUE)
  context <- posterior_metadata(mixed, "prior_context")
  context$formula_design <- NULL
  grid <- seq(0.02, 3, length.out = 30)
  density <- .prior_density_from_context(context, weights(context, setup$sd_names))
  expect_equal(random_block_density_ordinate(density, grid), stats::dgamma(grid, 2, 2), tolerance = 1e-8)
  expect_null(plot_transformed_prior(
    attr(setup$fit, "prior_list"), c("mu_intercept", "mu_x", setup$sd_names),
    attr(setup$fit, "formula_scale"), parameter = setup$sd_names
  ))
})

test_that("the law of a rescaled SD keeps spike-and-slab and mixture priors, truncation and its support", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  s <- stats::sd(df$x)
  grid <- c(0.1, 0.5, 1, 2)
  sd_density <- function(sd_prior){
    formula_result <- JAGS_formula(
      ~ 1 + x + diag(1 + x | id), "mu", df,
      prior_list = list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
      formula_scale = list(x = TRUE),
      prior_random = prior_random(id = random_block(sd = sd_prior, monitor = random_block_scale_monitor()))
    )
    columns <- unique(unlist(lapply(names(formula_result$prior_list), function(parameter){
      .prior_linear_prior_columns(parameter, formula_result$prior_list[[parameter]])
    })))
    context <- .prior_density_context(
      formula_result$prior_list, columns, list(mu = formula_result$formula_scale),
      formula_design = list(mu = .bt_formula_design_random_unscale_view(formula_result$formula_design))
    )
    weights <- stats::setNames(as.numeric(columns == "mu__xREx__id_x"), columns)
    list(
      density = .prior_density_from_context(context, weights),
      support = .posterior_support_from_prior_context_weights(context, weights),
      intercept = tryCatch(
        .prior_density_from_context(context, stats::setNames(as.numeric(columns == "mu__xREx__id_intercept"), columns)),
        BayesTools_formula_measure_unavailable = function(e) e
      )
    )
  }

  # a spike at zero stays an atom of the same mass; the continuous part is rescaled
  slab <- sd_density(prior_spike_and_slab(
    prior("gamma", list(2, 2)), prior_inclusion = prior("spike", list(0.5))
  ))
  expect_equal(prior_density_ordinate(slab$density, 0)$point_mass, 0.5, tolerance = 1e-12)
  expect_equal(
    random_block_density_ordinate(slab$density, grid),
    0.5 * stats::dgamma(grid * s, 2, 2) * s, tolerance = 1e-8
  )
  expect_s3_class(slab$intercept, "BayesTools_formula_measure_unavailable")

  mixture <- sd_density(prior_mixture(
    list(prior("gamma", list(2, 2)), prior("gamma", list(5, 3))), is_null = c(FALSE, FALSE)
  ))
  expect_equal(
    random_block_density_ordinate(mixture$density, grid),
    0.5 * stats::dgamma(grid * s, 2, 2) * s + 0.5 * stats::dgamma(grid * s, 5, 3) * s,
    tolerance = 1e-8
  )

  # a truncated SD prior keeps its normalization, and its support is rescaled
  truncated <- sd_density(prior("gamma", list(2, 2), list(0, 3)))
  normalizer <- stats::pgamma(3, 2, 2)
  expect_equal(
    random_block_density_ordinate(truncated$density, grid[grid * s < 3]),
    stats::dgamma(grid[grid * s < 3] * s, 2, 2) * s / normalizer, tolerance = 1e-8
  )
  expect_equal(random_block_density_ordinate(truncated$density, 3 / s + 0.1), 0)
  expect_equal(unname(truncated$support$bounds), c(0, 3 / s), tolerance = 1e-12)
})

test_that("a structured block and a block without a scaled column keep the fitted law of their SDs", {

  skip_if_not_installed("runjags")

  set.seed(5)
  df <- data.frame(
    x     = stats::rnorm(36, mean = 10, sd = 2),
    index = factor(rep(c("i1", "i2", "i3"), 12)),
    id    = factor(rep(paste0("g", 1:6), each = 6))
  )
  grid <- c(0.1, 0.5, 1, 2)
  cases <- list(
    list(~ 1 + x + cs(index | id), prior("normal", list(0, 0.5))),
    list(~ 1 + x + ar1(index | id), prior("normal", list(0, 0.5))),
    list(~ 1 + x + (1 | id), NULL)
  )
  for(case in cases){
    formula <- case[[1L]]
    formula_result <- JAGS_formula(
      formula, "mu", df,
      prior_list = list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
      formula_scale = list(x = TRUE),
      prior_random = prior_random(id = random_block(
        sd = prior("gamma", list(2, 2)), cor = case[[2L]],
        monitor = random_block_scale_monitor()
      ))
    )
    sd_name <- grep("__xREx__id_(sd|intercept)$", names(formula_result$prior_list), value = TRUE)
    expect_identical(length(sd_name), 1L)
    columns <- unique(unlist(lapply(names(formula_result$prior_list), function(parameter){
      .prior_linear_prior_columns(parameter, formula_result$prior_list[[parameter]])
    })))
    context <- .prior_density_context(
      formula_result$prior_list, columns, list(mu = formula_result$formula_scale),
      formula_design = list(mu = .bt_formula_design_random_unscale_view(formula_result$formula_design))
    )
    density <- .prior_density_from_context(context, stats::setNames(as.numeric(columns == sd_name), columns))
    expect_equal(
      random_block_density_ordinate(density, grid), stats::dgamma(grid, 2, 2), tolerance = 1e-8,
      info = deparse(formula)
    )
  }
})

test_that("the prior draws of the level comparisons unscale the random-effect SDs with the block design", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  setup <- random_block_prior_draw_fit(~ 1 + x + diag(1 + f + f:x | id), df)
  mixed <- as_mixed_posteriors(
    setup$fit, c("mu__xREx__id_intercept", "mu__xREx__id_f", "mu__xREx__id_f__xXx__x"),
    transform_scaled = TRUE
  )
  context <- posterior_metadata(mixed, "prior_context")
  design <- attr(setup$fit, "formula_design")

  set.seed(3)
  drawn <- .hypothesis_prior_samples_from_context(context, 40L)
  set.seed(3)
  by_design <- .generate_transformed_prior_samples(
    context$prior_list, context$column_names, 40L,
    formula_scale = context$formula_scale, formula_design = design
  )
  set.seed(3)
  by_names <- .generate_transformed_prior_samples(
    context$prior_list, context$column_names, 40L, formula_scale = context$formula_scale
  )
  expect_equal(drawn[, setup$sd_names], by_design[, setup$sd_names], tolerance = 1e-12)
  # pairing names would leave the intercept SD at its fitted value
  expect_gt(max(abs(by_names[, setup$sd_names[1L]] - by_design[, setup$sd_names[1L]])), 1)
  expect_gt(min(drawn[, setup$sd_names[1L]] - by_names[, setup$sd_names[1L]]), 0)
})

test_that("the formula design is an input of the memo of the transformed prior densities", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  setup <- random_block_prior_draw_fit(~ 1 + x + diag(1 + x | id), df)
  prior_list <- attr(setup$fit, "prior_list")
  formula_scale <- attr(setup$fit, "formula_scale")
  columns <- c("mu_intercept", "mu_x", setup$sd_names)
  memo <- .bt_prior_density_memo(setup$fit)
  with_design <- .generate_transformed_prior_densities(
    prior_list, columns, formula_scale, memo = memo,
    formula_design = .bt_formula_designs_random_unscale_view(attr(setup$fit, "formula_design"))
  )
  without_design <- .generate_transformed_prior_densities(
    prior_list, columns, formula_scale, memo = memo
  )
  # the memo does not hand the densities of one design to a call with another
  expect_false(is.null(with_design[["mu__xREx__id_x"]]))
  expect_null(without_design[["mu__xREx__id_x"]])
  expect_identical(
    posterior_metadata(without_design, "measure_unavailable")$cause,
    rep("random_effect_design_unavailable", 2L)
  )
  again <- .generate_transformed_prior_densities(
    prior_list, columns, formula_scale, memo = memo,
    formula_design = .bt_formula_designs_random_unscale_view(attr(setup$fit, "formula_design"))
  )
  expect_identical(again[["mu__xREx__id_x"]], with_design[["mu__xREx__id_x"]])
})

test_that("the density plot of a random-effect SD on the original scale takes no bounds from the fitted prior", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  setup <- random_block_prior_draw_fit(
    ~ 1 + x + diag(1 + x | id), df,
    sd_prior = prior("uniform", list(0, 5)), rsd = function(n) stats::runif(n, 0, 5)
  )
  composite <- "mu__xREx__id_intercept"
  pure <- "mu__xREx__id_x"
  mixed  <- as_mixed_posteriors(setup$fit, setup$sd_names, transform_scaled = TRUE)
  fitted <- as_mixed_posteriors(setup$fit, setup$sd_names, transform_scaled = FALSE)
  area <- function(density) sum(diff(density$x) * (utils::head(density$y, -1L) + utils::tail(density$y, -1L)) / 2)
  curve <- function(samples, parameter) .plot_data_samples.simple(samples, parameter, 256, NULL, NULL, FALSE)$density

  # the composite SD has no exact support: its original-scale draws exceed the
  # fitted bound 5 (they are tau_0 + tau_1 mean(x) combinations), and the density
  # covers them with unit area
  draws <- as.numeric(mixed[[composite]])
  expect_gt(mean(draws > 5), .5)
  expect_null(.posterior_support_bounds(mixed[[composite]], interval_only = TRUE))
  density <- curve(mixed, composite)
  expect_equal(min(density$x), 0)
  expect_gte(max(density$x), max(draws))
  expect_equal(area(density), 1, tolerance = .03)

  # both backends of the public plot
  figure <- plot_posterior(mixed, composite, prior = FALSE, plot_type = "ggplot")
  layer <- ggplot2::ggplot_build(figure)$data[[1L]]
  expect_gte(max(layer$x), max(draws))
  expect_equal(area(layer), 1, tolerance = .03)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_posterior(mixed, composite, prior = FALSE))
  expect_gte(graphics::par("usr")[2L], max(draws))

  # the pure rescaling keeps its rescaled support, and the fitted-scale element
  # of the same model keeps the bounds of its prior
  pure_density <- curve(mixed, pure)
  expect_equal(range(pure_density$x), c(0, 5 / stats::sd(df$x)))
  expect_equal(area(pure_density), 1, tolerance = .03)
  fitted_density <- curve(fitted, composite)
  expect_equal(range(fitted_density$x), c(0, 5))
  expect_equal(area(fitted_density), 1, tolerance = .03)
})

test_that("the posterior plot of a prior with multiply_by keeps the block design of the prior context", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  df$z <- stats::rnorm(nrow(df))
  s <- stats::sd(df$x)
  z_prior <- prior("normal", list(0, 1))
  attr(z_prior, "multiply_by") <- .5
  curves <- list()
  for(multiplied in c(FALSE, TRUE)){
    prior_list <- list(
      intercept = prior("normal", list(0, 1)),
      x = prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("spike", list(.5))),
      z = if(multiplied) z_prior else prior("normal", list(0, 1))
    )
    formula_result <- JAGS_formula(
      ~ 1 + x + z + diag(1 + x | id), "mu", df, prior_list = prior_list,
      formula_scale = list(x = TRUE),
      prior_random = prior_random(id = random_block(
        sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
      ))
    )
    sd_names <- unlist(lapply(formula_result$formula_design$random_effects, function(term){
      unique(term$sd_parameter_names)
    }))
    set.seed(2)
    n <- 2000L
    inclusion <- stats::rbinom(n, 1, .5)
    posterior <- cbind(
      mu_intercept = stats::rnorm(n), mu_x = inclusion * stats::rnorm(n), mu_x_indicator = inclusion,
      mu_z = stats::rnorm(n),
      matrix(stats::rgamma(n * length(sd_names), 2, 2), n, length(sd_names), dimnames = list(NULL, sd_names))
    )
    fit <- make_random_scale_table_fit(formula_result, posterior)
    mixed <- as_mixed_posteriors(fit, c("mu_x", sd_names), conditional = "mu_x", transform_scaled = TRUE)

    # the context rebuilt without the multiplier has the design of the stored context
    context <- .plot_data_prior_context_for_plot(
      posterior_metadata(mixed, "prior_context"), mixed, "mu__xREx__id_x", 64L
    )
    expect_false(is.null(context$formula_design$mu))
    expect_identical(.plot_data_prior_list_without_multiply_by(attr(mixed, "prior_list"))$changed, multiplied)

    # the prior curve of the pure SD is exact (the fitted gamma(2, 2) rescaled)
    figure <- plot_posterior(mixed, "mu__xREx__id_x", prior = TRUE, plot_type = "ggplot")
    prior_layer <- ggplot2::ggplot_build(figure)$data[[1L]]
    inside <- prior_layer$x > .05 & prior_layer$x < 1.2
    expect_gt(sum(inside), 10L)
    expect_equal(prior_layer$y[inside], stats::dgamma(prior_layer$x[inside] * s, 2, 2) * s, tolerance = 1e-8)
    curves[[as.character(multiplied)]] <- prior_layer[, c("x", "y")]
    marginal <- marginal_posterior(mixed, "mu__xREx__id_x", use_formula = FALSE, prior_samples = TRUE)
    expect_equal(
      exp(prior_density_ordinate(posterior_metadata(marginal, "prior_density"), .3)$log_density),
      stats::dgamma(.3 * s, 2, 2) * s, tolerance = 1e-8
    )
  }
  expect_equal(curves[["TRUE"]], curves[["FALSE"]])
})

test_that("the parameter map builds the unscale matrix of a scaled block once", {

  skip_if_not_installed("runjags")

  # the support of each composite SD reads one row of the block's unscale
  # matrix, which the map build computes once per block and not once per SD
  df <- random_block_scale_data()
  df$g <- factor(rep(c("p", "q"), length.out = 36))
  block <- prior_random(id = random_block(
    sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
  ))
  original <- BayesTools:::.bt_parameter_prior_density_random_unscale_matrix
  calls <- 0L
  setup <- testthat::with_mocked_bindings(
    random_block_scale_fit(~ 1 + x + diag(1 + f * g * x | id), df, block),
    .bt_parameter_prior_density_random_unscale_matrix = function(...){
      calls <<- calls + 1L
      original(...)
    },
    .package = "BayesTools"
  )

  expect_gt(length(setup$sd_names), 6L)
  expect_identical(calls, 1L)

  # every SD and variance of the block keeps its exact hull [0, Inf)
  catalog <- parameter_catalog(setup$fit)
  composite <- which(catalog$quantities$quantity %in% c("sd", "var"))
  expect_gt(length(composite), 12L)
  for(support in catalog$quantities$support[composite]){
    expect_true(support$exact)
    expect_identical(unname(support$bounds), c(0, Inf))
  }
})

test_that("a scaled random block with a linearly dependent design keeps the fit and refuses original-scale SDs", {

  skip_if_not_installed("runjags")

  # `x:f:f2` without its margins makes the block design linearly dependent: the
  # original-scale coefficients of the block are not identified. The map build
  # does not stop (the fit would be lost), the SDs of the block are composites
  # without an original-scale law, and only original-scale requests refuse.
  df <- random_block_scale_data()
  df$f2 <- factor(rep(c("u", "v"), each = 18))
  block <- prior_random(id = random_block(
    sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
  ))
  expect_not_identified <- function(expr, label){
    refusal <- expect_error(
      expr, class = "BayesTools_formula_transform_unavailable", info = label
    )
    expect_identical(refusal$reason, "original_scale_not_identified", info = label)
  }

  for(formula in list(
    ~ 1 + x + diag(1 + x + f + x:f:f2 | id),
    ~ 1 + x + diag(1 + x + x:f:f2 | id)
  )){
    label <- paste(deparse(formula), collapse = "")
    setup <- random_block_scale_fit(formula, df, block)
    fit <- setup$fit
    expect_gt(length(setup$sd_names), 6L)

    catalog <- parameter_catalog(fit)
    sd_rows <- which(catalog$quantities$quantity == "sd")
    expect_identical(length(sd_rows), length(setup$sd_names), info = label)
    expect_true(all(vapply(
      catalog$quantities$extraction_key[sd_rows],
      function(key) identical(key$source_type, "composite"), logical(1)
    )), info = label)

    # the fitted-scale tables work: the fixed effects and, in the raw summary,
    # the fitted SDs; the standard summary leaves the original-scale SDs out
    # and says so
    table <- JAGS_estimates_table(fit, remove_diagnostics = TRUE)
    expect_identical(rownames(table), c("(mu) intercept", "(mu) x"), info = label)
    expect_match(attr(table, "warnings"), "linearly dependent", fixed = TRUE, info = label)
    expect_match(attr(table, "warnings"), "random_effects_summary", fixed = TRUE, info = label)
    raw <- JAGS_estimates_table(fit, remove_diagnostics = TRUE, random_effects_summary = "raw")
    expect_equal(
      unname(as.data.frame(raw)[grep("sd(", rownames(raw), fixed = TRUE), "Mean"]),
      setup$tau, tolerance = 1e-12, info = label
    )

    # original-scale requests refuse with the class and the reason
    expect_not_identified(JAGS_estimates_table(fit, transform_scaled = TRUE), label)
    expect_not_identified(transform_scale_samples(fit), label)
    expect_not_identified(
      parameter_draws(fit, parameter_catalog_resolve(catalog, "(mu) sd(intercept)")), label
    )
    expect_not_identified(as_mixed_posteriors(fit, setup$sd_names, transform_scaled = TRUE), label)
    # and so does the prior law of an SD on the original scale
    context <- .prior_density_context(
      attr(fit, "prior_list"), c("mu_intercept", "mu_x", setup$sd_names),
      formula_scale = attr(fit, "formula_scale"),
      formula_design = list(mu = .bt_formula_design_random_unscale_view(attr(fit, "formula_design")$mu))
    )
    expect_not_identified(
      .prior_density_from_context(
        context, stats::setNames(as.numeric(context$column_names == setup$sd_names[2L]), context$column_names)
      ),
      label
    )
  }

  # a block with an independent design keeps its original-scale SDs
  setup <- random_block_scale_fit(~ 1 + x + diag(1 + x + f | id), df, block)
  table <- JAGS_estimates_table(setup$fit, remove_diagnostics = TRUE)
  expect_identical(sum(grepl("sd(", rownames(table), fixed = TRUE)), 4L)
  expect_length(attr(table, "warnings"), 0L)
})

test_that("the table of a block that drops induced terms says what its intercept SD is", {

  skip_if_not_installed("runjags")

  df <- random_block_scale_data()
  block <- prior_random(id = random_block(
    sd = prior("gamma", list(2, 2)), monitor = random_block_scale_monitor()
  ))
  note <- paste0(
    "(mu) sd(intercept): the SD of the group effect at the predictor mean; ",
    "the intercept shifts that the scaled slope induces by level are not ",
    "columns of the block."
  )

  # `f:x` induces the level intercepts, which no column of the block carries
  setup <- random_block_scale_fit(~ 1 + x + diag(1 + f:x | id), df, block)
  table <- JAGS_estimates_table(setup$fit, remove_diagnostics = TRUE)
  expect_identical(attr(table, "footnotes"), stats::setNames(note, "(mu) sd(intercept)"))
  expect_identical(
    attr(table[c("(mu) intercept", "(mu) sd(f[a]:x)"), ], "footnotes"), NULL
  )
  expect_identical(
    attr(table["(mu) sd(intercept)", ], "footnotes"), stats::setNames(note, "(mu) sd(intercept)")
  )
  # the fitted-scale rows of the raw summary need no remark
  raw <- JAGS_estimates_table(setup$fit, remove_diagnostics = TRUE, random_effects_summary = "raw")
  expect_length(attr(raw, "footnotes"), 0L)

  # blocks whose design is representable report the intercept of a level
  for(formula in list(~ 1 + x + diag(1 + x | id), ~ 1 + x + diag(1 + f + f:x | id))){
    setup <- random_block_scale_fit(formula, df, block)
    table <- JAGS_estimates_table(setup$fit, remove_diagnostics = TRUE)
    expect_length(attr(table, "footnotes"), 0L)
  }
})
