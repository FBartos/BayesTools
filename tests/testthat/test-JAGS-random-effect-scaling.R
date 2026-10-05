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
