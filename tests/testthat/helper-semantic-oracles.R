bayestools_oracle_gaussian_regression_data <- function(seed = 20260504L, n = 48L) {
  set.seed(seed)

  x1 <- rnorm(n, mean = 8, sd = 3)
  x2 <- rnorm(n, mean = -2, sd = 4)
  x1_scaled <- as.numeric(scale(x1))
  x2_scaled <- as.numeric(scale(x2))

  data.frame(
    y = 1.5 + 0.8 * x1_scaled - 0.5 * x2_scaled +
      0.35 * x1_scaled * x2_scaled + rnorm(n, mean = 0, sd = 0.6),
    x1 = x1,
    x2 = x2
  )
}

bayestools_manual_scaled_data <- function(data, variables) {
  out <- data
  scale_info <- list()

  for (variable in variables) {
    scale_info[[paste0("mu_", variable)]] <- list(
      mean = mean(data[[variable]]),
      sd = stats::sd(data[[variable]])
    )
    out[[variable]] <- (data[[variable]] - scale_info[[paste0("mu_", variable)]]$mean) /
      scale_info[[paste0("mu_", variable)]]$sd
  }

  attr(out, "manual_scale") <- scale_info
  out
}

bayestools_oracle_formula_design_data <- function() {
  data.frame(
    x = c(-3, -1, 0, 2, 4, 5),
    z = c(10, 9, 8, 7, 6, 4),
    x_alias = 2 * c(-3, -1, 0, 2, 4, 5),
    group = factor(c("b", "a", "c", "b", "a", "c"), levels = c("a", "b", "c", "unused")),
    group_reordered = factor(c("b", "a", "c", "b", "a", "c"), levels = c("c", "a", "b", "unused")),
    a = factor(c("a1", "a2", "a1", "a2", "a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b1", "b2", "b2", "b3", "b3"), levels = c("b1", "b2", "b3"))
  )
}

bayestools_oracle_nonsyntactic_formula_data <- function() {
  out <- data.frame(y = 1:4, check.names = FALSE)
  out[["x weird"]] <- c(-1, 0, 1, 2)
  out
}

bayestools_lm_coef_to_jags_names <- function(coef_names, prefix = "mu") {
  out <- coef_names
  out <- gsub("\\(Intercept\\)", "intercept", out)
  out <- gsub(":", "__xXx__", out, fixed = TRUE)
  paste0(prefix, "_", out)
}

expect_model_matrix_equal <- function(actual, formula, data, tolerance = 1e-12) {
  expected <- stats::model.matrix(formula, data = data)
  testthat::expect_equal(unname(actual), unname(expected), tolerance = tolerance)
  testthat::expect_equal(colnames(actual), colnames(expected))
  invisible(actual)
}

expect_jags_formula_design_matches_model_matrix <- function(design, formula, data,
                                                            tolerance = 1e-12) {
  expected <- stats::model.matrix(formula, data = data)
  expected_jags_names <- gsub(":", "__xXx__", colnames(expected), fixed = TRUE)

  testthat::expect_equal(unname(design$model_matrix), unname(expected), tolerance = tolerance)
  testthat::expect_equal(design$raw_column_names, colnames(expected))
  testthat::expect_equal(colnames(design$model_matrix), expected_jags_names)
  testthat::expect_equal(design$column_names, expected_jags_names)
  testthat::expect_equal(design$assign, attr(expected, "assign"))
  invisible(design)
}

expect_formula_scale_equal <- function(actual, expected, tolerance = 1e-12) {
  testthat::expect_setequal(names(actual), names(expected))
  for (term in names(expected)) {
    testthat::expect_equal(actual[[term]]$mean, expected[[term]]$mean, tolerance = tolerance)
    testthat::expect_equal(actual[[term]]$sd, expected[[term]]$sd, tolerance = tolerance)
  }
  invisible(actual)
}

expect_posterior_summary_close <- function(posterior, expected, tolerance, statistic = mean) {
  observed <- apply(posterior[, names(expected), drop = FALSE], 2, statistic)
  testthat::expect_equal(unname(observed), unname(expected), tolerance = tolerance)
  invisible(observed)
}

expect_lm_predictions_equal <- function(lm_fit, posterior, newdata, formula, prefix = "mu",
                                        tolerance = 0.25) {
  design <- stats::model.matrix(formula, data = newdata)
  posterior_names <- bayestools_lm_coef_to_jags_names(colnames(design), prefix = prefix)
  testthat::expect_true(all(posterior_names %in% colnames(posterior)))

  posterior_means <- colMeans(posterior[, posterior_names, drop = FALSE])
  observed <- as.numeric(design %*% posterior_means)
  expected <- as.numeric(stats::predict(lm_fit, newdata = newdata))

  testthat::expect_equal(observed, expected, tolerance = tolerance)
  invisible(observed)
}

expect_design_predictions_equal <- function(design, coefficients, expected, tolerance = 1e-12) {
  testthat::expect_true(all(colnames(design) %in% names(coefficients)))
  observed <- as.numeric(design %*% coefficients[colnames(design)])
  testthat::expect_equal(observed, as.numeric(expected), tolerance = tolerance)
  invisible(observed)
}

bayestools_gaussian_posterior_oracle <- function(X, y, sigma, prior_mean, prior_sd) {
  testthat::expect_equal(length(prior_mean), ncol(X))
  testthat::expect_equal(length(prior_sd), ncol(X))

  prior_precision <- diag(1 / prior_sd^2, nrow = ncol(X))
  posterior_cov <- solve(crossprod(X) / sigma^2 + prior_precision)
  posterior_mean <- posterior_cov %*% (crossprod(X, y) / sigma^2 + prior_precision %*% prior_mean)

  list(
    mean = as.numeric(posterior_mean),
    cov = posterior_cov,
    sd = sqrt(diag(posterior_cov))
  )
}

expect_probability_mass <- function(observed, expected = 1, tolerance = 1e-8) {
  testthat::expect_equal(sum(observed), expected, tolerance = tolerance)
  testthat::expect_true(all(is.finite(observed)))
  testthat::expect_true(all(observed >= 0))
  invisible(observed)
}

expect_table_semantics <- function(table, expected_rows = NULL, expected_columns = NULL) {
  if (!is.null(expected_rows)) {
    testthat::expect_equal(rownames(table), expected_rows)
  }
  if (!is.null(expected_columns)) {
    testthat::expect_true(all(expected_columns %in% colnames(table)))
  }
  invisible(table)
}

expect_plot_data_semantics <- function(plot, expected_layers = NULL) {
  testthat::expect_s3_class(plot, "ggplot")
  if (!is.null(expected_layers)) {
    testthat::expect_equal(length(plot$layers), expected_layers)
  }
  testthat::expect_false(is.null(plot$data))
  invisible(plot)
}
