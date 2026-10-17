skip_if_not_test_profile("unit")

.selected_replay_case <- function(formula, contrasts){
  data <- expand.grid(g = factor(c("c", "a", "b"), levels = c("c", "a", "b")),
    h = factor(c("v", "u"), levels = c("v", "u")), x = c(-2, 1, 3))
  labels <- attr(stats::terms(formula), "term.labels")
  priors <- list(intercept = prior("normal", list(0, 1)))
  for(label in labels){
    priors[[label]] <- if(any(strsplit(label, ":", fixed = TRUE)[[1]] %in% c("g", "h"))){
      prior_factor("normal", list(0, 1), contrast = contrasts)
    }else prior("normal", list(0, 1))
  }
  compiled <- JAGS_formula(formula, "mu", data, priors)
  design <- compiled$formula_design
  terms <- .bt_dnode_linear_predictor_design_terms(design, compiled$prior_list)
  names <- unlist(lapply(terms, .bt_dnode_linear_predictor_coefficient_names), use.names = FALSE)
  draws <- matrix(seq_along(names) * c(2, 7)[rep(1:2, length.out = length(names))], 1L,
    dimnames = list(NULL, names))
  list(data = data, priors = priors, design = design, draws = draws,
    carrier = JAGS_formula_draws(draws, formula, "mu", data, priors))
}

test_that("selected interactions retain fitted coding, columns, and variable order", {
  for(contrast in c("treatment", "independent")){
    for(formula in list(~ g * x, ~ g + g:x, ~ g * h, ~ x + g * h)){
      case <- .selected_replay_case(formula, contrast)
      interaction <- tail(attr(case$design$terms, "term.labels"), 1L)
      selected <- stats::as.formula(paste("~ 0 +", paste(rev(strsplit(interaction, ":", fixed = TRUE)[[1L]]), collapse = ":")))
      term_index <- match(interaction, attr(case$design$terms, "term.labels"))
      columns <- which(case$design$assign == term_index)
      expected <- case$design$model_matrix[, columns, drop = FALSE] %*% as.vector(case$draws[, columns])
      expect_equal(unname(JAGS_evaluate_formula(case$carrier, selected, "mu")), unname(expected))
      rows <- c(nrow(case$data), 2L, 2L, 1L)
      needed <- all.vars(selected)
      prediction <- JAGS_predict_formula(case$carrier, "mu", formula = selected,
        data = case$data[rows, needed, drop = FALSE], formula_target = "fixed")
      expect_equal(unname(prediction$value), unname(expected[rows, , drop = FALSE]))
    }
  }
})

test_that("reordered main effects and orthonormal factor subsets retain fitted replay", {
  data <- data.frame(g = factor(c("c", "a", "b", "c"), levels = c("c", "a", "b")), x = c(-2, -1, 1, 3))
  priors <- list(intercept = prior("point", list(2)),
    g = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"), x = prior("normal", list(0, 1)))
  draws <- cbind(mu_intercept = 2, "mu_g[1]" = 3, "mu_g[2]" = 5, mu_x = 7)
  carrier <- JAGS_formula_draws(draws, ~ g + x, "mu", data, priors)
  full <- JAGS_evaluate_formula(carrier, parameter = "mu")
  expect_equal(JAGS_evaluate_formula(carrier, ~ x + g, "mu"), full)
  expect_equal(unname(JAGS_evaluate_formula(carrier, ~ 0 + x, "mu", data = data["x"])), matrix(7 * data$x, 4L))
  expect_error(JAGS_evaluate_formula(carrier, ~ x:g, "mu"), "unavailable in the fitted formula")
  no_intercept <- JAGS_formula_draws(cbind(mu_x = c(3, 5)), ~ 0 + x, "mu", data["x"],
    list(x = prior("normal", list(0, 1))))
  expect_equal(unname(JAGS_evaluate_formula(no_intercept, parameter = "mu")), outer(data$x, c(3, 5)))
  scalar <- JAGS_formula_draws(cbind(mu_intercept = c(3, 5)), ~ 1, "mu", data["x"],
    list(intercept = prior("normal", list(0, 1))))
  expect_equal(unname(JAGS_evaluate_formula(scalar, parameter = "mu")), matrix(rep(c(3, 5), each = 4), 4, 2))
})

test_that("replay and selected terms follow the fitted variable order, not the order of the selected expression", {
  data <- data.frame(
    a = factor(c("no", "yes", "no", "yes", "yes", "no", "no", "yes"), levels = c("no", "yes")),
    b = factor(c("no", "no", "yes", "yes", "no", "no", "yes", "yes"), levels = c("no", "yes"))
  )
  normal <- prior("normal", list(0, 1))
  treatment <- prior_factor("normal", list(0, 1), contrast = "treatment")

  # `a:b + b` lists `b` before `a:b`, so the expression's first appearance of
  # variables (a, b) differs from the fitted variable order.
  draws <- cbind(mu_intercept = 1, mu_b = 3, "mu_a__xXx__b[1]" = 7, "mu_a__xXx__b[2]" = 11)
  carrier <- JAGS_formula_draws(draws, ~ a:b + b, "mu", data,
    list(intercept = normal, b = treatment, "a:b" = treatment))
  design <- attr(carrier, "formula_design")$mu
  fitted <- drop(design$model_matrix %*% c(1, 3, 7, 11))
  expect_equal(drop(JAGS_evaluate_formula(carrier, parameter = "mu")), fitted)
  rows <- c(8L, 1L, 6L, 3L)
  expect_equal(
    drop(JAGS_evaluate_formula(carrier, parameter = "mu", data = data[rows, ])),
    fitted[rows]
  )

  # A selection of a full-factorial fit keeps the fitted coding of each column.
  draws_full <- cbind(mu_intercept = 1, mu_a = 2, mu_b = 3, mu_a__xXx__b = 7)
  carrier_full <- JAGS_formula_draws(draws_full, ~ a * b, "mu", data,
    list(intercept = normal, a = treatment, b = treatment, "a:b" = treatment))
  design_full <- attr(carrier_full, "formula_design")$mu
  columns <- which(design_full$assign %in% c(2L, 3L))
  expect_equal(
    drop(JAGS_evaluate_formula(carrier_full, ~ 0 + b + a:b, "mu")),
    drop(design_full$model_matrix[, columns, drop = FALSE] %*% c(3, 7))
  )
})
