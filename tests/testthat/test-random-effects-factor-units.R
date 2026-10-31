skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Exchangeable Factor Units in Random-Effect Blocks
# ============================================================================ #
#
# PURPOSE:
#   A mean-difference or orthonormal factor term is one exchangeable unit of a
#   random-effect block: one shared SD, uncorrelated coordinates, the LKJ prior
#   restricted to an identity unit block under us(), and reported quantities
#   that do not depend on the arbitrary contrast basis.
#
# TAGS: @unit, @random-effects, @LKJ, @factor
# ============================================================================ #

.unit_data <- function(){
  set.seed(1)
  data.frame(
    x  = stats::rnorm(96),
    f  = factor(rep(c("a", "b", "c", "d"), 24)),
    f2 = factor(rep(c("u", "v"), each = 4, length.out = 96)),
    g  = factor(rep(c("p", "q", "r"), each = 8, length.out = 96)),
    h  = factor(rep(c("k", "l", "m"), each = 2, length.out = 96)),
    id = factor(rep(1:12, each = 8))
  )
}

.unit_sd_prior <- prior("normal", list(0, 1), list(0, Inf))

# `fixed`: the predictors of the fixed part of the formula
.unit_formula <- function(formula, fixed, prior_random = NULL, contrast = "meandif",
                          data = .unit_data()){
  priors <- list(
    intercept = prior("normal", list(0, 1)),
    x  = prior("normal", list(0, 1)),
    f  = prior_factor("mnormal", list(0, 1), contrast = contrast),
    f2 = prior_factor("mnormal", list(0, 1), contrast = contrast),
    g  = prior_factor("mnormal", list(0, 1), contrast = contrast),
    h  = prior_factor("normal", list(0, 1), contrast = "treatment")
  )
  JAGS_formula(
    formula,
    parameter = "mu",
    data = data,
    prior_list = priors[c("intercept", fixed)],
    prior_random = if(is.null(prior_random)){
      prior_random(id = random_block(sd = .unit_sd_prior))
    }else{
      prior_random
    }
  )
}

.unit_term <- function(result){
  result$formula_design$random_effects[[1L]]
}

test_that("a mean-difference or orthonormal factor term shares one SD", {

  for(contrast in c("meandif", "orthonormal")){
    for(formula in list(~ 1 + f + diag(1 + f | id), ~ 1 + f + us(1 + f | id))){
      result <- .unit_formula(formula, "f", contrast = contrast)
      term <- .unit_term(result)
      expect_identical(
        term$sd_parameter_names,
        c("mu__xREx__id_intercept", rep("mu__xREx__id_f", 3L))
      )
      expect_identical(unname(term$sd_leaves$leaf_terms), c("intercept", "f"))
      expect_identical(term$sd_leaves$unit_index_by_column, c(NA, 1L, 1L, 1L))
      expect_identical(term$sd_leaves$units[[1L]]$columns, 2:4)
      # one scalar SD prior of the term
      expect_true(is.prior.simple(result$prior_list[["mu__xREx__id_f"]]))
      expect_false(is.prior.factor(result$prior_list[["mu__xREx__id_f"]]))
      # the cell weights are the rows of the factor's contrast
      expect_equal(
        unname(term$sd_leaves$units[[1L]]$weights),
        unname(.factor_contrast_matrix(c("a", "b", "c", "d"), paste0("contr.", contrast))),
        tolerance = 1e-14
      )
      expect_identical(term$sd_leaves$units[[1L]]$cells, paste0("f[", c("a", "b", "c", "d"), "]"))
    }
  }

  # a homogeneous block keeps its one SD
  term <- .unit_term(.unit_formula(~ 1 + f + id(1 + f | id), "f"))
  expect_identical(term$sd_parameter_names, rep("mu__xREx__id_sd", 4L))
  expect_length(term$sd_leaves$units, 0L)
})

test_that("an interaction is one unit, or one unit per column of a partner coding", {

  term <- .unit_term(.unit_formula(~ 1 + x + f + diag(1 + x + x:f | id), c("x", "f")))
  expect_identical(
    term$sd_parameter_names,
    c("mu__xREx__id_intercept", "mu__xREx__id_x", rep("mu__xREx__id_x__xXx__f", 3L))
  )

  # f (meandif) crossed with h (treatment): one unit per treatment column of h
  term <- .unit_term(.unit_formula(~ 1 + f + h + diag(1 + f + h + f:h | id), c("f", "h")))
  units <- term$sd_leaves$units
  expect_length(units, 3L)
  expect_identical(
    vapply(units, `[[`, character(1), "term"),
    c("f", "f__xXx__h[l]", "f__xXx__h[m]")
  )
  columns <- term$column_names
  for(unit in units[-1L]){
    level <- sub("^f__xXx__h\\[(.*)\\]$", "\\1", unit$term)
    expect_true(all(endsWith(columns[unit$columns], paste0("__xXx__h", level))))
    expect_identical(length(unit$columns), 3L)
    expect_equal(
      unname(unit$weights),
      unname(contr.meandif(c("a", "b", "c", "d"))),
      tolerance = 1e-14
    )
  }
  expect_identical(
    unique(term$sd_parameter_names[unlist(lapply(units[-1L], `[[`, "columns"))]),
    paste0("mu__xREx__id_f__xXx__h[", 1:2, "]")
  )

  # an indicator-coded factor (its margin absent) has level cells, not a unit
  term <- .unit_term(.unit_formula(~ 1 + x + f + diag(1 + f:x | id), c("x", "f")))
  expect_length(term$sd_leaves$units, 0L)
})

test_that("us() orders a unit first in its LKJ vine with an identity unit block", {

  result <- .unit_formula(~ 1 + x + f + us(1 + x + f | id), c("x", "f"))
  term <- .unit_term(result)
  expect_identical(term$correlation$unit_columns, 3:5)
  layout <- .bt_random_effect_lkj_layout(term)
  expect_identical(layout$vine_order, c(3L, 4L, 5L, 1L, 2L))
  expect_identical(layout$position, c(4L, 5L, 1L, 2L, 3L))
  expect_identical(
    term$correlation$primitive_names,
    paste0("mu__xREx__id_xRE_CORx_lkj_u[", 4:10, "]")
  )
  for(p in 1:3){
    expect_match(result$formula_syntax, paste0("mu__xREx__id_xRE_CORx_lkj_u[", p, "] <- 0.5"), fixed = TRUE)
  }
  expect_match(
    result$formula_syntax,
    "mu__xREx__id_xRE_CORx_lkj_u[4:10] ~ dbt_lkj_cpc(mu__xREx__id_xRE_CORx_lkj_alpha[4:10])",
    fixed = TRUE
  )
  expect_match(result$formula_syntax, "bt_lkj_cholesky_ordered(", fixed = TRUE)

  # a two-level factor (one column) is an ordinary column in the design order
  result <- .unit_formula(~ 1 + f2 + us(1 + f2 | id), "f2")
  term <- .unit_term(result)
  expect_identical(term$correlation$unit_columns, integer())
  expect_identical(term$correlation$primitive_names, "mu__xREx__id_xRE_CORx_lkj_u[1]")
  expect_false(grepl("bt_lkj_cholesky_ordered", result$formula_syntax, fixed = TRUE))
  expect_identical(term$sd_parameter_names, c("mu__xREx__id_intercept", "mu__xREx__id_f2"))
  # with a unit of three columns, the two-level factor follows the unit
  term <- .unit_term(.unit_formula(~ 1 + f + f2 + us(1 + f2 + f | id), c("f", "f2")))
  expect_identical(term$correlation$unit_columns, 3:5)

  # a block that is one unit has no correlation parameter
  result <- .unit_formula(~ 1 + f + us(0 + f | id), "f")
  expect_null(.unit_term(result)$correlation)
  expect_false(grepl("xRE_CORx", result$formula_syntax, fixed = TRUE))
  expect_error(
    .unit_formula(~ 1 + f + us(0 + f | id), "f", prior_random = prior_random(id = random_block(
      sd = .unit_sd_prior,
      cor = prior_lkj(eta = 1)
    ))),
    paste0(
      "Random-effect structure 'us' of a block that is one meandif or orthonormal ",
      "factor term has no correlation parameter (its columns share one SD and are ",
      "uncorrelated); remove the 'cor' prior."
    ),
    fixed = TRUE
  )
})

test_that("models that a unit cannot represent are refused, not restricted", {

  condition <- tryCatch(
    .unit_formula(~ 1 + f + g + us(1 + f + g | id), c("f", "g")),
    error = function(error) error
  )
  expect_s3_class(condition, "BayesTools_random_effects_exchangeable_units_unavailable")
  expect_s3_class(condition, "BayesTools_random_effects_error")
  expect_identical(condition$reason, "several_correlated_units")
  expect_identical(condition$terms, c("f", "g"))
  expect_identical(
    conditionMessage(condition),
    paste0(
      "Correlations of random-effect block 'id' are unavailable: structure 'us' ",
      "correlates several meandif or orthonormal factor terms with three or more ",
      "levels ('f', 'g'), each of which is one exchangeable unit, and correlated ",
      "units are not implemented. Use structure 'diag', place the terms in ",
      "separate blocks, or code the factors with contr.treatment or contr.independent."
    )
  )
  # the same terms are available without correlations
  expect_no_error(.unit_formula(~ 1 + f + g + diag(1 + f + g | id), c("f", "g")))

  condition <- tryCatch(
    .unit_formula(~ 1 + f + diag(1 + f | id), "f", prior_random = prior_random(id = random_block(
      sd = .unit_sd_prior,
      terms = list(f = prior_factor("spike", list(0), contrast = "meandif"))
    ))),
    error = function(error) error
  )
  expect_s3_class(condition, "BayesTools_random_effects_exchangeable_units_unavailable")
  expect_identical(condition$reason, "unit_sd_prior")
  expect_identical(
    conditionMessage(condition),
    paste0(
      "Random-effect SD prior for 'f' must be an ordinary scalar prior: the columns ",
      "of a meandif or orthonormal factor term share one SD, because their coding ",
      "basis is arbitrary. Use a scalar prior (e.g., prior(\"point\", list(0)) for a ",
      "point mass), or contr.treatment or contr.independent coding for an SD per level."
    )
  )
  # the scalar point mass is the same model
  term <- .unit_term(.unit_formula(~ 1 + f + diag(1 + f | id), "f", prior_random = prior_random(id = random_block(
    sd = .unit_sd_prior,
    terms = list(f = prior("point", list(0)))
  ))))
  expect_identical(term$sd_parameter_names[2:4], rep("mu__xREx__id_f", 3L))
})

test_that("correlation rows report level cells and the multiple correlation", {

  term <- .unit_term(.unit_formula(~ 1 + f + f2 + us(1 + f + f2 | id), c("f", "f2")))
  rows <- .bt_random_effect_correlation_rows(term)
  expect_identical(
    vapply(rows, `[[`, character(1), "label"),
    c(paste0("intercept,f[", c("a", "b", "c", "d"), "]"), "intercept,f",
      paste0("intercept,f2[", c("u", "v"), "]"), "intercept,f2",
      paste0("f2[", rep(c("u", "v"), each = 4L), "],f[", c("a", "b", "c", "d"), "]"),
      "f2,f")
  )
  expect_identical(
    vapply(rows, `[[`, character(1), "quantity"),
    rep(c("cor", "cor_total", "cor", "cor_total", "cor", "cor_total"), c(4, 1, 2, 1, 8, 1))
  )

  # correlation matrices of the restricted LKJ prior
  layout <- .bt_random_effect_lkj_layout(term)
  set.seed(2)
  u <- matrix(stats::rbeta(50L * length(layout$free_pairs), 2, 2), ncol = length(layout$free_pairs))
  R <- .bt_lkj_layout_u_to_R(u, layout)
  values <- .bt_random_effect_correlation_row_values(rows, R)
  C <- contr.meandif(c("a", "b", "c", "d"))
  for(draw in seq_len(nrow(u))){
    R_draw <- R[draw, , ]
    # cor(intercept, f[a]) = C[a, ] r / |C[a, ]| for the correlations r of the
    # intercept with the unit's columns
    r <- R_draw[2:4, 1L]
    expect_equal(values[draw, 1:4], as.vector(C %*% r) / sqrt(rowSums(C^2)), tolerance = 1e-12)
    expect_equal(sum(values[draw, 1:4]), 0, tolerance = 1e-12)
    expect_equal(values[draw, 5L], sqrt(sum(r^2)), tolerance = 1e-12)
    # the two cells of a two-level factor are -1 times each other
    expect_equal(values[draw, 6L], -values[draw, 7L], tolerance = 1e-12)
    expect_equal(values[draw, 8L], abs(R_draw[5L, 1L]), tolerance = 1e-12)
    expect_equal(values[draw, 17L], sqrt(sum(R_draw[2:4, 5L]^2)), tolerance = 1e-12)
  }
})

test_that("the reported quantities do not depend on the contrast basis", {

  formula <- ~ 1 + x + f + us(1 + x + f | id)
  term <- .unit_term(.unit_formula(formula, c("x", "f")))
  # a rotation with a reflection of the meandif basis
  Q <- qr.Q(qr(matrix(c(0.3, -1.2, 0.5, 0.8, 0.1, -0.4, 1.1, 0.7, 0.2), 3L)))
  Q[, 1L] <- -Q[, 1L]
  original <- contr.meandif
  local_mocked_bindings(
    contr.meandif = function(n, contrasts = TRUE){
      out <- original(n, contrasts)
      if(isTRUE(contrasts) && ncol(out) == 3L) out %*% Q else out
    }
  )
  rotated <- .unit_term(.unit_formula(formula, c("x", "f")))
  expect_equal(
    unname(rotated$sd_leaves$units[[1L]]$weights),
    unname(term$sd_leaves$units[[1L]]$weights %*% Q),
    tolerance = 1e-14
  )
  expect_identical(rotated$sd_parameter_names, term$sd_parameter_names)
  expect_identical(rotated$correlation$primitive_names, term$correlation$primitive_names)
  rows <- .bt_random_effect_correlation_rows(term)
  rotated_rows <- .bt_random_effect_correlation_rows(rotated)
  expect_identical(
    vapply(rotated_rows, `[[`, character(1), "label"),
    vapply(rows, `[[`, character(1), "label")
  )

  # the same level-cell model in the rotated basis: b' = Q' b on the unit
  layout <- .bt_random_effect_lkj_layout(term)
  set.seed(3)
  u <- matrix(stats::rbeta(40L * length(layout$free_pairs), 1.5, 1.5), ncol = length(layout$free_pairs))
  R <- .bt_lkj_layout_u_to_R(u, layout)
  T <- diag(5)
  T[3:5, 3:5] <- t(Q)
  R_rotated <- R
  for(draw in seq_len(nrow(u))){
    R_rotated[draw, , ] <- T %*% R[draw, , ] %*% t(T)
  }
  expect_equal(
    .bt_random_effect_correlation_row_values(rotated_rows, R_rotated),
    .bt_random_effect_correlation_row_values(rows, R),
    tolerance = 1e-12
  )
})

test_that("level-cell correlations and cor_total have their prior laws", {

  # 1 + x + f under us(): K = 5 with a unit of k = 3 columns
  term <- .unit_term(.unit_formula(~ 1 + x + f + us(1 + x + f | id), c("x", "f")))
  layout <- .bt_random_effect_lkj_layout(term)
  eta <- 1.5
  K <- 5L
  k <- 3L
  alpha <- .bt_lkj_layout_alpha(layout, eta)
  set.seed(4)
  n <- 20000L
  u <- vapply(alpha, function(shape) stats::rbeta(n, shape, shape), numeric(n))
  R <- .bt_lkj_layout_u_to_R(u, layout)
  rows <- .bt_random_effect_correlation_rows(term)
  values <- .bt_random_effect_correlation_row_values(rows, R)
  labels <- vapply(rows, `[[`, character(1), "label")
  shape <- eta - 1 + K / 2
  for(label in c("intercept,x", "intercept,f[a]", "x,f[c]")){
    expect_gt(
      suppressWarnings(stats::ks.test((values[, labels == label] + 1) / 2, "pbeta", shape, shape)$p.value),
      1e-3
    )
  }
  for(label in c("intercept,f", "x,f")){
    expect_gt(
      suppressWarnings(stats::ks.test(values[, labels == label]^2, "pbeta", k / 2, eta + (K - k - 1) / 2)$p.value),
      1e-3
    )
  }
})

test_that("a shared SD of a correlated block is certified for a common rescaling", {

  leaves <- c("i", "f", "f", "f")
  # the unit's columns rescaled by one factor
  M <- diag(c(1, 0.5, 0.5, 0.5))
  expect_true(.random_sd_shared_leaf_is_invariant(M, leaves, 2:4, correlated = TRUE))
  # different factors, or a column that mixes fitted columns
  M_unequal <- diag(c(1, 0.5, 0.5, 0.25))
  expect_false(.random_sd_shared_leaf_is_invariant(M_unequal, leaves, 2:4, correlated = TRUE))
  M_mixed <- M
  M_mixed[3L, 1L] <- 0.2
  expect_false(.random_sd_shared_leaf_is_invariant(M_mixed, leaves, 2:4, correlated = TRUE))
})

test_that("a block that is one unit is one SD component of a block allocation", {

  allocation <- prior_random(random_variance_allocation(
    name = "tot", terms = c(id = "id", h = "h"),
    sd = prior("gamma", list(2, 2)),
    weights = prior("dirichlet", list(alpha = c(1, 1)))
  ))
  result <- .unit_formula(~ 1 + f + diag(0 + f | id) + (1 | h), "f", prior_random = allocation)
  term <- .unit_term(result)
  expect_identical(term$sd_parameter_names, rep("mu__xREx__id_f[1]", 3L))
  expect_match(result$formula_syntax, "mu__xREx__id_f[1] = mu__xRE_ALLOCx_tot__allocation_sd", fixed = TRUE)
  expect_error(
    .unit_formula(~ 1 + f + diag(1 + f | id) + (1 | h), "f", prior_random = allocation),
    "Block-level variance allocation supports random-effect blocks with one SD component.",
    fixed = TRUE
  )
})

test_that("the contrast coordinates of a unit are refused as selectors", {

  result <- .unit_formula(~ 1 + f + us(1 + f | id), "f", prior_random = prior_random(id = random_block(
    sd = .unit_sd_prior,
    monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
  )))
  term <- .unit_term(result)
  layout <- .bt_random_effect_lkj_layout(term)
  u <- matrix(c(0.3, 0.6, 0.8), nrow = 1L)
  R <- .bt_lkj_layout_u_to_R(u, layout)[1L, , ]
  L <- .bt_lkj_layout_u_to_L(u, layout)[1L, , ]
  cells <- function(name){
    outer(1:4, 1:4, function(row, column) paste0(term$parameter_stem, name, "[", row, ",", column, "]"))
  }
  values <- c(
    mu_intercept = 0, stats::setNames(rep(0, 3L), paste0("mu_f[", 1:3, "]")),
    stats::setNames(c(0.5, 1), unique(term$sd_parameter_names)),
    stats::setNames(as.vector(R), as.vector(cells("_xRE_CORx_R"))),
    stats::setNames(as.vector(L), as.vector(cells("_xRE_CORx_L"))),
    stats::setNames(as.vector(u), term$correlation$primitive_names)
  )
  posterior <- matrix(rep(values, 2L), nrow = 2L, byrow = TRUE, dimnames = list(NULL, names(values)))
  fit <- list(mcmc = coda::mcmc.list(coda::mcmc(posterior)), summary.pars = list(mutate = NULL),
              monitor = colnames(posterior), sample = nrow(posterior))
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "formula_scale") <- list(mu = result$formula_scale)
  fit <- attach_test_parameter_map(fit)
  catalog <- parameter_catalog(fit)

  expect_true(all(c("(mu) sd(f)", "(mu) cor(intercept,f[a])", "(mu) cor_total(intercept,f)") %in%
                    catalog$quantities$canonical_name))
  refusal <- expect_error(
    parameter_catalog_resolve(catalog, "(mu) sd(f{2})"),
    class = "BayesTools_selector_unavailable"
  )
  expect_identical(
    conditionMessage(refusal),
    paste0(
      "Selector '(mu) sd(f{2})' is unavailable: the columns of factor term 'f' are one ",
      "exchangeable unit, whose contrast coordinates have no own SD or correlation. ",
      "Select '(mu) sd(f)'."
    )
  )
  refusal <- expect_error(
    parameter_catalog_resolve(catalog, "(mu) cor(f{1},intercept)"),
    class = "BayesTools_selector_unavailable"
  )
  expect_match(conditionMessage(refusal), "Select '(mu) cor_total(intercept,f)'.", fixed = TRUE)
  # the unit's own quantities resolve, and their draws are the reported values
  draws <- parameter_draws(fit, parameter_catalog_resolve(catalog, "(mu) cor(intercept,f[b])"))
  C <- contr.meandif(c("a", "b", "c", "d"))
  expect_equal(
    as.numeric(as.matrix(draws)),
    rep(sum(C[2L, ] * R[2:4, 1L]) / sqrt(sum(C[2L, ]^2)), 2L),
    tolerance = 1e-12
  )
  draws <- parameter_draws(fit, parameter_catalog_resolve(catalog, "(mu) cor_total(intercept,f)"))
  expect_equal(as.numeric(as.matrix(draws)), rep(sqrt(sum(R[2:4, 1L]^2)), 2L), tolerance = 1e-12)
})

test_that("the marginal covariance of a unit block is that of its restricted correlation", {

  data <- .unit_data()[1:24, ]
  data$id <- droplevels(data$id)
  result <- .unit_formula(~ 1 + f + us(1 + f | id), "f", data = data, prior_random = prior_random(
    id = random_block(
      sd = prior("point", list(location = 0.4)),
      terms = list(f = prior("point", list(location = 0.7))),
      cor = prior_lkj(eta = 2, include_primitives = TRUE)
    )
  ))
  term <- .unit_term(result)
  u <- c(0.3, 0.65, 0.8)
  posterior <- matrix(u, nrow = 1L, dimnames = list(NULL, term$correlation$primitive_names))
  posterior <- .bt_random_effect_marginal_covariance_validate_posterior(posterior, allow_zero_columns = TRUE)
  covariance <- .bt_random_effect_marginal_covariance_samples(
    design = result$formula_design,
    posterior = posterior,
    prior_list = result$prior_list,
    blocks = "id"
  )$samples[1L, , ]

  # independent: the vine correlation of the full primitives with the unit
  # (columns 2:4) first, its within-unit pairs at 1/2, read in design order
  vine <- .bt_lkj_cholesky_cpc_u_to_L(c(0.5, 0.5, 0.5, u), K = 4)
  vine_R <- tcrossprod(vine)
  vine_order <- c(2L, 3L, 4L, 1L)
  R <- matrix(NA_real_, 4L, 4L)
  R[vine_order, vine_order] <- vine_R
  expect_equal(R[2:4, 2:4], diag(3), tolerance = 1e-15)
  S <- diag(c(0.4, 0.7, 0.7, 0.7))
  Z <- term$model_matrix
  same_group <- outer(data$id, data$id, "==")
  expected <- (Z %*% S %*% R %*% S %*% t(Z)) * same_group
  expect_equal(unname(covariance), unname(expected), tolerance = 1e-12)
})

# A carrier fit of the compiled formula `result` whose draws are the fixed
# SDs `sd` (by SD leaf) and the restricted correlation of the free primitives `u`.
.unit_carrier <- function(result, sd, u){
  term <- .unit_term(result)
  layout <- .bt_random_effect_lkj_layout(term)
  u <- matrix(u, nrow = 1L)
  R <- .bt_lkj_layout_u_to_R(u, layout)[1L, , ]
  L <- .bt_lkj_layout_u_to_L(u, layout)[1L, , ]
  K <- term$n_columns
  cells <- function(name){
    outer(seq_len(K), seq_len(K), function(row, column) paste0(term$parameter_stem, name, "[", row, ",", column, "]"))
  }
  fixed <- unlist(lapply(names(result$prior_list), function(name){
    if(!startsWith(name, "mu_") || startsWith(name, "mu__xREx__")) return(NULL)
    n <- attr(result$prior_list[[name]], "levels")
    if(is.null(n) || n == 1L) name else paste0(name, "[", seq_len(n - 1L), "]")
  }))
  values <- c(
    stats::setNames(rep(0, length(fixed)), fixed),
    stats::setNames(sd, unique(term$sd_parameter_names)),
    stats::setNames(as.vector(R), as.vector(cells("_xRE_CORx_R"))),
    stats::setNames(as.vector(L), as.vector(cells("_xRE_CORx_L"))),
    stats::setNames(as.vector(u), term$correlation$primitive_names)
  )
  posterior <- matrix(rep(values, 2L), nrow = 2L, byrow = TRUE, dimnames = list(NULL, names(values)))
  fit <- list(mcmc = coda::mcmc.list(coda::mcmc(posterior)), summary.pars = list(mutate = NULL),
              monitor = colnames(posterior), sample = nrow(posterior))
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "formula_scale") <- list(mu = result$formula_scale)
  list(fit = attach_test_parameter_map(fit), R = R, term = term)
}

test_that("original-scale level cells of a unit need its certified SD; cor_total does not", {

  skip_if_not_installed("runjags")

  data <- .unit_data()
  m <- mean(data$x)
  s <- stats::sd(data$x)
  levels <- c("a", "b", "c", "d")
  C <- contr.meandif(levels)
  # the original-scale map of the fitted coefficients of the block, from the
  # design on a grid of level cells
  unscale_map <- function(rhs){
    grid <- expand.grid(f = factor(levels, levels = levels), x = c(-1, 0, 1, 2))
    stats::contrasts(grid$f) <- C
    scaled <- grid
    scaled$x <- (grid$x - m) / s
    qr.solve(stats::model.matrix(rhs, grid), stats::model.matrix(rhs, scaled))
  }
  table_samples <- function(fit){
    JAGS_estimates_table(fit, transform_scaled = TRUE, random_effects_summary = "standard",
                         remove_diagnostics = TRUE, return_samples = TRUE)
  }
  scaled_formula <- function(formula, fixed){
    priors <- list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      f = prior_factor("mnormal", list(0, 1), contrast = "meandif")
    )
    JAGS_formula(formula, "mu", data, prior_list = priors[c("intercept", fixed)],
                 formula_scale = list(x = TRUE),
                 prior_random = prior_random(id = random_block(
                   sd = .unit_sd_prior,
                   monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
                 )))
  }

  # x:f is divided by sd(x): a common rescaling certifies the unit
  carrier <- .unit_carrier(scaled_formula(~ 1 + x + f + us(1 + x + x:f | id), c("x", "f")),
                           sd = c(0.4, 0.6, 0.9), u = c(0.4, 0.55, 0.7, 0.35, 0.6, 0.45, 0.65))
  A <- unscale_map(~ x + x:f)
  S <- diag(c(0.4, 0.6, 0.9, 0.9, 0.9))
  covariance <- A %*% S %*% carrier$R %*% S %*% t(A)
  samples <- table_samples(carrier$fit)
  for(a in seq_along(levels)){
    w <- c(0, 0, C[a, ])
    expected <- sum(w * covariance[, 2L]) / sqrt(sum(w * (covariance %*% w)) * covariance[2L, 2L])
    expect_equal(unname(samples[1L, paste0("(mu) cor(x,x:f[", levels[a], "])")]), expected,
                 tolerance = 1e-10, info = levels[a])
  }
  expect_equal(unname(samples[1L, "(mu) sd(x:f)"]), 0.9 / s, tolerance = 1e-12)

  # f next to the indicator-coded f:x: centering shifts the levels of f by
  # slopes with their own SDs, so the original-scale unit has no one SD
  carrier <- .unit_carrier(scaled_formula(~ 1 + x + f + us(1 + f + f:x | id), c("x", "f")),
                           sd = c(0.4, 0.6, 0.5, 0.7, 0.8, 0.9),
                           u = seq(0.3, 0.7, length.out = 25L))
  A <- unscale_map(~ f + f:x)
  S <- diag(c(0.4, 0.6, 0.6, 0.6, 0.5, 0.7, 0.8, 0.9))
  covariance <- A %*% S %*% carrier$R %*% S %*% t(A)
  correlation <- stats::cov2cor(covariance)
  table <- JAGS_estimates_table(carrier$fit, transform_scaled = TRUE, remove_diagnostics = TRUE)
  values <- stats::setNames(as.data.frame(table)[["Mean"]], rownames(table))
  level_rows <- paste0("(mu) cor(intercept,f[", levels, "])")
  expect_true(all(is.na(values[c("(mu) sd(f)", level_rows)])))
  notes <- attr(table, "footnotes")
  expect_true(all(level_rows %in% names(notes)))
  expect_match(notes[[level_rows[1L]]], "its shared design columns cannot be certified to have a single SD", fixed = TRUE)
  r <- correlation[2:4, 1L]
  expect_equal(unname(values[["(mu) cor_total(intercept,f)"]]),
               sqrt(sum(r * solve(correlation[2:4, 2:4], r))), tolerance = 1e-10)
})

test_that("the prior densities of level cells and cor_total follow their laws", {

  eta <- 2
  correlated <- prior_random(id = random_block(
    sd = .unit_sd_prior,
    cor = prior_lkj(eta = eta),
    monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
  ))
  log_density <- function(fit, name, value){
    density <- parameter_prior_density(fit, parameter_catalog_resolve(parameter_catalog(fit), name))
    prior_density_ordinate(density, value)$log_density
  }

  # K = 5 with a unit of k = 3 columns
  carrier <- .unit_carrier(.unit_formula(~ 1 + x + f + us(1 + x + f | id), c("x", "f"), prior_random = correlated),
                           sd = c(0.4, 0.6, 0.9), u = c(0.4, 0.55, 0.7, 0.35, 0.6, 0.45, 0.65))
  K <- 5L
  k <- 3L
  b <- eta + (K - k - 1) / 2
  for(t in c(0.2, 0.5, 0.8)){
    expect_equal(
      log_density(carrier$fit, "(mu) cor_total(intercept,f)", t),
      log(2) + (k - 1) * log(t) + (b - 1) * log(1 - t^2) - lbeta(k / 2, b),
      tolerance = 1e-8, info = t
    )
  }
  expect_identical(log_density(carrier$fit, "(mu) cor_total(intercept,f)", 0), -Inf)
  shape <- eta - 1 + K / 2
  for(r in c(-0.5, 0, 0.3)){
    expect_equal(
      log_density(carrier$fit, "(mu) cor(x,f[c])", r),
      stats::dbeta((r + 1) / 2, shape, shape, log = TRUE) - log(2),
      tolerance = 1e-8, info = r
    )
  }

  # K = 2 with a two-level factor: the multiple correlation is |r|, with a
  # finite ordinate at zero, and a level cell is +-(2u - 1)
  carrier <- .unit_carrier(.unit_formula(~ 1 + f2 + us(1 + f2 | id), "f2", prior_random = correlated),
                           sd = c(0.4, 0.6), u = 0.3)
  b <- eta
  expect_equal(log_density(carrier$fit, "(mu) cor_total(intercept,f2)", 0), log(2) - lbeta(1 / 2, b),
               tolerance = 1e-8)
  expect_equal(log_density(carrier$fit, "(mu) cor_total(intercept,f2)", 0.4),
               log(2) + (b - 1) * log(1 - 0.4^2) - lbeta(1 / 2, b), tolerance = 1e-8)
  C <- contr.meandif(c("u", "v"))
  catalog <- parameter_catalog(carrier$fit)
  for(level in 1:2){
    selection <- parameter_catalog_resolve(catalog, paste0("(mu) cor(intercept,f2[", c("u", "v")[level], "])"))
    expect_identical(selection$quantities$source_type, "one_to_one_transform")
    transform <- parameter_transform(carrier$fit, selection)
    expect_equal(parameter_transform_forward(0.3, transform), sign(C[level, 1L]) * (2 * 0.3 - 1),
                 tolerance = 1e-14)
    draws <- parameter_draws(carrier$fit, selection)
    expect_equal(as.numeric(as.matrix(draws)), rep(sign(C[level, 1L]) * (2 * 0.3 - 1), 2L),
                 tolerance = 1e-14)
  }
})

test_that("a us() block that is one unit is its diagonal block downstream", {

  result <- .unit_formula(~ 1 + f + us(0 + f | id), "f")
  term <- .unit_term(result)
  expect_identical(term$structure, "diag")
  expect_null(term$correlation)
  values <- c(mu_intercept = 0, stats::setNames(rep(0, 3L), paste0("mu_f[", 1:3, "]")), mu__xREx__id_f = 0.7)
  posterior <- matrix(rep(values, 2L), nrow = 2L, byrow = TRUE, dimnames = list(NULL, names(values)))
  fit <- list(mcmc = coda::mcmc.list(coda::mcmc(posterior)), summary.pars = list(mutate = NULL),
              monitor = colnames(posterior), sample = nrow(posterior))
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "formula_scale") <- list(mu = result$formula_scale)
  fit <- attach_test_parameter_map(fit)
  quantities <- parameter_catalog(fit)$quantities
  random <- quantities$canonical_name[quantities$owner_type == "random_block"]
  expect_identical(random, c("(mu) sd(f)", "(mu) var(f)"))
  table <- JAGS_estimates_table(fit, remove_diagnostics = TRUE)
  expect_equal(as.data.frame(table)["(mu) sd(f)", "Mean"], 0.7)
  # its SD prior is tagged with the stored structure, so structure filters
  # keep it in every table
  expect_identical(attr(result$prior_list[["mu__xREx__id_f"]], "random_structure"), "diag")
  raw <- JAGS_estimates_table(fit, remove_diagnostics = TRUE, random_effects_summary = "raw",
                              keep_random_structures = "diag")
  expect_true(any(grepl("sd(f", rownames(raw), fixed = TRUE)))
  # the declared us() formula matches the fitted block in predictions
  expect_identical(term$formula_structure, "us")
  expect_no_error(.bt_validate_random_effect_prediction_terms(
    requested = .bt_parse_random_effects(~ 1 + f + us(0 + f | id))$terms,
    fitted = .bt_formula_design_random_effects(result$formula_design),
    parameter = "mu"
  ))
  expect_error(
    .bt_validate_random_effect_prediction_terms(
      requested = .bt_parse_random_effects(~ 1 + f + id(0 + f | id))$terms,
      fitted = .bt_formula_design_random_effects(result$formula_design),
      parameter = "mu"
    ),
    "does not match the fitted formula", fixed = TRUE
  )
})

test_that("a multiple correlation of a singular unit block is undefined, not an error", {

  term <- .unit_term(.unit_formula(~ 1 + x + f + us(1 + x + f | id), c("x", "f")))
  rows <- .bt_random_effect_correlation_rows(term)
  # an original-scale correlation of rank one within the unit (columns 3:5)
  v <- c(1, -1, 0) / sqrt(2)
  R <- diag(5)
  R[3:5, 3:5] <- tcrossprod(c(1, 1, 1)) / 1
  R[1, 2] <- R[2, 1] <- 0.3
  R[1, 3:5] <- R[3:5, 1] <- 0.2
  R <- array(R, dim = c(1L, 5L, 5L))
  values <- .bt_random_effect_correlation_row_values(rows, R)
  labels <- vapply(rows, `[[`, character(1), "label")
  expect_true(is.na(values[1L, labels == "intercept,f"]))
  expect_equal(values[1L, labels == "intercept,x"], 0.3)

  # an ill-conditioned but nonsingular unit block keeps its value:
  # R_ff = [[1, 1 - e], [1 - e, 1]] and r = (1/2, 1/2) give r' R_ff^-1 r = 1 / (2 (2 - e))
  term <- .unit_term(.unit_formula(~ 1 + g + us(1 + g | id), "g"))
  rows <- .bt_random_effect_correlation_rows(term)
  e <- 1e-9
  R <- matrix(c(1, 0.5, 0.5, 0.5, 1, 1 - e, 0.5, 1 - e, 1), 3L)
  values <- .bt_random_effect_correlation_row_values(rows, array(R, dim = c(1L, 3L, 3L)))
  labels <- vapply(rows, `[[`, character(1), "label")
  expect_equal(values[1L, labels == "intercept,g"], sqrt(1 / (2 * (2 - e))), tolerance = 1e-6)

  # a nearly singular block with correlations at the boundary (exactly about
  # 0.999999) may round above one; it never leaves the support tolerance
  # of posterior summaries, and a fitted-scale identity block at the boundary
  # stays finite
  term <- .unit_term(.unit_formula(~ 1 + f + us(1 + f | id), "f"))
  rows <- .bt_random_effect_correlation_rows(term)
  labels <- vapply(rows, `[[`, character(1), "label")
  b <- 0.7071067811861939
  R <- diag(4)
  R[2:4, 2:4] <- matrix(c(1, 0, b, 0, 1, b, b, b, 1), 3L)
  R[1, 2:4] <- R[2:4, 1] <- c(0, 0, 1.0001620857883492e-6)
  value <- .bt_random_effect_correlation_row_values(rows, array(R, dim = c(1L, 4L, 4L)))[1L, labels == "intercept,f"]
  expect_true(is.na(value) || (value >= 0 && value <= 1 + sqrt(.Machine$double.eps)))
  R <- diag(4)
  R[1, 2:4] <- R[2:4, 1] <- c(0.6, 0.8, 0)
  value <- .bt_random_effect_correlation_row_values(rows, array(R, dim = c(1L, 4L, 4L)))[1L, labels == "intercept,f"]
  expect_equal(value, 1, tolerance = 1e-15)
})

test_that("convergence roles mark the constant cells of a unit block as structural", {

  roles <- function(result){
    K <- .unit_term(result)$n_columns
    cells <- expand.grid(row = seq_len(K), column = seq_len(K))
    names <- function(matrix_name){
      paste0("mu__xREx__id_xRE_CORx_", matrix_name, "[", cells$row, ",", cells$column, "]")
    }
    # the three fixed within-unit primitives lead the vine; the fourth is
    # free (all four monitored explicitly, by their base name)
    primitives <- paste0("mu__xREx__id_xRE_CORx_lkj_u[", 1:4, "]")
    coordinates <- .bt_build_parameter_coordinates(
      columns = c(names("L"), names("R"), primitives),
      prior_list = result$prior_list,
      formula_design = list(mu = result$formula_design),
      add_parameters = c(result$add_parameters, "mu__xREx__id_xRE_CORx_lkj_u"),
      model_syntax = paste0("model{\n", result$formula_syntax, "}"),
      data_names = names(result$data)
    )
    role <- function(matrix_name){
      matrix(coordinates$convergence_role[match(names(matrix_name), coordinates$coordinate_name)], K)
    }
    list(
      L = role("L"), R = role("R"),
      u = coordinates$convergence_role[match(primitives, coordinates$coordinate_name)]
    )
  }
  # derived cells, with the zero upper Cholesky triangle and the unit diagonal
  # and within-unit cells of the correlation matrix structural
  expected <- function(K, unit){
    L <- matrix("derived", K, K)
    L[upper.tri(L)] <- "structural"
    R <- matrix("derived", K, K)
    diag(R) <- "structural"
    R[unit, unit] <- "structural"
    list(L = L, R = R, u = c(rep("structural", 3L), "sampled"))
  }
  monitor <- prior_random(id = random_block(
    sd = .unit_sd_prior,
    monitor = random_monitor(correlation = TRUE)
  ))

  # a unit leading the design keeps the vine order: its Cholesky block is the
  # identity, so constant
  leading <- .unit_formula(~ 1 + f + x + us(0 + f + x | id), c("f", "x"), prior_random = monitor)
  expect_identical(.unit_term(leading)$correlation$unit_columns, 1:3)
  out <- expected(4L, 1:3)
  out$L[1:3, 1:3] <- "structural"
  expect_identical(roles(leading), out)

  # a reordered vine: the first Cholesky diagonal is a row norm (derived)
  reordered <- .unit_formula(~ 1 + x + f + us(1 + x + f | id), c("x", "f"), prior_random = monitor)
  expect_identical(.unit_term(reordered)$correlation$unit_columns, 3:5)
  expect_identical(roles(reordered), expected(5L, 3:5))
})

test_that("original-scale level cells need only the unit's own certified SD", {

  skip_if_not_installed("runjags")

  set.seed(4)
  data <- data.frame(
    x    = stats::rnorm(96),
    z    = stats::rnorm(96, mean = 3, sd = 2),
    flag = rep(c(TRUE, FALSE), 48L),
    f    = factor(rep(c("a", "b", "c", "d"), 24)),
    id   = factor(rep(1:12, each = 8))
  )
  levels <- c("a", "b", "c", "d")
  C <- contr.meandif(levels)
  result <- JAGS_formula(
    ~ 1 + f + us(1 + f + x:flag + x:z | id), "mu", data,
    prior_list = list(intercept = prior("normal", list(0, 1)),
                      f = prior_factor("mnormal", list(0, 1), contrast = "meandif")),
    formula_scale = list(z = TRUE),
    prior_random = prior_random(id = random_block(
      sd = .unit_sd_prior,
      monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
    ))
  )
  term <- .unit_term(result)
  expect_identical(term$sd_leaves$units[[1L]]$columns, 2:4)
  carrier <- .unit_carrier(result, sd = c(0.4, 0.6, 0.5, 0.8),
                           u = seq(0.35, 0.65, length.out = length(term$correlation$primitive_names)))
  # the shared SD of the x:flag columns is not certified (centering z shifts
  # both by the x:z slope); the unit f is
  table <- JAGS_estimates_table(carrier$fit, transform_scaled = TRUE, remove_diagnostics = TRUE)
  values <- stats::setNames(as.data.frame(table)[["Mean"]], rownames(table))
  shared <- grep("^[(]mu[)] sd[(]flag", names(values), value = TRUE)
  expect_true(all(is.na(values[shared])))
  level_rows <- grep("^[(]mu[)] cor[(]x:flag[[]TRUE[]]:?[^,]*,f[[]", names(values), value = TRUE)
  level_rows <- grep(",f[", grep("cor(", names(values), fixed = TRUE, value = TRUE), fixed = TRUE, value = TRUE)
  level_rows <- level_rows[grepl("flag", level_rows, fixed = TRUE)]
  expect_length(level_rows, 8L)
  expect_false(anyNA(values[level_rows]))

  # independent original-scale map of the block from its design on a grid
  grid <- expand.grid(f = factor(levels, levels = levels), x = c(-1, 0, 1, 2),
                      z = c(-1, 1, 4), flag = c(TRUE, FALSE))
  stats::contrasts(grid$f) <- C
  scaled <- grid
  scaled$z <- (grid$z - mean(data$z)) / stats::sd(data$z)
  rhs <- ~ f + x:flag + x:z
  Z <- stats::model.matrix(rhs, grid)
  A <- qr.solve(Z, stats::model.matrix(rhs, scaled))
  S <- diag(c(0.4, 0.6, 0.6, 0.6, 0.5, 0.5, 0.8))
  covariance <- A %*% S %*% carrier$R %*% S %*% t(A)
  o <- which(colnames(Z) == "x:flagTRUE")
  for(a in seq_along(levels)){
    w <- c(0, C[a, ], 0, 0, 0)
    expected <- sum(w * covariance[, o]) / sqrt(sum(w * (covariance %*% w)) * covariance[o, o])
    label <- level_rows[grepl("TRUE", level_rows, fixed = TRUE) & endsWith(level_rows, paste0("f[", levels[a], "])"))]
    expect_length(label, 1L)
    expect_equal(unname(values[[label]]), expected, tolerance = 1e-10, info = levels[a])
  }
})

test_that("tables count the defined draws of cor_total", {

  result <- .unit_formula(~ 1 + f + us(1 + f | id), "f", prior_random = prior_random(id = random_block(
    sd = .unit_sd_prior,
    monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
  )))
  term <- .unit_term(result)
  layout <- .bt_random_effect_lkj_layout(term)
  u <- matrix(c(0.3, 0.6, 0.8), nrow = 1L)
  R <- .bt_lkj_layout_u_to_R(u, layout)[1L, , ]
  L <- .bt_lkj_layout_u_to_L(u, layout)[1L, , ]
  cells <- function(name){
    outer(1:4, 1:4, function(row, column) paste0(term$parameter_stem, name, "[", row, ",", column, "]"))
  }
  draw <- function(sd_f){
    c(mu_intercept = 0, stats::setNames(rep(0, 3L), paste0("mu_f[", 1:3, "]")),
      stats::setNames(c(0.5, sd_f), unique(term$sd_parameter_names)),
      stats::setNames(as.vector(R), as.vector(cells("_xRE_CORx_R"))),
      stats::setNames(as.vector(L), as.vector(cells("_xRE_CORx_L"))),
      stats::setNames(as.vector(u), term$correlation$primitive_names))
  }
  posterior <- rbind(draw(0.7), draw(0), draw(0.9), draw(0.4))
  fit <- list(mcmc = coda::mcmc.list(coda::mcmc(posterior)), summary.pars = list(mutate = NULL),
              monitor = colnames(posterior), sample = nrow(posterior))
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "formula_scale") <- list(mu = result$formula_scale)
  fit <- attach_test_parameter_map(fit)
  table <- JAGS_estimates_table(fit, remove_diagnostics = TRUE)
  footnotes <- attr(table, "footnotes")
  expect_true(all(c("(mu) cor_total(intercept,f)", "(mu) cor(intercept,f[a])") %in% names(footnotes)))
  expect_match(footnotes[["(mu) cor_total(intercept,f)"]], "3 of 4 draws", fixed = TRUE)
  expect_match(footnotes[["(mu) cor_total(intercept,f)"]],
               "the correlation matrix of the factor term is not singular", fixed = TRUE)
  expect_match(footnotes[["(mu) cor(intercept,f[a])"]], "both SDs are positive", fixed = TRUE)
  quantities <- parameter_catalog(fit)$quantities
  expect_identical(
    quantities$definedness[quantities$canonical_name == "(mu) cor_total(intercept,f)"],
    "multiple_correlation"
  )
})
