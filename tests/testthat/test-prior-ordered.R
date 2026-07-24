skip_if_not_test_profile("unit")

test_that("prior_ordered() validates constructor inputs", {
  p <- prior_ordered(prior("normal", list(0, 1)))

  expect_true(is.prior(p))
  expect_true(is.prior.factor(p))
  expect_true(is.prior.ordered(p))
  expect_equal(p$contrast, "cumulative")
  expect_equal(p$allocation$type, "default_dirichlet")

  expect_error(
    prior_ordered(prior_factor("normal", list(0, 1), contrast = "treatment")),
    "scalar prior"
  )
  expect_error(
    prior_ordered(prior("normal", list(0, 1)), allocation = c(.2, .2)),
    "sum to one"
  )
  expect_error(
    prior_ordered(prior("normal", list(0, 1)), id = ""),
    "cannot contain empty strings"
  )
  expect_error(
    prior_ordered(prior("normal", list(0, 1)), id = "  "),
    "cannot contain empty strings"
  )
  expect_error(
    prior_ordered(prior("normal", list(0, 1)), allocation = prior("normal", list(0, 1))),
    "Dirichlet"
  )
  expect_error(
    prior_spike_and_slab(p),
    "inside prior_ordered"
  )
})

test_that("ordered cumulative contrasts encode level effects", {
  expect_equal(
    contr.ordered_cumulative(c("low", "mid", "high")),
    matrix(c(0, 0, 1, 0, 1, 1), nrow = 3, byrow = TRUE)
  )
  expect_equal(
    contr.ordered_cumulative_levels(c("low", "mid", "high")),
    matrix(c(1, 0, 0, 1, 1, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
  )
})

test_that("formula binding converts explicit ordered priors with current factor rules", {
  df <- data.frame(
    y  = seq_len(6),
    f  = factor(rep(c("mid", "low", "high"), 2), levels = c("mid", "low", "high")),
    ch = rep(c("b", "a", "c"), 2)
  )

  formula_info <- JAGS_formula(
    y ~ f + ch,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f  = prior_ordered(prior("normal", list(0, 1)), allocation = c(.2, .8)),
      ch = prior_ordered(prior("normal", list(0, 1)), allocation = c(.4, .6))
    )
  )

  expect_equal(formula_info$formula_design$xlevels$f, c("mid", "low", "high"))
  expect_equal(formula_info$formula_design$xlevels$ch, c("a", "b", "c"))
  expect_equal(formula_info$formula_design$contrasts$f, "contr.ordered_cumulative")
  expect_equal(formula_info$formula_design$contrasts$ch, "contr.ordered_cumulative")
  expect_equal(unname(formula_info$data$mu_data_f[1:3, ]), matrix(c(0, 0, 1, 0, 1, 1), nrow = 3, byrow = TRUE))

  expect_error(
    JAGS_formula(
      y ~ f,
      "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        f = prior_ordered(prior("normal", list(0, 1)), allocation = c(.2, .3, .5))
      )
    ),
    "has length 3, but 2 value"
  )
})

test_that("cumulative_levels uses one allocation piece per level", {
  p <- prior_ordered(
    prior("point", list(location = 10)),
    allocation = c(.2, .3, .5),
    contrast = "cumulative_levels"
  )
  attr(p, "levels") <- 3
  attr(p, "level_names") <- c("low", "mid", "high")

  samples <- rng(p, 2)
  expect_equal(unname(samples[1, ]), c(2, 5, 10))

  df <- data.frame(
    y = seq_len(6),
    f = factor(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )
  formula_info <- JAGS_formula(
    y ~ 0 + f,
    "mu",
    data = df,
    prior_list = list(
      f = prior_ordered(
        prior("normal", list(0, 1)),
        allocation = c(.2, .3, .5),
        contrast = "cumulative_levels"
      )
    )
  )

  expect_equal(unname(formula_info$data$mu_data_f[1:3, ]), contr.ordered_cumulative_levels(1:3))
})

test_that("JAGS syntax, inits, and monitors use latent total and allocation nodes", {
  df <- data.frame(
    y = seq_len(6),
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )
  formula_fixed <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), allocation = c(.25, .75))
    )
  )
  syntax_fixed <- JAGS_add_priors("model{}", formula_fixed$prior_list)

  expect_match(syntax_fixed, "mu_f_ordered_total ~ dnorm\\(0,1\\)")
  expect_match(syntax_fixed, "mu_f\\[1\\] <- mu_f_ordered_total \\* 0.25")
  expect_false(grepl("prior_par_eta_mu_f_ordered_alloc", syntax_fixed))

  formula_dirichlet <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)))
    )
  )
  syntax_dirichlet <- JAGS_add_priors("model{}", formula_dirichlet$prior_list)
  expect_match(syntax_dirichlet, "prior_par_eta_mu_f_ordered_alloc_f_1\\[1\\] ~ dgamma\\(1, 1\\)")
  expect_match(syntax_dirichlet, "mu_f\\[2\\] <- mu_f_ordered_total \\* mu_f_ordered_alloc_f_1\\[2\\]")

  monitors <- JAGS_to_monitor(formula_dirichlet$prior_list)
  expect_true(all(c("mu_f", "mu_f_ordered_total", "prior_par_eta_mu_f_ordered_alloc_f_1") %in% monitors))

  inits <- JAGS_get_inits(formula_dirichlet$prior_list, chains = 1, seed = 1)[[1]]
  expect_true("mu_f_ordered_total" %in% names(inits))
  expect_true("prior_par_eta_mu_f_ordered_alloc_f_1" %in% names(inits))
  expect_equal(length(inits$prior_par_eta_mu_f_ordered_alloc_f_1), 2L)

  formula_spike_slab <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(
        prior_spike_and_slab(prior("normal", list(0, 1))),
        allocation = c(.25, .75)
      )
    )
  )
  syntax_spike_slab <- JAGS_add_priors("model{}", formula_spike_slab$prior_list)
  expect_match(syntax_spike_slab, "mu_f_ordered_total_indicator ~ dbern")
  expect_match(syntax_spike_slab, "mu_f_ordered_total = mu_f_ordered_total_variable \\* mu_f_ordered_total_indicator")
})

test_that("ordered interactions expand through the formula binder", {
  df <- data.frame(
    y = seq_len(12),
    x = rep(c(-1, 1), 6),
    f = ordered(rep(c("low", "mid", "high"), 4), levels = c("low", "mid", "high")),
    g = factor(rep(c("A", "B"), each = 6))
  )

  formula_info <- JAGS_formula(
    y ~ f + x + g + f:x + f:g,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1))),
      x = prior("normal", list(0, 1)),
      g = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "f:x" = prior_ordered(prior("normal", list(0, 1))),
      "f:g" = prior_ordered(
        prior("normal", list(0, 1)),
        allocation = list(f = prior("dirichlet", list(alpha = c(2, 3))))
      )
    )
  )

  expect_equal(dim(attr(formula_info$prior_list$mu_f__xXx__x, "factor_design")), c(3L, 2L))
  expect_equal(dim(attr(formula_info$prior_list$mu_f__xXx__g, "factor_design")), c(6L, 2L))

  syntax <- JAGS_add_priors("model{}", formula_info$prior_list)
  expect_match(syntax, "mu_f__xXx__x\\[1\\] <- mu_f__xXx__x_ordered_total")
  expect_match(syntax, "prior_par_eta_mu_f__xXx__g_ordered_alloc_f_1\\[2\\] ~ dgamma\\(3, 1\\)")

  expect_error(
    JAGS_formula(
      y ~ f:g,
      "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "f:g" = prior_ordered(prior("normal", list(0, 1)))
      )
    ),
    "non-hierarchical factor interaction"
  )
})

test_that("ordered allocation id sharing emits one shared allocation", {
  df <- data.frame(
    y = seq_len(6),
    x = rep(c(-1, 1), 3),
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )

  formula_info <- JAGS_formula(
    y ~ f + x + f:x,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), id = "shape"),
      x = prior("normal", list(0, 1)),
      "f:x" = prior_ordered(prior("normal", list(0, 1)), id = "shape")
    )
  )

  syntax <- JAGS_add_priors("model{}", formula_info$prior_list)
  expect_equal(lengths(regmatches(syntax, gregexpr("prior_par_eta_ordered_alloc_shape_f[1] ~", syntax, fixed = TRUE))), 1L)
  expect_match(syntax, "mu_f__xXx__x\\[1\\] <- mu_f__xXx__x_ordered_total \\* ordered_alloc_shape_f\\[1\\]")

  monitors <- JAGS_to_monitor(formula_info$prior_list)
  expect_equal(sum(monitors == "prior_par_eta_ordered_alloc_shape_f"), 1L)

  formula_info_other <- JAGS_formula(
    y ~ f,
    "theta",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), id = "shape")
    )
  )
  combined_priors <- c(formula_info$prior_list, formula_info_other$prior_list)
  combined_syntax <- JAGS_add_priors("model{}", combined_priors)
  expect_equal(
    lengths(regmatches(
      combined_syntax,
      gregexpr("prior_par_eta_ordered_alloc_shape_f[1] ~", combined_syntax, fixed = TRUE)
    )),
    1L
  )

  expect_error(
    JAGS_formula(
      y ~ f + x + f:x,
      "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        f = prior_ordered(
          prior("normal", list(0, 1)),
          allocation = c(.25, .75),
          id = "shape"
        ),
        x = prior("normal", list(0, 1)),
        "f:x" = prior_ordered(
          prior("normal", list(0, 1)),
          allocation = c(.5, .5),
          id = "shape"
        )
      )
    ),
    "incompatible allocation specifications"
  )

  formula_info_incompatible <- JAGS_formula(
    y ~ f,
    "theta",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(
        prior("normal", list(0, 1)),
        allocation = prior("dirichlet", list(alpha = c(2, 3))),
        id = "shape"
      )
    )
  )
  incompatible_priors <- c(formula_info$prior_list, formula_info_incompatible$prior_list)
  expect_error(
    JAGS_add_priors("model{}", incompatible_priors),
    "incompatible allocation specifications"
  )
  expect_error(
    JAGS_get_inits(incompatible_priors, chains = 1, seed = 1),
    "incompatible allocation specifications"
  )
  expect_error(
    JAGS_to_monitor(incompatible_priors),
    "incompatible allocation specifications"
  )

  formula_info_collision_a <- JAGS_formula(
    y ~ f,
    "phi",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), id = "shape-a")
    )
  )
  formula_info_collision_b <- JAGS_formula(
    y ~ f,
    "theta",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), id = "shape a")
    )
  )
  colliding_priors <- c(
    formula_info_collision_a$prior_list,
    formula_info_collision_b$prior_list
  )
  expect_error(
    JAGS_add_priors("model{}", colliding_priors),
    "generate the same JAGS node"
  )
})

test_that("ordered hidden total nodes cannot collide with formula coefficients", {
  df <- data.frame(
    y = seq_len(6),
    f = ordered(
      rep(c("low", "mid", "high"), 2),
      levels = c("low", "mid", "high")
    ),
    f_ordered_total = rep(c(-1, 1), 3)
  )

  expect_error(
    JAGS_formula(
      y ~ f + f_ordered_total,
      "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        f = prior_ordered(prior("normal", list(0, 1))),
        f_ordered_total = prior("normal", list(0, 1))
      )
    ),
    "generates hidden total node 'mu_f_ordered_total'"
  )
})

test_that("multi-slice ordered expression totals omit initialization", {
  df <- expand.grid(
    f = ordered(
      c("low", "mid", "high"),
      levels = c("low", "mid", "high")
    ),
    g = factor(c("a", "b", "c"), levels = c("a", "b", "c"))
  )
  formula_info <- JAGS_formula(
    ~ f * g,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1))),
      g = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "f:g" = prior_ordered(
        prior("normal", list(0, expression(sigma)))
      )
    )
  )

  expect_equal(
    attr(formula_info$prior_list$mu_f__xXx__g, "ordered_metadata")$theta_dim,
    2L
  )
  inits <- JAGS_get_inits(
    formula_info$prior_list,
    chains = 1,
    seed = 1
  )[[1L]]
  expect_false("mu_f__xXx__g_ordered_total" %in% names(inits))
  expect_true(any(grepl(
    "prior_par_eta_mu_f__xXx__g_ordered_alloc",
    names(inits),
    fixed = TRUE
  )))
})

test_that("ordered fixed contrasts propagate to random slope designs", {
  df <- data.frame(
    y = seq_len(12),
    f = ordered(rep(c("low", "mid", "high"), 4), levels = c("low", "mid", "high")),
    id = factor(rep(seq_len(4), each = 3))
  )

  formula_info <- JAGS_formula(
    ~ 1 + f + (f || id),
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("normal", list(0, 1), list(0, Inf)))
    )
  )

  random_term <- formula_info$formula_design$random_effects[[1]]
  expect_equal(random_term$sd_leaves$leaf_names, c("mu__xREx__id_intercept", "mu__xREx__id_f"))
  expect_equal(
    random_term$sd_leaves$leaf_names_by_column,
    c("mu__xREx__id_intercept", "mu__xREx__id_f", "mu__xREx__id_f")
  )
  expect_equal(
    unname(formula_info$data$mu__xREx__id_xRE_DATAx[1:3, ]),
    matrix(c(1, 0, 0, 1, 1, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
  )
})

test_that("prior_ordered() can define ordered random slope SD components", {
  df <- data.frame(
    y = seq_len(12),
    f = ordered(rep(c("low", "mid", "high"), 4), levels = c("low", "mid", "high")),
    id = factor(rep(seq_len(4), each = 3))
  )

  formula_info <- JAGS_formula(
    ~ 1 + (0 + f || id),
    "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior_ordered(
          prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf)),
          allocation = c(.4, .6)
        )
      )
    )
  )

  random_term <- formula_info$formula_design$random_effects[[1]]
  expect_equal(random_term$sd_parameter_names, c("mu__xREx__id_f[1]", "mu__xREx__id_f[2]"))
  expect_s3_class(formula_info$prior_list$mu__xREx__id_f, "prior.ordered")

  syntax <- JAGS_add_priors("model{}", formula_info$prior_list)
  expect_match(syntax, "mu__xREx__id_f_ordered_total ~ dnorm\\(0,1\\)T\\(0,\\)")
  expect_match(syntax, "mu__xREx__id_f\\[1\\] <- mu__xREx__id_f_ordered_total \\* 0.4")
  expect_match(syntax, "mu__xREx__id_f\\[2\\] <- mu__xREx__id_f_ordered_total \\* 0.6")
})

test_that("ordered posterior extraction transforms coefficients to public levels", {
  df <- data.frame(
    y = seq_len(6),
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )
  formula_info <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), allocation = c(.25, .75))
    )
  )
  samples <- matrix(
    c(1, 2, 4, 8),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_f[1]", "mu_f[2]"))
  )

  transformed <- BayesTools:::.transform_factor_contrasts(
    samples,
    formula_info$prior_list,
    transform_factors = TRUE
  )

  expect_equal(colnames(transformed), c("mu_f[dif: low]", "mu_f[dif: mid]", "mu_f[dif: high]"))
  expect_equal(unname(transformed[1, ]), c(0, 1, 3))
  expect_equal(unname(transformed[2, ]), c(0, 4, 12))
})

test_that("ordered densities are direct when supported and bridge sampling stops for complex totals", {
  p <- prior_ordered(prior("normal", list(0, 1)), allocation = c(.25, .75))
  attr(p, "levels") <- 3
  density_fixed <- density(p, n_points = 21)
  expect_s3_class(density_fixed, "density.prior.ordered")
  expect_equal(attr(density_fixed, "method"), "direct")
  expect_equal(density_fixed[[2]]$y[11], stats::dnorm(density_fixed[[2]]$x[11] / .25) / .25)

  p_dir <- prior_ordered(prior("normal", list(0, 1)), allocation = prior("dirichlet", list(alpha = c(1, 1))))
  attr(p_dir, "levels") <- 3
  density_dir <- density(p_dir, n_points = 11)
  expect_equal(attr(density_dir, "method"), "direct")

  df <- data.frame(
    y = seq_len(6),
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
  )
  formula_info <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)))
    )
  )
  samples <- matrix(
    1,
    nrow = 2,
    ncol = 6,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_f[1]",
      "mu_f[2]",
      "mu_f_ordered_total",
      "prior_par_eta_mu_f_ordered_alloc_f_1[1]",
      "prior_par_eta_mu_f_ordered_alloc_f_1[2]"
    ))
  )
  bridge_samples <- JAGS_bridgesampling_posterior(samples, formula_info$prior_list)
  expect_true("mu_f_ordered_total" %in% colnames(bridge_samples))
  expect_false("mu_f[1]" %in% colnames(bridge_samples))

  formula_info$prior_list$mu_f$total <- prior_spike_and_slab(prior("normal", list(0, 1)))
  expect_error(
    JAGS_bridgesampling_posterior(samples, formula_info$prior_list),
    "only available when 'total' is a simple scalar prior"
  )

  formula_info$prior_list$mu_f$total <- prior("normal", list(0, expression(sigma)))
  expect_error(
    JAGS_bridgesampling_posterior(samples, formula_info$prior_list),
    "does not support parameter expressions in 'total'"
  )

  formula_info$prior_list$mu_f$total <- prior("bernoulli", list(.5))
  expect_error(
    JAGS_bridgesampling_posterior(samples, formula_info$prior_list),
    "requires a continuous or point-valued 'total' prior"
  )

  formula_info$prior_list$mu_f$total <- prior("point", list(1))
  expect_no_error(
    JAGS_bridgesampling_posterior(samples, formula_info$prior_list)
  )
})
