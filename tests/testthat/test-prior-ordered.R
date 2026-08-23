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

test_that("ordered random slope contrasts are specified independently", {
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
      id = random_block(
        sd = prior("normal", list(0, 1), list(0, Inf)),
        contrasts = c(f = "cumulative")
      )
    )
  )

  random_term <- formula_info$formula_design$random_effects[[1]]
  expect_equal(
    random_term$sd_leaves$leaf_names,
    c(
      "mu__xREx__id_intercept",
      "mu__xREx__id_f[1]",
      "mu__xREx__id_f[2]"
    )
  )
  expect_equal(
    random_term$sd_leaves$leaf_names_by_column,
    c(
      "mu__xREx__id_intercept",
      "mu__xREx__id_f[1]",
      "mu__xREx__id_f[2]"
    )
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

test_that("prior sample generation uses stored ordered-mixture dimensions", {
  make_prior <- function(total, prior_weights){
    p <- prior_ordered(
      prior("point", list(location = total)),
      allocation = c(.2, .3, .5),
      contrast = "cumulative_levels",
      prior_weights = prior_weights
    )
    attr(p, "levels") <- 3
    attr(p, "level_names") <- c("low", "mid", "high")
    p
  }

  p1 <- make_prior(10, 1)
  p2 <- make_prior(20, 3)

  single <- BayesTools:::.generate_prior_sample_matrix(
    prior_list = list(mu_f = p1),
    n_samples  = 4,
    seed       = 1
  )
  expect_equal(
    unname(single[, , drop = FALSE]),
    matrix(c(2, 3, 5), nrow = 4, ncol = 3, byrow = TRUE)
  )
  expect_equal(colnames(single), paste0("mu_f[", 1:3, "]"))

  mixture_prior <- prior_mixture(list(p1, p2))
  mixed <- BayesTools:::.generate_prior_sample_matrix(
    prior_list = list(mu_f = mixture_prior),
    n_samples  = 32,
    seed       = 2
  )
  row_values <- apply(mixed, 1L, paste, collapse = ",")
  expect_setequal(
    unique(row_values),
    c("2,3,5", "4,6,10")
  )
  expect_equal(colnames(mixed), paste0("mu_f[", 1:3, "]"))
})

test_that("fixed ordered allocations canonicalize only roundoff drift", {

  allocation <- c(.2, .3, .5 + .Machine$double.eps)
  spec <- BayesTools:::.prior_ordered_allocation_spec(allocation)

  expect_equal(sum(spec$weights), 1)
  expect_equal(
    spec$canonicalization$original_sum,
    sum(allocation)
  )
  expect_lte(
    spec$canonicalization$max_correction,
    spec$canonicalization$roundoff_bound
  )

  expect_error(
    BayesTools:::.prior_ordered_allocation_spec(c(.2, .3, .5 + 1e-10)),
    "exceeding the roundoff bound"
  )
})

test_that("public posterior mixing preserves ordered coefficient rows and metadata", {
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
  ordered_prior <- formula_info$prior_list$mu_f
  posterior <- matrix(
    c(1, 2, 4, 8),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_f[1]", "mu_f[2]"))
  )

  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_info$prior_list
  single <- as_mixed_posteriors(fit, parameters = "mu_f")

  expect_equal(unname(single$mu_f[, , drop = FALSE]), unname(posterior))
  expect_true(isTRUE(attr(single$mu_f, "ordered")))
  expect_equal(
    attr(single$mu_f, "ordered_metadata"),
    attr(ordered_prior, "ordered_metadata")
  )

  single_levels <- transform_factor_samples(single)$mu_f
  expect_equal(
    unname(single_levels[, , drop = FALSE]),
    matrix(c(0, 1, 3, 0, 4, 12), nrow = 2, byrow = TRUE)
  )
  expect_equal(
    BayesTools:::.prior_factor_level_weight_matrix(single_levels, "mu_f"),
    structure(
      attr(ordered_prior, "factor_design"),
      dimnames = list(
        colnames(single_levels),
        c("mu_f[1]", "mu_f[2]")
      )
    )
  )

  make_model <- function(samples, prior){
    samples <- coda::mcmc(samples)
    model_fit <- structure(
      list(
        mcmc = coda::mcmc.list(samples),
        sample = nrow(samples),
        summary.pars = list(mutate = NULL),
        monitor = colnames(samples)
      ),
      class = c("runjags", "BayesTools_fit", "list")
    )
    attr(model_fit, "prior_list") <- list(mu_f = prior)
    list(
      fit = model_fit,
      marglik = bridgesampling_object(0),
      prior_weights = 1
    )
  }

  posterior_1 <- matrix(
    c(1, 10, 2, 20, 3, 30),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_f[1]", "mu_f[2]"))
  )
  posterior_2 <- matrix(
    c(4, 40, 5, 50, 6, 60),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_f[1]", "mu_f[2]"))
  )
  ordered_prior_2 <- JAGS_formula(
    y ~ f,
    "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_ordered(prior("normal", list(0, 1)), allocation = c(.4, .6))
    )
  )$prior_list$mu_f

  mixed <- mix_posteriors(
    list(
      make_model(posterior_1, ordered_prior),
      make_model(posterior_2, ordered_prior_2)
    ),
    parameters = "mu_f",
    is_null_list = list(mu_f = c(FALSE, FALSE)),
    seed = 3,
    n_samples = 6
  )

  source_samples <- list(posterior_1, posterior_2)
  for(row_i in seq_len(nrow(mixed$mu_f))){
    model_i <- attr(mixed$mu_f, "models_ind")[[row_i]]
    sample_i <- attr(mixed$mu_f, "sample_ind")[[row_i]]
    expect_equal(
      unname(mixed$mu_f[row_i, ]),
      unname(source_samples[[model_i]][sample_i, ])
    )
  }
  expect_true(isTRUE(attr(mixed$mu_f, "ordered")))
  expect_false(is.null(attr(mixed$mu_f, "ordered_metadata")))
})

test_that("marginal posterior uses the stored full-rank ordered design", {
  df <- data.frame(
    y = seq_len(6),
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high"))
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
  posterior <- matrix(
    c(2, 3, 5, 4, 6, 10),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, paste0("mu_f[", 1:3, "]"))
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_info$prior_list
  mixed <- as_mixed_posteriors(fit, parameters = "mu_f")

  formula_marginal <- marginal_posterior(
    mixed,
    parameter = "mu_f",
    formula = ~ 0 + f,
    prior_samples = FALSE
  )
  direct_marginal <- marginal_posterior(
    mixed,
    parameter = "mu_f",
    use_formula = FALSE,
    prior_samples = FALSE
  )

  expected <- list(
    low = c(2, 4),
    mid = c(5, 10),
    high = c(10, 20)
  )
  expect_equal(names(formula_marginal), names(expected))
  expect_equal(names(direct_marginal), names(expected))
  for(level in names(expected)){
    expect_equal(as.numeric(formula_marginal[[level]]), expected[[level]])
    expect_equal(as.numeric(direct_marginal[[level]]), expected[[level]])
  }
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

test_that("ordered ranges cover zero and scaled total densities", {
  p_positive <- prior_ordered(
    prior("normal", list(10, 1)),
    allocation = c(.25, .75),
    contrast = "cumulative"
  )
  attr(p_positive, "levels") <- 3

  expected_range <- c(0, stats::qnorm(.995, mean = 10, sd = 1))
  expect_equal(range(p_positive), expected_range)
  expect_equal(
    range(p_positive, quantiles = .1),
    c(0, stats::qnorm(.9, mean = 10, sd = 1))
  )

  density_positive <- density(p_positive, n_points = 101)
  grid_step <- diff(expected_range) / 100
  expect_equal(attr(density_positive, "x_range"), expected_range)
  expect_lte(
    abs(density_positive[[2]]$x[which.max(density_positive[[2]]$y)] - 2.5),
    grid_step
  )
  expect_gt(max(density_positive[[2]]$y), 1)

  p_negative <- prior_ordered(
    prior("normal", list(-10, 1)),
    allocation = c(.25, .75),
    contrast = "cumulative"
  )
  expect_equal(
    range(p_negative),
    c(stats::qnorm(.005, mean = -10, sd = 1), 0)
  )
})

test_that("ordered continuous-mixture totals have deterministic density ranges", {
  total <- prior_mixture(
    list(
      prior("normal", list(-20, 1)),
      prior("normal", list(10, 1))
    ),
    components = c("lower", "upper")
  )
  p <- prior_ordered(
    total,
    allocation = c(.25, .75),
    contrast = "cumulative"
  )
  attr(p, "levels") <- 3

  expected_range <- c(
    stats::qnorm(.005, mean = -20, sd = 1),
    stats::qnorm(.995, mean = 10, sd = 1)
  )
  expect_equal(range(p), expected_range)

  set.seed(7403)
  density_mixture <- density(p, n_points = 101, n_samples = 2000)

  expect_s3_class(density_mixture, "density.prior.ordered")
  expect_equal(attr(density_mixture, "x_range"), expected_range)
  expect_s3_class(density_mixture[[1]], "density.prior.point")
  expect_true(all(vapply(density_mixture, function(component){
    all(is.finite(component$y))
  }, logical(1))))
  expect_equal(
    density_mixture[[2]]$samples * 4,
    density_mixture[[3]]$samples
  )
  expect_true(all(vapply(density_mixture[2:3], function(component){
    max(component$y) > 0
  }, logical(1))))
})

test_that("ordered sampled densities transform components and preserve cumulative atoms", {
  p <- prior_ordered(
    prior("normal", list(0, 1)),
    allocation = c(.25, .75),
    contrast = "cumulative"
  )
  attr(p, "levels") <- 3

  set.seed(7401)
  density_identity <- density(
    p,
    n_points = 31,
    n_samples = 2000,
    force_samples = TRUE
  )
  set.seed(7401)
  density_exp <- density(
    p,
    n_points = 31,
    n_samples = 2000,
    force_samples = TRUE,
    transformation = "exp"
  )

  expect_s3_class(density_exp[[1]], "density.prior.point")
  expect_equal(density_exp[[1]]$x, 1)
  expect_equal(unique(density_exp[[1]]$samples), 1)
  expect_equal(density_exp[[1]]$y, 1)
  expect_null(density_exp[[1]]$bw)

  expect_s3_class(density_exp[[2]], "density.prior.simple")
  expect_equal(density_exp[[2]]$x, exp(density_identity[[2]]$x))
  expect_equal(density_exp[[2]]$samples, exp(density_identity[[2]]$samples))
  expect_equal(
    density_exp[[2]]$y,
    density_identity[[2]]$y / density_exp[[2]]$x,
    tolerance = 1e-10
  )
  expect_true(all(density_exp[[2]]$x > 0))

  plot_exp <- plot(
    p,
    plot_type = "ggplot",
    show_figures = 1,
    n_points = 31,
    n_samples = 2000,
    force_samples = TRUE,
    transformation = "exp"
  )
  expect_s3_class(plot_exp, "ggplot")
  expect_s3_class(plot_exp$layers[[1]]$geom, "GeomSegment")
  expect_equal(NROW(plot_exp$layers[[1]]$data), 1)
  expect_equal(plot_exp$layers[[1]]$data$x, 1)
  expect_equal(plot_exp$layers[[1]]$data$xend, 1)
  expect_equal(plot_exp$layers[[1]]$data$yend, 1)
})

test_that("ordered cumulative atoms contribute to direct and transformed plot ranges", {
  p <- prior_ordered(
    prior("normal", list(10, 1)),
    allocation = c(.25, .75),
    contrast = "cumulative"
  )
  attr(p, "levels") <- 3

  density_direct <- density(p, n_points = 31)
  expect_s3_class(density_direct[[1]], "density.prior.point")
  expect_equal(density_direct[[1]]$x, 0)
  expect_equal(min(attr(density_direct[[1]], "x_range")), 0)
  expect_equal(min(attr(density_direct, "x_range")), 0)

  plot_direct <- plot(
    p,
    plot_type = "ggplot",
    show_figures = 1,
    n_points = 31
  )
  expect_s3_class(plot_direct$layers[[1]]$geom, "GeomSegment")
  expect_equal(plot_direct$layers[[1]]$data$x, 0)
  direct_limits <- plot_direct$scales$get_scales("x")$limits
  expect_lte(min(direct_limits), 0)
  expect_gte(max(direct_limits), 0)

  set.seed(7402)
  density_exp <- density(
    p,
    n_points = 31,
    n_samples = 2000,
    transformation = "exp"
  )
  expect_s3_class(density_exp[[1]], "density.prior.point")
  expect_equal(density_exp[[1]]$x, 1)
  expect_equal(min(attr(density_exp[[1]], "x_range")), 1)
  expect_equal(min(attr(density_exp, "x_range")), 1)

  set.seed(7402)
  plot_exp <- plot(
    p,
    plot_type = "ggplot",
    show_figures = 1,
    n_points = 31,
    n_samples = 2000,
    transformation = "exp"
  )
  expect_s3_class(plot_exp$layers[[1]]$geom, "GeomSegment")
  expect_equal(plot_exp$layers[[1]]$data$x, 1)
  exp_limits <- plot_exp$scales$get_scales("x")$limits
  expect_lte(min(exp_limits), 1)
  expect_gte(max(exp_limits), 1)
})

test_that("ordered mixed-measure densities preserve atoms and continuous mass", {
  total <- prior_spike_and_slab(
    prior("normal", list(0, 1)),
    prior_inclusion = prior("spike", list(.5))
  )
  fixed <- prior_ordered(
    total,
    allocation = c(.25, .75),
    contrast = "cumulative"
  )
  attr(fixed, "levels") <- 3

  density_fixed <- density(fixed, n_points = 201)
  expect_equal(attr(density_fixed, "method"), "analytic_mixed_measure")
  expect_equal(attr(density_fixed, "measure_schema_version"), 1L)
  expect_true(all(vapply(
    density_fixed,
    inherits,
    logical(1),
    what = "density.prior.mixed_measure"
  )))

  expect_equal(density_fixed[[1]]$atoms, data.frame(location = 0, mass = 1))
  expect_null(density_fixed[[1]]$continuous)
  for(i in 2:3){
    expect_equal(
      density_fixed[[i]]$atoms,
      data.frame(location = 0, mass = .5)
    )
    expect_equal(attr(density_fixed[[i]]$continuous, "mass"), .5)
    expect_equal(
      BayesTools:::.density.prior.ordered_curve_integral(
        density_fixed[[i]]$continuous$x,
        density_fixed[[i]]$continuous$density
      ),
      .5,
      tolerance = 1e-10
    )
    expect_equal(density_fixed[[i]]$diagnostics$atom_mass, .5)
    expect_equal(density_fixed[[i]]$diagnostics$continuous_mass, .5)
    expect_equal(
      density_fixed[[i]]$diagnostics$method,
      "analytic_components"
    )
  }

  density_exp <- density(
    fixed,
    n_points = 201,
    transformation = "exp"
  )
  expect_equal(density_exp[[2]]$atoms, data.frame(location = 1, mass = .5))
  expect_equal(attr(density_exp[[2]]$continuous, "mass"), .5)
  expect_equal(density_exp[[2]]$transformation$name, "exp")
  expect_true(all(density_exp[[2]]$continuous$x > 0))

  random_allocation <- prior_ordered(
    total,
    allocation = prior("dirichlet", list(alpha = c(2, 3))),
    contrast = "cumulative"
  )
  attr(random_allocation, "levels") <- 3
  density_random <- density(random_allocation, n_points = 201)
  expect_equal(density_random[[2]]$atoms$mass, .5)
  expect_equal(density_random[[3]]$atoms$mass, .5)
  expect_equal(attr(density_random[[2]]$continuous, "mass"), .5)
  expect_equal(attr(density_random[[3]]$continuous, "mass"), .5)

  plot_middle <- plot(
    fixed,
    plot_type = "ggplot",
    show_figures = 2,
    n_points = 201
  )
  geom_classes <- vapply(
    plot_middle$layers,
    function(layer) class(layer$geom)[1L],
    character(1)
  )
  expect_true("GeomLine" %in% geom_classes)
  expect_true("GeomSegment" %in% geom_classes)
  point_layer <- plot_middle$layers[[which(geom_classes == "GeomSegment")[1L]]]
  expect_equal(point_layer$data$x, 0)
  expect_equal(point_layer$data$yend, .5)

  layer_geoms <- geom_prior(fixed, show_parameter = 2, n_points = 201)
  expect_equal(
    vapply(layer_geoms, function(layer) class(layer$geom)[1L], character(1)),
    c("GeomLine", "GeomSegment")
  )
})

test_that("ordered mixed measures propagate through marginal inference", {
  df <- data.frame(
    y = seq_len(6),
    f = ordered(
      rep(c("low", "mid", "high"), 2),
      levels = c("low", "mid", "high")
    )
  )
  formula_info <- JAGS_formula(
    y ~ 0 + f,
    "mu",
    data = df,
    prior_list = list(
      f = prior_ordered(
        prior_spike_and_slab(
          prior("normal", list(0, 1)),
          prior_inclusion = prior("spike", list(.5))
        ),
        allocation = c(.25, .75),
        contrast = "cumulative"
      )
    )
  )
  posterior <- cbind(
    "mu_f[1]" = c(0, .1, .2, 0),
    "mu_f[2]" = c(0, .2, .3, 0),
    "mu_f_ordered_total_indicator" = c(0, 1, 1, 0)
  )
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_info$prior_list

  samples <- as_mixed_posteriors(fit, parameters = "mu_f")
  coefficient_atoms <- BayesTools:::.posterior_atoms_get(samples$mu_f)
  expect_equal(coefficient_atoms$locations, matrix(
    0,
    nrow = 1,
    ncol = 2,
    dimnames = list(NULL, c("mu_f[1]", "mu_f[2]"))
  ))
  expect_equal(coefficient_atoms$mass, .5)
  expect_equal(
    attr(samples$mu_f, "ordered_total_indicator"),
    c(0L, 1L, 1L, 0L)
  )

  for(use_formula in c(FALSE, TRUE)){
    marginal <- marginal_posterior(
      samples,
      parameter = "mu_f",
      formula = if(use_formula) ~ 0 + f else NULL,
      use_formula = use_formula,
      prior_samples = TRUE,
      n_samples = 128
    )
    expected_mass <- c(low = 1, mid = .5, high = .5)
    for(level in names(expected_mass)){
      prior_density <- attr(marginal[[level]], "prior_density")
      expect_equal(
        BayesTools:::.prior_linear_density_point_mass(prior_density, 0),
        expected_mass[[level]],
        tolerance = 1e-12
      )
      if(expected_mass[[level]] < 1){
        expect_equal(
          prior_density$density$mass,
          1 - expected_mass[[level]],
          tolerance = 1e-12
        )
      }
      posterior_atoms <- BayesTools:::.posterior_atoms_get(
        marginal[[level]]
      )
      expect_equal(posterior_atoms$mass, expected_mass[[level]])
      expect_equal(unname(posterior_atoms$locations[, 1L]), 0)
    }
  }

  prior_samples <- BayesTools:::.as_mixed_priors.factor(
    formula_info$prior_list$mu_f,
    parameter = "mu_f",
    seed = 991,
    n_samples = 50
  )
  expect_length(attr(prior_samples, "ordered_total_indicator"), 50)
  expect_equal(
    BayesTools:::.posterior_atoms_get(prior_samples)$mass,
    .5
  )

  expect_error(
    BayesTools:::.as_mixed_posteriors.factor(
      posterior[, c("mu_f[1]", "mu_f[2]"), drop = FALSE],
      formula_info$prior_list$mu_f,
      "mu_f"
    ),
    "required total-prior indicator"
  )
})

test_that("ordered Dirichlet marginals reject non-level combinations", {
  p <- prior_ordered(
    prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("spike", list(.5))
    ),
    allocation = prior("dirichlet", list(alpha = c(2, 3))),
    contrast = "cumulative"
  )
  attr(p, "levels") <- 3
  p <- BayesTools:::.prior_ordered_default_bound(p, "mu_f")

  expect_error(
    BayesTools:::.prior_ordered_linear_distribution(
      p,
      weights = c(1, 2),
      indices = 1:2,
      dx = .01,
      n_grid = 128
    ),
    "not a level or a single allocation subset"
  )
})
