skip_if_not_test_profile("certification")
skip_if_not_installed("rjags")

# Live assertions of the post-fit interface: input validation of JAGS_fit() and
# its helpers, convergence checks, formula and random-effect fits, the native
# LKJ modules, bridge-sampling and marginal-likelihood helpers, and the
# reproducibility of fits. They were the assertion blocks of
# test-00-model-fits.R, which runs only when its cache is stale; each block fits
# live, and the case is certified again when the code it reaches changes.

test_that("pre-fit sampler sources and the fitted catalog agree in their respective roles", {

  predictors <- data.frame(
    x = seq(1, 8.7, length.out = 9),
    z = c(-1, 0, 2, 3, -1, 1, 2, 0, 4),
    f = factor(rep(c("lo", "mid", "hi"), 3), levels = c("hi", "lo", "mid"))
  )
  bias <- prior_mixture(list(
    prior_none(),
    prior_weightfunction("one-sided", 0.05, wf_cumulative(c(1, 1))),
    prior_PET("normal", list(0, 1))
  ), is_null = c(TRUE, FALSE, FALSE))
  prefit <- new.env(parent = emptyenv())
  fit <- JAGS_fit(
    model_syntax = "model{ for(i in 1:N){ y[i] ~ dnorm(mu[i] + theta * s[i] + omega[2] + PET * s[i], 1) } }",
    data = list(N = 9L, y = c(-1, .4, .7, 1.1, .7, -.4, .6, .5, 1.5), s = seq(0, 1, length.out = 9)),
    prior_list = list(theta = prior("normal", list(0, 1)), bias = bias),
    formula_list = list(mu = ~ x + z + f),
    formula_data_list = list(mu = predictors),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      f = prior_factor("normal", list(0, 1), contrast = "treatment"),
      Signal = prior_R2(terms = "z", model = R2_model(law = "normal", sigma = 1))
    )),
    formula_scale_list = list(mu = list(x = TRUE)),
    chains = 1L, adapt = 100L, burnin = 50L, sample = 100L, seed = 123L, silent = TRUE,
    sampler_gates = TRUE,
    sampler_blocks = function(priors){
      prefit$block <- JAGS_sampler_sources(priors, c("theta", "mu_x", "mu_f", "mu_z", "mu_Signal"))
      prefit$bias <- JAGS_mixture_sources(priors, "bias")
      NULL
    }
  )
  expect_s3_class(fit, "BayesTools_fit")

  nodes <- fit$samplers$Node
  coordinates <- parameter_coordinates(fit)
  quantities <- parameter_catalog(fit, transform_scaled = FALSE)$quantities
  dependencies <- function(name){
    quantities$extraction_key[[match(name, quantities$canonical_name)]]$dependencies
  }
  eligible <- prefit$block$node[prefit$block$status == "eligible"]
  expect_setequal(eligible, c("theta", "mu_x", "mu_f[1]", "mu_f[2]", "mu_Signal"))
  expect_true(all(eligible %in% nodes))
  expect_true(all(eligible %in% coordinates$coordinate_name))
  expect_identical(dependencies("theta"), "theta")
  expect_identical(dependencies("mu_x"), "mu_x")
  expect_identical(quantities$display_scale[quantities$canonical_name == "mu_x"], "fitted_standardized")
  expect_identical(dependencies("mu_f[hi]"), character())
  expect_identical(dependencies("mu_f[lo]"), "mu_f[1]")
  expect_identical(dependencies("mu_f[mid]"), "mu_f[2]")

  expect_identical(prefit$block$status[prefit$block$prior == "mu_z"], "declined")
  expect_identical(coordinates$monitor_status[coordinates$coordinate_name == "mu_z"], "sampled")
  expect_identical(quantities$owner_type[quantities$canonical_name == "mu_z"], "r2_budget")
  expect_true("prior_par_z_mu_z" %in% nodes)
  expect_false("mu_z" %in% nodes)

  bias_gate <- prefit$bias$gate
  raw_pet <- prefit$bias$branches[[3L]]$source
  branch_weights <- prefit$bias$branches[[2L]]$omega
  expect_true(all(c(bias_gate, raw_pet) %in% nodes))
  expect_identical(dependencies(bias_gate), bias_gate)
  expect_identical(dependencies("PET"), "PET")
  expect_identical(dependencies("omega[2]"), "omega[2]")
  expect_false(any(c(raw_pet, branch_weights) %in% quantities$canonical_name))
})


# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-fit-edge-cases.R
# ============================================================================ #


# ============================================================================ #
# TEST FILE: JAGS Fit Edge Cases
# ============================================================================ #
#
# PURPOSE:
#   Edge case tests for JAGS fitting functions including input validation,
#   error handling, and boundary conditions.
#
# DEPENDENCIES:
#   - rjags: For JAGS model syntax generation and testing
#   - common-functions.R: REFERENCE_DIR, test_reference_text, skip_if_no_fits
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("rjags"): For all tests
#
# MODELS/FIXTURES:
#   - Some tests use pre-fitted models from test-00-model-fits.R
#
# TAGS: @edge-cases, @JAGS, @input-validation
# ============================================================================ #

# Reference directory for text output comparisons
REFERENCE_DIR <<- testthat::test_path("..", "results", "JAGS-fit-edge-cases")

source(testthat::test_path("common-functions.R"))


# ============================================================================ #
# SECTION 1: Input validation tests
# ============================================================================ #
test_that("JAGS_add_priors input validation works", {

  # Empty prior_list returns original syntax
  expect_equal(JAGS_add_priors("model{}", list()), "model{}")

  # prior_list must be a list of priors
  expect_error(JAGS_add_priors("model{}", list(x = 1)), "'prior_list' must be a list of priors.")
  expect_error(JAGS_add_priors("model{}", prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")

})


test_that("JAGS_get_inits input validation works", {

  # Empty prior_list returns only each chain's random-number generator entries
  expected_seeds <- local({
    set.seed(1)
    sample.int(.Machine$integer.max, 2)
  })
  expected_empty_inits <- lapply(expected_seeds, function(chain_seed){
    list(.RNG.seed = chain_seed, .RNG.name = "base::Super-Duper")
  })
  expect_identical(JAGS_get_inits(list(), chains = 2, seed = 1), expected_empty_inits)
  expect_identical(JAGS_get_inits(NULL, chains = 2, seed = 1), expected_empty_inits)

  # Input validation
  expect_error(JAGS_get_inits(list(x = 1), chains = 2, seed = 1), "'prior_list' must be a list of priors.")
  expect_error(JAGS_get_inits(prior("normal", list(0, 1)), chains = 2, seed = 1), "'prior_list' must be a list of priors.")

})


test_that("JAGS_to_monitor input validation works", {

  # Empty prior_list returns empty string
  expect_equal(JAGS_to_monitor(list()), "")

  # Input validation
  expect_error(JAGS_to_monitor(list(x = 1)), "'prior_list' must be a list of priors.")
  expect_error(JAGS_to_monitor(prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")

})


test_that(".check_JAGS_syntax validates syntax correctly", {

  # Test with valid syntax
  expect_silent(JAGS_add_priors("model{}", list(mu = prior("normal", list(0, 1)))))
  expect_equal(JAGS_add_priors(NULL, list()), "model{}")
  expect_match(JAGS_add_priors(NULL, list(mu = prior("normal", list(0, 1)))), "^model\\{")

  # Test with missing "model" keyword
  expect_error(
    JAGS_add_priors("invalid{}", list(mu = prior("normal", list(0, 1)))),
    "syntax must be a JAGS model syntax"
  )

  # Test with missing opening brace
  expect_error(
    JAGS_add_priors("model}", list(mu = prior("normal", list(0, 1)))),
    "syntax must be a JAGS model syntax"
  )

  # Test with missing closing brace
  expect_error(
    JAGS_add_priors("model{", list(mu = prior("normal", list(0, 1)))),
    "syntax must be a JAGS model syntax"
  )

  # Test with non-character input
  expect_error(
    JAGS_add_priors(123, list(mu = prior("normal", list(0, 1)))),
    "must be a character"
  )

})


test_that("JAGS_extend error handling", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  # Test error when fit is not a BayesTools_fit
  expect_error(
    JAGS_extend(list(), autofit_control = list()),
    "'fit' must be a fitted 'BayesTools_fit'"
  )

})


test_that("required packages are checked locally and on parallel workers", {

  expect_silent({
    local_loaded <- .JAGS_require_packages("stats")
  })
  expect_equal(unname(local_loaded), TRUE)
  expect_error(
    .JAGS_require_packages("BayesToolsMissingPackageForTest"),
    "Required packages are not available: 'BayesToolsMissingPackageForTest'.",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit("model{}", required_packages = NA_character_),
    "The 'required_packages' argument cannot contain NA/NaN values.",
    fixed = TRUE
  )

  cl <- parallel::makePSOCKcluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  expect_silent({
    worker_loaded <- .JAGS_require_packages("stats", cl)
  })
  expect_equal(unname(worker_loaded), TRUE)
  expect_error(
    .JAGS_require_packages(c("stats", "BayesToolsMissingPackageForTest"), cl),
    "Required packages are not available: 'BayesToolsMissingPackageForTest'.",
    fixed = TRUE
  )

})

test_that("JAGS_fit rejects unrecognized formula_scale_list entries", {


  data <- data.frame(x = c(-1, 0, 1))
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      formula_list = list(mu = ~ x),
      formula_data_list = list(mu = data),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      )),
      formula_scale_list = list(muu = list(x = TRUE))
    ),
    "not recognized"
  )
})

test_that("JAGS_fit stores formula design metadata on fitted formula models", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7),
    x = c(-1, 0, 1, 2)
  )

  fit <- JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 1)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ x),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1))
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 123
  )

  design <- JAGS_formula_design(fit, "mu")

  expect_s3_class(fit, "BayesTools_fit")
  expect_s3_class(design, "BayesTools_formula_design")
  expect_equal(colnames(design$model_matrix), c("(Intercept)", "x"))
  expect_equal(unname(design$model_matrix[, "x"]), df$x)
  expect_equal(design$jags_data_names$x, "mu_data_x")
  expect_identical(JAGS_formula_design(fit), list(mu = design))
})

test_that("JAGS_fit does not emit unused sampled data for marginalized random effects", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.3, -0.1),
    study = factor(c("s1", "s1", "s2", "s2")),
    estimate = factor(c("e1", "e2", "e3", "e4"))
  )
  sd_prior <- prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
  prior_random <- prior_random(
    study = random_block(sd = sd_prior),
    estimate = random_block(sd = sd_prior)
  )
  warnings <- character()

  fit <- withCallingHandlers(
    JAGS_fit(
      model_syntax = "model{
        for(i in 1:N_mu){
          y[i] ~ dnorm(mu[i], 1)
        }
      }",
      data = list(y = df$y),
      formula_list = list(
        mu = ~ 1 +
          random(1 | study, name = "study", covariance = "diag") +
          random(1 | estimate, name = "estimate", covariance = "diag")
      ),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(mean = 0, sd = 1))
      )),
      formula_random_prior_list = list(mu = prior_random),
      formula_random_effects_compile_list = list(
        mu = random_effects_compile(marginalized = "estimate")
      ),
      chains = 1,
      adapt = 100,
      burnin = 100,
      sample = 100,
      silent = FALSE,
      seed = 123
    ),
    warning = function(w){
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  design <- JAGS_formula_design(fit, "mu")
  estimate_term <- design$random_effects[[which(vapply(
    design$random_effects,
    function(random_term) identical(random_term$block_name, "estimate"),
    logical(1)
  ))]]

  expect_equal(warnings, character())
  expect_equal(estimate_term$compile_mode, "marginalized")
  expect_equal(estimate_term$jags_data_names, character())
  expect_false(grepl("mu__xREx__estimate_xRE_DATAx", attr(fit, "model"), fixed = TRUE))
  expect_false(grepl("mu__xREx__estimate_xRE_MAPx", attr(fit, "model"), fixed = TRUE))
})

test_that("JAGS_fit runs dummy structured random-effect formula models", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  sd_prior <- prior("normal", list(0, 0.5), truncation = list(lower = 0, upper = Inf))
  rho_prior <- prior("normal", list(0, 0.5))
  model_syntax <- "model{
    for(i in 1:N_mu){
      y[i] ~ dnorm(mu[i], 4)
    }
  }"

  df_cs <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3, 0.5, 0.9),
    x = c(-1, 0, 1, 2, -2, 0.5, 1.5, 2.5),
    idx = factor(rep(c("t1", "t2"), 4), levels = c("t1", "t2")),
    id = factor(rep(c("g1", "g2", "g3", "g4"), each = 2))
  )

  fit_cs <- suppressWarnings(JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = df_cs$y),
    formula_list = list(mu = ~ 1 + x + cs(idx | id)),
    formula_data_list = list(mu = df_cs),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = rho_prior,
        monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 101
  ))

  expect_s3_class(fit_cs, "BayesTools_fit")
  expect_equal(attr(fit_cs, "jags_modules"), character())
  expect_equal(JAGS_formula_design(fit_cs, "mu")$random_effects[[1]]$structure, "cs")
  expect_true("mu__xREx__id_rho" %in% colnames(as.matrix(fit_cs$mcmc)))

  df_ar1 <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3, 0.5, 0.9, 0.2, 0.6, 0.8, 1.0),
    f = factor(rep(c("a", "b", "c"), 4), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2", "g3", "g4"), each = 3))
  )

  fit_ar1 <- suppressWarnings(JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = df_ar1$y),
    formula_list = list(mu = ~ 1 + ar1(f | id)),
    formula_data_list = list(mu = df_ar1),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = rho_prior,
        monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 102
  ))

  expect_s3_class(fit_ar1, "BayesTools_fit")
  ar1_design <- JAGS_formula_design(fit_ar1, "mu")$random_effects[[1]]
  expect_equal(ar1_design$structure, "ar1")
  expect_equal(ar1_design$column_names, c("fa", "fb", "fc"))
  expect_true("mu__xREx__id_rho" %in% colnames(as.matrix(fit_ar1$mcmc)))

  df_us <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3),
    x = c(-1, 0, 1, 2, -2, 0.5),
    id = factor(c("g1", "g1", "g2", "g2", "g3", "g3"))
  )

  fit_us <- suppressWarnings(JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = df_us$y),
    formula_list = list(mu = ~ 1 + x + (1 + x | id)),
    formula_data_list = list(mu = df_us),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior_lkj(eta = 1, include_correlation = FALSE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 103
  ))

  expect_s3_class(fit_us, "BayesTools_fit")
  expect_equal(attr(fit_us, "jags_modules"), "BayesTools")
  expect_equal(JAGS_formula_design(fit_us, "mu")$random_effects[[1]]$structure, "us")
  expect_true(any(grepl("mu__xREx__id_xRE_CORx_L", colnames(as.matrix(fit_us$mcmc)), fixed = TRUE)))
})

test_that("emitted marginal random-effect covariance JAGS syntax reproduces random_effects_marginal_vcov()", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  # The syntax of JAGS_formula_random_marginal_covariance() is only text until
  # JAGS executes it: fit it for structures whose factor, distance, group
  # kernel and row-multiplier handling differ and compare the monitored nodes
  # of every draw with the covariance computed from the same draws in R.
  set.seed(93)
  data <- data.frame(
    g = factor(rep(c("a", "b", "c"), each = 4L)),
    x = round(stats::rnorm(12L), 2),
    z = round(stats::rnorm(12L), 2),
    f = factor(rep(c("p", "q", "r", "s"), 3L)),
    time = rep(c(0, 0.5, 2, 5), 3L),
    se = round(stats::runif(12L, 0.5, 1.5), 2)
  )
  y <- stats::rnorm(nrow(data))
  kernel <- matrix(
    c(1, .3, .1, .3, 1, .4, .1, .4, 1), 3L,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  sd_prior <- prior("gamma", list(4, 6))
  by_group <- unname(split(seq_len(nrow(data)), data$g))
  row_sd_syntax <- paste0(
    "for(i in 1:N){ tau[i] <- se[i] * s_tau }\n",
    "s_tau ~ dgamma(4, 6)\n"
  )
  cases <- list(
    us = list(
      formula = ~ 1 + us(1 + x + z | g),
      prior_random = prior_random(sd = sd_prior),
      row_blocks = by_group
    ),
    hcs = list(
      formula = ~ 1 + hcs(f | g),
      prior_random = prior_random(sd = sd_prior),
      row_blocks = by_group
    ),
    ar1 = list(
      formula = ~ 1 + ar1(f | g),
      prior_random = prior_random(sd = sd_prior),
      row_blocks = by_group
    ),
    car = list(
      formula = ~ 1 + car(time | g),
      prior_random = prior_random(sd = sd_prior),
      row_blocks = by_group
    ),
    known_us = list(
      formula = random_effects_formula(
        ~ 1 + us(1 + x | g),
        group_covariance = random_group_covariance(kernel)
      ),
      prior_random = prior_random(sd = sd_prior),
      row_blocks = list(seq_len(nrow(data)))
    ),
    row_indexed_sd = list(
      formula = ~ 1 + us(1 + x | g),
      prior_random = prior_random(
        g = random_block(sd_source = random_sd_source("tau", shape = "row"))
      ),
      row_blocks = by_group,
      syntax = row_sd_syntax,
      monitor = "tau"
    )
  )

  max_error <- function(compiled, draws, covariance, row_blocks){
    error <- 0
    for(block in seq_along(row_blocks)){
      rows <- row_blocks[[block]]
      if(identical(compiled$representation, "dense")){
        pairs <- BayesTools:::.bt_lower_triangle_pairs(rows)
        columns <- paste0(compiled$lower_names[block], "[", seq_len(nrow(pairs)), "]")
        reference <- t(vapply(seq_len(dim(covariance)[1L]), function(draw){
          matrix(covariance[draw, , ], dim(covariance)[2L])[cbind(pairs$row_1, pairs$row_2)]
        }, numeric(nrow(pairs))))
        if(nrow(pairs) == 1L){
          reference <- t(reference)
        }
        node <- draws[, columns, drop = FALSE]
        error <- max(error, abs(node - reference) / pmax(1, abs(reference)))
      }else{
        rank <- compiled$loading_ranks[block]
        for(draw in seq_len(dim(covariance)[1L])){
          diagonal <- draws[draw, paste0(compiled$diagonal_names[block], "[", seq_along(rows), "]")]
          loading <- if(rank > 0L){
            matrix(
              draws[draw, paste0(
                compiled$loading_names[block], "[",
                rep(seq_along(rows), rank), ",", rep(seq_len(rank), each = length(rows)), "]"
              )],
              length(rows), rank
            )
          }else{
            matrix(0, length(rows), 0L)
          }
          reference <- covariance[draw, rows, rows]
          error <- max(
            error,
            abs(diag(diagonal, length(rows)) + tcrossprod(loading) - reference) /
              pmax(1, abs(reference))
          )
        }
      }
    }

    error
  }

  for(case_name in names(cases)){
    case <- cases[[case_name]]
    compile <- random_effects_compile(marginalized = "g")
    specification <- JAGS_formula(
      case$formula, "mu", data,
      list(intercept = prior("normal", list(0, 1))),
      prior_random = case$prior_random,
      random_effects_compile = compile
    )
    for(representation in c("dense", "auto")){
      info <- paste(case_name, representation)
      compiled <- JAGS_formula_random_marginal_covariance(
        specification$formula_design, case$row_blocks,
        prefix = "rmc", representation = representation
      )
      expect_identical(
        compiled$representation,
        if(identical(representation, "dense")) "dense" else "diagonal_factor",
        info = info
      )
      monitors <- if(identical(compiled$representation, "dense")){
        compiled$lower_names
      }else{
        c(compiled$diagonal_names, compiled$loading_names[!is.na(compiled$loading_names)])
      }
      fit <- suppressWarnings(JAGS_fit(
        model_syntax = paste0(
          "model{\n for(i in 1:N){ y[i] ~ dnorm(mu[i], 1) }\n",
          case$syntax, compiled$syntax, "\n}\n"
        ),
        data = c(list(y = y, N = length(y), se = data$se), compiled$data),
        formula_list = list(mu = case$formula),
        formula_data_list = list(mu = data),
        formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)))),
        formula_random_prior_list = list(mu = case$prior_random),
        formula_random_effects_compile_list = list(mu = compile),
        add_parameters = c(monitors, case$monitor),
        chains = 1, adapt = 100, burnin = 50, sample = 100, thin = 5,
        silent = TRUE, seed = 1
      ))
      expect_s3_class(fit, "BayesTools_fit")
      if(!inherits(fit, "BayesTools_fit")){
        next
      }
      draws <- as.matrix(coda::as.mcmc.list(fit$mcmc))
      covariance <- random_effects_marginal_vcov(fit, parameter = "mu")$samples

      expect_identical(dim(covariance)[1L], nrow(draws), info = info)
      error <- max_error(compiled, draws, covariance, case$row_blocks)
      expect_lte(error, 1e-12, label = paste(info, "maximum relative error"))
    }
  }
})

test_that("JAGS_fit monitors random coefficients for observed-level prediction", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3),
    id = factor(c("g1", "g1", "g2", "g2", "g3", "g3"))
  )

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ 1 + diag(1 | id)),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(coefficients = TRUE, correlation = FALSE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 103
  ))

  posterior <- as.matrix(fit$mcmc)
  expect_true(any(grepl("mu__xREx__id_xRE_COEFx", colnames(posterior), fixed = TRUE)))

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = attr(fit, "prior_list"),
    formula_target = "conditional"
  )

  expect_equal(dim(prediction), c(nrow(df), nrow(posterior)))
  expect_true(all(is.finite(prediction)))
})

test_that("JAGS_fit predicts observed random effects from latent monitors", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3),
    id = factor(c("g1", "g1", "g2", "g2", "g3", "g3"))
  )

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ 1 + diag(1 | id)),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(coefficients = FALSE, correlation = FALSE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 104
  ))

  posterior <- as.matrix(fit$mcmc)
  expect_true(any(grepl("mu__xREx__id_xRE_Zx", colnames(posterior), fixed = TRUE)))
  expect_false(any(grepl("mu__xREx__id_xRE_COEFx", colnames(posterior), fixed = TRUE)))

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = attr(fit, "prior_list"),
    formula_target = "conditional"
  )

  expect_equal(dim(prediction), c(nrow(df), nrow(posterior)))
  expect_true(all(is.finite(prediction)))
})

test_that("mean-centered fits read the semantic predictor through their group locations", {

  skip_if_not_installed("rjags")

  # four groups observed twice on each level of a three-level factor; the
  # fitted predictor (from the group locations) equals the evaluator's
  # X c + Z b (from the standardized latents) at every draw, for correlated
  # (LKJ, compound-symmetry, and AR(1)) and nested blocks, and the fits use
  # all their JAGS data (no unused-data warning)
  set.seed(5)
  data <- expand.grid(rep = 1:2, f = factor(c("a", "b", "c")), g = factor(1:4))
  data$y <- c(0, 0.6, 1.1)[as.integer(data$f)] + stats::rnorm(4L, 0, 0.5)[as.integer(data$g)] +
    stats::rnorm(nrow(data), 0, 0.3)
  sd_prior <- prior("normal", list(0, 1), list(0, Inf))
  centered <- function(...) random_block(sd = sd_prior, parameterization = "mean_centered", ...)
  models <- list(
    us     = list(~ f + us(0 + f | g), prior_random(g = centered(cor = prior_lkj(eta = 1), contrasts = c(f = "independent")))),
    cs     = list(~ f + cs(f | g), prior_random(g = centered())),
    ar1    = list(~ f + ar1(f | g), prior_random(g = centered())),
    nested = list(~ f + (1 | g/f), prior_random(sd = sd_prior, g = centered(), f_g = centered()))
  )
  for(model in models){
    fit <- JAGS_fit(
      model_syntax = "model{ for(i in 1:N_mu){ y[i] ~ dnorm(mu[i], pow(0.3, -2)) } }",
      data = list(y = data$y),
      formula_list = list(mu = model[[1L]]), formula_data_list = list(mu = data),
      formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)),
                                          f = prior_factor("normal", list(0, 1), contrast = "treatment"))),
      formula_random_prior_list = list(mu = model[[2L]]),
      add_parameters = "mu", chains = 1, adapt = 100, burnin = 200, sample = 300, silent = TRUE, seed = 3
    )
    draws <- as.matrix(fit$mcmc)
    evaluated <- JAGS_evaluate_formula(fit, parameter = "mu", formula_target = "conditional")
    expect_lt(max(abs(t(evaluated) - draws[, paste0("mu[", seq_len(nrow(data)), "]")])), 1e-12)
  }
})

test_that("centered and noncentered fits preserve substantive output schemas", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.4, -0.1, 0.2, 0.5, 0.8, 0.1, 0.3, 0.6, 0.9, 1.1,
          -0.3, 0, 0.1, 0.4, 0.7, 0.2, 0.5, 0.7, 1, 1.2),
    index = factor(rep(c("i1", "i2"), 10L)),
    id = factor(rep(paste0("g", seq_len(4L)), each = 5L))
  )
  # Integrate the Gaussian intercept and the two correlated effects per group.
  # The marginal observation covariance is I / 4 + 11' + Z Sigma Z'.
  # Regular quadrature on log SD and Fisher z then gives the posterior means
  # independently of either JAGS parameterization.
  log_sd <- seq(-16, 4, length.out = 101)
  fisher_z <- seq(-4.5, 4.5, length.out = 81)
  grid <- expand.grid(log_sd = log_sd, z = fisher_z)
  grid$sd <- exp(grid$log_sd)
  grid$rho <- tanh(grid$z)
  group <- outer(as.integer(df$id), as.integer(df$id), "==")
  level <- outer(as.integer(df$index), as.integer(df$index), "==")
  fixed_cov <- 1 + diag(0.25, nrow(df))
  moments <- vapply(seq_len(nrow(grid)), function(i) {
    covariance <- fixed_cov + group * grid$sd[i]^2 * ifelse(level, 1, grid$rho[i])
    upper <- chol(covariance)
    residual <- forwardsolve(t(upper), df$y)
    c(loglik = -sum(log(diag(upper))) - sum(residual^2) / 2,
      intercept = sum(backsolve(upper, residual)))
  }, numeric(2L))
  log_weight <- moments["loglik", ] + stats::dgamma(grid$sd, 2, rate = 2, log = TRUE) +
    grid$log_sd + stats::dnorm(grid$z, 0, 0.5, log = TRUE)
  weight <- exp(log_weight - max(log_weight))
  reference <- vapply(list(moments["intercept", ], grid$sd, grid$rho),
                      stats::weighted.mean, numeric(1L), w = weight)

  fit_parameterization <- function(parameterization, seed){

    suppressWarnings(JAGS_fit(
      model_syntax = "model{
        for(i in 1:N_mu){
          y[i] ~ dnorm(mu[i], 4)
        }
      }",
      data = list(y = df$y),
      formula_list = list(mu = ~ 1 + cs(index | id)),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1))
      )),
      formula_random_prior_list = list(mu = prior_random(
        id = random_block(
          sd = prior("gamma", list(2, 2)),
          cor = prior("normal", list(0, 0.5)),
          parameterization = parameterization
        )
      )),
      chains = 4,
      adapt = 1000,
      burnin = 1800,
      sample = 5000,
      silent = TRUE,
      seed = seed
    ))
  }

  noncentered <- fit_parameterization("noncentered", 207)
  centered    <- fit_parameterization("centered", 208)
  noncentered_draws <- as.matrix(noncentered$mcmc)
  centered_draws    <- as.matrix(centered$mcmc)
  centered_term <- JAGS_formula_design(centered, "mu")$random_effects[[1L]]
  sd_name       <- centered_term$sd_parameter_names[[1L]]
  rho_name      <- centered_term$correlation$rho_name

  expect_setequal(colnames(centered_draws), colnames(noncentered_draws))
  for (fit in list(centered, noncentered)) {
    columns <- c("mu_intercept", sd_name, rho_name)
    chains <- coda::as.mcmc.list(fit)[, columns, drop = FALSE]
    draws <- as.matrix(chains)
    ess <- coda::effectiveSize(chains)
    mcse <- apply(draws, 2L, stats::sd) / sqrt(ess)
    expect_gt(min(ess), 150)
    expect_lt(max(abs(colMeans(draws) - reference) / mcse), 4)
  }
  expect_identical(centered_term$parameterization_resolved, "centered")
})

test_that("group-local structured fits monitor only active latent cells", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  K <- 40L
  df <- data.frame(
    y = seq(-0.5, 0.5, length.out = K),
    index = factor(
      paste0("level_", seq_len(K)),
      levels = paste0("level_", seq_len(K))
    ),
    id = factor(rep(paste0("group_", seq_len(8L)), length.out = K))
  )
  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ 1 + cs(index | id)),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = prior("point", list(location = 1)),
        cor = prior("point", list(location = 0.2))
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 209
  ))
  term <- JAGS_formula_design(fit, "mu")$random_effects[[1L]]
  posterior <- as.matrix(fit$mcmc)
  z_names <- grep("_xRE_Zx[", colnames(posterior), fixed = TRUE, value = TRUE)

  expect_s3_class(
    term$latent_layout,
    "BayesTools_random_effect_structured_local_layout"
  )
  expect_length(z_names, K)
  expect_lt(length(z_names), term$n_groups * term$n_columns)
  expect_false(any(grepl("_xRE_CORx_L[", colnames(posterior), fixed = TRUE)))
  expect_false(any(grepl("_xRE_CORx_R[", colnames(posterior), fixed = TRUE)))

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + cs(index | id),
    parameter = "mu",
    data = df,
    prior_list = attr(fit, "prior_list"),
    formula_target = "conditional"
  )
  expect_equal(dim(prediction), c(K, nrow(posterior)))
  expect_true(all(is.finite(prediction)))
})

test_that("JAGS_fit predicts row-indexed external SD random effects from latent monitors", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7),
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b")),
    tau_data = c(0.5, 0.75, 1.0, 1.25)
  )
  formula <- ~ 1 +
    random(1 | study, name = "study", covariance = "diag") +
    random(1 | drug, name = "drug", covariance = "diag")

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
        tau[i] <- tau_data[i]
      }
    }",
    data = list(y = df$y, tau_data = df$tau_data),
    formula_list = list(mu = formula),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )),
    add_parameters = c("tau", "mu"),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 105
  ))

  posterior <- as.matrix(fit$mcmc)
  tau_names <- paste0("tau[", seq_len(nrow(df)), "]")
  mu_names <- paste0("mu[", seq_len(nrow(df)), "]")

  expect_true(all(tau_names %in% colnames(posterior)))
  expect_true(all(mu_names %in% colnames(posterior)))
  expect_true(any(grepl("mu__xREx__study_xRE_Zx", colnames(posterior), fixed = TRUE)))
  expect_true(any(grepl("mu__xREx__drug_xRE_Zx", colnames(posterior), fixed = TRUE)))
  expect_false(any(grepl("_xRE_COEFx", colnames(posterior), fixed = TRUE)))
  expect_equal(
    unname(posterior[, tau_names, drop = FALSE]),
    unname(matrix(
      rep(df$tau_data, each = nrow(posterior)),
      nrow = nrow(posterior),
      ncol = nrow(df)
    )),
    tolerance = 1e-12
  )

  # The monitored mu includes the random-effect contributions.
  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = formula,
    parameter = "mu",
    data = df,
    fitted_rows = seq_len(nrow(df)),
    prior_list = attr(fit, "prior_list"),
    formula_target = "conditional"
  )

  expect_equal(dim(prediction), c(nrow(df), nrow(posterior)))
  expect_true(all(is.finite(prediction)))
  expect_equal(
    unname(t(prediction)),
    unname(posterior[, mu_names, drop = FALSE]),
    tolerance = 1e-8
  )
})

test_that("JAGS_fit stores formula design metadata on failed sampling objects", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(x = c(-1, 0, 1))

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      broken_node <- missing_node
    }",
    formula_list = list(mu = ~ x),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1))
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 123
  ))

  design <- JAGS_formula_design(fit, "mu")

  expect_s3_class(fit, "error")
  expect_false(inherits(fit, "BayesTools_fit"))
  expect_match(conditionMessage(fit), "missing_node", fixed = TRUE)
  expect_s3_class(design, "BayesTools_formula_design")
  expect_equal(unname(design$model_matrix[, "x"]), df$x)
})


# ============================================================================ #
# SECTION 2: Convergence edge cases
# ============================================================================ #
test_that("autofit settings keep indicator checks off by default", {

  settings <- JAGS_check_and_list_autofit_settings(list(
    max_Rhat = 1.05,
    min_ESS = 500,
    max_error = 0.01,
    max_SD_error = 0.05,
    max_time = list(time = 60, unit = "mins"),
    sample_extend = 1000,
    restarts = 10,
    max_extend = 10
  ))
  expect_false(settings$check_indicators)

  settings <- JAGS_check_and_list_autofit_settings(list(
    max_Rhat = 1.05,
    min_ESS = 500,
    max_error = 0.01,
    max_SD_error = 0.05,
    max_time = list(time = 60, unit = "mins"),
    sample_extend = 1000,
    restarts = 10,
    max_extend = 10,
    check_indicators = TRUE
  ))
  expect_true(settings$check_indicators)

})


test_that("JAGS_check_convergence ignores indicator variables unless requested", {

  set.seed(1)
  chain_1 <- cbind(mu = rnorm(100), mu_indicator = rep(1, 100))
  chain_2 <- cbind(mu = rnorm(100), mu_indicator = rep(2, 100))
  fit <- list(
    mcmc         = coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2))
  )
  class(fit) <- "BayesTools_fit"

  # The mixture prior declares 'mu_indicator' as its component indicator.
  prior_list <- list(mu = prior_mixture(list(
    prior("normal", list(0, 1)),
    prior("normal", list(1, 1))
  )))
  attr(fit, "prior_list") <- prior_list
  fit <- attach_test_parameter_map(fit)

  expect_true(JAGS_check_convergence(
    fit,
    prior_list       = prior_list,
    max_Rhat        = 1.05,
    min_ESS         = NULL,
    max_error       = NULL,
    max_SD_error    = NULL,
    check_indicators = FALSE
  ))

  with_indicators <- JAGS_check_convergence(
    fit,
    prior_list       = prior_list,
    max_Rhat        = 1.05,
    min_ESS         = NULL,
    max_error       = NULL,
    max_SD_error    = NULL,
    check_indicators = TRUE
  )
  expect_false(with_indicators)
  expect_match(attr(with_indicators, "errors"), "R-hat")

})

test_that("JAGS_check_convergence ignores add_parameters without priors", {

  set.seed(2)
  chain_1 <- cbind(mu = rnorm(100), "aux[1]" = rep(0, 100))
  chain_2 <- cbind(mu = rnorm(100), "aux[1]" = rep(1, 100))
  fit <- list(
    mcmc         = coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2))
  )
  class(fit) <- "BayesTools_fit"

  prior_list <- list(mu = prior("normal", list(0, 1)))
  attr(fit, "prior_list") <- prior_list
  fit <- attach_test_parameter_map(fit)

  without_aux <- JAGS_check_convergence(
    fit,
    prior_list    = prior_list,
    max_Rhat     = 1.05,
    min_ESS      = NULL,
    max_error    = NULL,
    max_SD_error = NULL
  )
  expect_false(without_aux)

  expect_true(JAGS_check_convergence(
    fit,
    prior_list     = prior_list,
    max_Rhat      = 1.05,
    min_ESS       = NULL,
    max_error     = NULL,
    max_SD_error  = NULL,
    add_parameters = "aux"
  ))
})


test_that("single-chain fits skip R-hat and retain the remaining criteria", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  prior_list <- list(mu = prior("normal", list(0, 1)))

  fit <- suppressWarnings(JAGS_fit(
    "model{}", prior_list = prior_list,
    chains = 1,  # Single chain - R-hat cannot be computed
    adapt = 50,
    burnin = 50,
    sample = 100,
    seed = 1
  ))

  expect_warning(
    convergence <- JAGS_check_convergence(
      fit,
      prior_list = prior_list,
      max_Rhat = 1.05,
      min_ESS = NULL,
      max_error = NULL,
      max_SD_error = NULL
    ),
    "Only one chain was run", fixed = TRUE
  )
  expect_true(convergence)
  expect_true(is.na(attr(convergence, "diagnostics")$Rhat))

  expect_warning(autofit <- JAGS_fit(
    "model{}", prior_list = prior_list, chains = 1,
    adapt = 50, burnin = 50, sample = 1000, seed = 1, silent = TRUE,
    autofit = TRUE,
    autofit_control = list(max_Rhat = 1.05, min_ESS = 10,
      max_error = NULL, max_SD_error = NULL, max_extend = 1,
      max_time = NULL, sample_extend = 100, restarts = 1,
      check_indicators = FALSE)
  ), "Only one chain was run", fixed = TRUE)
  expect_s3_class(autofit, "BayesTools_fit")
  expect_equal(nrow(.extract_posterior_samples(autofit, as_list = TRUE)[[1L]]), 1000)
  expect_false(any(grepl("max_extend", attr(autofit, "warnings"), fixed = TRUE)))

})


test_that("JAGS_check_convergence handles ESS and error checks", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  prior_list <- list(mu = prior("normal", list(0, 1)))

  fit <- suppressWarnings(JAGS_fit(
    "model{}", prior_list = prior_list,
    chains = 2,
    adapt = 50,
    burnin = 50,
    sample = 100,  # Small sample for testing convergence failures
    seed = 1
  ))

  # Test with very strict ESS requirement (should fail)
  result_ess <- JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = NULL, min_ESS = 10000, max_error = NULL, max_SD_error = NULL, fail_fast = FALSE)
  expect_false(result_ess)
  expect_true(!is.null(attr(result_ess, "errors")))

  # Test with very strict error requirement
  result_err <- JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = NULL, min_ESS = NULL, max_error = 0.00001, max_SD_error = NULL, fail_fast = FALSE)
  expect_false(result_err)

  # Test with very strict SD error requirement
  result_sd <- JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = 0.00001, fail_fast = FALSE)
  expect_false(result_sd)

})


# ============================================================================ #
# SECTION 3: JAGS_fit with is_JASP mode
# ============================================================================ #
test_that("JAGS_fit works with is_JASP mode", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  # Simple model for testing is_JASP mode
  set.seed(1)
  data <- list(
    y = rnorm(20, 0.5, 1),
    N = 20
  )

  prior_list <- list(
    mu    = prior("normal", list(0, 1)),
    sigma = prior("normal", list(0, 1), list(0, Inf))
  )

  model_syntax <- "model{
    for(i in 1:N){
      y[i] ~ dnorm(mu, 1/pow(sigma, 2))
    }
  }"

  # Mock JASP progress bar functions (they should be skipped if not available)
  # The is_JASP mode should work but simply skip progress bars if functions don't exist
  fit_jasp <- capture.output(tryCatch({
    suppressWarnings(JAGS_fit(
      model_syntax = model_syntax,
      data = data,
      prior_list = prior_list,
      chains = 1,
      adapt = 50,
      burnin = 50,
      sample = 100,
      seed = 1,
      silent = TRUE,
      is_JASP = TRUE,
      is_JASP_prefix = "Test"
    ))
  }, error = function(e) {
    # If JASP functions don't exist, this should still produce a fit
    # or fail gracefully
    if (grepl("JASP", e$message)) {
      skip("JASP progress bar functions not available")
    }
    stop(e)
  }))

  test_reference_text(paste0(fit_jasp, collapse = ","), "fit_jasp.txt")

})

test_that("runtime reports come from the process of a serial fit and from every parallel worker", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  set.seed(1)
  data <- list(y = rnorm(20, 0.5, 1), N = 20)
  prior_list <- list(mu = prior("normal", list(0, 1)))
  model_syntax <- "model{
    for(i in 1:N){
      y[i] ~ dnorm(mu, 1)
    }
  }"
  report <- function(context){

    list(pid = Sys.getpid(), role = context$role, process_id = context$process_id)
  }
  environment(report) <- baseenv()
  fit_with <- function(...){
    JAGS_fit(model_syntax = model_syntax, data = data, prior_list = prior_list,
      chains = 2, adapt = 100, burnin = 100, sample = 200, seed = 1, silent = TRUE,
      runtime_report = report, ...)
  }

  serial <- fit_with()
  expect_identical(attr(serial, "runtime_report_values"),
    list(list(list(pid = Sys.getpid(), role = "local", process_id = 1L))))
  expect_identical(attr(serial, "runtime_report"), report)

  parallel_fit <- fit_with(parallel = TRUE, cores = 2)
  runs <- attr(parallel_fit, "runtime_report_values")
  expect_length(runs, 1L)
  expect_length(runs[[1L]], 2L)
  pids <- vapply(runs[[1L]], `[[`, integer(1L), "pid")
  expect_identical(anyDuplicated(pids), 0L)
  expect_false(Sys.getpid() %in% pids)
  expect_identical(vapply(runs[[1L]], `[[`, character(1L), "role"), rep("worker", 2L))
  expect_identical(vapply(runs[[1L]], `[[`, integer(1L), "process_id"), 1:2)
  expect_null(attr(parallel_fit, "warnings"))
})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-fit-lm-oracles.R
# ============================================================================ #


# ============================================================================ #
# TEST FILE: JAGS Fit LM Oracles
# ============================================================================ #
#
# PURPOSE:
#   Fit-profile semantic oracle tests for Gaussian JAGS formula models. These
#   compare automatic formula scaling and manual standardization against lm()
#   predictions on the original data scale, and compare a known-sigma Gaussian
#   model against its closed-form posterior.
#
# TAGS: @fit, @JAGS, @formula, @standardization
# ============================================================================ #


.fit_gaussian_formula_oracle <- function(data, formula_data, formula_scale_list = NULL, seed = 1L) {
  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], pow(sigma, -2))\n",
    "}\n",
    "}"
  )

  JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = data$y, N = nrow(data)),
    prior_list = list(
      sigma = prior("lognormal", list(0, 1))
    ),
    formula_list = list(mu = ~ x1 * x2),
    formula_data_list = list(mu = formula_data),
    formula_prior_list = list(
      mu = list(
        "intercept" = prior("normal", list(0, 10)),
        "x1" = prior("normal", list(0, 5)),
        "x2" = prior("normal", list(0, 5)),
        "x1:x2" = prior("normal", list(0, 5))
      )
    ),
    formula_scale_list = formula_scale_list,
    chains = 2,
    adapt = 250,
    burnin = 500,
    sample = 1200,
    seed = seed,
    silent = TRUE
  )
}

.fit_known_sigma_gaussian_formula_oracle <- function(data, sigma_known, seed = 1L) {
  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], pow(sigma_known, -2))\n",
    "}\n",
    "}"
  )

  JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = data$y, N = nrow(data), sigma_known = sigma_known),
    prior_list = NULL,
    formula_list = list(mu = ~ x1 * x2),
    formula_data_list = list(mu = data[c("x1", "x2")]),
    formula_prior_list = list(
      mu = list(
        "intercept" = prior("normal", list(0, 10)),
        "x1" = prior("normal", list(0, 5)),
        "x2" = prior("normal", list(0, 5)),
        "x1:x2" = prior("normal", list(0, 5))
      )
    ),
    formula_scale_list = list(mu = list(x1 = TRUE, x2 = TRUE)),
    chains = 2,
    adapt = 500,
    burnin = 1000,
    sample = 3000,
    seed = seed,
    silent = TRUE
  )
}

test_that("Gaussian JAGS formula fit agrees with lm oracle after scaling", {
  data <- bayestools_oracle_gaussian_regression_data()
  manual_formula_data <- bayestools_manual_scaled_data(data[c("x1", "x2")], c("x1", "x2"))
  manual_scale <- attr(manual_formula_data, "manual_scale")

  fit_manual <- .fit_gaussian_formula_oracle(
    data = data,
    formula_data = manual_formula_data,
    formula_scale_list = NULL,
    seed = 7701L
  )
  attr(fit_manual, "manual_scale") <- manual_scale

  fit_auto <- .fit_gaussian_formula_oracle(
    data = data,
    formula_data = data[c("x1", "x2")],
    formula_scale_list = list(mu = list(x1 = TRUE, x2 = TRUE)),
    seed = 7701L
  )

  expect_formula_scale_equal(attr(fit_auto, "formula_scale")$mu, manual_scale)

  scaled_parameters <- c("mu_intercept", "mu_x1", "mu_x2", "mu_x1__xXx__x2")
  posterior_manual <- as.matrix(suppressWarnings(coda::as.mcmc(fit_manual)))
  posterior_auto <- as.matrix(suppressWarnings(coda::as.mcmc(fit_auto)))

  expect_posterior_summary_close(
    posterior_auto[, scaled_parameters, drop = FALSE],
    colMeans(posterior_manual[, scaled_parameters, drop = FALSE]),
    tolerance = 0.05
  )

  posterior_original <- BayesTools:::.bt_transform_scale_posterior(
    posterior_auto[, scaled_parameters, drop = FALSE],
    attr(fit_auto, "formula_scale")
  )
  # the manual standardization with the fitted design of ~ x1 * x2, which
  # original-scale transforms require
  posterior_manual_original <- BayesTools:::.bt_transform_scale_posterior(
    posterior_manual[, scaled_parameters, drop = FALSE],
    list(mu = formula_scale_for_test(
      ~ x1 * x2,
      list(x1 = manual_scale$mu_x1, x2 = manual_scale$mu_x2),
      data = data[c("x1", "x2")]
    ))
  )

  lm_fit <- stats::lm(y ~ x1 * x2, data = data)
  lm_coef <- stats::coef(lm_fit)
  names(lm_coef) <- bayestools_lm_coef_to_jags_names(names(lm_coef))

  expect_posterior_summary_close(
    posterior_original[, names(lm_coef), drop = FALSE],
    lm_coef,
    tolerance = 0.30
  )

  newdata <- data.frame(
    x1 = stats::quantile(data$x1, probs = c(0.10, 0.50, 0.90), names = FALSE),
    x2 = stats::quantile(data$x2, probs = c(0.20, 0.60, 0.80), names = FALSE)
  )
  prediction_design <- stats::model.matrix(~ x1 * x2, data = newdata)
  prediction_names <- bayestools_lm_coef_to_jags_names(colnames(prediction_design))
  prediction_auto <- as.numeric(prediction_design %*% colMeans(
    posterior_original[, prediction_names, drop = FALSE]
  ))
  prediction_manual <- as.numeric(prediction_design %*% colMeans(
    posterior_manual_original[, prediction_names, drop = FALSE]
  ))
  expect_equal(prediction_auto, prediction_manual, tolerance = 0.15)

  expect_lm_predictions_equal(
    lm_fit = lm_fit,
    posterior = posterior_original,
    newdata = newdata,
    formula = ~ x1 * x2,
    tolerance = 0.30
  )

  sigma_samples <- posterior_auto[, "sigma"]
  expect_equal(mean(sigma_samples), summary(lm_fit)$sigma, tolerance = 0.25)

  posterior_intervals <- apply(
    posterior_original[, names(lm_coef), drop = FALSE],
    2,
    stats::quantile,
    probs = c(0.005, 0.995)
  )
  expect_true(all(lm_coef >= posterior_intervals[1, ] & lm_coef <= posterior_intervals[2, ]))
})

test_that("known-sigma Gaussian JAGS formula fit matches closed-form posterior", {
  data <- bayestools_oracle_gaussian_regression_data()
  sigma_known <- 0.6
  parameters <- c("mu_intercept", "mu_x1", "mu_x2", "mu_x1__xXx__x2")

  fit <- .fit_known_sigma_gaussian_formula_oracle(
    data = data,
    sigma_known = sigma_known,
    seed = 7711L
  )

  posterior <- as.matrix(suppressWarnings(coda::as.mcmc(fit)))
  posterior <- posterior[, parameters, drop = FALSE]

  scale_info <- attr(fit, "formula_scale")$mu
  expected_scale <- bayestools_manual_scaled_data(data[c("x1", "x2")], c("x1", "x2"))
  expect_formula_scale_equal(scale_info, attr(expected_scale, "manual_scale"))

  design <- stats::model.matrix(~ x1 * x2, data = expected_scale)
  jags_design <- JAGS_formula_design(fit, "mu")$model_matrix
  expect_equal(colnames(design), c("(Intercept)", "x1", "x2", "x1:x2"))
  expect_equal(colnames(jags_design), c("(Intercept)", "x1", "x2", "x1__xXx__x2"))
  expect_equal(unname(jags_design), unname(design), tolerance = 1e-12)
  expect_setequal(colnames(posterior), parameters)

  oracle <- bayestools_gaussian_posterior_oracle(
    X = design,
    y = data$y,
    sigma = sigma_known,
    prior_mean = c(0, 0, 0, 0),
    prior_sd = c(10, 5, 5, 5)
  )
  names(oracle$mean) <- parameters
  names(oracle$sd) <- parameters

  # Each quantity is compared at its own scale. The variances of the
  # coefficients are of order 0.01, below any absolute tolerance that suits
  # the means, so the covariance is compared as correlations and variances.
  # The Monte Carlo error follows from the effective sample size of each
  # coefficient. A comparison with a NaN fails (e.g. the correlations of a
  # zero covariance matrix).
  ess <- coda::effectiveSize(suppressWarnings(coda::as.mcmc.list(fit))[, parameters])
  posterior_sd <- apply(posterior, 2, stats::sd)
  # means: within three Monte Carlo standard errors of the closed-form means
  mean_error <- abs(colMeans(posterior) - oracle$mean) / (posterior_sd / sqrt(ess))
  expect_true(all(mean_error < 3), info = paste(round(mean_error, 2), collapse = " "))
  expect_equal(posterior_sd, oracle$sd, tolerance = 0.04)
  # variances: within 10% of the closed-form variances
  variance_ratio <- diag(stats::cov(posterior)) / diag(oracle$cov)
  expect_true(all(abs(variance_ratio - 1) < 0.1), info = paste(round(variance_ratio, 3), collapse = " "))
  # correlations: within three Monte Carlo standard errors of a sample
  # correlation at the smallest effective sample size
  correlation_error <- abs(stats::cov2cor(stats::cov(posterior)) - stats::cov2cor(oracle$cov))
  expect_true(all(correlation_error < 3 / sqrt(min(ess))), info = paste(round(max(correlation_error), 3)))

  posterior_intervals <- apply(posterior, 2, stats::quantile, probs = c(0.01, 0.99))
  expect_true(all(oracle$mean >= posterior_intervals[1, ] & oracle$mean <= posterior_intervals[2, ]))
})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-lkj-cholesky-fit.R
# ============================================================================ #


# ============================================================================ #
# TEST FILE: JAGS LKJ-Cholesky Fit Oracles
# ============================================================================ #
#
# PURPOSE:
#   Fit-profile MCMC checks for the package-shipped JAGS LKJ-Cholesky module.
#
# TAGS: @fit, @JAGS, @Stan, @LKJ, @Cholesky
# ============================================================================ #


.fit_jags_lkj_cholesky_prior <- function(K, eta, sample = 3000, seed = 1L) {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  module <- JAGS_lkj_corr_cholesky(
    name = "Omega",
    K = K,
    eta = eta,
    include_correlation = TRUE,
    include_primitives = TRUE
  )


  fit <- suppressWarnings(.bt_rjags_fit(
    model = paste0("model{\n", module$syntax, "\n}"),
    data = list(),
    monitor = module$monitor,
    adapt = 500,
    burnin = 500,
    sample = sample,
    thin = 1,
    inits = lapply(seq_len(2), function(i){
      list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed + i)
    })
  ))

  list(
    fit = fit,
    module = module,
    samples = as.matrix(fit$mcmc)
  )
}

.fit_jags_lkj_cholesky_prior_cached <- local({
  cache <- new.env(parent = emptyenv())

  function(K, eta, sample = 3000, seed = 1L) {
    key <- paste(K, eta, sample, seed, sep = "|")
    if(!exists(key, envir = cache, inherits = FALSE)){
      assign(
        key,
        .fit_jags_lkj_cholesky_prior(
          K = K,
          eta = eta,
          sample = sample,
          seed = seed
        ),
        envir = cache
      )
    }

    get(key, envir = cache, inherits = FALSE)
  }
})

.eval_jags_lkj_geometry <- function(function_name, u, K, declared = FALSE){

  target <- if(declared) "out[1:9]" else "out"
  data <- list(u = u, K = K)
  if(function_name == "bt_lkj_cholesky_ordered") data$position <- 1:2
  syntax <- sprintf("model{ %s <- %s(%s) }", target, function_name, paste(names(data), collapse = ", "))
  con <- textConnection(syntax)
  on.exit(close(con), add = TRUE)
  model <- rjags::jags.model(
    file = con, data = data, n.chains = 1L,
    n.adapt = 0L, quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model, variable.names = "out", n.iter = 1L,
    quiet = TRUE, progress.bar = "none"
  )
  as.numeric(as.matrix(samples)[1L, ])
}

test_that("native JAGS LKJ transforms reject invalid geometry before evaluation", {

  skip_if_not_installed("rjags")
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )
  for(function_name in c("bt_lkj_cholesky", "bt_lkj_cholesky_ordered", "bt_lkj_corr")){
    for(setting in list(
      list(K = 65537, u = rep(.5, 32768L)),
      list(K = 2^32, u = .5),
      list(K = 1, u = .5),
      list(K = 3, u = .5)
    )){
      expect_error(
        .eval_jags_lkj_geometry(function_name, setting$u, setting$K),
        "Invalid LKJ geometry", fixed = TRUE
      )
    }
    expect_error(
      .eval_jags_lkj_geometry(function_name, .5, 1.5),
      "Failed check for discrete-valued parameters", fixed = TRUE
    )
    expect_error(
      .eval_jags_lkj_geometry(function_name, .5, 3, declared = TRUE),
      "Invalid LKJ geometry", fixed = TRUE
    )

    expected <- if(function_name != "bt_lkj_corr"){
      c(1, 0, -.5, sqrt(3) / 2)
    }else{
      c(1, -.5, -.5, 1)
    }
    expect_equal(
      .eval_jags_lkj_geometry(function_name, .25, 2), expected,
      tolerance = 1e-12
    )
    expected_extension <- if(function_name != "bt_lkj_corr"){
      c(1, 0, -1, 0)
    }else{
      c(1, -1, -1, 1)
    }
    expect_equal(
      .eval_jags_lkj_geometry(function_name, -.25, 2), expected_extension,
      tolerance = 1e-12
    )
  }
})

test_that("native JAGS ordered LKJ factor is the R kernel and keeps the finite extension", {

  skip_if_not_installed("rjags")
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )
  evaluate <- function(syntax, data){
    con <- textConnection(syntax)
    on.exit(close(con), add = TRUE)
    model <- rjags::jags.model(
      file = con, data = data, n.chains = 1L, n.adapt = 0L, quiet = TRUE
    )
    samples <- rjags::coda.samples(
      model = model, variable.names = "out", n.iter = 1L,
      quiet = TRUE, progress.bar = "none"
    )
    as.numeric(as.matrix(samples)[1L, ])
  }
  ordered_syntax <- "model{ out[1:9] <- bt_lkj_cholesky_ordered(u, 3, position) }"
  corr_syntax    <- "model{ out[1:9] <- bt_lkj_corr(u, 3) }"
  for(position in list(c(1, 2, 3), c(3, 1, 2), c(2, 3, 1))){
    for(u in list(c(0.3, 0.6, 0.8), c(1e-40, 0.5, 1 - 1e-16))){
      # JAGS flattens the factor row by row
      expect_identical(
        evaluate(ordered_syntax, list(u = u, position = position)),
        as.vector(t(.bt_lkj_cholesky_cpc_u_to_L_ordered(u, K = 3, position = position)))
      )
    }
    # an out-of-support proposal keeps the finite extension of the clamped
    # primitives: the factor of their vine correlation matrix, reordered
    u <- c(-0.25, 1.3, 0.5)
    L <- matrix(evaluate(ordered_syntax, list(u = u, position = position)), 3L, byrow = TRUE)
    R <- matrix(evaluate(corr_syntax, list(u = u)), 3L, byrow = TRUE)
    expect_true(all(is.finite(L)))
    expect_true(all(L[upper.tri(L)] == 0))
    expect_true(all(diag(L) >= 0))
    expect_lt(max(abs(tcrossprod(L) - R[position, position])), 1e-14)
  }
  for(position in list(c(1, 1, 2), c(1, 2, 4))){
    expect_error(evaluate(ordered_syntax, list(u = c(0.3, 0.6, 0.8), position = position)))
  }

  # primitives that are partly fixed and partly sampled
  con <- textConnection(paste0(
    "model{ u[1] <- 0.5; u[2:3] ~ dbt_lkj_cpc(alpha[1:2]); ",
    "out[1:9] <- bt_lkj_cholesky_ordered(u[1:3], 3, position[1:3]) }"
  ))
  on.exit(close(con), add = TRUE)
  model <- rjags::jags.model(
    file = con, data = list(alpha = c(1.5, 1), position = c(2, 3, 1)),
    inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 1L),
    n.chains = 1L, n.adapt = 0L, quiet = TRUE
  )
  samples <- as.matrix(rjags::coda.samples(
    model = model, variable.names = c("u", "out"), n.iter = 50L,
    quiet = TRUE, progress.bar = "none"
  ))
  u <- samples[, paste0("u[", 1:3, "]"), drop = FALSE]
  expect_true(all(u[, 1L] == 0.5))
  expect_gt(length(unique(u[, 2L])), 1L)
  expected <- .bt_lkj_cholesky_cpc_u_to_L_ordered(u, K = 3, position = c(2, 3, 1))
  expect_identical(
    unname(samples[, paste0("out[", 1:9, "]"), drop = FALSE]),
    matrix(aperm(expected, c(1L, 3L, 2L)), nrow = nrow(u))
  )
})

test_that("native JAGS LKJ module of a unit layout matches its R evaluator", {

  skip_if_not_installed("rjags")
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )
  module <- .bt_JAGS_lkj_corr_cholesky_module(
    name = "Omega", K = 4, eta = 2, include_correlation = TRUE,
    include_primitives = TRUE, unit_columns = c(2L, 4L)
  )
  con <- textConnection(paste0("model{\n", module$syntax, "\n}"))
  on.exit(close(con), add = TRUE)
  model <- rjags::jags.model(
    file = con, n.chains = 1L, n.adapt = 100L, quiet = TRUE,
    inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 1L)
  )
  samples <- as.matrix(.bt_rjags_samples(model, module$monitor, 200L, 1L, "none"))
  node <- .bt_dnode_lkj(
    "Omega", K = 4, include_correlation = TRUE, include_primitives = TRUE,
    unit_columns = c(2L, 4L)
  )
  values <- .bt_deterministic_node_evaluate(
    node,
    .bt_deterministic_lookup(samples[, module$primitive_names, drop = FALSE])
  )
  expect_identical(unname(samples[, node$coordinates]), unname(values))
  expect_true(all(samples[, "Omega_R[2,4]"] == 0))
  expect_gt(length(unique(samples[, "Omega_R[1,2]"])), 1L)
})

test_that("native JAGS LKJ blocks with an exchangeable unit extend their fits", {

  skip_if_not_installed("rjags")
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  set.seed(3)
  data_formula <- data.frame(
    f  = factor(rep(c("a", "b", "c", "d"), 12L)),
    id = factor(rep(LETTERS[1:6], each = 8L))
  )
  data <- list(y = stats::rnorm(48), N = 48L)
  unit_fit <- function(seed){
    JAGS_fit(
      model_syntax = "model{\nfor(i in 1:N){\n  y[i] ~ dnorm(mu[i], 1)\n}\n}",
      data = data,
      formula_list = list(mu = ~ 1 + f + us(1 + f | id)),
      formula_data_list = list(mu = data_formula),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1)),
        f = prior_factor("mnormal", list(0, 1), contrast = "meandif")
      )),
      formula_random_prior_list = list(mu = prior_random(id = random_block(
        sd = prior("normal", list(0, 1), list(0, Inf)),
        cor = prior_lkj(eta = 2),
        monitor = random_monitor(correlation = TRUE)
      ))),
      chains = 2, adapt = 100, burnin = 100, sample = 200, seed = seed, silent = TRUE
    )
  }
  # the unit (columns 2:4) has exact zero within-unit correlations
  within <- paste0("mu__xREx__id_xRE_CORx_R[", c(2, 2, 3, 3, 4, 4), ",", c(3, 4, 2, 4, 2, 3), "]")

  fit <- unit_fit(22L)
  samples <- as.matrix(fit$mcmc)
  expect_true(all(samples[, within] == 0))
  expect_gt(length(unique(samples[, "mu__xREx__id_xRE_CORx_R[1,2]"])), 1L)
  extended <- JAGS_extend(
    fit,
    autofit_control = list(
      max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = NULL,
      max_time = list(time = 60, unit = "secs"),
      sample_extend = 100, restarts = 1, max_extend = 1
    ),
    silent = TRUE
  )
  expect_identical(JAGS_draw_geometry(extended)$total_draws, 600L)
  extended_samples <- as.matrix(extended$mcmc)
  expect_true(all(extended_samples[, within] == 0))
  expect_identical(parameter_catalog(extended), parameter_catalog(fit))

  # a unit leading the design keeps the vine order: its Cholesky block is the
  # constant identity, which the enabled convergence checks treat as structural
  data_formula$x <- stats::rnorm(48)
  leading <- JAGS_fit(
    model_syntax = "model{\nfor(i in 1:N){\n  y[i] ~ dnorm(mu[i], 1)\n}\n}",
    data = data,
    formula_list = list(mu = ~ 1 + f + x + us(0 + f + x | id)),
    formula_data_list = list(mu = data_formula),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      f = prior_factor("mnormal", list(0, 1), contrast = "meandif")
    )),
    formula_random_prior_list = list(mu = prior_random(id = random_block(
      sd = prior("normal", list(0, 1), list(0, Inf)),
      monitor = random_monitor(correlation = TRUE)
    ))),
    chains = 2, adapt = 100, burnin = 100, sample = 200, seed = 23L, silent = TRUE
  )
  unit_cholesky <- paste0("mu__xREx__id_xRE_CORx_L[", c(1, 2, 3, 2, 3, 3), ",", c(1, 1, 1, 2, 2, 3), "]")
  leading_samples <- as.matrix(leading$mcmc)[, unit_cholesky]
  # one value per cell over all draws: the identity up to rounding of the
  # diagonal row norms
  expect_identical(apply(leading_samples, 2L, function(cell) length(unique(cell))), stats::setNames(rep(1L, 6L), unit_cholesky))
  expect_equal(unname(leading_samples[1L, ]), c(1, 0, 0, 1, 0, 1), tolerance = 1e-14)
  check <- JAGS_check_convergence(leading, prior_list = attr(leading, "prior_list"),
    max_Rhat = 10, min_ESS = 1, max_error = NULL, max_SD_error = NULL)
  diagnostics <- attr(check, "diagnostics")
  expect_identical(
    diagnostics$state[match(unit_cholesky, diagnostics$parameter)],
    rep("structural_constant", length(unit_cholesky))
  )
  expect_false(any(diagnostics$state == "not_assessable"))
})

test_that("native JAGS LKJ primitive sampling refuses a draw rounded to the boundary of (0, 1)", {

  skip_if_not_installed("rjags")
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  # Beta(eta, eta) puts about 1% of its mass within 1e-16 of 1 at eta = 0.1:
  # a prior draw there rounds to u = 1 (log density -Inf, Cholesky diagonal 0).
  # Such a draw is refused, not clamped or redrawn.
  sample_prior <- function(eta, seed, n_iter){
    con <- textConnection(paste0(
      "model{ u[1:1] ~ dbt_lkj_cpc(alpha[1:1]); ",
      "L[1:4] <- bt_lkj_cholesky(u[1:1], 2); dummy ~ dnorm(0, 1) }"
    ))
    on.exit(close(con), add = TRUE)
    model <- rjags::jags.model(
      file = con, data = list(alpha = eta),
      inits = list(u = .5, dummy = 0, .RNG.name = "base::Wichmann-Hill", .RNG.seed = seed),
      n.chains = 1L, n.adapt = 0L, quiet = TRUE
    )
    as.matrix(rjags::coda.samples(
      model = model, variable.names = c("u", "L"), n.iter = n_iter,
      quiet = TRUE, progress.bar = "none"
    ))
  }

  expect_error(
    sample_prior(.1, 301L, 2000L),
    "rounded to the boundary of (0, 1)", fixed = TRUE
  )
  expect_error(
    sample_prior(.1, 301L, 2000L),
    "Supply initial values or use a larger 'eta'.", fixed = TRUE
  )

  draws <- sample_prior(2, 301L, 2000L)
  expect_true(all(draws[, "u"] > 0 & draws[, "u"] < 1))
  expect_true(all(draws[, "L[4]"] > 0))
})

.eval_jags_lkj_cholesky_transform <- function(u, K) {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  K2 <- K * K
  syntax <- c(
    paste0("  Omega_L_flat[1:", K2, "] <- bt_lkj_cholesky(u, ", K, ")"),
    paste0("  Omega_R_flat[1:", K2, "] <- bt_lkj_corr(u, ", K, ")")
  )
  for(i in seq_len(K)){
    for(j in seq_len(K)){
      flat_index <- BayesTools:::.bt_lkj_cholesky_flat_index(i, j, K)
      syntax <- c(
        syntax,
        paste0("  Omega_L[", i, ",", j, "] <- Omega_L_flat[", flat_index, "]"),
        paste0("  Omega_R[", i, ",", j, "] <- Omega_R_flat[", flat_index, "]")
      )
    }
  }

  con <- textConnection(paste0("model{\n", paste(syntax, collapse = "\n"), "\n}\n"))
  on.exit(close(con), add = TRUE)

  model <- rjags::jags.model(
    file = con,
    data = list(u = u),
    n.chains = 1,
    n.adapt = 0,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = c("Omega_L", "Omega_R"),
    n.iter = 1,
    quiet = TRUE,
    progress.bar = "none"
  )

  as.matrix(samples)
}

.jags_lkj_vector_column <- function(samples, base_name, index) {
  indexed_name <- paste0(base_name, "[", index, "]")
  if(indexed_name %in% colnames(samples)){
    return(indexed_name)
  }
  if(index == 1L && base_name %in% colnames(samples)){
    return(base_name)
  }
  stop("Missing monitored JAGS column: ", indexed_name, call. = FALSE)
}

.jags_lkj_matrix_draw <- function(samples, row, prefix, K) {
  out <- matrix(NA_real_, K, K)
  for(i in seq_len(K)){
    for(j in seq_len(K)){
      out[i, j] <- samples[row, paste0(prefix, "[", i, ",", j, "]")]
    }
  }
  out
}

.jags_lkj_offdiag_columns <- function(samples, prefix, K) {
  out <- character(0)
  for(row in 2:K){
    for(column in seq_len(row - 1L)){
      out <- c(out, paste0(prefix, "[", column, ",", row, "]"))
    }
  }
  if(!all(out %in% colnames(samples))){
    stop("Missing monitored JAGS correlation columns.", call. = FALSE)
  }

  out
}

.eval_jags_lkj_cpc_deviance <- function(u, alpha) {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  con <- textConnection(paste0(
    "model{\n",
    "  u[1:", length(u), "] ~ dbt_lkj_cpc(alpha)\n",
    "}\n"
  ))
  on.exit(close(con), add = TRUE)

  model <- rjags::jags.model(
    file = con,
    data = list(u = u, alpha = alpha),
    n.chains = 2,
    n.adapt = 0,
    quiet = TRUE
  )
  dic <- rjags::dic.samples(
    model = model,
    n.iter = 1,
    type = "pD",
    quiet = TRUE,
    progress.bar = "none"
  )

  list(
    deviance = as.numeric(dic[["deviance"]]),
    penalty = as.numeric(dic[["penalty"]])
  )
}

.lkj_bivariate_data <- function(N = 28L, rho = 0.55, seed = 4711L) {
  set.seed(seed)
  z_1 <- stats::rnorm(N)
  z_2 <- rho * z_1 + sqrt(1 - rho^2) * stats::rnorm(N)
  cbind(z_1, z_2)
}

.lkj_bivariate_loglik <- function(rho, y) {
  denom <- 1 - rho^2
  quad <- (y[, 1]^2 - 2 * rho * y[, 1] * y[, 2] + y[, 2]^2) / denom
  sum(-0.5 * log(denom) - 0.5 * quad)
}

.lkj_bivariate_grid_oracle <- function(y, eta, grid_size = 20001L) {
  rho <- seq(-.999, .999, length.out = grid_size)
  log_density <- vapply(rho, .lkj_bivariate_loglik, numeric(1), y = y) +
    (eta - 1) * log1p(-rho^2)
  log_density <- log_density - max(log_density)
  weights <- exp(log_density)
  weights <- weights / sum(weights)
  cdf <- cumsum(weights)
  posterior_mean <- sum(weights * rho)

  list(
    mean = posterior_mean,
    sd = sqrt(sum(weights * (rho - posterior_mean)^2)),
    quantiles = stats::approx(cdf, rho, xout = c(.1, .5, .9), ties = "ordered")$y
  )
}

.fit_jags_lkj_bivariate_posterior <- function(y, eta, sample = 5000,
                                              seed = 812L) {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  model <- paste0(
    "model{\n",
    "  alpha[1] <- eta\n",
    "  u[1:1] ~ dbt_lkj_cpc(alpha)\n",
    "  R_flat[1:4] <- bt_lkj_corr(u, 2)\n",
    "  rho <- R_flat[2]\n",
    "  for(i in 1:N){\n",
    "    loglik[i] <- -0.5 * log(1 - pow(rho, 2)) - ",
    "(pow(y[i,1], 2) - 2 * rho * y[i,1] * y[i,2] + pow(y[i,2], 2)) / ",
    "(2 * (1 - pow(rho, 2)))\n",
    "    zeros[i] ~ dpois(C - loglik[i])\n",
    "  }\n",
    "}\n"
  )


  fit <- suppressWarnings(.bt_rjags_fit(
    model = model,
    data = list(y = y, eta = eta, N = nrow(y), zeros = rep(0, nrow(y)), C = 100),
    monitor = "rho",
    adapt = 1000,
    burnin = 1000,
    sample = sample,
    thin = 1,
    inits = lapply(seq_len(2), function(i){
      list(u = c(0.5), .RNG.name = "base::Wichmann-Hill", .RNG.seed = seed + i)
    })
  ))

  as.numeric(as.matrix(fit$mcmc)[, "rho"])
}

test_that("JAGS LKJ-Cholesky module matches primitive and correlation marginals", {
  settings <- data.frame(
    K = c(2L, 3L, 4L, 5L),
    eta = c(0.35, 1, 2.5, 0.75),
    sample = c(3000L, 3000L, 3000L, 4000L)
  )

  for(setting_i in seq_len(nrow(settings))){
    K <- settings$K[setting_i]
    eta <- settings$eta[setting_i]
    result <- .fit_jags_lkj_cholesky_prior_cached(
      K = K,
      eta = eta,
      sample = settings$sample[setting_i],
      seed = 100 + setting_i
    )
    samples <- result$samples
    pairs <- result$module$pairs

    u_columns <- character(nrow(pairs))
    for(p in seq_len(nrow(pairs))){
      u_name <- .jags_lkj_vector_column(samples, "Omega_lkj_u", p)
      u_columns[p] <- u_name
      alpha <- pairs$alpha[p]
      expect_equal(mean(samples[, u_name]), 0.5, tolerance = 0.05)
      expect_true(abs(stats::var(samples[, u_name]) - 1 / (4 * (2 * alpha + 1))) < 0.012)
      expect_equal(
        stats::quantile(samples[, u_name], c(.1, .5, .9), names = FALSE),
        stats::qbeta(c(.1, .5, .9), alpha, alpha),
        tolerance = 0.06
      )
    }

    if(length(u_columns) >= 2L){
      primitive_cor <- stats::cor(samples[, u_columns, drop = FALSE])
      expect_lt(max(abs(primitive_cor[upper.tri(primitive_cor)])), 0.10)
    }

    rho_columns <- .jags_lkj_offdiag_columns(samples, "Omega_R", K)
    alpha_rho <- eta - 1 + K / 2
    rho_quantiles <- 2 * stats::qbeta(c(.1, .5, .9), alpha_rho, alpha_rho) - 1
    rho_means <- vapply(rho_columns, function(name) mean(samples[, name]), numeric(1))
    rho_sds <- vapply(rho_columns, function(name) stats::sd(samples[, name]), numeric(1))
    for(rho_name in rho_columns){
      expect_equal(mean(samples[, rho_name]), 0, tolerance = 0.045)
      expect_true(abs(stats::var(samples[, rho_name]) - 1 / (2 * alpha_rho + 1)) < 0.025)
      expect_equal(
        stats::quantile(samples[, rho_name], c(.1, .5, .9), names = FALSE),
        rho_quantiles,
        tolerance = 0.08
      )
    }
    expect_lt(max(rho_means) - min(rho_means), 0.08)
    expect_lt(max(rho_sds) - min(rho_sds), 0.08)

    check_rows <- unique(round(seq(1, nrow(samples), length.out = min(25, nrow(samples)))))
    for(row in check_rows){
      L <- .jags_lkj_matrix_draw(samples, row, "Omega_L", K)
      R <- .jags_lkj_matrix_draw(samples, row, "Omega_R", K)

      expect_equal(L[upper.tri(L)], rep(0, K * (K - 1) / 2), tolerance = 1e-10)
      expect_true(all(diag(L) > 0))
      expect_equal(rowSums(L^2), rep(1, K), tolerance = 1e-8)
      expect_equal(R, L %*% t(L), tolerance = 1e-8)
      expect_equal(diag(R), rep(1, K), tolerance = 1e-8)
      expect_true(all(eigen(R, symmetric = TRUE, only.values = TRUE)$values > 0))
    }
  }
})

test_that("compiled LKJ-Cholesky functions match hand-coded transform oracles", {
  settings <- list(
    list(
      K = 2L,
      u = c(0.31),
      expected_L = matrix(c(
        1, 0,
        -0.38, sqrt(1 - .38^2)
      ), nrow = 2, byrow = TRUE)
    ),
    list(
      K = 3L,
      u = c(0.61, 0.22, 0.74),
      expected_L = matrix(c(
        1, 0, 0,
        .22, sqrt(1 - .22^2), 0,
        -.56, .48 * sqrt(1 - .56^2), sqrt(1 - .56^2) * sqrt(1 - .48^2)
      ), nrow = 3, byrow = TRUE)
    )
  )

  for(setting in settings){
    samples <- .eval_jags_lkj_cholesky_transform(setting$u, setting$K)
    expected_L <- setting$expected_L
    expected_R <- expected_L %*% t(expected_L)
    observed_L <- .jags_lkj_matrix_draw(samples, 1, "Omega_L", setting$K)
    observed_R <- .jags_lkj_matrix_draw(samples, 1, "Omega_R", setting$K)

    expect_equal(observed_L, expected_L, tolerance = 1e-12)
    expect_equal(observed_R, expected_R, tolerance = 1e-12)
    expect_equal(observed_L[upper.tri(observed_L)], rep(0, setting$K * (setting$K - 1L) / 2L), tolerance = 1e-12)
    expect_true(all(diag(observed_L) > 0))
    expect_equal(rowSums(observed_L^2), rep(1, setting$K), tolerance = 1e-12)
    expect_equal(diag(observed_R), rep(1, setting$K), tolerance = 1e-12)
  }
})

test_that("compiled LKJ transforms evaluate boundary proposals", {
  settings <- list(
    list(K = 2L, u = -1e-8),
    list(K = 2L, u = 0),
    list(K = 2L, u = 1),
    list(K = 2L, u = 1 + 1e-8),
    list(K = 3L, u = c(0, 1, 0.5))
  )

  for(setting in settings){
    samples <- .eval_jags_lkj_cholesky_transform(setting$u, setting$K)
    L <- .jags_lkj_matrix_draw(samples, 1, "Omega_L", setting$K)
    R <- .jags_lkj_matrix_draw(samples, 1, "Omega_R", setting$K)

    expect_true(all(is.finite(L)))
    expect_true(all(is.finite(R)))
    expect_equal(R, L %*% t(L), tolerance = 1e-12)
    expect_equal(diag(R), rep(1, setting$K), tolerance = 1e-12)
  }
})

test_that("compiled LKJ CPC distribution logDensity matches beta density for observed nodes", {
  K <- 4L
  eta <- 0.75
  u <- c(0.08, 0.91, 0.47, 0.12, 0.63, 0.88)
  alpha <- BayesTools:::.bt_lkj_cholesky_cpc_pairs(K = K, eta = eta)$alpha
  deviance <- .eval_jags_lkj_cpc_deviance(u = u, alpha = alpha)
  expected_deviance <- -2 * sum(stats::dbeta(u, alpha, alpha, log = TRUE))

  expect_equal(deviance$deviance, expected_deviance, tolerance = 1e-12)
  expect_equal(deviance$penalty, 0, tolerance = 1e-12)
})

test_that("generated LKJ primitives are scalar nodes with bounded univariate slice samplers", {
  skip_if_not(
    isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)),
    "BayesTools JAGS module is unavailable."
  )

  # a bounded univariate slice sampler evaluates the model only inside (0, 1),
  # so the multivariate-normal likelihood never sees a singular correlation
  # matrix at a proposal outside the support
  for(unit_columns in list(integer(), c(1L, 3L))){
    module <- BayesTools:::.bt_JAGS_lkj_corr_cholesky_module(
      name = "Omega", K = 4L, eta = 0.8, include_correlation = TRUE,
      include_primitives = FALSE, unit_columns = unit_columns
    )
    y <- matrix(stats::rnorm(4 * 30), ncol = 4L) %*% chol(matrix(0.5, 4, 4) + diag(0.5, 4))
    model <- rjags::jags.model(
      textConnection(paste0(
        "model{\n", module$syntax,
        "for(i in 1:N){ y[i,1:4] ~ dmnorm.vcov(zero[1:4], Omega_R[1:4,1:4]) }\n}\n"
      )),
      data = list(y = y, N = nrow(y), zero = rep(0, 4)),
      n.chains = 1, n.adapt = 0, quiet = TRUE,
      inits = list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 11)
    )
    samplers <- rjags::list.samplers(model)
    primitive <- vapply(samplers, function(nodes) any(grepl("Omega_lkj_u", nodes, fixed = TRUE)), logical(1))
    expect_setequal(unlist(samplers[primitive]), module$primitive_names)
    expect_true(all(lengths(samplers[primitive]) == 1L))
    expect_true(all(names(samplers)[primitive] == "base::RealSlicer"))
    expect_no_error(stats::update(model, 500, progress.bar = "none"))
  }
})

test_that("JAGS LKJ posterior for bivariate normal correlation matches grid oracle", {
  eta <- 0.8
  y <- .lkj_bivariate_data()
  oracle <- .lkj_bivariate_grid_oracle(y = y, eta = eta)
  rho <- .fit_jags_lkj_bivariate_posterior(y = y, eta = eta, sample = 5000)

  expect_equal(mean(rho), oracle$mean, tolerance = 0.06)
  expect_equal(stats::sd(rho), oracle$sd, tolerance = 0.06)
  expect_equal(
    stats::quantile(rho, c(.1, .5, .9), names = FALSE),
    oracle$quantiles,
    tolerance = 0.09
  )
})

test_that("JAGS_fit loads BayesTools module from generated LKJ metadata", {
  module <- JAGS_lkj_corr_cholesky(
    name = "Omega",
    K = 2,
    eta = 1.1,
    include_correlation = TRUE
  )

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = paste0("model{\n", module$syntax, "\n}"),
    prior_list = NULL,
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    add_parameters = module$monitor,
    required_packages = module$required_packages,
    jags_modules = module$jags_module,
    silent = TRUE,
    seed = 10
  ))

  expect_s3_class(fit, "BayesTools_fit")
  expect_equal(attr(fit, "jags_modules"), "BayesTools")
  samples <- as.matrix(fit$mcmc)
  expect_true(all(c("Omega_L[1,1]", "Omega_R[1,2]") %in% colnames(samples)))
})

test_that("one-dimensional generated LKJ module runs without native LKJ calls", {
  module <- JAGS_lkj_corr_cholesky(
    name = "One",
    K = 1,
    eta = 3,
    include_correlation = TRUE
  )

  expect_false(grepl("dbt_lkj_cpc", module$syntax, fixed = TRUE))
  expect_false(grepl("bt_lkj_cholesky", module$syntax, fixed = TRUE))
  expect_false(grepl("bt_lkj_corr", module$syntax, fixed = TRUE))

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = paste0("model{\n", module$syntax, "\n}"),
    prior_list = NULL,
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    add_parameters = module$monitor,
    required_packages = module$required_packages,
    jags_modules = module$jags_module,
    silent = TRUE,
    seed = 11
  ))

  expect_s3_class(fit, "BayesTools_fit")
  samples <- as.matrix(fit$mcmc)
  expect_equal(colnames(samples), c("One_L", "One_R"))
  expect_equal(unname(samples[, "One_L"]), rep(1, nrow(samples)))
  expect_equal(unname(samples[, "One_R"]), rep(1, nrow(samples)))
})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-marglik.R
# ============================================================================ #


# ============================================================================ #
# TEST FILE: JAGS Marginal Likelihood Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for JAGS marginal likelihood computation functions.
#   Uses simple models where the log marginal likelihood is known to be 0
#   (for prior samples, the marginal likelihood for any proper prior is 1).
#
# DEPENDENCIES:
#   - rjags: For JAGS model fitting
#   - bridgesampling: For marginal likelihood computation
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("rjags")
#   - Note: Creates fresh models, does not need pre-fitted models
#
# MODELS/FIXTURES:
#   - Creates models with known analytical marginal likelihoods for validation
#
# TAGS: @evaluation, @JAGS, @marginal-likelihood
# ============================================================================ #

# Load common test helpers

make_bridge_random_fixture <- function(formula = ~ 1 + x + us(1 + x | id),
                                       data = NULL,
                                       prior_list = NULL,
                                       prior_random_list = NULL){

  if(is.null(data)){
    data <- data.frame(
      x = c(-1, 0, 1, 2),
      id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
    )
  }
  if(is.null(prior_list)){
    prior_list <- list(intercept = prior("normal", list(0, 1)))
    if("x" %in% all.vars(formula)){
      prior_list$x <- prior("normal", list(0, 1))
    }
  }
  if(is.null(prior_random_list)){
    prior_random_list <- prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1, include_correlation = FALSE)
      )
    )
  }

  list(
    data = data,
    formula = formula,
    prior_list = prior_list,
    prior_random_list = prior_random_list,
    result = JAGS_formula(
      formula = formula,
      parameter = "mu",
      data = data,
      prior_list = prior_list,
      prior_random = prior_random_list
    )
  )
}

expect_formula_random_prior_only_bridge <- function(formula, data, prior_list,
                                                    prior_random_list,
                                                    n_iter = 8000,
                                                    tolerance = 0.08,
                                                    seed = 1,
                                                    maxiter = 2000){

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(seed)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = seed),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- .bt_rjags_samples(model, monitor, n_iter, 1L, "none")
  completed_scale <- .bt_formula_scale_finalize(formula_result$formula_scale,
    formula_result$formula_design, formula_result$prior_list, formula_result$data, owner_scope = "fit")
  formula_result$formula_design$formula_scale <- completed_scale
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)
  attr(samples, "formula_scale") <- list(mu = completed_scale)

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(),
    prior_list = NULL,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = data),
    formula_prior_list = list(mu = prior_list),
    formula_random_prior_list = list(mu = prior_random_list),
    maxiter = maxiter,
    seed = 1
  )

  expect_s3_class(marglik, "BayesTools_marglik")
  expect_equal(marglik$logml, 0, tolerance = tolerance)

  invisible(list(
    marglik = marglik,
    formula_result = formula_result,
    samples = samples
  ))
}

test_that("direct composed bias priors support bridge-sampling helpers", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 2)))
  phacking  <- prior_phacking(form = "linear", alpha = prior("beta", list(2, 3)))
  bias      <- prior_bias(selection = selection, phacking = phacking)

  posterior <- matrix(
    c(
      .625, .2,
      .65625, .4
    ),
    ncol = 2,
    byrow = TRUE
  )
  colnames(posterior) <- c("omega[2]", "alpha")

  prepared <- JAGS_bridgesampling_posterior(posterior, list(bias = bias))
  expect_equal(colnames(prepared), c("omega[2]", "alpha"))
  expect_equal(attr(prepared, "lb"), c("omega[2]" = 0, "alpha" = 0))
  expect_equal(attr(prepared, "ub"), c("omega[2]" = 1, "alpha" = 1))

  samples <- posterior[1, ]
  expected_prior_density <-
    stats::dbeta(samples[["omega[2]"]], 2, 1, log = TRUE) +
    stats::dbeta(samples[["alpha"]], 2, 3, log = TRUE)
  expect_equal(
    JAGS_marglik_priors(samples, list(bias = bias)),
    expected_prior_density,
    tolerance = 1e-12
  )

  parameters <- JAGS_marglik_parameters(samples, list(bias = bias))
  constants <- phack_backend_constants(phacking$form, phacking$source, phacking$destination, target = phacking$target)
  expect_equal(parameters$omega, c(1, samples[["omega[2]"]]))
  expect_equal(parameters$alpha, samples[["alpha"]])
  expect_equal(
    parameters$pi_null,
    samples[["alpha"]] * constants$pi_null_per_alpha,
    tolerance = 1e-12
  )
  expect_equal(
    parameters$beta_null,
    samples[["alpha"]] * constants$beta_null_per_alpha,
    tolerance = 1e-12
  )
})

test_that("p-hacking bridge helpers support point and inverse-gamma alpha priors", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 2)))

  point_alpha <- prior_phacking(
    form  = "linear",
    alpha = prior("point", list(.25))
  )
  point_bias <- prior_bias(selection = selection, phacking = point_alpha)
  point_posterior <- matrix(
    c(
      .625,
      .65625
    ),
    ncol = 1,
    byrow = TRUE
  )
  colnames(point_posterior) <- "omega[2]"

  point_prepared <- JAGS_bridgesampling_posterior(point_posterior, list(bias = point_bias))
  expect_equal(colnames(point_prepared), "omega[2]")
  expect_equal(attr(point_prepared, "lb"), c("omega[2]" = 0))
  expect_equal(attr(point_prepared, "ub"), c("omega[2]" = 1))

  point_samples <- point_posterior[1, ]
  expect_equal(
    JAGS_marglik_priors(point_samples, list(bias = point_bias)),
    stats::dbeta(point_samples[["omega[2]"]], 2, 1, log = TRUE),
    tolerance = 1e-12
  )

  point_parameters <- JAGS_marglik_parameters(point_samples, list(bias = point_bias))
  point_constants <- phack_backend_constants(
    point_alpha$form, point_alpha$source, point_alpha$destination,
    target = point_alpha$target
  )
  expect_equal(point_parameters$alpha, .25)
  expect_equal(
    point_parameters$pi_null,
    .25 * point_constants$pi_null_per_alpha,
    tolerance = 1e-12
  )
  expect_equal(
    point_parameters$beta_null,
    .25 * point_constants$beta_null_per_alpha,
    tolerance = 1e-12
  )

  invgamma_alpha_prior <- prior("invgamma", list(shape = 3, scale = .4), list(0, 1))
  invgamma_alpha <- prior_phacking(
    form  = "linear",
    alpha = invgamma_alpha_prior
  )
  invgamma_bias <- prior_bias(selection = selection, phacking = invgamma_alpha)
  invgamma_posterior <- matrix(
    c(
      .625, 0.4,
      .65625, 0.5
    ),
    ncol = 2,
    byrow = TRUE
  )
  colnames(invgamma_posterior) <- c("omega[2]", "alpha")

  invgamma_prepared <- JAGS_bridgesampling_posterior(invgamma_posterior, list(bias = invgamma_bias))
  expect_equal(colnames(invgamma_prepared), c("omega[2]", "alpha"))
  expect_equal(attr(invgamma_prepared, "lb"), c("omega[2]" = 0, "alpha" = 0))
  expect_equal(attr(invgamma_prepared, "ub"), c("omega[2]" = 1, "alpha" = 1))

  invgamma_samples <- invgamma_posterior[1, ]
  expected_invgamma_prior_density <-
    stats::dbeta(invgamma_samples[["omega[2]"]], 2, 1, log = TRUE) +
    lpdf(invgamma_alpha_prior, invgamma_samples[["alpha"]])
  expect_equal(
    JAGS_marglik_priors(invgamma_samples, list(bias = invgamma_bias)),
    expected_invgamma_prior_density,
    tolerance = 1e-12
  )

  invgamma_parameters <- JAGS_marglik_parameters(invgamma_samples, list(bias = invgamma_bias))
  invgamma_constants <- phack_backend_constants(
    invgamma_alpha$form, invgamma_alpha$source, invgamma_alpha$destination,
    target = invgamma_alpha$target
  )
  expect_equal(invgamma_parameters$alpha, invgamma_samples[["alpha"]])
  expect_equal(
    invgamma_parameters$pi_null,
    invgamma_samples[["alpha"]] * invgamma_constants$pi_null_per_alpha,
    tolerance = 1e-12
  )
  expect_equal(
    invgamma_parameters$beta_null,
    invgamma_samples[["alpha"]] * invgamma_constants$beta_null_per_alpha,
    tolerance = 1e-12
  )

  # BayesTools 0.3.0 'inv_' precision coordinates are not read
  legacy_invgamma_posterior <- invgamma_posterior
  colnames(legacy_invgamma_posterior)[2] <- "inv_alpha"
  legacy_invgamma_posterior[, "inv_alpha"] <- 1 / legacy_invgamma_posterior[, "inv_alpha"]
  expect_error(
    JAGS_bridgesampling_posterior(legacy_invgamma_posterior, list(bias = invgamma_bias)),
    "'posterior' does not contain all of the parameters",
    fixed = TRUE
  )
})

test_that("bias mixtures fail explicitly in bridge sampling and give the active branch's parameters", {

  bias <- prior_mixture(list(
    prior_none(),
    prior_phacking(form = "linear")
  ))
  samples <- c("bias_indicator" = 1, "alpha" = .2)
  posterior <- matrix(samples, nrow = 1)

  # Bridge sampling cannot sample the discrete bias indicator.
  expect_error(
    JAGS_bridgesampling_posterior(posterior, list(bias = bias)),
    "bias mixture priors"
  )
  expect_error(
    JAGS_marglik_priors(samples, list(bias = bias)),
    "bias mixture priors"
  )

  # JAGS_marglik_parameters() returns the parameters of the branch that
  # 'bias_indicator' selects: the weights on the global bins (unit weights
  # outside the weight-function branch), the PET term (0 outside its branch),
  # and the p-hacking parameters (0 outside the p-hacking branch). The
  # expected values are computed by hand from the branch definitions.
  phacking <- prior_phacking(form = "linear")
  mixture <- prior_mixture(list(
    prior_none(),
    prior_PET("normal", list(0, 1)),
    prior_weightfunction("one-sided", c(0.025, 0.05), wf_cumulative(c(1, 1, 1))),
    prior_bias(phacking = phacking)
  ))
  constants <- phack_backend_constants(
    phacking$form, phacking$source, phacking$destination,
    target = phacking$target
  )
  draw <- c(
    "PET_1" = .37,
    "eta_component_3[1]" = 1, "eta_component_3[2]" = 2, "eta_component_3[3]" = 5,
    "alpha_component_4" = .2
  )
  expected <- list(
    list(omega = c(1, 1, 1), PET = 0, alpha = 0, pi_null = 0, beta_null = 0),
    list(omega = c(1, 1, 1), PET = .37, alpha = 0, pi_null = 0, beta_null = 0),
    # eta / sum(eta) = (1, 2, 5) / 8; omega[j] = sum of the shares j..3
    list(omega = c(1, 7 / 8, 5 / 8), PET = 0, alpha = 0, pi_null = 0, beta_null = 0),
    list(omega = c(1, 1, 1), PET = 0, alpha = .2,
         pi_null = .2 * constants$pi_null_per_alpha,
         beta_null = .2 * constants$beta_null_per_alpha)
  )
  for(branch in seq_along(expected)){
    parameters <- JAGS_marglik_parameters(
      c("bias_indicator" = branch, draw),
      list(bias = mixture)
    )
    expect_equal(parameters[names(expected[[branch]])], expected[[branch]], tolerance = 1e-15)
  }
  expect_error(
    JAGS_marglik_parameters(draw, list(bias = mixture)),
    "'samples' does not contain all monitored bias-mixture parameters of 'bias'.",
    fixed = TRUE
  )
})

# This file tests the JAGS marginal likelihood computation functions
# It uses simple models where the log marginal likelihood is known to be 0
# (for prior samples, the marginal likelihood for any proper prior is 1, log(1) = 0)
# More complex consistency tests (e.g., including formulas etc part of `test-00-model-fits.R`)

test_that("JAGS model functions work (simple)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("normal", list(0, 1), list(1, Inf)),
    p3  = prior("lognormal", list(0, .5)),
    p4  = prior("t", list(0, .5, 5)),
    p5  = prior("Cauchy", list(1, 0.1), list(-10, 0)),
    p6  = prior("gamma", list(2, 1)),
    p7  = prior("invgamma", list(3, 2), list(1, 3)),
    p8  = prior("exp", list(1.5)),
    p9  = prior("beta", list(3, 2)),
    p10 = prior("uniform", list(1, 5)),
    PET = prior_PET("normal", list(0, 1)),
    PEESE = prior_PEESE("gamma", list(1, 1))
    #p12 = prior("bernoulli", list(0.75)) discrete priors are not supported with bridgesampling
  )
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

# skip the rest as it takes too long

test_that("JAGS model functions work (vector)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    p1  = prior("mnormal", list(mean = 0, sd = 1, K = 3),),
    p2  = prior("mcauchy", list(location = 0, scale = 1.5, K = 2)),
    p3  = prior("mt",      list(location = 2, scale = 0.5, df = 5, K = 2))
  )
  log_posterior <- STANDARD_LOG_POSTERIOR

  # Absolute tolerances around the analytic logml = 0. The mcauchy (p2)
  # estimate shows no bias but a larger Monte Carlo spread, because the fitted
  # normal bridge proposal takes the sample covariance of Cauchy draws. With
  # 2 chains x 10000 draws its logml had SD 0.043 over independent proposal
  # streams on fixed JAGS draws (mean -0.004, SE 0.004, n = 120) and SD 0.050
  # over iid draws (mean -0.001, SE 0.004, n = 200), so 0.2 is ~4 SD. Bridge
  # sampling is seeded per prior: JAGS does not draw from R's stream, so an
  # unseeded bridge after 'set.seed(1)' used the same proposal variates in every
  # run, and the runs shared one realization of the proposal error (mean
  # -0.026 over the same 120 JAGS draws) that an earlier review had read as a
  # bias. BayesToolsVerse logs: .work/logs/rscript-20260924-075719-9bf9f112
  # (fixed JAGS draws) and .work/logs/rscript-20260924-080645-95104665 (iid
  # draws; BayesTools and upstream bridgesampling agree to 3.6e-15).
  logml_tolerance <- c(p1 = 5e-2, p2 = 0.2, p3 = 5e-2)

  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 10000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior, seed = i)
    expect_equal(
      marglik$logml, 0,
      tolerance = logml_tolerance[[names(prior_list)]],
      info = names(prior_list)
    )
  }

})

test_that("JAGS model functions work (factor)", {

  skip_if_not_installed("rjags")
  all_priors   <- list(
    p1  = prior_factor("mnorm", list(mean = 0, sd = 1),    contrast = "orthonormal"),
    p2  = prior_factor("beta",  list(alpha = 1, beta = 1), contrast = "treatment"),
    p3  = prior_factor("beta",  list(alpha = 2, beta = 2), contrast = "treatment"),
    p4  = prior_factor("gamma",   list(shape = 2, rate = 3), contrast = "independent"),
    p5  = prior_factor("uniform", list(a = -0.5, b = 1.5),   contrast = "independent"),
    p6  = prior_factor("mnorm", list(mean = 0, sd = 1),     contrast = "meandif")
  )

  # add levels
  attr(all_priors[[1]], "levels") <- 3
  attr(all_priors[[2]], "levels") <- 2
  attr(all_priors[[3]], "levels") <- 3
  attr(all_priors[[4]], "levels") <- 1
  attr(all_priors[[5]], "levels") <- 3
  attr(all_priors[[6]], "levels") <- 3
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 10000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("JAGS marginal-likelihood helpers reject spike-and-slab priors explicitly", {

  prior_list <- list(
    theta = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("beta", list(1, 1))
    )
  )
  posterior <- matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "theta"))
  samples <- c(theta = 0, theta_inclusion = .5)

  expect_error(
    JAGS_bridgesampling_posterior(posterior, prior_list),
    "spike and slab priors is not implemented"
  )
  expect_error(
    JAGS_marglik_priors(samples, prior_list),
    "prior mixture priors is not implemented"
  )
  # JAGS_marglik_parameters() evaluates the spike-and-slab node itself (see
  # test-JAGS-deterministic-nodes.R); these samples lack its components.
  expect_error(
    JAGS_marglik_parameters(samples, prior_list),
    "'samples' does not contain all monitored spike-and-slab parameters of 'theta'.",
    fixed = TRUE
  )
})

test_that("JAGS model functions work (weightfunctions)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    prior_weightfunction("one-sided", c(.05), wf_cumulative(c(1, 1))),
    prior_weightfunction("one-sided", c(.05, 0.10), wf_cumulative(c(1, 2, 3))),
    prior_weightfunction("one-sided", c(.05, 0.60), wf_independent(prior("beta", list(1, 1)))),
    prior_weightfunction("two-sided", c(.05), wf_cumulative(c(1, 1)))
  )
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    names(prior_list) <- "omega"
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("JAGS model functions work (spikes)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    p1    = prior("spike", list(1)),
    p2.2  = prior_factor("spike", list(location = 2), contrast = "treatment"),
    p3.2  = prior_factor("spike", list(location = 3), contrast = "independent"),
    p4.2  = prior_factor("spike", list(location = 0), contrast = "orthonormal"),
    p5.2  = prior_factor("spike", list(location = 0), contrast = "meandif"),
    p2.5  = prior_factor("spike", list(location = 2), contrast = "treatment"),
    p3.5  = prior_factor("spike", list(location = 3), contrast = "independent"),
    p4.5  = prior_factor("spike", list(location = 0), contrast = "orthonormal"),
    p5.5  = prior_factor("spike", list(location = 0), contrast = "meandif")
  )
  attr(all_priors$p2.2, "levels") <- 2
  attr(all_priors$p3.2, "levels") <- 2
  attr(all_priors$p4.2, "levels") <- 2
  attr(all_priors$p5.2, "levels") <- 2
  attr(all_priors$p2.5, "levels") <- 2
  attr(all_priors$p3.5, "levels") <- 2
  attr(all_priors$p4.5, "levels") <- 2
  attr(all_priors$p5.5, "levels") <- 2
  nuisance_prior <- list(sigma = prior("normal", list(0, 1)))
  log_posterior <- STANDARD_LOG_POSTERIOR


  for(i in seq_along(all_priors)){
    prior_list   <- c(all_priors[i], nuisance_prior)
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("bridge sampling object function works",{

  marglik0 <- bridgesampling_object()
  marglik1 <- bridgesampling_object(1)

  expect_equal(marglik0$logml, -Inf)
  expect_equal(marglik1$logml, 1)
  expect_s3_class(marglik0, "BayesTools_marglik")

})

test_that("JAGS marglik with formula works", {

  # Test marginal likelihood computation with formula interface
  # Uses intercept-only formula with various priors
  # When sampling from prior and computing marglik, the result should be ~0 (log(1))

  skip_if_not_installed("rjags")

  # Simple data for the formula
  set.seed(1)
  df_test <- data.frame(x = rnorm(10))
  log_posterior <- STANDARD_LOG_POSTERIOR

  # Create formula prior list with intercept only
  prior_list <- list(
    "intercept" = prior("gamma",  list(2, 2)),
    "x"         = prior("normal", list(0, 1))
  )

  # Process formula to get JAGS syntax
  formula_result <- JAGS_formula(~ 1 + x, parameter = "mu", data = df_test, prior_list = prior_list)

  # Build JAGS model with formula priors
  model_syntax <- JAGS_add_priors("model{}", formula_result$prior_list)
  monitor      <- JAGS_to_monitor(formula_result$prior_list)
  inits        <- JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1)

  # Sample from prior using JAGS
  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
  completed_scale <- .bt_formula_scale_finalize(formula_result$formula_scale,
    formula_result$formula_design, formula_result$prior_list, formula_result$data, owner_scope = "fit")
  formula_result$formula_design$formula_scale <- completed_scale
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)
  attr(samples, "formula_scale") <- list(mu = completed_scale)

  # Compute marginal likelihood using formula interface
  marglik <- JAGS_bridgesampling(
    fit                = samples,
    log_posterior      = log_posterior,
    data               = list(),
    prior_list         = NULL,
    formula_list       = list(mu = ~ 1 + x),
    formula_data_list  = list(mu = df_test),
    formula_prior_list = list(mu = prior_list)
  )

  expect_equal(marglik$logml, 0, tolerance = 1e-3)
})

test_that("JAGS formula expressions replay sampled indexed parameters", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  set.seed(41)
  mapping_id <- rep(1:2, each = 4)
  model_data <- list(
    y = stats::rnorm(8, c(-0.5, 0.75)[mapping_id], 1),
    mapping_id = mapping_id,
    N = 8L
  )
  formula_data <- data.frame(row_id = seq_len(model_data$N))
  prior_list <- list(
    mu_id = prior_factor(
      "normal",
      list(0, 1),
      contrast = "independent"
    )
  )
  prior_list$mu_id <- prior_factor_levels(prior_list$mu_id, 2L)
  formula <- ~ expression(mu_id[mapping_id[i]])
  fit <- JAGS_fit(
    model_syntax = paste0(
      "model{\n",
      "for(i in 1:N){\n",
      "  y[i] ~ dnorm(mu[i], 1)\n",
      "}\n",
      "}"
    ),
    data = model_data,
    prior_list = prior_list,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = formula_data),
    formula_prior_list = list(
      mu = list(intercept = prior("point", list(0)))
    ),
    chains = 2,
    adapt = 100,
    burnin = 100,
    sample = 300,
    seed = 41
  )
  design <- JAGS_formula_design(fit, "mu")
  expect_equal(design$expression_data$mapping_id, mapping_id)
  expect_identical(
    design$expression_specs[[1L]]$parameter_dependencies,
    "mu_id"
  )

  posterior <- as.matrix(BayesTools:::.fit_to_posterior(fit))
  expected <- t(posterior[, c("mu_id[1]", "mu_id[2]")][, mapping_id])
  fitted_prediction <- JAGS_evaluate_formula(
    fit,
    formula = NULL,
    parameter = "mu"
  )
  expect_equal(unname(fitted_prediction), unname(expected))

  new_mapping <- c(2L, 1L, 2L)
  new_prediction <- JAGS_evaluate_formula(
    fit,
    formula = formula,
    parameter = "mu",
    data = data.frame(
      row_id = seq_along(new_mapping),
      mapping_id = new_mapping
    )
  )
  expect_equal(
    unname(new_prediction),
    unname(t(posterior[, c("mu_id[1]", "mu_id[2]")][, new_mapping]))
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters$mu, 1, log = TRUE))
  }
  marglik <- JAGS_bridgesampling(
    fit,
    log_posterior = log_posterior,
    data = model_data,
    maxiter = 1000,
    seed = 1
  )
  expect_true(is.finite(marglik$logml))
  changed_data <- model_data
  changed_data$mapping_id <- rev(changed_data$mapping_id)
  expect_error(
    JAGS_bridgesampling(
      fit,
      log_posterior = log_posterior,
      data = changed_data,
      maxiter = 100,
      seed = 1
    ),
    "conflict with the fitted source snapshot",
    fixed = TRUE
  )

  expect_error(
    JAGS_fit(
      model_syntax = "model{ theta <- 0 }",
      data = list(N = 2L),
      prior_list = list(),
      formula_list = list(mu = ~ expression(theta)),
      formula_data_list = list(mu = data.frame(row_id = 1:2)),
      formula_prior_list = list(
        mu = list(intercept = prior("normal", list(0, 1)))
      ),
      chains = 1,
      adapt = 50,
      burnin = 50,
      sample = 100
    ),
    "unknown replay dependency 'theta'",
    fixed = TRUE
  )

  sparse_fit <- JAGS_fit(
    model_syntax = paste0(
      "model{\n",
      "theta[1] ~ dnorm(0, 1)\n",
      "theta[2] ~ dnorm(0, 1)\n",
      "for(i in 1:N){\n",
      "  y[i] ~ dnorm(sparse_mu[i], 1)\n",
      "}\n",
      "}"
    ),
    data = list(y = model_data$y, N = model_data$N),
    formula_list = list(sparse_mu = ~ expression(theta[2])),
    formula_data_list = list(sparse_mu = formula_data),
    formula_prior_list = list(
      sparse_mu = list(intercept = prior("normal", list(0, 1)))
    ),
    add_parameters = "theta[2]",
    chains = 1,
    adapt = 100,
    burnin = 100,
    sample = 100,
    seed = 42
  )
  sparse_posterior <- as.matrix(BayesTools:::.fit_to_posterior(sparse_fit))
  expect_identical(
    colnames(sparse_posterior),
    c("sparse_mu_intercept", "theta[2]")
  )
  sparse_prediction <- JAGS_evaluate_formula(
    sparse_fit,
    parameter = "sparse_mu"
  )
  expect_equal(
    unname(sparse_prediction),
    matrix(
      sparse_posterior[, "sparse_mu_intercept"] +
        sparse_posterior[, "theta[2]"],
      nrow = model_data$N,
      ncol = nrow(sparse_posterior),
      byrow = TRUE
    )
  )
})

test_that("JAGS marglik with exp(intercept) formula works", {

  # Test marginal likelihood computation with formula interface
  # Uses intercept-only formula with various priors
  # When sampling from prior and computing marglik, the result should be ~0 (log(1))

  skip_if_not_installed("rjags")

  # Simple data for the formula
  set.seed(1)
  df_test <- data.frame(x = rnorm(10))
  log_posterior <- STANDARD_LOG_POSTERIOR

  # Create formula prior list with intercept only
  prior_list <- list(
    "intercept" = prior("gamma",  list(2, 2)),
    "x"         = prior("normal", list(0, 1))
  )

  # Process formula to get JAGS syntax
  formula <- ~ 1 + x
  attr(formula, "log(intercept)") <- TRUE
  formula_result <- JAGS_formula(formula, parameter = "mu", data = df_test, prior_list = prior_list)
  expect_equal(formula_result$formula_syntax, "for(i in 1:N_mu){\n  mu[i] = log(mu_intercept) + mu_x * mu_data_x[i]\n}\n")

  # Build JAGS model with formula priors
  model_syntax <- JAGS_add_priors("model{}", formula_result$prior_list)
  monitor      <- JAGS_to_monitor(formula_result$prior_list)
  inits        <- JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1)

  # Sample from prior using JAGS
  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
  completed_scale <- .bt_formula_scale_finalize(formula_result$formula_scale,
    formula_result$formula_design, formula_result$prior_list, formula_result$data, owner_scope = "fit")
  formula_result$formula_design$formula_scale <- completed_scale
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)
  attr(samples, "formula_scale") <- list(mu = completed_scale)

  # Compute marginal likelihood using formula interface
  marglik <- JAGS_bridgesampling(
    fit                = samples,
    log_posterior      = log_posterior,
    data               = list(),
    prior_list         = NULL,
    formula_list       = list(mu = formula),
    formula_data_list  = list(mu = df_test),
    formula_prior_list = list(mu = prior_list)
  )

  expect_equal(marglik$logml, 0, tolerance = 1e-3)
})

test_that("JAGS bridgesampling infers formula scaling metadata from fits", {

  df_test <- data.frame(x = c(10, 20, 30))
  prior_list <- list(
    "intercept" = prior("normal", list(0, 1)),
    "x"         = prior("normal", list(0, 1))
  )

  scaled_formula <- JAGS_formula(
    ~ 1 + x,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    formula_scale = list(x = TRUE)
  )
  fit <- matrix(c(0, 1), ncol = 2)
  attr(fit, "formula_scale") <- list(mu = scaled_formula$formula_scale)

  inferred_scale <- BayesTools:::.JAGS_formula_scale_list_from_fit(fit, "mu")
  expect_equal(inferred_scale, list(mu = list(x = TRUE)))

  bridge_formula <- JAGS_formula(
    ~ 1 + x,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    formula_scale = inferred_scale$mu
  )

  expect_equal(bridge_formula$data$mu_data_x, scaled_formula$data$mu_data_x)
})

test_that("JAGS bridgesampling posterior supports add-only parameters", {

  posterior <- matrix(
    c(
      -1, 0.2,
       0, 0.5,
       1, 0.8
    ),
    ncol = 2,
    byrow = TRUE
  )
  colnames(posterior) <- c("x", "prob")

  info <- BayesTools:::.JAGS_bridgesampling_posterior_info(NULL)
  expect_type(info, "character")
  expect_equal(length(info), 0L)
  expect_equal(attr(info, "lb"), numeric())
  expect_equal(attr(info, "ub"), numeric())

  result_null <- JAGS_bridgesampling_posterior(
    posterior = posterior,
    prior_list = NULL,
    add_parameters = c("x", "prob"),
    add_bounds = list(
      lb = c(x = -Inf, prob = 0),
      ub = c(x = Inf, prob = 1)
    )
  )
  expect_equal(colnames(result_null), c("x", "prob"))
  expect_equal(attr(result_null, "lb"), c(x = -Inf, prob = 0))
  expect_equal(attr(result_null, "ub"), c(x = Inf, prob = 1))

  result_empty <- JAGS_bridgesampling_posterior(
    posterior = posterior,
    prior_list = list(),
    add_parameters = "x",
    add_bounds = list(lb = c(x = -Inf), ub = c(x = Inf))
  )
  expect_equal(colnames(result_empty), "x")
  expect_equal(attr(result_empty, "lb"), c(x = -Inf))
  expect_equal(attr(result_empty, "ub"), c(x = Inf))
})

test_that("JAGS bridgesampling passes requested bridge context to callback", {

  skip_if_not_installed("bridgesampling")
  skip_if_not_installed("coda")

  set.seed(1)
  posterior <- matrix(
    rnorm(2000),
    ncol = 1,
    dimnames = list(NULL, "mu")
  )
  posterior <- coda::as.mcmc(posterior)
  seen <- new.env(parent = emptyenv())
  seen$calls <- 0L
  seen$problems <- character()
  # The callback runs for every draw of every bridge iteration. It checks the
  # context on each call and records what fails; the expectations are made once
  # after the run, not one per call.
  log_posterior <- function(parameters, data, bridge_context){
    seen$calls <- seen$calls + 1L
    if(!exists("context", envir = seen, inherits = FALSE)){
      seen$context <- bridge_context
      seen$mu <- parameters$mu
    }
    if(!inherits(bridge_context, "BayesTools_bridge_context")){
      seen$problems <- c(seen$problems, "the context is not a BayesTools_bridge_context")
    }
    if(!"mu" %in% names(bridge_context$state)){
      seen$problems <- c(seen$problems, "the context state has no 'mu'")
    }
    differences <- waldo::compare(
      bridge_context$nodes[["mu"]], parameters$mu,
      tolerance = testthat::testthat_tolerance()
    )
    if(length(differences) > 0L){
      seen$problems <- c(seen$problems, paste("node 'mu' differs from the parameter:", differences))
    }
    0
  }

  marglik <- JAGS_bridgesampling(
    fit = posterior,
    log_posterior = log_posterior,
    data = list(),
    prior_list = list(mu = prior("normal", list(0, 1))),
    bridge_context = TRUE,
    maxiter = 1000
  )

  expect_s3_class(marglik, "BayesTools_marglik")
  expect_s3_class(seen$context, "BayesTools_bridge_context")
  expect_equal(seen$context$nodes[["mu"]], seen$mu)
  # the callback ran for the draws of the bridge, and every call passed its checks
  expect_gte(seen$calls, nrow(posterior))
  expect_identical(seen$problems, character())
})

test_that("JAGS bridgesampling validates rebuilt formula random design metadata", {

  fixture <- make_bridge_random_fixture()
  fitted <- list(mu = fixture$result$formula_design)
  rebuilt <- list(mu = fixture$result$formula_design)

  expect_silent(BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = fitted[["mu"]],
      rebuilt = rebuilt[["mu"]]
    ))

  changed <- rebuilt
  changed$mu$random_effects[[1]]$structure <- "diag"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = fitted[["mu"]],
      rebuilt = changed[["mu"]]
    ),
    "covariance structure",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$structure <- NULL
  changed$mu$random_effects[[1]]$covariance <- "us"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = fitted[["mu"]],
      rebuilt = changed[["mu"]]
    ),
    "missing canonical 'random_term\\$structure'"
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$homogeneous_sd <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = fitted[["mu"]],
      rebuilt = changed[["mu"]]
    ),
    "missing canonical 'random_term\\$homogeneous_sd'"
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$correlation <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = fitted[["mu"]],
      rebuilt = changed[["mu"]]
    ),
    "missing canonical 'random_term\\$correlation'"
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$block_name <- "other"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = fitted[["mu"]],
      rebuilt = changed[["mu"]]
    ),
    "block names",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$group_levels <- rev(changed$mu$random_effects[[1]]$group_levels)
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = fitted[["mu"]],
      rebuilt = changed[["mu"]]
    ),
    "group levels",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$column_names[1] <- "changed"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = fitted[["mu"]],
      rebuilt = changed[["mu"]]
    ),
    "column names",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$n_columns <- changed$mu$random_effects[[1]]$n_columns + 1L
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = fitted[["mu"]],
      rebuilt = changed[["mu"]]
    ),
    "group or column counts",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$model_matrix <- changed$mu$random_effects[[1]]$model_matrix[-1, , drop = FALSE]
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = fitted[["mu"]],
      rebuilt = changed[["mu"]]
    ),
    "model matrix shape or columns",
    fixed = TRUE
  )

  changed <- rebuilt
  changed$mu$random_effects[[1]]$group_map <- rev(changed$mu$random_effects[[1]]$group_map)
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = fitted[["mu"]],
      rebuilt = changed[["mu"]]
    ),
    "group map",
    fixed = TRUE
  )

  structured_fixture <- make_bridge_random_fixture(
    formula = ~ 1 + ar1(idx | id),
    data = data.frame(
      idx = factor(c("t1", "t2", "t1", "t2"), levels = c("t1", "t2")),
      id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior("normal", list(0, 0.5))
      )
    )
  )
  structured_fitted <- list(mu = structured_fixture$result$formula_design)
  structured_changed <- structured_fitted
  structured_changed$mu$random_effects[[1]]$structured_index$label <- "changed"
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = structured_fitted[["mu"]],
      rebuilt = structured_changed[["mu"]]
    ),
    "structured random-effect index metadata",
    fixed = TRUE
  )

  allocation_data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b", "a", "b"), levels = c("a", "b"))
  )
  allocation_formula <- ~ 1 +
    random(1 | study, name = "study", covariance = "diag") +
    random(1 | drug, name = "drug", covariance = "diag")
  allocation_prior_list <- list(intercept = prior("normal", list(0, 1)))
  allocation_fixture <- make_bridge_random_fixture(
    formula = allocation_formula,
    data = allocation_data,
    prior_list = allocation_prior_list,
    prior_random_list = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  independent_fixture <- make_bridge_random_fixture(
    formula = allocation_formula,
    data = allocation_data,
    prior_list = allocation_prior_list,
    prior_random_list = prior_random(
      study = random_block(sd = prior("gamma", list(2, 2))),
      drug = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_design(
      parameter = "mu",
      fitted = allocation_fixture$result$formula_design,
      rebuilt = independent_fixture$result$formula_design
    ),
    "scale/allocation",
    fixed = TRUE
  )
})

test_that("JAGS bridgesampling errors on fitted/rebuilt random design mismatches", {

  fixture <- make_bridge_random_fixture()
  fit <- coda::mcmc(matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "dummy")))
  completed_scale <- .bt_formula_scale_finalize(fixture$result$formula_scale,
    fixture$result$formula_design, fixture$result$prior_list, fixture$result$data, owner_scope = "fit")
  fixture$result$formula_design$formula_scale <- completed_scale
  attr(fit, "formula_design") <- list(mu = fixture$result$formula_design)
  attr(fit, "formula_scale") <- list(mu = completed_scale)

  mismatch_data <- fixture$data
  mismatch_data$id <- factor(as.character(mismatch_data$id), levels = c("b", "a"))

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = STANDARD_LOG_POSTERIOR,
      data = list(),
      prior_list = NULL,
      formula_list = list(mu = fixture$formula),
      formula_data_list = list(mu = mismatch_data),
      formula_prior_list = list(mu = fixture$prior_list),
      formula_random_prior_list = list(mu = fixture$prior_random_list),
      maxiter = 10,
      seed = 1
    ),
    "original formula source data differ",
    fixed = TRUE
  )
})

test_that("JAGS bridgesampling supports formula random effects through prior_random", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  if(!isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE))){
    skip("BayesTools JAGS module is not available")
  }

  df_test <- data.frame(
    x = c(-0.5, 0.5),
    id = factor(c("a", "b"), levels = c("a", "b"))
  )
  formula <- ~ 1 + x + us(1 + x | id)
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x         = prior("normal", list(0, 1))
  )
  prior_random_list <- prior_random(
    id = random_block(
      sd = prior("gamma", list(2, 2)),
      cor = prior_lkj(eta = 1, include_correlation = FALSE)
    )
  )

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  expect_true("mu__xREx__id_xRE_CORx_lkj_u[1]" %in% formula_result$add_parameters)

  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(1)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- .bt_rjags_samples(model, monitor, 5000L, 1L, "none")
  completed_scale <- .bt_formula_scale_finalize(formula_result$formula_scale,
    formula_result$formula_design, formula_result$prior_list, formula_result$data, owner_scope = "fit")
  formula_result$formula_design$formula_scale <- completed_scale
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)
  attr(samples, "formula_scale") <- list(mu = completed_scale)

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(),
    prior_list = NULL,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = df_test),
    formula_prior_list = list(mu = prior_list),
    formula_random_prior_list = list(mu = prior_random_list),
    maxiter = 1000,
    seed = 1
  )

  expect_s3_class(marglik, "BayesTools_marglik")
  expect_equal(marglik$logml, 0, tolerance = 0.08)
})

test_that("JAGS bridgesampling supports continuous-time CAR formula random effects", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  df_test <- data.frame(
    time = c(0, 0.5, 2, 0, 0.5, 2),
    id = factor(c("a", "a", "a", "b", "b", "b"), levels = c("a", "b"))
  )
  formula <- ~ 1 + car(0 + time | id)
  prior_list <- list(
    intercept = prior("normal", list(0, 1))
  )
  prior_random_list <- prior_random(
    id = random_block(
      sd = prior("gamma", list(2, 2)),
      cor = prior("normal", list(0, 0.5))
    )
  )

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  expect_equal(formula_result$formula_design$random_effects[[1]]$structure, "car")
  expect_equal(formula_result$formula_design$random_effects[[1]]$correlation$bounds, c(lower = 0, upper = 1))
  expect_true("mu__xREx__id_rho" %in% formula_result$add_parameters)

  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(1)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = 8000,
    quiet = TRUE,
    progress.bar = "none"
  )
  completed_scale <- .bt_formula_scale_finalize(formula_result$formula_scale,
    formula_result$formula_design, formula_result$prior_list, formula_result$data, owner_scope = "fit")
  formula_result$formula_design$formula_scale <- completed_scale
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)
  attr(samples, "formula_scale") <- list(mu = completed_scale)

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(),
    prior_list = NULL,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = df_test),
    formula_prior_list = list(mu = prior_list),
    formula_random_prior_list = list(mu = prior_random_list),
    maxiter = 1000,
    seed = 1
  )

  expect_s3_class(marglik, "BayesTools_marglik")
  expect_equal(marglik$logml, 0, tolerance = 0.08)
})

test_that("centered continuous-time CAR syntax samples sequential conditionals", {

  skip_if_not_installed("rjags")

  df_test <- data.frame(
    time = c(0, 0.5, 2, 0, 0.5, 2),
    id = factor(c("a", "a", "a", "b", "b", "b"), levels = c("a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = df_test,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("point", list(location = 1)),
        cor = prior("normal", list(0, 0.5)),
        monitor = random_monitor(coefficients = TRUE),
        parameterization = "centered"
      )
    )
  )

  expect_false(grepl(
    "dmnorm.vcov",
    formula_result$formula_syntax,
    fixed = TRUE
  ))
  expect_match(
    formula_result$formula_syntax,
    "pexp(-2 * mu__xREx__id_xRE_CAR_LOG_PHIX[2], 1)",
    fixed = TRUE
  )

  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters,
    "mu__xREx__id_xRE_CAR_PHIX",
    "mu__xREx__id_xRE_CAR_INNOV_VARx"
  ))
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 1, seed = 31),
    n.chains = 1,
    n.adapt = 100,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = 250,
    quiet = TRUE,
    progress.bar = "none"
  )
  draws <- as.matrix(samples[[1L]])
  coefficient_columns <- grep(
    "mu__xREx__id_xRE_COEFx[",
    colnames(draws),
    fixed = TRUE
  )
  transition_columns <- c(
    grep(
      "mu__xREx__id_xRE_CAR_PHIX[",
      colnames(draws),
      fixed = TRUE
    ),
    grep(
      "mu__xREx__id_xRE_CAR_INNOV_VARx[",
      colnames(draws),
      fixed = TRUE
    )
  )

  expect_length(coefficient_columns, 6L)
  expect_length(transition_columns, 4L)
  expect_true(all(is.finite(draws[, coefficient_columns, drop = FALSE])))
  expect_true(all(is.finite(draws[, transition_columns, drop = FALSE])))
  expect_true(all(apply(
    draws[, coefficient_columns, drop = FALSE],
    2L,
    stats::sd
  ) > 0))
  expect_true(all(apply(
    draws[, transition_columns, drop = FALSE],
    2L,
    stats::sd
  ) > 0))
})

test_that("JAGS bridgesampling supports Dirichlet variance-allocation random effects", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  df_test <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b", "a", "b"), levels = c("a", "b"))
  )
  formula <- ~ 1 +
    random(1 | study, name = "study", covariance = "diag") +
    random(1 | drug, name = "drug", covariance = "diag")
  prior_list <- list(
    intercept = prior("normal", list(0, 1))
  )
  prior_random_list <- prior_random(
    allocation = random_variance_allocation(name = "allocation",
      sd = prior("gamma", list(2, 2)),
      weights = prior("dirichlet", list(alpha = c(2, 3)))
    )
  )

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  expect_true("mu__xRE_ALLOCx_allocation__weight" %in% names(formula_result$prior_list))
  expect_true("mu__xREx__study_xRE_Zx" %in% formula_result$add_parameters)
  expect_true("mu__xREx__drug_xRE_Zx" %in% formula_result$add_parameters)

  model_syntax <- JAGS_add_priors(
    paste0("model{\n", formula_result$formula_syntax, "\n}"),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(1)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = formula_result$data,
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 1),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = 8000,
    quiet = TRUE,
    progress.bar = "none"
  )
  completed_scale <- .bt_formula_scale_finalize(formula_result$formula_scale,
    formula_result$formula_design, formula_result$prior_list, formula_result$data, owner_scope = "fit")
  formula_result$formula_design$formula_scale <- completed_scale
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)
  attr(samples, "formula_scale") <- list(mu = completed_scale)

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(),
    prior_list = NULL,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = df_test),
    formula_prior_list = list(mu = prior_list),
    formula_random_prior_list = list(mu = prior_random_list),
    maxiter = 1000,
    seed = 1
  )

  expect_s3_class(marglik, "BayesTools_marglik")
  expect_equal(marglik$logml, 0, tolerance = 0.08)
})

test_that("JAGS bridgesampling reconstructs row-indexed external SD sources from values", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  df_test <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b", "a", "b"), levels = c("a", "b")),
    tau_factor = c(0.5, 0.75, 1.0, 1.25)
  )
  formula <- ~ 1 +
    random(1 | study, name = "study", covariance = "diag") +
    random(1 | drug, name = "drug", covariance = "diag")
  prior_list <- list(
    intercept = prior("normal", list(0, 1))
  )
  tau_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      data$tau_factor[seq_len(n_rows)]
    },
    inputs = character()
  )
  prior_random_list <- prior_random(
    allocation = random_variance_allocation(name = "allocation",
      sd_source = random_sd_source(tau_source),
      weights = prior("dirichlet", list(alpha = c(2, 3)))
    )
  )
  no_values_prior_random_list <- prior_random(
    allocation = random_variance_allocation(name = "allocation",
      sd_source = random_sd_source("tau", shape = "row"),
      weights = prior("dirichlet", list(alpha = c(2, 3)))
    )
  )

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = prior_random_list
  )
  no_values_formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = df_test,
    prior_list = prior_list,
    prior_random = no_values_prior_random_list
  )
  model_syntax <- JAGS_add_priors(
    paste0(
      "model{\n",
      "for(i in 1:N_mu){\n",
      "  tau[i] = tau_factor[i]\n",
      "}\n",
      formula_result$formula_syntax,
      "\n}"
    ),
    formula_result$prior_list
  )
  monitor <- unique(c(
    JAGS_to_monitor(formula_result$prior_list),
    formula_result$add_parameters
  ))

  set.seed(11)
  model <- rjags::jags.model(
    file = textConnection(model_syntax),
    data = c(formula_result$data, list(tau_factor = df_test$tau_factor)),
    inits = JAGS_get_inits(formula_result$prior_list, chains = 2, seed = 11),
    n.chains = 2,
    quiet = TRUE
  )
  samples <- rjags::coda.samples(
    model = model,
    variable.names = monitor,
    n.iter = 8000,
    quiet = TRUE,
    progress.bar = "none"
  )
  completed_scale <- .bt_formula_scale_finalize(formula_result$formula_scale,
    formula_result$formula_design, formula_result$prior_list, formula_result$data, owner_scope = "fit")
  formula_result$formula_design$formula_scale <- completed_scale
  attr(samples, "formula_design") <- list(mu = formula_result$formula_design)
  attr(samples, "formula_scale") <- list(mu = completed_scale)
  expect_false(any(grepl("^tau\\[", colnames(as.matrix(samples)))))

  marglik <- JAGS_bridgesampling(
    fit = samples,
    log_posterior = STANDARD_LOG_POSTERIOR,
    data = list(tau_factor = df_test$tau_factor),
    prior_list = NULL,
    maxiter = 1000,
    seed = 1
  )

  expect_s3_class(marglik, "BayesTools_marglik")
  expect_equal(marglik$logml, 0, tolerance = 0.08)

  completed_scale <- .bt_formula_scale_finalize(no_values_formula_result$formula_scale,
    no_values_formula_result$formula_design, no_values_formula_result$prior_list,
    no_values_formula_result$data, owner_scope = "fit")
  no_values_formula_result$formula_design$formula_scale <- completed_scale
  graft_samples <- samples
  attr(graft_samples, "formula_design") <- list(mu = no_values_formula_result$formula_design)
  expect_error(
    JAGS_bridgesampling(
      fit = graft_samples,
      log_posterior = STANDARD_LOG_POSTERIOR,
      data = list(),
      prior_list = NULL,
      formula_list = list(mu = formula),
      formula_data_list = list(mu = df_test),
      formula_prior_list = list(mu = prior_list),
      formula_random_prior_list = list(mu = prior_random_list),
      maxiter = 1000,
      seed = 1
    ),
    "scale/allocation metadata differ",
    fixed = TRUE
  )
})

test_that("JAGS bridgesampling gives unit marglik for prior-only random-effect settings", {

  sd_prior <- prior("gamma", list(2, 2))
  fixed_priors <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  continuous_data <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("a", "a", "b", "b", "c", "c"), levels = c("a", "b", "c"))
  )
  factor_data <- data.frame(
    f = factor(rep(c("a", "b", "c"), 3), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2", "g3"), each = 3), levels = c("g1", "g2", "g3"))
  )

  diag_heterogeneous <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + x + diag(1 + x | id, hom = FALSE),
    data = continuous_data,
    prior_list = fixed_priors,
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        terms = list(x = prior("gamma", list(3, 1)))
      )
    ),
    seed = 10
  )
  expect_false(diag_heterogeneous$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    diag_heterogeneous$formula_result$formula_design$random_effects[[1]]$sd_parameter_names,
    c("mu__xREx__id_intercept", "mu__xREx__id_x")
  )

  id_homogeneous <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + x + id(1 + x | id),
    data = continuous_data,
    prior_list = fixed_priors,
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        monitor = random_monitor(
          latent = TRUE,
          coefficients = TRUE,
          correlation = FALSE
        )
      )
    ),
    seed = 11
  )
  expect_true(id_homogeneous$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    unique(id_homogeneous$formula_result$formula_design$random_effects[[1]]$sd_parameter_names),
    "mu__xREx__id_sd"
  )

  lkj_module <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + x + us(1 + x | id),
    data = continuous_data,
    prior_list = fixed_priors,
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior_lkj(
          eta = 2,
          include_correlation = FALSE
        )
      )
    ),
    seed = 12
  )
  expect_equal(lkj_module$formula_result$jags_modules, "BayesTools")
  expect_true(any(grepl("_xRE_CORx_lkj_u", lkj_module$formula_result$add_parameters, fixed = TRUE)))

  # A three-level mean-difference factor is one exchangeable unit: the LKJ
  # prior restricted to its identity block is normalized over the free
  # primitives, and a diagonal unit has one SD.
  unit_priors <- list(
    intercept = prior("normal", list(0, 1)),
    f = prior_factor("mnormal", list(0, 1), contrast = "meandif")
  )
  lkj_unit <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + f + us(1 + f | id),
    data = factor_data,
    prior_list = unit_priors,
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior_lkj(eta = 1.5, include_correlation = FALSE)
      )
    ),
    seed = 13
  )
  lkj_unit_term <- lkj_unit$formula_result$formula_design$random_effects[[1]]
  expect_identical(lkj_unit_term$correlation$unit_columns, 2:3)
  expect_identical(lkj_unit_term$correlation$primitive_names, paste0("mu__xREx__id_xRE_CORx_lkj_u[", 2:3, "]"))
  diag_unit <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + f + diag(1 + f | id),
    data = factor_data,
    prior_list = unit_priors,
    prior_random_list = prior_random(id = random_block(sd = sd_prior)),
    seed = 14
  )
  expect_identical(
    unique(diag_unit$formula_result$formula_design$random_effects[[1]]$sd_parameter_names),
    c("mu__xREx__id_intercept", "mu__xREx__id_f")
  )

  cs_fisher_z <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + cs(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior("normal", list(0, 0.5))
      )
    ),
    seed = 13
  )
  expect_true(cs_fisher_z$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    cs_fisher_z$formula_result$formula_design$random_effects[[1]]$correlation$rho_scale,
    "fisher_z"
  )
  expect_true("mu__xREx__id_rho_z" %in% names(cs_fisher_z$formula_result$prior_list))

  hcs_logit <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + hcs(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          cor = prior("normal", list(0, 0.5)),
          cor_scale = "logit"
        )
      )
    ),
    seed = 14
  )
  expect_false(hcs_logit$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_true("mu__xREx__id_rho_logit" %in% names(hcs_logit$formula_result$prior_list))
  expect_equal(
    hcs_logit$formula_result$formula_design$random_effects[[1]]$correlation$rho_scale,
    "logit"
  )

  ar1_fisher_z <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + ar1(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior("normal", list(0, 0.5))
      )
    ),
    seed = 15
  )
  expect_equal(ar1_fisher_z$formula_result$formula_design$random_effects[[1]]$structure, "ar1")
  expect_true("mu__xREx__id_rho_z" %in% names(ar1_fisher_z$formula_result$prior_list))

  har_raw <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + har(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          cor = prior("normal", list(0, 0.5), truncation = list(lower = -1, upper = 1)),
          cor_scale = "cor"
        )
      )
    ),
    seed = 16
  )
  expect_equal(har_raw$formula_result$formula_design$random_effects[[1]]$structure, "har")
  expect_false(har_raw$formula_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_true("mu__xREx__id_rho" %in% names(har_raw$formula_result$prior_list))
  expect_false("mu__xREx__id_rho_z" %in% names(har_raw$formula_result$prior_list))

  nested_allocation <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    data = data.frame(
      study = factor(c("s1", "s1", "s2", "s2"), levels = c("s1", "s2")),
      paper = factor(c("p1", "p2", "p1", "p2"), levels = c("p1", "p2")),
      drug = factor(c("d1", "d1", "d2", "d2"), levels = c("d1", "d2"))
    ),
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(nested = "nested", drug = "drug"),
        sd = sd_prior,
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      ),
      random_variance_allocation(
        name = "nested_split",
        parent = allocation_ref("total_re", "nested"),
        terms = c(study = "study", paper = "paper"),
        weights = prior("dirichlet", list(alpha = c(3, 2)))
      )
    ),
    n_iter = 10000,
    seed = 17
  )
  expect_true("mu__xRE_ALLOCx_total_re__weight" %in% names(nested_allocation$formula_result$prior_list))
  expect_true("mu__xRE_ALLOCx_nested_split__weight" %in% names(nested_allocation$formula_result$prior_list))
  expect_equal(
    length(nested_allocation$formula_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$factors),
    2L
  )

  sd_leaf_allocation <- expect_formula_random_prior_only_bridge(
    formula = ~ 1 + hcs(f | id),
    data = factor_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random_list = prior_random(
      allocation = random_variance_allocation(
        name = "leaf_alloc",
        terms = "id",
        target = "sd_component",
        scale = "mean_variance",
        sd = sd_prior,
        weights = prior("dirichlet", list(alpha = c(2, 3, 4)))
      ),
      id = random_block(cor = prior("normal", list(0, 0.5)))
    ),
    n_iter = 10000,
    seed = 18
  )
  expect_true("mu__xRE_ALLOCx_leaf_alloc__weight" %in% names(sd_leaf_allocation$formula_result$prior_list))
  expect_equal(
    sd_leaf_allocation$formula_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$target,
    "sd_component"
  )
  expect_equal(
    sd_leaf_allocation$formula_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$scale,
    "mean_variance"
  )
})

test_that("JAGS formula marglik reconstructs inverse-gamma terms on natural scale", {

  samples <- c(
    "mu_intercept" = 0.5,
    "mu_x"         = 0.25
  )
  formula_output <- JAGS_formula(
    ~ 1 + x, "mu", data.frame(x = c(10, 20)),
    list(
      intercept = prior("invgamma", list(2, 1)),
      x         = prior("invgamma", list(2, 1))
    )
  )
  formula_data_list <- list(mu = formula_output$data)
  formula_prior_list <- list(mu = formula_output$prior_list)
  formula_design_list <- list(mu = formula_output$formula_design)

  parameters <- JAGS_marglik_parameters_formula(
    samples            = samples,
    formula_list       = list(mu = ~ 1 + x),
    formula_data_list  = formula_data_list,
    formula_prior_list = formula_prior_list,
    prior_list_parameters = list(),
    formula_design_list = formula_design_list
  )

  expect_equal(parameters$mu, c(0.5 + 0.25 * 10, 0.5 + 0.25 * 20))

  formula_log_intercept <- ~ 1 + x
  attr(formula_log_intercept, "log(intercept)") <- TRUE
  parameters_log <- JAGS_marglik_parameters_formula(
    samples            = samples,
    formula_list       = list(mu = formula_log_intercept),
    formula_data_list  = formula_data_list,
    formula_prior_list = formula_prior_list,
    prior_list_parameters = list(),
    formula_design_list = formula_design_list
  )

  expect_equal(parameters_log$mu, c(log(0.5) + 0.25 * 10, log(0.5) + 0.25 * 20))

  # BayesTools 0.3.0 'inv_' precision coordinates are not read
  legacy_samples <- c(
    "inv_mu_intercept" = 2,
    "inv_mu_x"         = 4
  )
  expect_error(
    JAGS_marglik_parameters_formula(
      samples            = legacy_samples,
      formula_list       = list(mu = ~ 1 + x),
      formula_data_list  = formula_data_list,
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list(),
      formula_design_list = formula_design_list
    ),
    "'samples' does not contain all monitored formula prior parameters.",
    fixed = TRUE
  )

  # the design-less reconstruction from formula data alone is removed
  expect_error(
    JAGS_marglik_parameters_formula(
      samples            = samples,
      formula_list       = list(mu = ~ 1 + x),
      formula_data_list  = formula_data_list,
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list()
    ),
    "requires the formula design of parameter 'mu' in 'formula_design_list'",
    fixed = TRUE
  )
})


test_that("JAGS marglik reconstructs indexed factor inverse-gamma parameters", {
  theta_prior <- prior_factor("invgamma", list(2, 1), contrast = "independent")
  attr(theta_prior, "levels") <- 2

  parameters <- BayesTools:::.JAGS_marglik_parameters.factor(
    samples = c("theta[1]" = 0.5, "theta[2]" = 0.25),
    prior = theta_prior,
    parameter_name = "theta"
  )

  expect_equal(parameters$theta, c(0.5, 0.25))

  # BayesTools 0.3.0 'inv_' precision coordinates are not read
  expect_error(
    BayesTools:::.JAGS_marglik_parameters.factor(
      samples = c("inv_theta[1]" = 2, "inv_theta[2]" = 4),
      prior = theta_prior,
      parameter_name = "theta"
    ),
    "'samples' does not contain all monitored",
    fixed = TRUE
  )
})


test_that("JAGS formula marglik preserves predictor names containing _data", {
  samples <- c(
    "mu_intercept" = 1,
    "mu_x_data"   = 2
  )
  formula_output <- JAGS_formula(
    ~ 1 + x_data, "mu", data.frame(x_data = c(10, 20)),
    list(
      intercept = prior("normal", list(0, 1)),
      x_data    = prior("normal", list(0, 1))
    )
  )

  parameters <- JAGS_marglik_parameters_formula(
    samples            = samples,
    formula_list       = list(mu = ~ 1 + x_data),
    formula_data_list  = list(mu = formula_output$data),
    formula_prior_list = list(mu = formula_output$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = formula_output$formula_design)
  )

  expect_equal(parameters$mu, c(1 + 2 * 10, 1 + 2 * 20))
})


# Targeted tests for uncovered code paths in JAGS-marglik.R

test_that("JAGS_bridgesampling_posterior input validation works", {

  posterior <- matrix(rnorm(30), nrow = 10, ncol = 3)
  colnames(posterior) <- c("mu", "sigma", "x")

  # Input validation errors

  expect_error(JAGS_bridgesampling_posterior(data.frame(x = 1), prior_list = NULL), "'posterior' must be a matrix")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = "x"), "'prior_list' must be a list.")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = list(x = 1)), "'prior_list' must be a list of priors.")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = 1), "'add_parameters' must be a character")
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = NULL,
      add_bounds = list(lb = -Inf, ub = Inf)
    ),
    "requires at least one 'add_parameters'",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = NULL,
      add_parameters = character(),
      add_bounds = list(lb = numeric(), ub = numeric())
    ),
    "requires at least one 'add_parameters'",
    fixed = TRUE
  )
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = "x", add_bounds = "x"), "'add_bounds' must be a list")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = "x", add_bounds = list(a = 1)), "'add_bounds' must contain lower and upper bounds")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = c("x", "y"), add_bounds = list(lb = 0, ub = 1)), "'lb' and 'ub' must have the same length")
  expect_error(JAGS_bridgesampling_posterior(posterior, prior_list = NULL, add_parameters = "x", add_bounds = list(lb = "a", ub = "b")), "'lb' and 'ub' must be numeric")
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = NULL,
      add_parameters = "x",
      add_bounds = list(lb = -Inf, ub = Inf)
    ),
    "names must be unique and match",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = NULL,
      add_parameters = "x",
      add_bounds = list(lb = stats::setNames(-Inf, "wrong"), ub = stats::setNames(Inf, "x"))
    ),
    "names must match",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = NULL,
      add_parameters = "x",
      add_bounds = list(lb = c(x = 1), ub = c(x = 0))
    ),
    "smaller than upper",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = list(mu = prior("normal", list(0, 1))),
      add_parameters = "mu",
      add_bounds = list(lb = c(mu = -Inf), ub = c(mu = Inf))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      posterior,
      prior_list = list(sigma = prior("invgamma", list(2, 1))),
      add_parameters = "sigma",
      add_bounds = list(lb = c(sigma = 0), ub = c(sigma = Inf))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_add_parameters_not_formula(
      add_parameters = "mu",
      formula_design_list = list(mu = list()),
      formula_prior_list = list()
    ),
    "BayesTools-owned formula parameter",
    fixed = TRUE
  )

  dirichlet_posterior <- matrix(
    c(1, 3, 0.25),
    nrow = 1,
    dimnames = list(NULL, c("prior_par_eta_w[1]", "prior_par_eta_w[2]", "w[1]"))
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      dirichlet_posterior,
      prior_list = list(w = prior("dirichlet", list(alpha = c(1, 1)))),
      add_parameters = "w[1]",
      add_bounds = list(lb = c("w[1]" = 0), ub = c("w[1]" = 1))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )

  # Unsupported prior types
  expect_error(
    JAGS_bridgesampling_posterior(posterior, prior_list = list(p1 = prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("beta", list(1, 1))))),
    "spike and slab"
  )
  expect_error(
    JAGS_bridgesampling_posterior(posterior, prior_list = list(p1 = prior_mixture(list(prior("normal", list(0, 1)), prior("normal", list(1, 1))), is_null = c(TRUE, FALSE)))),
    "prior mixture"
  )

  # Missing parameters
  posterior_small <- matrix(rnorm(20), nrow = 10, ncol = 2)
  colnames(posterior_small) <- c("a", "b")
  expect_error(JAGS_bridgesampling_posterior(posterior_small, prior_list = list(x = prior("normal", list(0, 1)))), "'posterior' does not contain all")

  # Successful case with add_parameters
  result <- JAGS_bridgesampling_posterior(posterior, prior_list = list(mu = prior("normal", list(0, 1))), add_parameters = "x", add_bounds = list(lb = c(x = -Inf), ub = c(x = Inf)))
  expect_true(is.matrix(result))
  expect_true("x" %in% colnames(result))

})

test_that("JAGS_marglik_priors input validation and edge cases work", {

  # Empty prior_list contributes zero log prior density

  expect_equal(JAGS_marglik_priors(list(), prior_list = list()), 0)
  expect_equal(JAGS_marglik_priors_formula(list(), formula_prior_list = list(mu = list())), 0)

  # Input validation
  expect_error(JAGS_marglik_priors(list(), prior_list = "x"), "'prior_list' must be a list.")
  expect_error(JAGS_marglik_priors(list(), prior_list = prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")
  expect_error(JAGS_marglik_priors(list(), prior_list = list(x = 1)), "'prior_list' must be a list of priors.")

})

test_that("JAGS_marglik_parameters input validation and edge cases work", {

  # Test: empty prior_list returns empty list
  result <- JAGS_marglik_parameters(list(), prior_list = list())
  expect_equal(result, list())

  # Test: prior_list must be a list
  expect_error(
    JAGS_marglik_parameters(list(), prior_list = "not_a_list"),
    "'prior_list' must be a list."
  )

  # Test: prior_list must be a list of priors (single prior passed)
  expect_error(
    JAGS_marglik_parameters(list(), prior_list = prior("normal", list(0, 1))),
    "'prior_list' must be a list of priors."
  )

  # Test: prior_list must be a list of priors (non-prior elements)
  expect_error(
    JAGS_marglik_parameters(list(), prior_list = list(x = 1)),
    "'prior_list' must be a list of priors."
  )

})

test_that(".fit_to_posterior handles different input types", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("coda")

  prior_list <- list(mu = prior("normal", list(0, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  monitor <- JAGS_to_monitor(prior_list)
  inits <- JAGS_get_inits(prior_list, chains = 2, seed = 1)
  log_posterior <- STANDARD_LOG_POSTERIOR

  set.seed(1)
  model <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)

  # mcmc.list (rjags::coda.samples)
  samples_mcmc_list <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 100, quiet = TRUE, progress.bar = "none")
  marglik <- JAGS_bridgesampling(samples_mcmc_list, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik, "BayesTools_marglik")

  # mcmc (coda::as.mcmc)
  samples_mcmc <- coda::as.mcmc(samples_mcmc_list[[1]])
  marglik_mcmc <- JAGS_bridgesampling(samples_mcmc, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik_mcmc, "BayesTools_marglik")

  # Error for unsupported input
  expect_error(JAGS_bridgesampling("bad_input", prior_list = prior_list, data = list(), log_posterior = log_posterior), "not implemented")

})

test_that(".fit_to_posterior handles jags.samples output", {

  skip_if_not_installed("rjags")

  # Scalar parameter
  prior_list <- list(mu = prior("normal", list(0, 1)), sigma = prior("gamma", list(1, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  monitor <- JAGS_to_monitor(prior_list)
  inits <- JAGS_get_inits(prior_list, chains = 2, seed = 1)
  log_posterior <- STANDARD_LOG_POSTERIOR

  set.seed(1)
  model <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples_jags <- rjags::jags.samples(model = model, variable.names = monitor, n.iter = 100, progress.bar = "none")
  marglik_jags <- JAGS_bridgesampling(samples_jags, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik_jags, "BayesTools_marglik")

})

test_that(".fit_to_posterior handles vector parameters in jags.samples", {

  skip_if_not_installed("rjags")

  # Vector parameter (K > 1)
  prior_list <- list(p = prior("mnormal", list(mean = 0, sd = 1, K = 3)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  monitor <- JAGS_to_monitor(prior_list)
  inits <- JAGS_get_inits(prior_list, chains = 2, seed = 1)
  log_posterior <- STANDARD_LOG_POSTERIOR

  set.seed(1)
  model <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples_jags <- rjags::jags.samples(model = model, variable.names = monitor, n.iter = 100, progress.bar = "none")
  marglik_jags <- JAGS_bridgesampling(samples_jags, prior_list = prior_list, data = list(), log_posterior = log_posterior)
  expect_s3_class(marglik_jags, "BayesTools_marglik")

})

test_that("JAGS_bridgesampling handles a BayesTools fit", {

  skip_if_not_installed("rjags")

  prior_list <- list(mu = prior("normal", list(0, 1)))
  log_posterior <- STANDARD_LOG_POSTERIOR
  fit <- JAGS_fit("model{}", prior_list = prior_list, chains = 2, adapt = 100,
                  burnin = 100, sample = 500, silent = TRUE, seed = 1,
                  jags_modules = "glm")

  marglik <- JAGS_bridgesampling(
    fit, prior_list = prior_list, data = list(),
    log_posterior = log_posterior, seed = 11
  )
  expect_s3_class(marglik, "BayesTools_marglik")
  expect_equal(marglik$logml, 0, tolerance = 0.1)

})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-nonlocal-fit.R
# ============================================================================ #


expect_nonlocal_prior_only_samples <- function(prior, samples, tolerance = .08) {
  expect_true(all(is.finite(samples)))
  expect_true(all(samples >= prior$truncation[["lower"]]))
  expect_true(all(samples <= prior$truncation[["upper"]]))

  probs <- c(.1, .25, .5, .75, .9)
  quantiles <- quant(prior, probs)
  sample_cdf <- vapply(quantiles, function(x) mean(samples <= x), numeric(1))
  expect_equal(sample_cdf, probs, tolerance = tolerance)

  expect_warning(prior_mean <- mean(prior), NA)
  expect_warning(prior_sd <- sd(prior), NA)
  if (is.finite(prior_mean) && is.finite(prior_sd)) {
    expect_equal(mean(samples), prior_mean, tolerance = tolerance)
    expect_equal(stats::sd(samples), prior_sd, tolerance = tolerance)
  }
}

test_that("BayesTools JAGS module initializes truncated nonlocal priors", {
  skip_if_not_installed("rjags")
  skip_on_cran()

  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  priors <- list(
    moment = prior(
      "moment",
      list(mode = .5),
      truncation = list(lower = -.1, upper = .1)
    ),
    invmoment = prior(
      "invmoment",
      list(mode = .5, df = 6),
      truncation = list(lower = -.1, upper = .1)
    )
  )

  for(prior_i in priors){
    syntax <- JAGS_add_priors("model{}", list(theta = prior_i))
    expect_silent(local({
      con <- textConnection(syntax)
      on.exit(close(con), add = TRUE)
      rjags::jags.model(
        file     = con,
        data     = list(),
        n.chains = 1,
        n.adapt  = 0,
        quiet    = TRUE
      )
    }))
  }
})

test_that("BayesTools JAGS module samples nonlocal priors", {
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  skip_on_cran()

  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  priors <- list(
    moment = prior("moment", list(mode = .5)),
    invmoment = prior("invmoment", list(mode = .5, df = 6)),
    moment_truncated = prior("moment", list(mode = .5), truncation = list(lower = -Inf, upper = 0)),
    invmoment_truncated = prior("invmoment", list(mode = .5, df = 6), truncation = list(lower = -Inf, upper = 0))
  )

  for(prior_name in names(priors)){
    prior_list <- list(theta = priors[[prior_name]])
    fit <- suppressWarnings(JAGS_fit(
      model_syntax = "model{}",
      data = NULL,
      prior_list = prior_list,
      chains = 2,
      adapt = 250,
      burnin = 250,
      sample = 4000,
      silent = TRUE,
      seed = 1
    ))

    expect_s3_class(fit, "BayesTools_fit")
    expect_true("BayesTools" %in% attr(fit, "jags_modules"))
    samples <- as.matrix(fit$mcmc)
    expect_true("theta" %in% colnames(samples))
    expect_nonlocal_prior_only_samples(priors[[prior_name]], samples[, "theta"])

    marglik <- JAGS_bridgesampling(
      fit = fit,
      log_posterior = STANDARD_LOG_POSTERIOR,
      data = list(),
      prior_list = prior_list,
      maxiter = 2000,
      seed = 1
    )
    expect_s3_class(marglik, "BayesTools_marglik")
    expect_equal(marglik$logml, 0, tolerance = .08)
  }
})

test_that("Gamma-coordinate priors compile and sample with automatic and manual bounded starts", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  priors <- list(
    prior("invgamma", list(shape = 3, scale = 1e-310)),
    prior("moment", list(mode = .5), list(lower = 14, upper = Inf)),
    prior("moment", list(mode = .5), list(lower = -1, upper = 1)),
    prior("moment", list(mode = .5), list(lower = -1, upper = -.1)),
    prior("moment", list(tau = .125), list(lower = -1e-8, upper = 1e-8)),
    prior("invmoment", list(tau = 1, df = 3), list(lower = -.01, upper = .01))
  )
  for(p in priors){
    syntax <- JAGS_add_priors("model{}", list(theta = p))
    for(manual in c(FALSE, TRUE)){
      inits <- list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 620L)
      if(manual) inits$theta <- quant(p, .25)
      model <- local({
        con <- textConnection(syntax)
        on.exit(close(con), add = TRUE)
        rjags::jags.model(con, data = list(), inits = inits, n.chains = 1,
                          n.adapt = 0, quiet = TRUE)
      })
      if(manual){
        expect_equal(model$state()[[1L]]$theta, inits$theta)
      }
      draws <- as.matrix(rjags::coda.samples(model, "theta", n.iter = 100, progress.bar = "none"))[, "theta"]
      expect_true(all(is.finite(draws)))
      expect_true(all(draws >= p$truncation$lower & draws <= p$truncation$upper))
      expect_true(all(draws != if(p$distribution %in% c("moment", "invmoment")) p$parameters$location else 0))
    }
  }
})

test_that("JAGS numerical refusals leave the real backend session usable", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  compile <- function(p, value = NULL, observed = FALSE){
    syntax <- if(observed) "model{ y ~ dnorm(theta, 1) }" else "model{}"
    data <- if(observed) list(y = 0) else list()
    con <- textConnection(JAGS_add_priors(syntax, list(theta = p)))
    on.exit(close(con), add = TRUE)
    inits <- list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 621L)
    if(!is.null(value)) inits$theta <- value
    rjags::jags.model(con, data = data, inits = inits, n.chains = 1,
                      n.adapt = 0, quiet = TRUE)
  }
  unresolved <- prior("invgamma", list(shape = .Machine$double.xmin * .Machine$double.eps, scale = 1),
                       list(lower = .4, upper = .5))
  # An observed child forces JAGS to initialize its parent at compilation.
  expect_error(compile(unresolved, observed = TRUE), "numerically unavailable", fixed = TRUE)
  manual_unresolved <- compile(unresolved, .45)
  expect_error(rjags::coda.samples(manual_unresolved, "theta", n.iter = 2, progress.bar = "none"),
               "numerically unavailable", fixed = TRUE)
  central <- prior("moment", list(mode = .5), list(lower = -.1, upper = .1))
  # A prior-only model can accept an external start without evaluating its density.
  expect_identical(lpdf(central, 0), -Inf)
  manual_central <- compile(central, 0)
  expect_identical(manual_central$state()[[1L]]$theta, 0)
  automatic_central <- compile(central, observed = TRUE)
  automatic_value <- automatic_central$state()[[1L]]$theta
  expect_true(automatic_value >= central$truncation$lower && automatic_value <= central$truncation$upper)
  expect_true(is.finite(lpdf(central, automatic_value)))
  model <- compile(prior("invmoment", list(tau = 1, df = 1e-310)))
  expect_error(rjags::coda.samples(model, "theta", n.iter = 2, progress.bar = "none"),
               "numerically unavailable", fixed = TRUE)
  expect_silent(compile(prior("moment", list(mode = .5))))
})

test_that("fully structural fits retain deterministic draw geometry", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  fit <- NULL
  expect_silent(
    fit <- JAGS_fit(
      model_syntax = "model{}",
      prior_list = list(theta = prior("point", list(0))),
      chains = 2,
      adapt = 100,
      burnin = 100,
      sample = 100,
      silent = TRUE,
      seed = 1
    )
  )

  coordinates <- parameter_coordinates(fit)
  expect_false(.bt_backend_anchor_name %in% coordinates$coordinate_name)
  expect_identical(
    coordinates$monitor_status[coordinates$coordinate_name == "theta"],
    "structural"
  )
  expect_identical(
    coordinates$fixed_value[coordinates$coordinate_name == "theta"],
    0
  )

  catalog <- parameter_catalog(fit)
  expect_false(.bt_backend_anchor_name %in% catalog$quantities$canonical_name)
  theta <- catalog$quantities[catalog$quantities$canonical_name == "theta", ]
  expect_identical(theta$status, "structural")
  expect_identical(theta$fixed_value, 0)
  expect_identical(JAGS_fit_contract(fit)$parameter_map_version,
                   BayesTools:::.bt_parameter_map_version)

  geometry <- JAGS_draw_geometry(fit)
  expect_identical(geometry$chains$iterations, c(100L, 100L))
  expect_identical(geometry$total_draws, 200L)

  draws <- JAGS_materialize_draws(fit)
  expect_identical(colnames(draws[[1L]]), "theta")
  expect_identical(as.numeric(draws[[1L]][, "theta"]), rep(0, 100))
  expect_false(.bt_backend_anchor_name %in% colnames(draws[[1L]]))

  empty_draws <- JAGS_materialize_draws(fit, character())
  expect_identical(dim(empty_draws[[1L]]), c(100L, 0L))
  expect_identical(attr(empty_draws[[1L]], "mcpar"), c(201, 300, 1))

  extended <- JAGS_extend(
    fit,
    autofit_control = list(
      max_Rhat = NULL,
      min_ESS = NULL,
      max_error = NULL,
      max_SD_error = NULL,
      max_time = list(time = 30, unit = "secs"),
      sample_extend = 100,
      restarts = 1,
      max_extend = 1
    ),
    silent = TRUE
  )
  extended_geometry <- JAGS_draw_geometry(extended)
  expect_identical(extended_geometry$chains$iterations, c(200L, 200L))
  expect_identical(extended_geometry$total_draws, 400L)
  expect_identical(extended_geometry$chains$end, c(400L, 400L))
  expect_identical(parameter_map(extended), parameter_map(fit))
  expect_identical(parameter_catalog(extended), catalog)
})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-fit-settings.R
# ============================================================================ #
#
# PURPOSE:
#   Seeded fits and extensions: distinct chain seeds, reproducible draws, and
#   the caller's random-number state around every seeded public function that
#   fits or reads a fit.
#
# TAGS: @fit, @JAGS, @seed
# ============================================================================ #

test_that("prior-only JAGS draws differ between chains of adjacent seeds", {

  skip_if_not_installed("rjags")

  fit_prior_only <- function(seed){
    JAGS_fit(
      model_syntax = "model{}",
      prior_list   = list(mu = prior("normal", list(0, 1))),
      chains       = 2,
      adapt        = 50,
      burnin       = 50,
      sample       = 100,
      seed         = seed
    )
  }
  chain_draws <- function(fit, chain) as.numeric(fit$mcmc[[chain]][, "mu"])

  fit_1       <- fit_prior_only(1)
  fit_1_again <- fit_prior_only(1)
  fit_2       <- fit_prior_only(2)

  # With '.RNG.seed = seed + chain' these two chains were identical, because a
  # prior-only node is sampled from the prior regardless of its initial value.
  expect_false(identical(chain_draws(fit_1, 2), chain_draws(fit_2, 1)))
  expect_false(identical(chain_draws(fit_1, 1), chain_draws(fit_1, 2)))
  for(chain in 1:2){
    expect_identical(chain_draws(fit_1, chain), chain_draws(fit_1_again, chain))
  }
})

test_that("JAGS fits with priors only in the model syntax are reproducible for a seed", {

  skip_if_not_installed("rjags")

  fit_syntax_prior <- function(seed){
    JAGS_fit(
      model_syntax   = "model{ mu ~ dnorm(0, 1) }",
      prior_list     = NULL,
      add_parameters = "mu",
      chains         = 2,
      adapt          = 50,
      burnin         = 50,
      sample         = 100,
      seed           = seed
    )
  }
  chain_draws <- function(fit, chain) as.numeric(fit$mcmc[[chain]][, "mu"])

  fit_1       <- fit_syntax_prior(1)
  fit_1_again <- fit_syntax_prior(1)
  fit_2       <- fit_syntax_prior(2)

  # Without '.RNG.seed' entries, the backend seeded these chains itself.
  for(chain in 1:2){
    expect_identical(chain_draws(fit_1, chain), chain_draws(fit_1_again, chain))
    expect_false(identical(chain_draws(fit_1, chain), chain_draws(fit_2, chain)))
  }
})

# The JAGS-free part of these checks, JAGS_get_inits(), is a unit test
# (test-JAGS-fit-settings.R); '.expect_scoped_rng()' is in common-functions.R.
test_that("seeded public functions leave the caller's random-number state unchanged", {

  skip_if_not_installed("rjags")
  withr::local_preserve_seed()
  withr::defer(RNGkind("default", "default", "default"))

  set.seed(3)
  data <- list(x = stats::rnorm(20, 0.3), N = 20L)
  syntax <- "model{ for(i in 1:N){ x[i] ~ dnorm(mu, pow(s, -2)) } }"
  priors <- list(
    mu = prior("normal", list(0, 1)),
    s = prior("normal", list(0, 1), list(0, Inf))
  )
  priors_null <- list(
    mu = prior("spike", list(0)),
    s = prior("normal", list(0, 1), list(0, Inf))
  )
  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$x, parameters$mu, parameters$s, log = TRUE))
  }
  fit <- function(prior_list, seed, ...){
    JAGS_fit(syntax, data, prior_list, chains = 2, adapt = 1000, burnin = 50,
             sample = 100, seed = seed, ...)
  }
  fit_alternative <- fit(priors, seed = 7)
  fit_null <- fit(priors_null, seed = 8)
  models <- list(
    list(
      fit = fit_alternative,
      marglik = JAGS_bridgesampling(fit_alternative, log_posterior, data, priors, seed = 1),
      prior_weights = 1
    ),
    list(
      fit = fit_null,
      marglik = JAGS_bridgesampling(fit_null, log_posterior, data, priors_null, seed = 1),
      prior_weights = 1
    )
  )
  is_null_list <- list(mu = c(FALSE, TRUE), s = c(FALSE, FALSE))

  .expect_scoped_rng("JAGS_fit", function(){
    as.matrix(fit(priors, seed = 7)$mcmc)
  })
  .expect_scoped_rng("JAGS_fit (autofit)", function(){
    as.matrix(fit(
      priors, seed = 7, autofit = TRUE,
      autofit_control = list(min_ESS = 1e6, max_extend = 1, sample_extend = 50)
    )$mcmc)
  })
  .expect_scoped_rng("JAGS_fit (parallel)", function(){
    as.matrix(fit(priors, seed = 7, parallel = TRUE, cores = 2)$mcmc)
  })
  .expect_scoped_rng("JAGS_bridgesampling", function(){
    JAGS_bridgesampling(fit_alternative, log_posterior, data, priors, seed = 2)$logml
  })
  .expect_scoped_rng("mix_posteriors", function(){
    lapply(
      mix_posteriors(models, c("mu", "s"), is_null_list, seed = 3, n_samples = 200),
      as.numeric
    )
  })
  .expect_scoped_rng("marginal_inference", function(){
    inference <- marginal_inference(
      models, marginal_parameters = "mu", parameters = c("mu", "s"),
      is_null_list = is_null_list, formula = NULL, n_samples = 200,
      seed = 4, silent = TRUE
    )
    list(as.numeric(inference$averaged$mu), as.numeric(inference$conditional$mu))
  })
  .expect_scoped_rng("transform_prior_samples", function(){
    transform_prior_samples(fit_alternative, n_samples = 100, seed = 5)
  })

  # Unseeded calls take their seed from the caller's stream: one draw.
  set.seed(9)
  mix_posteriors(models, c("mu", "s"), is_null_list, seed = NULL, n_samples = 200)
  after_mix <- .Random.seed
  set.seed(9)
  sample(.Machine$integer.max, 1)
  expect_identical(after_mix, .Random.seed)
})

test_that("generated JAGS names retain requested monitor cells and their order", {

  model <- paste(
    "model {",
    "p1 ~ dnorm(0,1); p2 ~ dnorm(0,1); p10 ~ dnorm(0,1)",
    "mu__xREx__id_scalar[1] <- p1",
    "for (i in 1:4) { mu__xREx__id_vector[i] <- p2 + i }",
    "}", sep = "\n"
  )
  fit <- .bt_rjags_fit(model, list(),
    inits = list(list(.RNG.name = "base::Super-Duper", .RNG.seed = 91L)),
    monitor = c("p2", "mu__xREx__id_vector[2:3]", "p10", "mu__xREx__id_vector[3]",
                "p1", "mu__xREx__id_scalar[1]"),
    adapt = 50L, burnin = 50L, sample = 20L, thin = 1L)
  draws <- as.matrix(fit$mcmc)
  expect_identical(colnames(draws), c("p2", "mu__xREx__id_vector[2]",
    "mu__xREx__id_vector[3]", "p10", "p1", "mu__xREx__id_scalar[1]"))
  expect_equal(draws[, "mu__xREx__id_vector[2]"], draws[, "p2"] + 2)
  expect_equal(draws[, "mu__xREx__id_vector[3]"], draws[, "p2"] + 3)
  expect_identical(draws[, "mu__xREx__id_scalar[1]"], draws[, "p1"])

  sparse <- .bt_rjags_fit("model { x[1] <- 1; x[3] <- 3 }", list(),
    inits = list(list(.RNG.name = "base::Super-Duper", .RNG.seed = 92L)),
    monitor = "x[1:3]", adapt = 50L, burnin = 50L, sample = 20L, thin = 1L)
  expect_identical(colnames(as.matrix(sparse$mcmc)), c("x[1]", "x[3]"))
})

test_that("data-block-only monitors and sampled cells survive model reconstruction", {

  fit <- JAGS_fit(
    model_syntax = paste(
      "data { y[1] ~ dnorm(0, 1); y[3] <- -1; generated[1] <- 1; generated[2] ~ dnorm(0, 1) }",
      "model { mu ~ dnorm(0, 1); for (i in 1:3) { y[i] ~ dnorm(mu, 1) } }", sep = "\n"
    ),
    add_parameters = c("mu", "generated[2]", "y"), chains = 2,
    adapt = 50, burnin = 50, sample = 100, seed = 492, silent = TRUE
  )
  coordinates <- parameter_coordinates(fit)
  roles <- setNames(coordinates$convergence_role, coordinates$coordinate_name)
  expect_identical(roles[["y[2]"]], "sampled")
  expect_identical(roles[["generated[2]"]], "structural")
  expect_true(is.na(fit$data$y[2L]))
  for(state in fit$end_state){
    expect_identical(is.na(as.numeric(state$y)), c(TRUE, FALSE, TRUE))
  }

  extended <- JAGS_extend(unserialize(serialize(fit, NULL)),
    autofit_control = list(max_Rhat = NULL, min_ESS = NULL, max_error = NULL,
      max_SD_error = NULL, max_extend = 1, sample_extend = 25))
  expect_identical(extended$data, fit$data)
  for(chain in extended$mcmc){
    expect_equal(nrow(chain), 125L)
    expect_identical(as.numeric(chain[, "generated[2]"]), rep(fit$data$generated[2L], 125L))
    expect_identical(as.numeric(chain[, "y[1]"]), rep(fit$data$y[1L], 125L))
    expect_true(length(unique(chain[101:125, "y[2]"])) > 1L)
  }
})

test_that("PSOCK workers and a reloaded extension share generated data and preserve thinned chains", {

  fit_once <- function(){
    JAGS_fit(
      model_syntax = paste("data { y_raw ~ dnorm(0, 1) }",
                           "model { y_raw ~ dnorm(mu, 1); y <- y_raw }", sep = "\n"),
      prior_list = list(mu = prior("normal", list(0, 1))),
      add_parameters = "y", chains = 3, cores = 2, parallel = TRUE,
      adapt = 50, burnin = 50, sample = 100, thin = 3, seed = 49, silent = TRUE
    )
  }
  fit <- fit_once()
  expect_identical(class(fit), "BayesTools_fit")
  expect_identical(coda::as.mcmc.list(fit), fit$mcmc)
  expect_true(all(vapply(fit$end_state, function(chain) length(chain$.RNG.state) > 0L, logical(1))))
  generated <- fit$data$y_raw
  expect_true(length(generated) == 1L && is.finite(generated))
  expect_identical(fit_once()$data$y_raw, generated)

  saved <- unserialize(serialize(fit, NULL))
  control <- list(max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = NULL,
                  max_extend = 1, sample_extend = 25)
  extended <- JAGS_extend(saved, autofit_control = control, parallel = TRUE, cores = 2)
  expect_identical(extended$data$y_raw, generated)
  expect_length(extended$mcmc, 3L)
  for(i in seq_along(fit$mcmc)){
    before <- fit$mcmc[[i]]
    after <- extended$mcmc[[i]]
    expect_equal(nrow(before), 100)
    expect_equal(nrow(after), 125)
    expect_equal(coda::thin(after), 3)
    expect_equal(stats::start(after), stats::start(before))
    expect_equal(stats::end(after), stats::end(before) + 3 * 25)
    expect_identical(as.matrix(after)[seq_len(100), , drop = FALSE], as.matrix(before))
    expect_identical(as.numeric(after[, "y"]), rep(as.numeric(generated), 125L))
  }
})

test_that("JAGS_extend draws depend only on the stored fit", {

  skip_if_not_installed("rjags")
  withr::local_preserve_seed()

  set.seed(3)
  data <- list(x = stats::rnorm(20, 0.3), N = 20L)
  fit <- JAGS_fit(
    "model{ for(i in 1:N){ x[i] ~ dnorm(mu, pow(s, -2)) } }", data,
    list(mu = prior("normal", list(0, 1)), s = prior("normal", list(0, 1), list(0, Inf))),
    chains = 2, adapt = 1000, burnin = 50, sample = 100, seed = 1
  )
  control <- list(
    max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = NULL,
    sample_extend = 50, max_extend = 1
  )
  draws <- function(extended) as.matrix(extended$mcmc)

  # The session's compiled model is not continued: extending the same object
  # twice, a saved and reloaded copy, or under another caller state gives the
  # same draws.
  first    <- JAGS_extend(fit, autofit_control = control)
  second   <- JAGS_extend(fit, autofit_control = control)
  reloaded <- JAGS_extend(unserialize(serialize(fit, NULL)), autofit_control = control)
  set.seed(99)
  other_caller <- JAGS_extend(fit, autofit_control = control)
  expect_identical(draws(second), draws(first))
  expect_identical(draws(reloaded), draws(first))
  expect_identical(draws(other_caller), draws(first))
  expect_false(identical(draws(first), as.matrix(fit$mcmc)))

  # The public fit drops the compiled model and rebuilds from its saved state.
  expect_null(fit$live_model)
  expect_null(first$live_model)
  with_model <- fit
  model_text <- textConnection(fit$model)
  with_model$live_model <- rjags::jags.model(model_text, fit$data, fit$end_state,
                                            n.chains = length(fit$mcmc),
                                            n.threads = 1, n.adapt = 0, quiet = TRUE)
  close(model_text)
  model_size <- length(serialize(with_model$live_model, NULL))
  expect_gt(
    length(serialize(with_model, NULL)) - length(serialize(fit, NULL)),
    0.9 * model_size
  )
  expect_identical(draws(JAGS_extend(with_model, autofit_control = control)), draws(first))
})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-convergence.R
# ============================================================================ #

test_that("the backend anchor of a model without monitors is auxiliary", {

  skip_if_not_installed("rjags")
  fit <- suppressWarnings(JAGS_fit(
    "model{ x ~ dnorm(0, 1) }", prior_list = NULL,
    chains = 2, adapt = 50, burnin = 50, sample = 100,
    autofit = TRUE,
    autofit_control = list(max_extend = 2, sample_extend = 50),
    seed = 1
  ))
  # Autofit does not extend: the anchor is not assessable but not selected.
  expect_null(attr(fit, "warnings"))
  expect_equal(nrow(summary(fit)), 0L)
  coordinates <- parameter_coordinates(fit)
  expect_identical(
    stats::setNames(coordinates$convergence_role, coordinates$coordinate_name),
    c(BayesTools_backend_anchor = "auxiliary")
  )
  result <- JAGS_check_convergence(fit)
  expect_true(result)
  diagnostics <- attr(result, "diagnostics")
  expect_identical(
    stats::setNames(diagnostics$state, diagnostics$parameter),
    c(BayesTools_backend_anchor = "not_requested")
  )
})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS FROM test-JAGS-selection-inits.R
# ============================================================================ #

test_that("heterogeneous weightfunction mixtures compile and adapt in JAGS", {

  # The end-to-end guard: a wrong init node name is only fatal once JAGS sees
  # the model, and only when the expansion changes the array length.
  mixture <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)), prior_weights = 1),
    prior_weightfunction("one-sided", c(.05, .10), wf_independent(prior("gamma", list(shape = 9, rate = 3))), prior_weights = 1),
    prior_weightfunction("one-sided", c(.025), wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega"), prior_weights = 1),
    prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .4)), prior_weights = 1)
  ))

  fit <- suppressWarnings(JAGS_fit(
    "model{}",
    data       = NULL,
    prior_list = list(bias = mixture),
    chains     = 1,
    adapt      = 50,
    burnin     = 50,
    sample     = 100,
    seed       = 14
  ))

  expect_s3_class(fit, "BayesTools_fit")
  expect_false(inherits(fit, "condition"))
})

# ============================================================================ #
# CENTRALIZED LIVE-FIT TESTS: STORED ENVIRONMENTS
# ============================================================================ #
#
# PURPOSE:
#   A fit keeps no environment of the code that created it, so a saved fit
#   does not carry the calling workspace.
#
# TAGS: @fit, @JAGS, @formula, @random
# ============================================================================ #

test_that("a fit created in a function does not keep that function's frame", {

  skip_if_not_installed("rjags")
  withr::local_preserve_seed()

  # The same random-effect model, fitted in a function that holds a large
  # local object and in one that does not.
  fit_random <- function(local_object_size){
    local_object <- stats::runif(local_object_size)
    data <- data.frame(
      x = c(-1.2, 0.4, 0.9, -0.3, 1.5, -0.8, 0.2, 1.1, -0.6, 0.7, 1.3, -1.0),
      g = factor(rep(c("g1", "g2", "g3"), each = 4))
    )
    set.seed(1)
    y <- stats::rnorm(nrow(data), 0.2 * data$x, 1)
    suppressWarnings(JAGS_fit(
      model_syntax = "model{\n  for(i in 1:N_mu){\n    y[i] ~ dnorm(mu[i], 1)\n  }\n}",
      data = list(y = y),
      formula_list = list(mu = ~ 1 + x + (1 + x | g)),
      formula_data_list = list(mu = data),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      )),
      formula_random_prior_list = list(mu = prior_random(
        g = random_block(sd = prior("normal", list(0, 1), list(0, Inf)), cor = prior_lkj(eta = 1))
      )),
      chains = 1, adapt = 50, burnin = 50, sample = 100, seed = 1, silent = TRUE
    ))
  }
  without_id <- function(parameter_map){
    attr(parameter_map, "runtime_cache_id") <- NULL
    parameter_map
  }
  # The size of 'small' is taken before 'large' exists: a kept frame reaches
  # this test's environment, and through it both fits.
  small <- fit_random(0)
  small_size <- length(serialize(small, NULL))
  large <- fit_random(1e6)

  expect_identical(stored_environment_paths(large, "fit"), character())
  # A kept frame adds the 8 MB local object; the parameter map's runtime cache
  # id, unique per fit, may differ by a few bytes.
  expect_lt(abs(length(serialize(large, NULL)) - small_size), 1024)
  expect_identical(as.matrix(large$mcmc), as.matrix(small$mcmc))
  expect_identical(JAGS_formula_design(large), JAGS_formula_design(small))
  expect_identical(
    without_id(attr(large, "parameter_map")),
    without_id(attr(small, "parameter_map"))
  )
})

test_that("the JAGS module registers its six distributions and two functions", {

  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))
  # JAGS compiles (and initializes, without sampling) syntax that uses each of
  # the names only if the module registered all of them
  syntax <- JAGS_add_priors(
    paste0("model{
", JAGS_lkj_corr_cholesky("Sigma", K = 3, eta = 1)$syntax, "
}"),
    list(
      a = prior("moment", list(location = 0, tau = 1, order = 2)),
      b = prior("invmoment", list(location = 0, tau = 1, order = 2, df = 3)),
      c = prior("invgamma", list(2, 1))
    )
  )
  # the sampler markers are emitted by JAGS_fit() only; they are added by hand
  syntax <- sub("\\}[[:space:]]*$", paste0(
    "BayesTools_gate_1 ~ dbt_gate(a, 1)\n",
    "BayesTools_block_1 ~ dbt_block(b, 1)\n}"
  ), syntax)
  for(name in c("dbt_lkj_cpc", "dbt_invgamma", "dbt_moment", "dbt_invmoment", "dbt_gate", "dbt_block",
                "bt_lkj_cholesky", "bt_lkj_corr")){
    expect_true(grepl(name, syntax, fixed = TRUE), info = name)
  }
  model <- rjags::jags.model(
    textConnection(syntax), data = list(BayesTools_gate_1 = 0, BayesTools_block_1 = 0),
    n.chains = 1, n.adapt = 0, quiet = TRUE
  )
  expect_s3_class(model, "jags")
})

# ============================================================================ #
# Gated and factor slice samplers of the BayesTools JAGS module
# ============================================================================ #
#
# Declarations are observed marker nodes, 'BayesTools_gate_<k> ~
# dbt_gate(<node>, <active>)' and 'BayesTools_block_<k> ~ dbt_block(<node>, <id>)'
# with data value 0. These blocks compile hand-written models with rjags.

.sampler_test_node_samplers <- function(model){
  samplers <- rjags::list.samplers(model)
  out <- unlist(lapply(seq_along(samplers), function(i){
    stats::setNames(rep(names(samplers)[i], length(samplers[[i]])), samplers[[i]])
  }))
  out[order(names(out))]
}

.sampler_test_compile <- function(syntax, data, chains = 2L, seed = 1L, n_adapt = 0L, inits = list()){
  rjags::jags.model(
    textConnection(syntax), data = data, n.chains = chains, n.adapt = n_adapt, quiet = TRUE,
    inits = lapply(seq_len(chains), function(chain){
      c(list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = seed + chain), inits)
    })
  )
}

.sampler_test_markers_zero <- function(names){
  stats::setNames(as.list(rep(0, length(names))), names)
}

test_that("sampler factories claim exactly the nodes their markers declare", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  set.seed(1)
  N <- 20L
  x1 <- stats::rnorm(N)
  x2 <- x1 + stats::rnorm(N, 0, 0.3)
  data <- list(N = N, x1 = x1, x2 = x2, y = 0.5 * x1 - 0.2 * x2 + stats::rnorm(N))
  syntax <- function(markers) paste0("
model{
  ind ~ dcat(c(1, 1))
  x ~ dnorm(0, 1)T(0,)
  b1 ~ dnorm(0, 1)
  b2 ~ dnorm(0, 1)
  sigma ~ dunif(0.1, 5)
  m <- x * equals(ind, 2)
  for(i in 1:N){
    y[i] ~ dnorm(m + b1 * x1[i] + b2 * x2[i], 1 / sigma^2)
  }
  ", paste(markers, collapse = "\n  "), "
}")
  gate <- "BayesTools::GatedSlice"
  block <- "BayesTools::FactorSlice"

  # declared: nodes whose sampler may change; claimed: the expected new samplers
  cases <- list(
    list(label = "gate", markers = "BayesTools_gate_1 ~ dbt_gate(x, equals(ind, 2))",
         observed = "BayesTools_gate_1", declared = "x", claimed = c(x = gate)),
    list(label = "block of two", markers = c("BayesTools_block_1 ~ dbt_block(b1, 1)", "BayesTools_block_2 ~ dbt_block(b2, 1)"),
         observed = paste0("BayesTools_block_", 1:2), declared = c("b1", "b2"), claimed = c(b1 = block, b2 = block)),
    list(label = "block with a doubly bounded member",
         markers = c("BayesTools_block_1 ~ dbt_block(b1, 1)", "BayesTools_block_2 ~ dbt_block(b2, 1)",
                     "BayesTools_block_3 ~ dbt_block(sigma, 1)"),
         observed = paste0("BayesTools_block_", 1:3), declared = c("b1", "b2", "sigma"),
         claimed = c(b1 = block, b2 = block, sigma = block)),
    list(label = "two blocks", markers = c("BayesTools_block_1 ~ dbt_block(b1, 1)", "BayesTools_block_2 ~ dbt_block(b2, 1)",
                                           "BayesTools_block_3 ~ dbt_block(x, 2)", "BayesTools_block_4 ~ dbt_block(sigma, 2)"),
         observed = paste0("BayesTools_block_", 1:4), declared = c("b1", "b2", "x", "sigma"),
         claimed = c(b1 = block, b2 = block, x = block, sigma = block), blocks = 2L),
    list(label = "gate wins over a block",
         markers = c("BayesTools_gate_1 ~ dbt_gate(x, equals(ind, 2))", "BayesTools_block_1 ~ dbt_block(x, 1)",
                     "BayesTools_block_2 ~ dbt_block(b1, 1)", "BayesTools_block_3 ~ dbt_block(b2, 1)"),
         observed = c("BayesTools_gate_1", paste0("BayesTools_block_", 1:3)), declared = c("x", "b1", "b2"),
         claimed = c(x = gate, b1 = block, b2 = block)),
    list(label = "block of one", markers = "BayesTools_block_1 ~ dbt_block(b1, 1)",
         observed = "BayesTools_block_1", declared = "b1", claimed = character()),
    list(label = "gate leaves a block of one", markers = c("BayesTools_gate_1 ~ dbt_gate(x, equals(ind, 2))",
                                                             "BayesTools_block_1 ~ dbt_block(x, 1)", "BayesTools_block_2 ~ dbt_block(b1, 1)"),
         observed = c("BayesTools_gate_1", paste0("BayesTools_block_", 1:2)), declared = c("x", "b1"), claimed = c(x = gate)),
    list(label = "markers on deterministic expressions",
         markers = c("BayesTools_block_1 ~ dbt_block(b1 + 0, 1)", "BayesTools_block_2 ~ dbt_block(b2 + 0, 1)",
                     "BayesTools_gate_1 ~ dbt_gate(x * 1, equals(ind, 2))"),
         observed = c(paste0("BayesTools_block_", 1:2), "BayesTools_gate_1"), declared = c("b1", "b2", "x"), claimed = character()),
    list(label = "unobserved markers",
         markers = c("BayesTools_block_1 ~ dbt_block(b1, 1)", "BayesTools_block_2 ~ dbt_block(b2, 1)",
                     "BayesTools_gate_1 ~ dbt_gate(x, equals(ind, 2))"),
         observed = character(), declared = c("b1", "b2", "x"), claimed = character()),
    list(label = "non-integer block id", markers = c("BayesTools_block_1 ~ dbt_block(b1, 1.5)", "BayesTools_block_2 ~ dbt_block(b2, 1.5)"),
         observed = paste0("BayesTools_block_", 1:2), declared = c("b1", "b2"), claimed = character()),
    list(label = "zero block id", markers = c("BayesTools_block_1 ~ dbt_block(b1, 0)", "BayesTools_block_2 ~ dbt_block(b2, 0)"),
         observed = paste0("BayesTools_block_", 1:2), declared = c("b1", "b2"), claimed = character()),
    list(label = "block id that is not fixed", markers = c("BayesTools_block_1 ~ dbt_block(b1, ind)", "BayesTools_block_2 ~ dbt_block(b2, ind)"),
         observed = paste0("BayesTools_block_", 1:2), declared = c("b1", "b2"), claimed = character()),
    list(label = "gate that depends on its node", markers = "BayesTools_gate_1 ~ dbt_gate(x, step(x - 1))",
         observed = "BayesTools_gate_1", declared = "x", claimed = character()),
    list(label = "two gates on one node",
         markers = c("BayesTools_gate_1 ~ dbt_gate(x, equals(ind, 2))", "BayesTools_gate_2 ~ dbt_gate(x, equals(ind, 1))"),
         observed = paste0("BayesTools_gate_", 1:2), declared = "x", claimed = character()),
    list(label = "marker naming the node as its second parent", markers = "BayesTools_gate_1 ~ dbt_gate(b1, x)",
         observed = "BayesTools_gate_1", declared = c("b1", "x"), claimed = c(b1 = gate))
  )

  # The block samplers of the glm module ('glm::Generic') update the
  # coefficients of a linear predictor jointly and refuse a node with a marker
  # child, so with glm loaded (earlier blocks of this case load it) an
  # undeclared coefficient that glm sampled jointly with a declared one gets
  # another stock sampler. The comparison is exact without glm; with glm,
  # undeclared nodes keep stock samplers. The session's module state is
  # restored afterwards.
  glm_loaded <- "glm" %in% rjags::list.modules()
  withr::defer({
    if(glm_loaded){
      rjags::load.module("glm", quiet = TRUE)
    }else if("glm" %in% rjags::list.modules()){
      rjags::unload.module("glm", quiet = TRUE)
    }
  })
  for(with_glm in c(FALSE, TRUE)){
    if(with_glm){
      rjags::load.module("glm", quiet = TRUE)
    }else if(glm_loaded){
      rjags::unload.module("glm", quiet = TRUE)
    }
    reference <- .sampler_test_node_samplers(.sampler_test_compile(syntax(character()), data))
    for(case in cases){
      label <- paste0(case$label, if(with_glm) " (glm loaded)")
      model <- .sampler_test_compile(syntax(case$markers), c(data, .sampler_test_markers_zero(case$observed)))
      current <- .sampler_test_node_samplers(model)
      expect_identical(names(current), names(reference), info = label)
      undeclared <- setdiff(names(reference), case$declared)
      if(with_glm){
        expect_false(any(current[undeclared] %in% c(gate, block)), info = label)
      }else{
        expect_identical(current[undeclared], reference[undeclared], info = label)
      }
      expect_identical(unname(current[names(case$claimed)]), unname(case$claimed), info = label)
      unclaimed <- setdiff(case$declared, names(case$claimed))
      expect_false(any(current[unclaimed] %in% c(gate, block)), info = label)
      if(!is.null(case$blocks)){
        samplers <- rjags::list.samplers(model)
        expect_identical(sum(names(samplers) == block), case$blocks, info = label)
      }
      expect_no_error(rjags::adapt(model, 200, end.adaptation = TRUE, progress.bar = "none"))
      expect_no_error(stats::update(model, 200, progress.bar = "none"))
    }
  }

  # a marker is an observed node with value 0
  expect_error(
    .sampler_test_compile(syntax("BayesTools_gate_1 ~ dbt_gate(x, equals(ind, 2))"), c(data, list(BayesTools_gate_1 = 1))),
    "inconsistent"
  )
})

# Monte Carlo checks of a fit: each quantity is the posterior mean of a series
# of the draws, compared with its exact value in units of its Monte Carlo
# standard error (sd / sqrt(ESS), with the ESS of coda summed over chains).
.sampler_test_mcse_checks <- function(samples, quantities){

  rows <- lapply(names(quantities), function(name){
    quantity <- quantities[[name]]
    chains <- coda::mcmc.list(lapply(samples, function(chain){
      coda::mcmc(quantity$series(as.matrix(chain)))
    }))
    values <- unlist(lapply(chains, as.numeric), use.names = FALSE)
    ess <- sum(coda::effectiveSize(chains))
    estimate <- mean(values)
    mcse <- stats::sd(values) / sqrt(ess)
    # a chain that never moves has no MCSE (sd 0, ESS 0): it fails the checks
    z <- if(isTRUE(mcse > 0)) (estimate - quantity$target) / mcse else if(estimate == quantity$target) 0 else Inf
    data.frame(quantity = name, target = quantity$target, estimate = estimate,
               mcse = mcse, ess = ess, z = z, stringsAsFactors = FALSE)
  })
  checks <- do.call(rbind, rows)
  if(isTRUE(getOption("BayesTools.sampler_test_print"))) print(checks, digits = 4)
  checks
}

.sampler_test_expect_checks <- function(checks, label, min_ess = 100){
  info <- paste(c(label, utils::capture.output(print(checks, digits = 4))), collapse = "\n")
  expect_true(all(abs(checks$z) <= 4), info = info)
  expect_true(all(checks$ess >= min_ess), info = info)
}

# Moments, a covariance and two tail probabilities of a bivariate normal
# posterior of the columns 'names' (exact mean 'm' and covariance 'V').
.sampler_test_gaussian_quantities <- function(names, m, V){

  sd <- sqrt(diag(V))
  list(
    mean_1 = list(series = function(x) x[, names[1]], target = m[1]),
    mean_2 = list(series = function(x) x[, names[2]], target = m[2]),
    var_1  = list(series = function(x) (x[, names[1]] - m[1])^2, target = V[1, 1]),
    var_2  = list(series = function(x) (x[, names[2]] - m[2])^2, target = V[2, 2]),
    cov    = list(series = function(x) (x[, names[1]] - m[1]) * (x[, names[2]] - m[2]), target = V[1, 2]),
    tail_1 = list(series = function(x) as.numeric(x[, names[1]] > m[1] + sd[1]),
                  target = stats::pnorm(1, lower.tail = FALSE)),
    tail_2 = list(series = function(x) as.numeric(x[, names[2]] < m[2] - 1.5 * sd[2]),
                  target = stats::pnorm(-1.5))
  )
}

.sampler_test_scalar_quantities <- function(name, mean, variance, threshold, tail){
  stats::setNames(list(
    list(series = function(x) x[, name], target = mean),
    list(series = function(x) (x[, name] - mean)^2, target = variance),
    list(series = function(x) as.numeric(x[, name] > threshold), target = tail)
  ), paste0(name, c("_mean", "_var", "_tail")))
}

# Draws of a hand-written model whose factor slice block adapts for
# 'n_adapt' iterations (0: coordinate axes and unit widths throughout).
.sampler_test_block_samples <- function(syntax, data, variables, n_adapt, n_iter, chains = 4L, seed = 1L, inits = list()){

  model <- .sampler_test_compile(syntax, data, chains = chains, seed = seed, inits = inits)
  samplers <- rjags::list.samplers(model)
  expect_identical(sum(names(samplers) == "BayesTools::FactorSlice"), 1L)
  expect_setequal(samplers[[which(names(samplers) == "BayesTools::FactorSlice")]], variables)
  # the block reports incomplete adaptation until its first adaptive update,
  # so zero adaptation leaves rjags' adaptation check FALSE (the stock state)
  adapted <- rjags::adapt(model, n.iter = n_adapt, end.adaptation = TRUE, progress.bar = "none")
  expect_identical(adapted, n_adapt > 0L)
  stats::update(model, 1000, progress.bar = "none")
  rjags::coda.samples(model, variables, n.iter = n_iter, progress.bar = "none")
}

.sampler_test_block_markers <- function(nodes, id = 1L){
  list(
    syntax = paste0("BayesTools_block_", seq_along(nodes), " ~ dbt_block(", nodes, ", ", id, ")", collapse = "\n"),
    data = .sampler_test_markers_zero(paste0("BayesTools_block_", seq_along(nodes)))
  )
}

test_that("a closed gate draws its node directly from the node's distribution", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # every node reaches its likelihood child only through the closed gate 'off'
  syntax <- "
model{
  x1 ~ dnorm(0.5, 1)T(0,)
  x2 ~ dt(1, 0.25, 3)
  x3 ~ dbeta(2, 3)
  x4 ~ dgamma(2, 1.5)
  x5 ~ dbt_invgamma(3, 2)
  y[1] ~ dnorm(x1 * off, 1)
  y[2] ~ dnorm(x2 * off, 1)
  y[3] ~ dnorm(x3 * off, 1)
  y[4] ~ dnorm(x4 * off, 1)
  y[5] ~ dnorm(x5 * off, 1)
  BayesTools_gate_1 ~ dbt_gate(x1, off)
  BayesTools_gate_2 ~ dbt_gate(x2, off)
  BayesTools_gate_3 ~ dbt_gate(x3, off)
  BayesTools_gate_4 ~ dbt_gate(x4, off)
  BayesTools_gate_5 ~ dbt_gate(x5, off)
}"
  data <- c(list(y = c(0.3, -0.2, 0.1, 0.4, -0.5), off = 0), .sampler_test_markers_zero(paste0("BayesTools_gate_", 1:5)))
  model <- .sampler_test_compile(syntax, data, chains = 2L, seed = 3L)
  samplers <- .sampler_test_node_samplers(model)
  expect_identical(unname(samplers[paste0("x", 1:5)]), rep("BayesTools::GatedSlice", 5L))
  rjags::adapt(model, n.iter = 200, end.adaptation = TRUE, progress.bar = "none")
  draws <- as.matrix(rjags::coda.samples(model, paste0("x", 1:5), n.iter = 5000, progress.bar = "none"))

  cdf <- list(
    x1 = function(q) (stats::pnorm(q, 0.5, 1) - stats::pnorm(0, 0.5, 1)) / stats::pnorm(0, 0.5, 1, lower.tail = FALSE),
    x2 = function(q) stats::pt((q - 1) * sqrt(0.25), 3),
    x3 = function(q) stats::pbeta(q, 2, 3),
    x4 = function(q) stats::pgamma(q, 2, 1.5),
    x5 = function(q) stats::pgamma(2 / q, 3, lower.tail = FALSE)
  )
  # the draws of a closed gate are independent: the Kolmogorov-Smirnov test of
  # each distribution at level 0.001
  tests <- do.call(rbind, lapply(names(cdf), function(name){
    test <- suppressWarnings(stats::ks.test(draws[, name], cdf[[name]]))
    data.frame(node = name, D = unname(test$statistic), p = test$p.value)
  }))
  if(isTRUE(getOption("BayesTools.sampler_test_print"))) print(tests, digits = 4)
  expect_true(all(tests$p > 0.001), info = paste(utils::capture.output(print(tests)), collapse = "\n"))
})

test_that("a closed gate keeps its node where its deterministic descendants are finite", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # exp(lw) overflows above log(.Machine$double.xmax) = 709.78, which the prior
  # N(700, 9.44^2) exceeds with probability 0.15; infinity times the closed gate
  # would be NaN in the likelihood
  sd <- 9.44
  bound <- log(.Machine$double.xmax)
  expect_equal(stats::pnorm(bound, 700, sd, lower.tail = FALSE), 0.15, tolerance = 0.01)
  syntax <- sprintf("
model{
  lw ~ dnorm(700, %s)
  w <- exp(lw)
  y ~ dnorm(w * off, 1)
  BayesTools_gate_1 ~ dbt_gate(lw, off)
}", format(1 / sd^2, digits = 17))
  data <- list(y = 0.2, off = 0, BayesTools_gate_1 = 0)
  model <- rjags::jags.model(
    textConnection(syntax), data = data, n.chains = 2L, n.adapt = 0, quiet = TRUE,
    inits = lapply(1:2, function(chain) list(lw = 700, .RNG.name = "base::Mersenne-Twister", .RNG.seed = 20 + chain))
  )
  expect_identical(unname(.sampler_test_node_samplers(model)["lw"]), "BayesTools::GatedSlice")
  expect_no_error(rjags::adapt(model, n.iter = 500, end.adaptation = TRUE, progress.bar = "none"))
  # a rejected draw repeats the state (acceptance 0.85): thinned by 5
  draws <- as.matrix(rjags::coda.samples(model, c("lw", "w"), n.iter = 25000, thin = 5, progress.bar = "none"))
  expect_true(all(is.finite(draws[, "w"])))
  expect_true(all(draws[, "lw"] <= bound))
  truncated_cdf <- function(q) stats::pnorm(q, 700, sd) / stats::pnorm(bound, 700, sd)
  test <- suppressWarnings(stats::ks.test(draws[, "lw"], truncated_cdf))
  if(isTRUE(getOption("BayesTools.sampler_test_print"))) print(c(D = unname(test$statistic), p = test$p.value), digits = 4)
  expect_gt(test$p.value, 0.001)
})

test_that("a gate on a node that reaches the likelihood outside its switch fails during adaptation", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  syntax <- "
model{
  x ~ dnorm(0, 1)
  y ~ dnorm(x + x * off, 1)
  BayesTools_gate_1 ~ dbt_gate(x, off)
}"
  data <- list(y = 0.5, off = 0, BayesTools_gate_1 = 0)
  model <- .sampler_test_compile(syntax, data, chains = 1L)
  expect_identical(unname(.sampler_test_node_samplers(model)["x"]), "BayesTools::GatedSlice")
  expect_error(
    rjags::adapt(model, n.iter = 10, end.adaptation = TRUE, progress.bar = "none"),
    "gated node changes the likelihood while its gate is closed", fixed = TRUE
  )
})

test_that("gated fits of a bias mixture reproduce the quadrature posterior", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # y_i ~ N(mu + PET se_i, se_i^2), mu ~ N(0, 1), PET ~ N(0, 1)T(0,), and the
  # branches 'none' and 'PET' with equal prior weights
  K <- 12L
  se <- seq(0.05, 0.5, length.out = K)
  set.seed(11)
  y <- 0.15 + 0.6 * se + stats::rnorm(K, 0, se)
  data <- list(y = y, se = se, K = K)
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = prior_mixture(list(prior_none(), prior_PET("normal", list(0, 1))), is_null = c(TRUE, FALSE))
  )

  # oracle: mu integrated analytically (y | PET ~ N(PET se, S + 11')), PET by
  # one-dimensional quadrature, checked against two-dimensional quadrature
  S_inv <- diag(1 / se^2)
  Sigma <- diag(se^2) + 1
  log_normal <- function(x, mean, covariance){
    chol_cov <- chol(covariance)
    z <- backsolve(chol_cov, x - mean, transpose = TRUE)
    -0.5 * length(x) * log(2 * pi) - sum(log(diag(chol_cov))) - 0.5 * sum(z^2)
  }
  log_m1 <- log_normal(y, rep(0, K), Sigma)
  log_m2_PET <- function(PET) vapply(PET, function(p) log(2) + stats::dnorm(p, log = TRUE) + log_normal(y, p * se, Sigma), numeric(1))
  shift <- log_m1
  m2 <- stats::integrate(function(p) exp(log_m2_PET(p) - shift), 0, Inf, rel.tol = 1e-10)$value
  precision_mu <- 1 + sum(1 / se^2)
  mean_mu <- function(p) vapply(p, function(pp) sum((y - pp * se) / se^2) / precision_mu, numeric(1))
  e2 <- stats::integrate(function(p) mean_mu(p) * exp(log_m2_PET(p) - shift), 0, Inf, rel.tol = 1e-10)$value / m2
  e1 <- sum(y / se^2) / precision_mu
  probability <- m2 / (exp(log_m1 - shift) + m2)
  model_averaged_mu <- (1 - probability) * e1 + probability * e2

  log_likelihood <- function(mu, p) sum(stats::dnorm(y, mu + p * se, se, log = TRUE))
  m2_2d <- stats::integrate(function(p) vapply(p, function(pp){
    stats::integrate(function(mu) vapply(mu, function(m){
      exp(log_likelihood(m, pp) + stats::dnorm(m, log = TRUE) + log(2) + stats::dnorm(pp, log = TRUE) - shift)
    }, numeric(1)), -Inf, Inf, rel.tol = 1e-10)$value
  }, numeric(1)), 0, Inf, rel.tol = 1e-8)$value
  m1_1d <- stats::integrate(function(mu) vapply(mu, function(m){
    exp(log_likelihood(m, 0) + stats::dnorm(m, log = TRUE) - shift)
  }, numeric(1)), -Inf, Inf, rel.tol = 1e-10)$value
  expect_equal(m2_2d, m2, tolerance = 1e-6)
  expect_equal(m1_1d, 1, tolerance = 1e-6)
  expect_gt(probability, 0.2)
  expect_lt(probability, 0.8)

  for(gates in c(FALSE, TRUE)){
    fit <- JAGS_fit(
      model_syntax = "model{\n  for(i in 1:K){\n    y[i] ~ dnorm(mu + PET * se[i], 1 / pow(se[i], 2))\n  }\n}",
      data = data, prior_list = prior_list,
      chains = 4, adapt = 500, burnin = 1000, sample = 10000, seed = 5, silent = TRUE,
      sampler_gates = gates
    )
    declarations <- attr(fit, "sampler_declarations")
    if(gates){
      expect_identical(declarations$gates$node, "PET_1")
      expect_true(grepl("BayesTools_gate_1 ~ dbt_gate(PET_1, equals(bias_indicator, 2))", fit$model, fixed = TRUE))
    }else{
      expect_null(declarations)
      expect_false(grepl("dbt_gate", fit$model, fixed = TRUE))
    }
    checks <- .sampler_test_mcse_checks(fit$mcmc, list(
      probability_PET = list(series = function(x) as.numeric(x[, "bias_indicator"] == 2), target = probability),
      model_averaged_mu = list(series = function(x) x[, "mu"], target = model_averaged_mu)
    ))
    .sampler_test_expect_checks(checks, paste0("gates = ", gates))
  }
})

test_that("factor slice blocks reproduce Gaussian and transformed targets for every adaptation length", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # (1) a linear model with two correlated predictors, known variance and
  # normal priors: conjugate posterior N(Q^-1 X'y, Q^-1), Q = X'X + I/100
  set.seed(21)
  N <- 30L
  x1 <- stats::rnorm(N)
  x2 <- 0.95 * x1 + sqrt(1 - 0.95^2) * stats::rnorm(N)
  y <- 0.5 * x1 - 0.3 * x2 + stats::rnorm(N)
  X <- cbind(x1, x2)
  Q <- crossprod(X) + diag(0.01, 2)
  linear <- list(
    syntax = "
model{
  b1 ~ dnorm(0, 0.01)
  b2 ~ dnorm(0, 0.01)
  for(i in 1:N){
    y[i] ~ dnorm(b1 * x1[i] + b2 * x2[i], 1)
  }
  %s
}",
    data = list(N = N, x1 = x1, x2 = x2, y = y),
    nodes = c("b1", "b2"),
    quantities = .sampler_test_gaussian_quantities(c("b1", "b2"), drop(solve(Q, crossprod(X, y))), solve(Q))
  )

  # (2) Gaussian targets N(m, V) from an observation m of N(b, V) with nearly
  # flat priors: posterior N((P + I/10^4)^-1 P m, (P + I/10^4)^-1), P = V^-1.
  # With equal variances the likelihood is that of the independent sum and
  # difference, b1 + b2 ~ N(m1 + m2, 2 sd^2 (1 + r)) and b1 - b2 ~ N(m1 - m2, 2 sd^2 (1 - r)).
  # 'narrow' (sd 0.05, correlation 0.995) is much narrower than the unit width
  # of an unadapted block; 'wide' (sd 9, correlation 0.9, conditional sd 3.9)
  # is wider than the stepping-out limit of 10 unit widths reaches.
  gaussian <- function(sd, correlation, m){
    P <- solve(sd^2 * matrix(c(1, correlation, correlation, 1), 2))
    posterior_precision <- P + diag(1e-4, 2)
    list(
      syntax = "
model{
  b1 ~ dnorm(0, 1e-4)
  b2 ~ dnorm(0, 1e-4)
  z_sum ~ dnorm(b1 + b2, precision_sum)
  z_difference ~ dnorm(b1 - b2, precision_difference)
  %s
}",
      data = list(z_sum = m[1] + m[2], z_difference = m[1] - m[2],
                  precision_sum = 1 / (2 * sd^2 * (1 + correlation)),
                  precision_difference = 1 / (2 * sd^2 * (1 - correlation))),
      nodes = c("b1", "b2"),
      quantities = .sampler_test_gaussian_quantities(
        c("b1", "b2"), drop(solve(posterior_precision, P %*% m)), solve(posterior_precision))
    )
  }

  # (3) every transform: unbounded (mu, theta; theta is a child of mu in the
  # block), lower-bounded (s), doubly bounded (p, and t truncated to [-1, 3]),
  # upper-bounded (u truncated to (-Inf, 2])
  y_theta <- c(0.8, 1.4, 1.1)
  Q_theta <- matrix(c(5, -4, -4, 4 + length(y_theta)), 2)
  m_theta <- drop(solve(Q_theta, c(0, sum(y_theta))))
  V_theta <- solve(Q_theta)
  ratio_u <- stats::dnorm(1) / stats::pnorm(1)
  Z_t <- stats::pnorm(3) - stats::pnorm(-1)
  mean_t <- (stats::dnorm(-1) - stats::dnorm(3)) / Z_t
  transforms <- list(
    syntax = "
model{
  mu ~ dnorm(0, 1)
  theta ~ dnorm(mu, 4)
  for(i in 1:3){
    y[i] ~ dnorm(theta, 1)
  }
  s ~ dgamma(3, 2)
  p ~ dbeta(2, 5)
  u ~ dnorm(1, 1)T(,2)
  t ~ dnorm(0, 1)T(-1, 3)
  %s
}",
    data = list(y = y_theta),
    nodes = c("mu", "theta", "s", "p", "u", "t"),
    quantities = c(
      .sampler_test_gaussian_quantities(c("mu", "theta"), m_theta, V_theta),
      .sampler_test_scalar_quantities("s", 1.5, 0.75, 3, stats::pgamma(3, 3, 2, lower.tail = FALSE)),
      .sampler_test_scalar_quantities("p", 2 / 7, 10 / (49 * 8), 0.5, stats::pbeta(0.5, 2, 5, lower.tail = FALSE)),
      .sampler_test_scalar_quantities("u", 1 - ratio_u, 1 - ratio_u - ratio_u^2, 0,
                                      (stats::pnorm(1) - stats::pnorm(-1)) / stats::pnorm(1)),
      .sampler_test_scalar_quantities("t", mean_t,
                                      1 + (-1 * stats::dnorm(-1) - 3 * stats::dnorm(3)) / Z_t - mean_t^2, 1,
                                      (stats::pnorm(3) - stats::pnorm(1)) / Z_t)
    )
  )

  # (4) a member whose scale is another member: s ~ Gamma(3, 2), v | s ~ N(0, 1/s)
  # (E[v^2] = E[1/s] = 2 / (3 - 1)); wrong when the density of v omits its
  # normalising constant in s
  scale_member <- list(
    syntax = "
model{
  s ~ dgamma(3, 2)
  v ~ dnorm(0, s)
  %s
}",
    data = list(),
    nodes = c("s", "v"),
    quantities = c(
      .sampler_test_scalar_quantities("s", 1.5, 0.75, 3, stats::pgamma(3, 3, 2, lower.tail = FALSE)),
      list(v_square = list(series = function(x) x[, "v"]^2, target = 1))
    )
  )

  # (5) a member whose truncation normaliser depends on another member:
  # m ~ N(0, 1), v | m ~ N(m, 1) T(0,); the normaliser pnorm(m) of v cancels
  # in the joint density only when the summed density keeps it, so m stays N(0, 1)
  truncation_member <- list(
    syntax = "
model{
  m ~ dnorm(0, 1)
  v ~ dnorm(m, 1)T(0,)
  %s
}",
    data = list(),
    nodes = c("m", "v"),
    quantities = list(
      m_mean = list(series = function(x) x[, "m"], target = 0),
      m_square = list(series = function(x) x[, "m"]^2, target = 1),
      m_tail = list(series = function(x) as.numeric(x[, "m"] > 1), target = stats::pnorm(-1))
    )
  )

  targets <- list(
    linear = linear,
    narrow = gaussian(0.05, 0.995, c(0.3, -0.2)),
    wide = gaussian(9, 0.9, c(4, -3)),
    transforms = transforms,
    scale_member = scale_member,
    truncation_member = truncation_member
  )
  for(target_name in names(targets)){
    target <- targets[[target_name]]
    markers <- .sampler_test_block_markers(target$nodes)
    for(n_adapt in c(0L, 150L, 1000L)){
      samples <- .sampler_test_block_samples(
        sprintf(target$syntax, markers$syntax), c(target$data, markers$data),
        variables = target$nodes, n_adapt = n_adapt, n_iter = 15000L, seed = 31L
      )
      checks <- .sampler_test_mcse_checks(samples, target$quantities)
      .sampler_test_expect_checks(checks, paste0(target_name, ", adapt = ", n_adapt))
    }
  }
})

test_that("factor slice blocks sample members whose transformed coordinate is beyond the range of exp()", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # A member with two bounds is sampled on u = log(x - lower) - log(upper - x);
  # its back-transform must hold wherever x is representable strictly inside the
  # support, also where |u| exceeds log(.Machine$double.xmax) = 709.78.
  # (1) Wide intervals: a normal with sd 0.1 truncated to (0, 1e308), a
  # half-normal, and its mirror on (-1e308, 0); |u| is 711.5 at |x| = 0.1.
  expect_lt(log(0.1) - log(1e308), -log(.Machine$double.xmax))
  half_normal_mean <- 0.1 * sqrt(2 / pi)
  half_normal_variance <- 0.01 * (1 - 2 / pi)
  normal_quantities <- function(sign) list(
    x_mean = list(series = function(x) x[, "x"], target = sign * half_normal_mean),
    x_var  = list(series = function(x) (x[, "x"] - sign * half_normal_mean)^2, target = half_normal_variance),
    x_tail = list(series = function(x) as.numeric(abs(x[, "x"]) > 0.1), target = 2 * stats::pnorm(-1))
  )
  # (2) Unit intervals at subnormal distances from a bound: x ~ U(0, 1) (or
  # U(-1, 0)) and an observation z = -720 of N(log|x|, 2^2). The prior density
  # of w = log|x| is exp(w) on w < 0, so w has the posterior N(-716, 2^2)
  # restricted to w < 0 (358 SDs away); |u| = |w| exceeds 709.78 for 99.9% of
  # the mass, and |x| stays far above the smallest subnormal, exp(-744.4)
  # (14 SDs away), so no draw rounds to a bound.
  log_quantities <- list(
    w_mean = list(series = function(x) log(abs(x[, "x"])), target = -716),
    w_var  = list(series = function(x) (log(abs(x[, "x"])) + 716)^2, target = 4),
    w_tail = list(series = function(x) as.numeric(log(abs(x[, "x"])) > -714), target = stats::pnorm(-1))
  )
  cases <- list(
    # At this scale JAGS's truncated-normal RNG overflows while standardizing
    # the far bound. Supply interior starts to exercise the block sampler.
    lower_wide = list(member = "x ~ dnorm(0, 100)T(0, 1e308)", data = list(), inits = list(x = 0.1),
                      quantities = normal_quantities(1)),
    upper_wide = list(member = "x ~ dnorm(0, 100)T(-1e308, 0)", data = list(), inits = list(x = -0.1),
                      quantities = normal_quantities(-1)),
    lower_unit = list(member = "x ~ dunif(0, 1)\n  z ~ dnorm(log(x), 0.25)", data = list(z = -720),
                      inits = list(x = exp(-714)), quantities = log_quantities),
    upper_unit = list(member = "x ~ dunif(-1, 0)\n  z ~ dnorm(log(-x), 0.25)", data = list(z = -720),
                      inits = list(x = -exp(-714)), quantities = log_quantities)
  )
  # each member shares a block with an ordinary normal
  markers <- .sampler_test_block_markers(c("x", "b"))
  for(case_name in names(cases)){
    case <- cases[[case_name]]
    syntax <- paste0("model{\n  ", case$member, "\n  b ~ dnorm(0, 1)\n", markers$syntax, "\n}")
    for(n_adapt in c(0L, 500L)){
      samples <- .sampler_test_block_samples(
        syntax, c(case$data, markers$data), variables = c("x", "b"),
        n_adapt = n_adapt, n_iter = 5000L, seed = 41L, inits = case$inits
      )
      checks <- .sampler_test_mcse_checks(samples, case$quantities)
      .sampler_test_expect_checks(checks, paste0(case_name, ", adapt = ", n_adapt))
    }
  }
})

test_that("factor slice blocks sample members far from distant bounds", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # A standard normal truncated to (-1e9, 1e9) or (-1e13, 1e13), or bounded
  # below by -1e6 or -1e12, in a block with an ordinary normal. The sampler
  # holds each member through its distance from the anchoring bound, to about
  # 16 significant digits of that distance (here a grid of up to about 0.02
  # posterior SDs), so the central mass, the moments and the partner's moments
  # are those of the target; their draws equal development's.
  # Seed 1: with seed 47 the partner's variance gave z = -4.02 under adaptation
  # 0 with the lower bound -1e6. Those draws are bit-identical on development
  # (90fee58), so the transform is not involved. Independent repetitions (200
  # chains with disjoint seeds; neighbouring seeds of this helper share chains,
  # .RNG.seed = seed + chain) give E[b^2] = 1.0023 +- 0.0021 against 1, which is
  # consistent with a Monte Carlo tail event at seed 47; it does not exclude a
  # small bias.
  quantities <- list(
    x_mean    = list(series = function(x) x[, "x"], target = 0),
    x_var     = list(series = function(x) x[, "x"]^2, target = 1),
    x_central = list(series = function(x) as.numeric(abs(x[, "x"]) <= 1), target = 2 * stats::pnorm(1) - 1),
    b_mean    = list(series = function(x) x[, "b"], target = 0),
    b_var     = list(series = function(x) x[, "b"]^2, target = 1)
  )
  markers <- .sampler_test_block_markers(c("x", "b"))
  members <- c("x ~ dnorm(0, 1)T(-1e9, 1e9)", "x ~ dnorm(0, 1)T(-1e13, 1e13)",
               "x ~ dnorm(0, 1)T(-1e6, )", "x ~ dnorm(0, 1)T(-1e12, )")
  for(member in members){
    syntax <- paste0("model{\n  ", member, "\n  b ~ dnorm(0, 1)\n", markers$syntax, "\n}")
    for(n_adapt in c(0L, 500L)){
      samples <- .sampler_test_block_samples(syntax, markers$data, variables = c("x", "b"),
                                             n_adapt = n_adapt, n_iter = 5000L, seed = 1L)
      checks <- .sampler_test_mcse_checks(samples, quantities)
      .sampler_test_expect_checks(checks, paste0(member, ", adapt = ", n_adapt))
    }
  }
})

test_that("a two-sided member whose support span is not finite stops the fit at the start", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # upper - lower overflows for (-1e308, 1e308): the member cannot be held
  # through its distance from a bound, which the bounds alone decide; a finite
  # span such as (0, 1e308) samples (the far-tail test above)
  markers <- .sampler_test_block_markers(c("x", "b"))
  syntax <- paste0("model{\n  x ~ dnorm(0, 1)T(-1e308, 1e308)\n  b ~ dnorm(0, 1)\n", markers$syntax, "\n}")
  expect_error(
    .sampler_test_compile(syntax, markers$data, chains = 2L, seed = 1L),
    paste0("Factor slice sampler: the support bounds of this node are too wide for the sampler's scale ",
           "(their difference exceeds the largest double); declare the block without this node"),
    fixed = TRUE
  )
})

test_that("generated declarations change the samplers of the declared nodes only", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # every branch type of a bias mixture, and a block of a simple prior, a
  # treatment factor, and the allocation primitives of a variance allocation
  allocation_sd <- prior("gamma", list(2, 2))
  attr(allocation_sd, "random_allocation_sd") <- TRUE
  attr(allocation_sd, "random_allocation") <- "tot"
  weights <- prior("dirichlet", list(alpha = c(1.5, 2.5)))
  attr(weights, "random_allocation") <- "tot"
  phacking <- prior_phacking(form = "linear")
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    f = prior_factor_levels(prior_factor("normal", list(0, 1), contrast = "treatment"), 3),
    sd_alloc = allocation_sd,
    w_alloc = weights,
    bias = prior_mixture(list(
      prior_none(),
      prior_weightfunction("one-sided", c(0.025, 0.05), wf_cumulative(c(1, 1, 1))),
      prior_weightfunction("one-sided", c(0.05), wf_independent(prior("normal", list(0, 1)), "log_omega")),
      phacking,
      prior_bias(prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1))), phacking),
      prior_PEESE("normal", list(0, 1)),
      prior_PET("normal", list(0, 1))
    ), is_null = c(TRUE, rep(FALSE, 6)))
  )
  syntax <- "model{
  for(i in 1:N){
    y[i] ~ dnorm(mu + f[1] * d1[i] + f[2] * d2[i] + PET * se[i] + PEESE * se[i]^2 + omega[2] + omega[3] + alpha,
                 1 / (se[i]^2 + sd_alloc^2 * w_alloc[1]))
  }
}"
  set.seed(4)
  N <- 20L
  data <- list(N = N, y = stats::rnorm(N), se = stats::runif(N, 0.1, 0.4),
               d1 = rep(c(1, 0), 10), d2 = rep(c(0, 1), 10))
  text <- JAGS_add_priors(syntax, prior_list)
  declared <- BayesTools:::.bt_sampler_declarations(
    text, data, prior_list, NULL, sampler_gates = TRUE,
    sampler_blocks = list(c("mu", "f", "sd_alloc", "w_alloc"))
  )
  expect_identical(nrow(declared$declarations$declined), 0L)
  gated <- declared$declarations$gates$node
  blocked <- declared$declarations$blocks$node
  expect_setequal(gated, c("eta_component_2[1]", "eta_component_2[2]", "eta_component_2[3]",
                           "log_omega_component_3[2]", "alpha_component_4",
                           "omega_ratio_component_5", "alpha_component_5", "PEESE_1", "PET_1"))
  expect_setequal(blocked, c("mu", "f[1]", "f[2]", "sd_alloc", "prior_par_eta_w_alloc[1]", "prior_par_eta_w_alloc[2]"))

  # JAGS may report a node under the name of a deterministic alias (the binary
  # ratio as omega_local_component_5[2]), so the nodes are matched by counts:
  # every declared node gets its sampler and no other node changes
  reference <- .sampler_test_node_samplers(.sampler_test_compile(text, data, seed = 5L))
  current_model <- .sampler_test_compile(declared$model, declared$data, seed = 5L)
  current <- .sampler_test_node_samplers(current_model)
  expect_identical(names(current), names(reference))
  changed <- names(current)[current != reference]
  expect_true(all(current[changed] %in% c("BayesTools::GatedSlice", "BayesTools::FactorSlice")))
  expect_identical(sum(current == "BayesTools::GatedSlice"), length(gated))
  expect_identical(sum(current == "BayesTools::FactorSlice"), length(blocked))
  expect_true(all(current[intersect(blocked, names(current))] == "BayesTools::FactorSlice"))
  expect_identical(sum(names(rjags::list.samplers(current_model)) == "BayesTools::FactorSlice"), 1L)
  expect_no_error(rjags::adapt(current_model, 200, end.adaptation = TRUE, progress.bar = "none"))
  expect_no_error(stats::update(current_model, 200, progress.bar = "none"))
})

test_that("markers that no factory claims leave draws, deviance, and DIC unchanged", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))
  rjags::load.module("dic", quiet = TRUE)

  # a block of one leaves 'sigma' with the stock slice sampler
  set.seed(8)
  N <- 15L
  x1 <- stats::rnorm(N)
  data <- list(N = N, x1 = x1, y = 0.4 * x1 + stats::rnorm(N, 0, 0.7))
  syntax <- function(marker) paste0("
model{
  b1 ~ dnorm(0, 0.01)
  sigma ~ dunif(0.1, 5)
  for(i in 1:N){
    y[i] ~ dnorm(b1 * x1[i], 1 / sigma^2)
  }
  ", marker, "
}")
  plain <- .sampler_test_compile(syntax(""), data, seed = 9L)
  marked <- .sampler_test_compile(syntax("BayesTools_block_1 ~ dbt_block(sigma, 1)"),
                                  c(data, list(BayesTools_block_1 = 0)), seed = 9L)
  expect_identical(.sampler_test_node_samplers(marked), .sampler_test_node_samplers(plain))
  expect_identical(unname(.sampler_test_node_samplers(marked)["sigma"]), "base::RealSlicer")
  expect_identical(plain$state(internal = TRUE), marked$state(internal = TRUE))
  for(model in list(plain, marked)){
    rjags::adapt(model, 100, end.adaptation = TRUE, progress.bar = "none")
    stats::update(model, 100, progress.bar = "none")
  }
  expect_identical(
    as.matrix(rjags::coda.samples(marked, c("b1", "sigma", "deviance"), n.iter = 500, progress.bar = "none")),
    as.matrix(rjags::coda.samples(plain, c("b1", "sigma", "deviance"), n.iter = 500, progress.bar = "none"))
  )
  dic_marked <- rjags::dic.samples(marked, n.iter = 200, type = "pD", progress.bar = "none")
  dic_plain <- rjags::dic.samples(plain, n.iter = 200, type = "pD", progress.bar = "none")
  # the marker is an observed node of the deviance with contribution 0
  values <- function(x) stats::setNames(as.numeric(x), names(x))
  for(component in c("deviance", "penalty")){
    marked <- values(dic_marked[[component]])
    plain_values <- values(dic_plain[[component]])
    expect_identical(marked[names(plain_values)], plain_values)
    expect_identical(unname(marked["BayesTools_block_1"]), 0)
    expect_identical(sum(marked), sum(plain_values))
  }
})

test_that("a rarely selected gated branch fits without adaptation warnings", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  K <- 10L
  se <- seq(0.1, 0.4, length.out = K)
  set.seed(12)
  data <- list(y = stats::rnorm(K, 0.2, se), se = se, K = K)
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = prior_mixture(list(prior_none(prior_weights = 200), prior_PET("normal", list(0, 1), prior_weights = 1)),
                         is_null = c(TRUE, FALSE))
  )
  warnings <- character()
  fit <- withCallingHandlers(
    JAGS_fit(
      model_syntax = "model{\n  for(i in 1:K){\n    y[i] ~ dnorm(mu + PET * se[i], 1 / pow(se[i], 2))\n  }\n}",
      data = data, prior_list = prior_list,
      chains = 2, adapt = 100, burnin = 200, sample = 1000, seed = 2, silent = FALSE,
      sampler_gates = TRUE
    ),
    warning = function(w){
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(attr(fit, "sampler_declarations")$gates$node, "PET_1")
  expect_lt(mean(as.matrix(fit$mcmc)[, "bias_indicator"] == 2), 0.05)
  expect_identical(warnings, character())
})

test_that("declared fits run serially, in parallel, with autofit, and extend after reloading", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  K <- 10L
  se <- seq(0.1, 0.4, length.out = K)
  set.seed(14)
  data <- list(y = stats::rnorm(K, 0.2, sqrt(se^2 + 0.1^2)), se = se, K = K)
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(0, 1), list(0, Inf)),
    bias = prior_mixture(list(prior_none(), prior_PET("normal", list(0, 1))), is_null = c(TRUE, FALSE))
  )
  syntax <- "model{\n  for(i in 1:K){\n    y[i] ~ dnorm(mu + PET * se[i], 1 / (pow(se[i], 2) + pow(tau, 2)))\n  }\n}"
  fit_with <- function(...){
    JAGS_fit(model_syntax = syntax, data = data, prior_list = prior_list,
             chains = 2, adapt = 100, burnin = 100, sample = 500, seed = 3, silent = TRUE,
             sampler_gates = TRUE, sampler_blocks = list(c("mu", "tau")), ...)
  }
  expect_declared <- function(fit, draws){
    expect_s3_class(fit, "BayesTools_fit")
    declarations <- attr(fit, "sampler_declarations")
    expect_identical(declarations$gates$node, "PET_1")
    expect_identical(declarations$blocks$node, c("mu", "tau"))
    expect_true("BayesTools" %in% attr(fit, "jags_modules"))
    expect_true("BayesTools" %in% attr(fit, "required_packages"))
    expect_true(grepl("BayesTools_block_2 ~ dbt_block(tau, 1)", fit$model, fixed = TRUE))
    expect_false(grepl("BayesTools_", attr(fit, "model_syntax"), fixed = TRUE))
    expect_identical(JAGS_draw_geometry(fit)$total_draws, as.integer(draws))
    samples <- as.matrix(fit$mcmc)
    expect_false(any(startsWith(colnames(samples), "BayesTools_")))
    expect_true(all(is.finite(samples[, c("mu", "tau")])))
  }

  serial <- fit_with()
  expect_declared(serial, 1000L)
  plain <- JAGS_fit(model_syntax = syntax, data = data, prior_list = prior_list,
                    chains = 2, adapt = 100, burnin = 100, sample = 500, seed = 3, silent = TRUE)
  expect_identical(colnames(as.matrix(serial$mcmc)), colnames(as.matrix(plain$mcmc)))
  expect_identical(parameter_catalog(serial), parameter_catalog(plain))

  # 'sampler_blocks' as a function of the combined prior list (here: the
  # simple continuous scalar priors) gives the fit of the equivalent list:
  # the same declarations, model text, sampler assignment and draws
  calls <- 0L
  select_scalar <- function(priors){
    calls <<- calls + 1L
    list(names(priors)[vapply(priors, function(prior){
      is.prior.simple(prior) && !is.prior.point(prior) && !is.prior.discrete(prior)
    }, logical(1))])
  }
  from_function <- JAGS_fit(model_syntax = syntax, data = data, prior_list = prior_list,
                            chains = 2, adapt = 100, burnin = 100, sample = 500, seed = 3, silent = TRUE,
                            sampler_gates = TRUE, sampler_blocks = select_scalar)
  expect_identical(calls, 1L)
  expect_declared(from_function, 1000L)
  expect_identical(attr(from_function, "sampler_declarations"), attr(serial, "sampler_declarations"))
  expect_identical(from_function$model, serial$model)
  expect_identical(as.matrix(from_function$mcmc), as.matrix(serial$mcmc))
  compiled_samplers <- function(fit){
    .sampler_test_node_samplers(rjags::jags.model(
      textConnection(fit$model), data = fit$data, n.chains = 1, n.adapt = 0, quiet = TRUE
    ))
  }
  samplers <- compiled_samplers(from_function)
  expect_identical(samplers, compiled_samplers(serial))
  expect_identical(unname(samplers[c("mu", "tau")]), rep("BayesTools::FactorSlice", 2L))
  expect_identical(unname(samplers["PET_1"]), "BayesTools::GatedSlice")

  parallel_fit <- fit_with(parallel = TRUE, cores = 2)
  expect_declared(parallel_fit, 1000L)

  # autofit extends at least once (no criterion can be met by 'min_ESS'); the
  # only warning is the 'max_extend' stop, which a silent fit stores on the fit
  # instead of raising it; no adaptation warning is raised or stored
  autofit_warnings <- character()
  autofit <- withCallingHandlers(
    fit_with(autofit = TRUE, autofit_control = list(
      max_Rhat = NULL, min_ESS = 1e6, max_error = NULL, max_SD_error = NULL,
      max_time = list(time = 60, unit = "secs"), sample_extend = 200, restarts = 1, max_extend = 1
    )),
    warning = function(w){
      autofit_warnings <<- c(autofit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(autofit_warnings, character())
  expect_identical(attr(autofit, "warnings"),
                   "The automatic model fitting was terminated due to the 'max_extend' constraint.")
  expect_declared(autofit, 1400L)

  control <- list(
    max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = NULL,
    max_time = list(time = 60, unit = "secs"), sample_extend = 200, restarts = 1, max_extend = 1
  )
  extended <- JAGS_extend(serial, autofit_control = control, silent = TRUE)
  expect_declared(extended, 1400L)
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(serial, file)
  reloaded <- JAGS_extend(readRDS(file), autofit_control = control, silent = TRUE)
  expect_declared(reloaded, 1400L)
  expect_identical(as.matrix(reloaded$mcmc), as.matrix(extended$mcmc))
  extended_parallel <- JAGS_extend(serial, autofit_control = control, silent = TRUE, parallel = TRUE, cores = 2)
  expect_declared(extended_parallel, 1400L)
})

test_that("bridge sampling of a fit with a factor slice block matches the exact marginal likelihood", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))
  skip_if_not_installed("bridgesampling")

  set.seed(21)
  N <- 30L
  x1 <- stats::rnorm(N)
  x2 <- 0.95 * x1 + sqrt(1 - 0.95^2) * stats::rnorm(N)
  y <- 0.5 * x1 - 0.3 * x2 + stats::rnorm(N)
  data <- list(N = N, x1 = x1, x2 = x2, y = y)
  prior_list <- list(b1 = prior("normal", list(0, 10)), b2 = prior("normal", list(0, 10)))
  fit <- JAGS_fit(
    model_syntax = "model{\n  for(i in 1:N){\n    y[i] ~ dnorm(b1 * x1[i] + b2 * x2[i], 1)\n  }\n}",
    data = data, prior_list = prior_list,
    chains = 4, adapt = 500, burnin = 500, sample = 5000, seed = 6, silent = TRUE,
    sampler_blocks = list(c("b1", "b2"))
  )
  expect_identical(attr(fit, "sampler_declarations")$blocks$node, c("b1", "b2"))
  marglik <- suppressWarnings(JAGS_bridgesampling(
    fit = fit,
    log_posterior = function(parameters, data){
      sum(stats::dnorm(data$y, parameters[["b1"]] * data$x1 + parameters[["b2"]] * data$x2, 1, log = TRUE))
    },
    data = data, prior_list = prior_list, seed = 1
  ))
  # y ~ N(0, I + 100 X X')
  X <- cbind(x1, x2)
  covariance <- diag(N) + 100 * tcrossprod(X)
  chol_cov <- chol(covariance)
  z <- backsolve(chol_cov, y, transpose = TRUE)
  exact <- -0.5 * N * log(2 * pi) - sum(log(diag(chol_cov))) - 0.5 * sum(z^2)
  error <- sqrt(bridgesampling::error_measures(marglik)$re2)
  if(isTRUE(getOption("BayesTools.sampler_test_print"))) print(c(logml = marglik$logml, exact = exact, error = error), digits = 8)
  expect_lt(abs(marglik$logml - exact), max(4 * error, 0.01))
})

test_that("a block of a random-effect SD and a scalar correlation keeps the posterior of the fit without it", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # y_ij = mu + u_ij + e_ij with a compound-symmetry random effect of the four
  # time points within eight groups (sd 0.6, correlation 0.5)
  set.seed(23)
  n_groups <- 8L
  n_times <- 4L
  formula_data <- data.frame(g = factor(rep(seq_len(n_groups), each = n_times)), t = rep(seq_len(n_times), n_groups))
  covariance <- 0.6^2 * (0.5 + 0.5 * diag(n_times))
  effects <- as.vector(t(matrix(stats::rnorm(n_groups * n_times), n_groups) %*% chol(covariance)))
  data <- list(y = 0.3 + effects + stats::rnorm(n_groups * n_times, 0, 0.3))
  random_names <- c("mu__xREx__g_sd", "mu__xREx__g_rho")
  random_fit <- function(...){
    JAGS_fit(
      model_syntax = "model{\n  for(i in 1:N_mu){\n    y[i] ~ dnorm(mu[i], 1 / pow(0.3, 2))\n  }\n}",
      data = data, prior_list = list(),
      formula_list = list(mu = ~ 1 + cs(t | g)), formula_data_list = list(mu = formula_data),
      formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)))),
      formula_random_prior_list = list(mu = prior_random(g = random_block(
        sd = prior("normal", list(0, 1), list(0, Inf)),
        covariance = random_covariance(cor = prior("uniform", list(-0.2, 0.95)), cor_scale = "cor")
      ))),
      chains = 4, adapt = 500, burnin = 1000, sample = 10000, seed = 7, silent = TRUE, ...
    )
  }
  plain <- random_fit()
  declared <- random_fit(sampler_blocks = list(c("mu_intercept", random_names)))

  # the block holds the random-effect SD and correlation under their own names
  declarations <- attr(declared, "sampler_declarations")
  expect_identical(declarations$blocks$node, c("mu_intercept", random_names))
  expect_identical(nrow(declarations$declined), 0L)
  expect_true(grepl("BayesTools_block_3 ~ dbt_block(mu__xREx__g_rho, 1)", declared$model, fixed = TRUE))
  model <- rjags::jags.model(textConnection(declared$model), data = declared$data,
                             n.chains = 1, n.adapt = 0, quiet = TRUE)
  samplers <- rjags::list.samplers(model)
  expect_identical(sum(names(samplers) == "BayesTools::FactorSlice"), 1L)
  expect_setequal(samplers[[which(names(samplers) == "BayesTools::FactorSlice")]], c("mu_intercept", random_names))
  expect_identical(colnames(as.matrix(declared$mcmc)), colnames(as.matrix(plain$mcmc)))

  # posterior means, variances about the plain fit's mean and probabilities
  # above its median agree with the fit without the block within Monte Carlo
  # error: each quantity of the declared fit is compared with the plain fit's
  # estimate in units of the combined Monte Carlo standard error of both
  plain_samples <- as.matrix(plain$mcmc)
  quantities <- unlist(lapply(c("mu_intercept", random_names), function(name){
    center <- mean(plain_samples[, name])
    .sampler_test_scalar_quantities(name, center, 0, stats::median(plain_samples[, name]), 0)
  }), recursive = FALSE)
  plain_checks <- .sampler_test_mcse_checks(plain$mcmc, quantities)
  declared_checks <- .sampler_test_mcse_checks(declared$mcmc, quantities)
  checks <- data.frame(
    quantity = plain_checks$quantity, plain = plain_checks$estimate, declared = declared_checks$estimate,
    ess_plain = plain_checks$ess, ess = declared_checks$ess,
    z = (declared_checks$estimate - plain_checks$estimate) / sqrt(plain_checks$mcse^2 + declared_checks$mcse^2)
  )
  info <- paste(utils::capture.output(print(checks, digits = 4)), collapse = "\n")
  if(isTRUE(getOption("BayesTools.sampler_test_print"))) print(checks, digits = 4)
  expect_true(all(abs(checks$z) <= 4), info = info)
  expect_true(all(checks$ess >= 100), info = info)
  expect_true(all(checks$ess_plain >= 100), info = info)
})

test_that("guarded spike-and-slab priors keep stock samplers and refuse blocks", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  expect_stock <- function(fit, guarded){
    data <- fit$data
    model <- rjags::jags.model(textConnection(fit$model), data = data, n.chains = 1, n.adapt = 0, quiet = TRUE)
    samplers <- .sampler_test_node_samplers(model)
    guarded_nodes <- names(samplers)[startsWith(names(samplers), guarded)]
    expect_gt(length(guarded_nodes), 0L)
    expect_false(any(samplers[guarded_nodes] %in% c("BayesTools::GatedSlice", "BayesTools::FactorSlice")))
    markers <- regmatches(fit$model, gregexpr("dbt_(gate|block)\\([^,]+", fit$model))[[1L]]
    expect_false(any(grepl(guarded, markers, fixed = TRUE)))
  }

  set.seed(1)
  formula_data <- data.frame(x = stats::rnorm(40), g = factor(rep(1:5, 8)))
  data <- list(y = stats::rnorm(40), se = rep(0.3, 40))
  bias <- prior_mixture(list(prior_none(), prior_PET("normal", list(0, 1))), is_null = c(TRUE, FALSE))
  random_fit <- function(...){
    JAGS_fit(
      model_syntax = "model{\n  for(i in 1:N_mu){\n    y[i] ~ dnorm(mu[i] + PET * se[i], 1)\n  }\n}",
      data = data, prior_list = list(bias = bias),
      formula_list = list(mu = ~ 1 + x + (1 | g)), formula_data_list = list(mu = formula_data),
      formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1)))),
      formula_random_prior_list = list(mu = prior_random(g = random_block(
        sd = prior_spike_and_slab(prior("normal", list(0, 1), list(0, Inf)))))),
      chains = 2, adapt = 100, burnin = 100, sample = 300, seed = 4, silent = TRUE, ...
    )
  }
  fit <- random_fit(sampler_gates = TRUE, sampler_blocks = list(c("mu_intercept", "mu_x")))
  expect_identical(attr(fit, "sampler_declarations")$gates$node, "PET_1")
  expect_stock(fit, "mu__xREx__g_intercept")
  expect_error(random_fit(sampler_blocks = list(c("mu_intercept", "mu__xREx__g_intercept"))),
               "it is a spike-and-slab prior", class = "BayesTools_sampler_declaration_error")

  ordered_data <- data.frame(f = ordered(rep(c("lo", "mid", "hi"), 10L), levels = c("lo", "mid", "hi")))
  ordered_fit <- function(...){
    JAGS_fit(
      model_syntax = "model{\n  for(i in 1:N_mu){\n    y[i] ~ dnorm(mu[i] + PET * 0.3, 1)\n  }\n}",
      data = list(y = stats::rnorm(30)), prior_list = list(bias = bias),
      formula_list = list(mu = ~ 1 + f), formula_data_list = list(mu = ordered_data),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1)),
        f = prior_ordered(prior_spike_and_slab(prior("normal", list(0, 1)), prior("bernoulli", list(0.5))))
      )),
      chains = 2, adapt = 100, burnin = 100, sample = 300, seed = 4, silent = TRUE, ...
    )
  }
  fit <- ordered_fit(sampler_gates = TRUE)
  expect_identical(attr(fit, "sampler_declarations")$gates$node, "PET_1")
  expect_stock(fit, "mu_f_ordered_total")
  expect_error(ordered_fit(sampler_blocks = list(c("mu_intercept", "mu_f"))),
               "it is an ordered prior", class = "BayesTools_sampler_declaration_error")
})

test_that("the new samplers report incomplete adaptation until their first adaptive update", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # models whose only adaptive sampler is a factor slice block or a gated node:
  # Both sampler factories initially report that adaptation is pending.
  block <- paste0("
model{
  b1 ~ dnorm(0, 1e-4)
  b2 ~ dnorm(0, 1e-4)
  z_sum ~ dnorm(b1 + b2, 100)
  z_difference ~ dnorm(b1 - b2, 1)
  ", .sampler_test_block_markers(c("b1", "b2"))$syntax, "
}")
  block_data <- c(list(z_sum = 0.5, z_difference = 0.1), .sampler_test_block_markers(c("b1", "b2"))$data)
  gated <- "
model{
  x ~ dnorm(0, 1)
  y ~ dnorm(x * off, 1)
  BayesTools_gate_1 ~ dbt_gate(x, off)
}"
  gated_data <- list(y = 0.5, off = 0, BayesTools_gate_1 = 0)
  for(case in list(list(syntax = block, data = block_data, sampler = "BayesTools::FactorSlice"),
                   list(syntax = gated, data = gated_data, sampler = "BayesTools::GatedSlice"))){
    model <- .sampler_test_compile(case$syntax, case$data, chains = 2L)
    expect_identical(names(rjags::list.samplers(model)), case$sampler)
    expect_false(rjags::adapt(model, 0, end.adaptation = FALSE, progress.bar = "none"))
    expect_true(rjags::adapt(model, 1, end.adaptation = FALSE, progress.bar = "none"))
    expect_true(rjags::adapt(model, 0, end.adaptation = TRUE, progress.bar = "none"))
  }
})

test_that("a serial extension adapts a factor slice block again", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # the narrow target (sd 0.05, correlation 0.995) needs adapted axes: along the
  # unit coordinate axes the block moves by about the conditional sd 0.005
  sd <- 0.05
  correlation <- 0.995
  data <- list(z_sum = 0.1, z_difference = 0.5,
               precision_sum = 1 / (2 * sd^2 * (1 + correlation)),
               precision_difference = 1 / (2 * sd^2 * (1 - correlation)))
  fit <- JAGS_fit(
    model_syntax = "model{\n  z_sum ~ dnorm(b1 + b2, precision_sum)\n  z_difference ~ dnorm(b1 - b2, precision_difference)\n}",
    data = data,
    prior_list = list(b1 = prior("normal", list(0, 100)), b2 = prior("normal", list(0, 100))),
    chains = 2, adapt = 500, burnin = 500, sample = 2000, seed = 8, silent = TRUE,
    sampler_blocks = list(c("b1", "b2"))
  )
  extended <- JAGS_extend(fit, autofit_control = list(
    max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = NULL,
    max_time = list(time = 60, unit = "secs"), sample_extend = 2000, restarts = 1, max_extend = 1
  ), silent = TRUE)
  expect_identical(JAGS_draw_geometry(extended)$total_draws, 8000L)
  # the effective sample size of b1 in the fit and in the 2 x 2000 extension draws
  ess <- function(fit, rows){
    sum(vapply(fit$mcmc, function(chain){
      coda::effectiveSize(coda::mcmc(as.matrix(chain)[rows, "b1"]))
    }, numeric(1)))
  }
  ess_fit <- ess(fit, 1:2000)
  ess_extension <- ess(extended, 2001:4000)
  # adapted axes give about 3,300 (both 3,419 and 3,268 measured); the unit axes
  # of an unadapted block give about 20 to 30 per 4,000 draws (the 'narrow'
  # target at adaptation 0 above, and the review's measurement)
  if(isTRUE(getOption("BayesTools.sampler_test_print"))) print(c(fit = ess_fit, extension = ess_extension))
  expect_gt(ess_fit, 1000)
  expect_gt(ess_extension, 1000)
})

test_that("a NaN or infinite density at a factor slice proposal stops the fit", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # A NaN (or +Inf) density is a failure to calculate, not zero density:
  # rejecting such proposals would silently truncate the target (in model 2
  # the density grows without bound where its computation overflows). The
  # block stops with a node error, as the stock samplers do without the block.
  run <- function(model, block){
    markers <- if(block) .sampler_test_block_markers(model$variables) else list(syntax = "", data = list())
    compiled <- rjags::jags.model(
      textConnection(sprintf(model$syntax, markers$syntax)), data = c(model$data, markers$data),
      n.chains = 4L, n.adapt = 0, quiet = TRUE,
      inits = lapply(1:4, function(chain){
        c(model$values, list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 40 + chain))
      })
    )
    if(block){
      expect_identical(names(rjags::list.samplers(compiled)), "BayesTools::FactorSlice")
    }
    rjags::adapt(compiled, 500, end.adaptation = TRUE, progress.bar = "none")
    stats::update(compiled, 1000, progress.bar = "none")
    rjags::coda.samples(compiled, model$variables, n.iter = 10000, progress.bar = "none")
  }
  models <- list(
    # (1) exp(1000 b1) - exp(1000 b1) is NaN where exp overflows (b1 > 0.7098)
    exp_difference = list(syntax = "
model{
  b1 ~ dnorm(0, 1)
  b2 ~ dnorm(0, 1)
  m <- exp(1000 * b1) - exp(1000 * b1)
  y ~ dnorm(m + b2, 1)
  %s
}", data = list(y = 0.8), values = list(b1 = 0, b2 = 0), variables = c("b1", "b2"),
      # without the block: base::RealSlicer on b1
      stock = "Failure to calculate log density"),
    # (2) y = 0.5 ~ Beta(s, s), s = exp(1000 x): the density grows like
    # exp(500 x) until its calculation overflows
    beta_overflow = list(syntax = "
model{
  x ~ dnorm(0, 1)
  b ~ dnorm(0, 1)
  shape <- exp(1000 * x)
  y ~ dbeta(shape, shape)
  %s
}", data = list(y = 0.5), values = list(x = 0, b = 0), variables = c("x", "b"),
      # without the block: base::RealSlicer on x
      stock = "Slicer stuck at value with infinite density")
  )
  for(name in names(models)){
    expect_error(run(models[[name]], block = TRUE),
                 "Factor slice sampler: failure to calculate the log density of the block at a proposal",
                 fixed = TRUE, info = name)
    expect_error(run(models[[name]], block = FALSE), models[[name]]$stock, fixed = TRUE, info = name)
  }
})

test_that("an infinite density of the block at initialization stops the fit", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # At b2 = 0.8 the precision overflows. With b1 = y, the observed node's
  # density is +Inf and JAGS 5 rejects the supplied starting state.
  markers <- .sampler_test_block_markers(c("b1", "b2"))
  con <- textConnection(sprintf("
model{
  b1 ~ dnorm(0, 1)
  b2 ~ dnorm(0, 1)
  precision <- exp(1000 * b2)
  y ~ dnorm(b1, precision)
  %s
}", markers$syntax))
  on.exit(close(con))
  expect_error(
    rjags::jags.model(con, data = c(list(y = 0.8), markers$data), n.chains = 1L, n.adapt = 0L,
      quiet = TRUE, inits = list(b1 = 0.8, b2 = 0.8, .RNG.name = "base::Mersenne-Twister", .RNG.seed = 3L)),
    "Error in node y\nInfinite log density", fixed = TRUE
  )
})

test_that("a gate whose activity depends on its node through a stochastic child is not claimed", {
  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  # p(x | active = 0) is proportional to dnorm(x) pnorm(-x), the skew-normal
  # with shape -1, whose mean is -1 / sqrt(pi) = -0.5642
  syntax <- "
model{
  x ~ dnorm(0, 1)
  active ~ dbern(phi(x))
  BayesTools_gate_1 ~ dbt_gate(x, active)
}"
  model <- .sampler_test_compile(syntax, list(active = 0, BayesTools_gate_1 = 0), chains = 4L, seed = 6L)
  samplers <- .sampler_test_node_samplers(model)
  expect_false(samplers[["x"]] %in% c("BayesTools::GatedSlice", "BayesTools::FactorSlice"))
  rjags::adapt(model, 500, end.adaptation = TRUE, progress.bar = "none")
  stats::update(model, 1000, progress.bar = "none")
  samples <- rjags::coda.samples(model, "x", n.iter = 10000, progress.bar = "none")
  checks <- .sampler_test_mcse_checks(samples, list(
    x_mean = list(series = function(x) x[, "x"], target = -1 / sqrt(pi))
  ))
  .sampler_test_expect_checks(checks, "stochastic gate ancestry")
})
