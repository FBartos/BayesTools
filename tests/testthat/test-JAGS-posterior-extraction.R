skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Posterior Density Extraction Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for posterior extraction helper functions including
#   .extract_posterior_samples and .remove_auxiliary_parameters.
#
# DEPENDENCIES:
#   - rjags, runjags, coda: For JAGS model handling
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("rjags"), skip_if_not_installed("runjags")
#   - Note: Creates mock objects, does not need pre-fitted models
#
# MODELS/FIXTURES:
#   - Creates mock runjags objects for testing
#
# TAGS: @evaluation, @JAGS, @posterior-extraction
# ============================================================================ #

.make_mock_mcarray <- function(values, dim, varname = NULL, iterations = NULL){
  parameter <- array(values, dim = dim)
  class(parameter) <- "mcarray"
  if(!is.null(varname)){
    attr(parameter, "varname") <- varname
  }
  if(!is.null(iterations)){
    attr(parameter, "iterations") <- iterations
  }
  parameter
}

# Tests for posterior extraction helper functions
test_that(".extract_posterior_samples extracts samples correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  
  # Load runjags to ensure S3 methods are registered
  suppressPackageStartupMessages(library(runjags))

  # Create a proper runjags object structure for testing
  # The runjags package has an as.mcmc method that handles mcmc.list objects
  set.seed(123)
  mcmc1 <- coda::mcmc(matrix(rnorm(100), ncol = 1, dimnames = list(NULL, "mu")),
                      start = 1, end = 100, thin = 1)
  mcmc2 <- coda::mcmc(matrix(rnorm(100), ncol = 1, dimnames = list(NULL, "mu")),
                      start = 1, end = 100, thin = 1)
  mcmc_list <- coda::mcmc.list(mcmc1, mcmc2)
  
  # Create a minimal runjags object
  fit <- structure(
    list(mcmc = mcmc_list),
    class = c("runjags", "list")
  )

  # Test matrix extraction (as_list = FALSE)
  # This calls coda::as.mcmc on the runjags object which returns an mcmc object
  samples_matrix <- BayesTools:::.extract_posterior_samples(fit, as_list = FALSE)
  # mcmc objects inherit from matrix
  expect_true(inherits(samples_matrix, "mcmc"))
  expect_equal(ncol(samples_matrix), 1)
  expect_true("mu" %in% colnames(samples_matrix))
  expect_equal(nrow(samples_matrix), 200)  # 100 samples x 2 chains merged

  # Test list extraction (as_list = TRUE)
  samples_list <- BayesTools:::.extract_posterior_samples(fit, as_list = TRUE)
  expect_true(inherits(samples_list, "mcmc.list"))
  expect_equal(length(samples_list), 2) # 2 chains
})

test_that(".fit_to_posterior flattens multidimensional mcarrays in coda order", {

  iterations <- c(start = 1, end = 3, thin = 1)
  theta <- .make_mock_mcarray(
    values = seq_len(2 * 3 * 3 * 2),
    dim = c(2, 3, 3, 2),
    varname = "theta",
    iterations = iterations
  )
  psi <- .make_mock_mcarray(
    values = 1000 + seq_len(2 * 2 * 2 * 3 * 2),
    dim = c(2, 2, 2, 3, 2),
    varname = "psi",
    iterations = iterations
  )
  alpha <- .make_mock_mcarray(
    values = 2000 + seq_len(1 * 3 * 2),
    dim = c(1, 3, 2),
    varname = "alpha",
    iterations = iterations
  )

  posterior <- BayesTools:::.fit_to_posterior(list(
    theta = theta,
    psi = psi,
    alpha = alpha
  ))

  expected_theta <- do.call(rbind, lapply(seq_len(2), function(chain){
    do.call(rbind, lapply(seq_len(3), function(iteration){
      as.vector(theta[, , iteration, chain])
    }))
  }))
  colnames(expected_theta) <- c(
    "theta[1,1]", "theta[2,1]",
    "theta[1,2]", "theta[2,2]",
    "theta[1,3]", "theta[2,3]"
  )

  expected_psi <- do.call(rbind, lapply(seq_len(2), function(chain){
    do.call(rbind, lapply(seq_len(3), function(iteration){
      as.vector(psi[, , , iteration, chain])
    }))
  }))
  colnames(expected_psi) <- c(
    "psi[1,1,1]", "psi[2,1,1]",
    "psi[1,2,1]", "psi[2,2,1]",
    "psi[1,1,2]", "psi[2,1,2]",
    "psi[1,2,2]", "psi[2,2,2]"
  )

  expected_alpha <- do.call(rbind, lapply(seq_len(2), function(chain){
    do.call(rbind, lapply(seq_len(3), function(iteration){
      alpha[1, iteration, chain]
    }))
  }))
  colnames(expected_alpha) <- "alpha"

  expect_equal(
    posterior,
    cbind(expected_theta, expected_psi, expected_alpha)
  )
})

test_that(".fit_to_posterior preserves monitored mcarray indices", {

  iterations <- c(start = 10, end = 11, thin = 1)
  theta <- .make_mock_mcarray(
    values = seq_len(4),
    dim = c(2, 2, 1),
    varname = "theta[2,4:5]",
    iterations = iterations
  )
  phi <- .make_mock_mcarray(
    values = 10 + seq_len(4),
    dim = c(2, 2, 1),
    varname = "phi[2:3,4]",
    iterations = iterations
  )
  eta <- .make_mock_mcarray(
    values = 20 + seq_len(4),
    dim = c(2, 2, 1),
    iterations = iterations
  )

  posterior <- BayesTools:::.fit_to_posterior(list(
    theta_subset = theta,
    phi_subset = phi,
    eta = eta
  ))

  expect_identical(
    colnames(posterior),
    c(
      "theta[2,4]", "theta[2,5]",
      "phi[2,4]", "phi[3,4]",
      "eta[1]", "eta[2]"
    )
  )
})

test_that(".fit_to_posterior rejects misaligned mcarray draws", {

  by_chain <- .make_mock_mcarray(
    values = seq_len(6),
    dim = c(1, 3, 2),
    varname = "x"
  )
  different_layout <- .make_mock_mcarray(
    values = seq_len(6),
    dim = c(1, 2, 3),
    varname = "y"
  )
  expect_error(
    BayesTools:::.fit_to_posterior(list(
      x = by_chain,
      y = different_layout
    )),
    "matching iteration and chain dimensions",
    fixed = TRUE
  )

  first_iterations <- .make_mock_mcarray(
    values = seq_len(6),
    dim = c(1, 3, 2),
    varname = "x",
    iterations = c(start = 1, end = 3, thin = 1)
  )
  shifted_iterations <- .make_mock_mcarray(
    values = seq_len(6),
    dim = c(1, 3, 2),
    varname = "y",
    iterations = c(start = 2, end = 4, thin = 1)
  )
  expect_error(
    BayesTools:::.fit_to_posterior(list(
      x = first_iterations,
      y = shifted_iterations
    )),
    "matching iteration metadata",
    fixed = TRUE
  )

  malformed <- matrix(seq_len(6), nrow = 2)
  class(malformed) <- "mcarray"
  attr(malformed, "varname") <- "z"
  expect_error(
    BayesTools:::.fit_to_posterior(list(z = malformed)),
    "at least one parameter dimension",
    fixed = TRUE
  )
})


test_that(".remove_auxiliary_parameters keeps user 'inv_' columns of inverse-gamma priors", {

  # Inverse-gamma priors are monitored by the parameter itself, so an
  # 'inv_<parameter>' column is a node of the user's model syntax (e.g. from
  # add_parameters), neither removed nor treated as a BayesTools 0.3.0
  # precision coordinate. Fits of 0.3.0 stop earlier at the fit-contract checks.
  model_samples <- matrix(rnorm(100), ncol = 2)
  colnames(model_samples) <- c("sigma", "inv_sigma")
  prior_list <- list(sigma = prior("invgamma", list(1, 1)))
  result <- BayesTools:::.remove_auxiliary_parameters(model_samples, prior_list, NULL)
  expect_identical(result$model_samples, model_samples)

  model_samples <- matrix(rnorm(400), ncol = 4)
  colnames(model_samples) <- c("theta[1]", "theta[2]", "inv_theta[1]", "inv_theta[2]")
  theta_prior <- prior_factor("invgamma", list(2, 1), contrast = "independent")
  theta_prior$parameters$K <- 2
  result <- BayesTools:::.remove_auxiliary_parameters(model_samples, list(theta = theta_prior), NULL)
  expect_identical(result$model_samples, model_samples)
})

test_that(".remove_auxiliary_parameters removes vector prior columns by base name", {

  model_samples <- matrix(
    seq_len(6),
    nrow = 1,
    dimnames = list(
      NULL,
      c("w[1]", "w[2]", "w[3]", "prior_par_eta_w[1]", "prior_par_eta_w[2]", "prior_par_eta_w[3]")
    )
  )
  prior_list <- list(w = prior("dirichlet", list(alpha = c(1, 2, 3))))

  result <- BayesTools:::.remove_auxiliary_parameters(
    model_samples,
    prior_list,
    remove_parameters = "w"
  )

  expect_equal(ncol(result$model_samples), 0L)
  expect_equal(names(result$prior_list), character())
})


test_that(".remove_auxiliary_parameters drops p-hacking internals by report scale and omega request", {

  model_samples <- matrix(seq_len(50), ncol = 5)
  colnames(model_samples) <- c("omega[1]", "omega[2]", "alpha", "pi_null", "phack_kind")

  phacking_prior <- prior_phacking(form = "linear", report_scale = "alpha")

  result <- BayesTools:::.remove_auxiliary_parameters(
    model_samples,
    list(ph = phacking_prior),
    remove_parameters = "omega"
  )

  expect_equal(colnames(result$model_samples), c("alpha", "phack_kind"))
  expect_true("ph" %in% names(result$prior_list))
})


test_that(".remove_auxiliary_parameters renames selection omegas and drops unreported p-hacking columns", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking <- prior_phacking(form = "linear", report_scale = "alpha")
  bias <- prior_bias(selection, phacking)

  model_samples <- matrix(seq_len(50), ncol = 5)
  colnames(model_samples) <- c("omega[1]", "omega[2]", "alpha", "pi_null", "phack_kind")

  result <- BayesTools:::.remove_auxiliary_parameters(
    model_samples,
    list(pub_bias = bias),
    remove_parameters = NULL
  )

  expect_equal(colnames(result$model_samples), c("omega[0,0.025]", "omega[0.025,1]", "alpha", "phack_kind"))
  expect_true("pub_bias" %in% names(result$prior_list))
})


test_that(".remove_auxiliary_parameters drops cumulative-weight auxiliaries from old and new fits", {

  selection <- prior_weightfunction(
    "one-sided",
    .025,
    wf_cumulative(c(1, 2))
  )
  model_samples <- matrix(seq_len(50), ncol = 5)
  colnames(model_samples) <- c(
    "omega[1]", "omega[2]", "eta[1]", "eta[2]", "omega_ratio"
  )

  result <- BayesTools:::.remove_auxiliary_parameters(
    model_samples,
    list(pub_bias = selection),
    remove_parameters = NULL
  )

  expect_equal(
    colnames(result$model_samples),
    c("omega[0,0.025]", "omega[0.025,1]")
  )
})


test_that(".process_spike_and_slab handles conditional samples", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples with spike and slab
  model_samples <- matrix(c(
    rnorm(50, 0, 1),  # mu values
    rep(1, 50)        # indicator (all in slab)
  ), ncol = 2)
  colnames(model_samples) <- c("mu", "mu_indicator")

  prior_list <- list(
    mu = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("spike", list(0.5))
    )
  )

  result <- BayesTools:::.process_spike_and_slab(
    model_samples, prior_list, "mu",
    conditional = TRUE, remove_inclusion = FALSE, warnings = NULL
  )

  expect_true("mu (inclusion)" %in% colnames(result$model_samples))
  expect_false("mu_indicator" %in% colnames(result$model_samples))
  expect_true(is.prior.simple(result$prior_list$mu))
})


test_that(".apply_parameter_transformations applies transformations", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples
  model_samples <- matrix(rnorm(100, 0, 1), ncol = 1)
  colnames(model_samples) <- "mu"

  prior_list <- list(
    mu = prior("normal", list(0, 1))
  )

  # Apply exp transformation
  transformations <- list(
    mu = list(fun = exp, arg = list())
  )

  result <- BayesTools:::.apply_parameter_transformations(
    model_samples, transformations, prior_list
  )

  expect_true(all(result[, "mu"] > 0))  # exp makes all values positive
  expect_equal(ncol(result), 1)
})

test_that("requested transformations reach independent and ordered factor coefficients", {

  data <- data.frame(
    g = factor(c("A", "B", "C"), levels = c("A", "B", "C")),
    o = ordered(c("lo", "mid", "hi"), levels = c("lo", "mid", "hi"))
  )
  prior_list <- JAGS_formula(
    ~ g + o, "mu", data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      g         = prior_factor("normal", list(0, 1), contrast = "independent"),
      o         = prior_ordered(prior("normal", list(0, 1)), allocation = c(.4, .6))
    )
  )$prior_list
  set.seed(8)
  model_samples <- matrix(rnorm(5 * 20), ncol = 5)
  colnames(model_samples) <- c("mu_g[1]", "mu_g[2]", "mu_g[3]", "mu_o[1]", "mu_o[2]")
  transformations <- list(
    mu_g = list(fun = exp, arg = list()),
    mu_o = list(fun = exp, arg = list())
  )

  coefficients <- BayesTools:::.apply_parameter_transformations(
    model_samples, transformations, prior_list, transform_factors = FALSE
  )
  expect_equal(coefficients, exp(model_samples))

  # with transform_factors = TRUE, ordered coefficients are transformed only
  # after the contrast transformation to level effects
  levels <- BayesTools:::.apply_parameter_transformations(
    model_samples, transformations, prior_list, transform_factors = TRUE
  )
  expect_equal(levels[, 1:3], exp(model_samples[, 1:3]))
  expect_equal(levels[, 4:5], model_samples[, 4:5])
  level_effects <- BayesTools:::.transform_factor_contrasts(
    levels, prior_list, transform_factors = TRUE, transformations = transformations
  )
  expect_equal(
    unname(level_effects[, ncol(level_effects)]),
    exp(model_samples[, "mu_o[1]"] + model_samples[, "mu_o[2]"])
  )
})


test_that("treatment factor coordinates are named by their level cells", {

  prior_obj <- prior_factor_levels(
    prior_factor("normal", list(0, 1), contrast = "treatment"),
    c("A", "B", "C", "D")
  )

  expect_identical(
    BayesTools:::.bt_label_prior_column_names("group", prior_obj),
    c("group[B]", "group[C]", "group[D]")
  )
})


test_that("contrast coefficients are named with braces", {

  # With numeric level labels, `[j]` would read as the level labelled j.
  data <- data.frame(
    m = factor(rep(c(5, 10, 20), 2), levels = c(5, 10, 20)),
    o = ordered(rep(c(5, 10, 20), 2), levels = c(5, 10, 20))
  )
  formula_result <- JAGS_formula(~ m + o, "mu", data, list(
    intercept = prior("normal", list(0, 1)),
    m = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
    o = prior_ordered(prior("normal", list(0, 1)))
  ))
  parts <- unlist(lapply(names(formula_result$prior_list), function(name){
    BayesTools:::.bt_label_parts_prior(name, formula_result$prior_list[[name]])$parts
  }), recursive = FALSE)

  # Mean-difference coordinates are coefficients; the first cumulative
  # ordered coordinate is level "10" and the second is an increment.
  expect_identical(
    BayesTools:::.bt_label(parts, "selector"),
    c("mu_intercept", "mu_m{1}", "mu_m{2}", "mu_o[10]", "mu_o{2}")
  )
  expect_identical(
    BayesTools:::.bt_label(parts, "table"),
    c("(mu) intercept", "(mu) m{1}", "(mu) m{2}", "(mu) o[10]", "(mu) o{2}")
  )

  # Ordinary factor priors label their levels 1..K by construction.
  ordinary <- prior_factor_levels(
    prior_factor("mnormal", list(0, 1), contrast = "meandif"), 3L
  )
  expect_identical(
    BayesTools:::.bt_label_prior_column_names("p1", ordinary),
    c("p1{1}", "p1{2}")
  )
  # A two-level contrast has one coefficient.
  two_level <- prior_factor_levels(
    prior_factor("mnormal", list(0, 1), contrast = "orthonormal"), 2L
  )
  expect_identical(
    BayesTools:::.bt_label_prior_column_names("p2", two_level),
    "p2{1}"
  )
})


test_that("treatment interactions with ordered factors label increments as coefficients", {

  levels <- c("5", "10", "20")
  data <- data.frame(
    f = factor(rep(c("a", "b"), 6)),
    o = ordered(rep(levels, each = 4), levels = levels)
  )
  formula_result <- JAGS_formula(~ f * o, "mu", data, list(
    intercept = prior("normal", list(0, 1)),
    f = prior_factor("normal", list(0, 1), contrast = "treatment"),
    o = prior_ordered(prior("normal", list(0, 1))),
    "f:o" = prior_factor("normal", list(0, 1), contrast = "treatment")
  ))
  columns <- unlist(lapply(names(formula_result$prior_list), function(name){
    prior <- formula_result$prior_list[[name]]
    if(BayesTools:::.bt_prior_is_factor_family(prior)){
      BayesTools:::.JAGS_prior_factor_names(name, prior)
    }else{
      name
    }
  }), use.names = FALSE)
  set.seed(1)
  samples <- matrix(stats::rnorm(20L * length(columns)), nrow = 20L,
                    dimnames = list(NULL, columns))

  # The interaction design is treatment by cumulative ordered coding:
  # coordinate 1 is the cell (b, 10) and coordinate 2 is the increment from
  # (b, 10) to (b, 20), so it is coefficient 2, never the cell "[20]".
  parts <- BayesTools:::.bt_label_parts_prior(
    "mu_f__xXx__o",
    formula_result$prior_list[["mu_f__xXx__o"]]
  )$parts
  expect_identical(
    BayesTools:::.bt_label(parts, "table"),
    c("(mu) f[b]:o[10]", "(mu) f:o{2}")
  )

  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = 20L),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- attach_test_parameter_map(fit)
  mixed <- as_mixed_posteriors(fit, parameters = "mu_f__xXx__o")
  expect_identical(colnames(mixed$mu_f__xXx__o),
                   c("mu_f__xXx__o[f=b, o=10]", "mu_f__xXx__o{2}"))

  # Every displayed interaction row selects the quantity that produced it
  # (the table summarizes the same draws: agreement up to rounding).
  catalog <- parameter_catalog(fit)
  table <- JAGS_estimates_table(fit)
  rows <- grep("f[^ ]*:o", rownames(table), value = TRUE)
  expect_identical(rows, c("(mu) f[b]:o[10]", "(mu) f:o{2}"))
  for(row in rows){
    selection <- parameter_catalog_resolve(catalog, row)
    expect_equal(mean(as.matrix(parameter_draws(fit, selection))),
                 table[row, "Mean"], tolerance = 1e-10, info = row)
  }
  # The cell (b, 20) itself is the sum of both coordinates.
  cell <- parameter_catalog_resolve(catalog, "mu_f__xXx__o[f=b, o=20]")
  expect_equal(
    as.numeric(as.matrix(parameter_draws(fit, cell))),
    unname(samples[, "mu_f__xXx__o[1]"] + samples[, "mu_f__xXx__o[2]"]),
    tolerance = 1e-14
  )
})


test_that("factor interaction labels keep each level on its factor component", {

  # labels of the term coefficients (a formula without an intercept keeps
  # its intercept fixed at zero)
  label_rows <- function(formula, data, prior_list){
    formula_result <- JAGS_formula(formula, "mu", data, prior_list)
    terms <- setdiff(names(formula_result$prior_list), "mu_intercept")
    unlist(lapply(terms, function(name){
      BayesTools:::.bt_label(
        BayesTools:::.bt_label_parts_prior(name, formula_result$prior_list[[name]])$parts,
        "table"
      )
    }), use.names = FALSE)
  }

  data <- data.frame(
    alloc = factor(rep(c("random", "systematic"), 4)),
    year  = seq(-1, 1, length.out = 8)
  )
  expect_identical(
    label_rows(~ 0 + alloc + alloc:year, data, list(
      alloc        = prior_factor("normal", list(0, 1), contrast = "independent"),
      "alloc:year" = prior_factor("normal", list(0, 1), contrast = "independent")
    )),
    c(
      "(mu) alloc[random]",
      "(mu) alloc[systematic]",
      "(mu) alloc[random]:year",
      "(mu) alloc[systematic]:year"
    )
  )

  # multi-factor interactions with a covariate between the factors
  data <- data.frame(
    a    = factor(rep(c("a1", "a2"), 6)),
    b    = factor(rep(c("b1", "b2", "b3"), each = 4)),
    year = seq(-1, 1, length.out = 12)
  )
  expect_identical(
    label_rows(~ 0 + a:year:b, data, list(
      "a:year:b" = prior_factor("normal", list(0, 1), contrast = "independent")
    )),
    c(
      "(mu) a[a1]:year:b[b1]",
      "(mu) a[a2]:year:b[b1]",
      "(mu) a[a1]:year:b[b2]",
      "(mu) a[a2]:year:b[b2]",
      "(mu) a[a1]:year:b[b3]",
      "(mu) a[a2]:year:b[b3]"
    )
  )
  treatment <- function() prior_factor("normal", list(0, 1), contrast = "treatment")
  rows <- label_rows(~ a * year * b, data, list(
    intercept  = prior("normal", list(0, 1)),
    a          = treatment(),
    b          = treatment(),
    year       = prior("normal", list(0, 1)),
    "a:year"   = treatment(),
    "a:b"      = treatment(),
    "year:b"   = treatment(),
    "a:year:b" = treatment()
  ))
  expect_identical(
    grep("^\\(mu\\) a\\[[^]]*\\]:year:b", rows, value = TRUE),
    c("(mu) a[a2]:year:b[b2]", "(mu) a[a2]:year:b[b3]")
  )
})

test_that("random factor interaction SDs are named by their level cells", {

  data <- data.frame(
    f = factor(rep(c("a", "b", "c"), 4), levels = c("a", "b", "c")),
    g = factor(rep(c("u", "v"), each = 6), levels = c("u", "v")),
    id = factor(rep(c("s1", "s2", "s3", "s4"), each = 3))
  )
  sd_labels <- function(random_formula){
    result <- JAGS_formula(
      formula = random_formula,
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
    )
    sd_names <- result$formula_design$random_effects[[1L]]$sd_parameter_names
    coordinates <- build_test_parameter_coordinates(
      columns = c("mu_intercept", sd_names),
      prior_list = result$prior_list,
      formula_design = list(mu = result$formula_design)
    )
    coordinates$display_label[match(sd_names, coordinates$coordinate_name)]
  }

  expect_identical(
    sd_labels(~ 1 + diag(0 + f:g | id)),
    c(
      "(mu) id: sd(f[a]:g[u])",
      "(mu) id: sd(f[b]:g[u])",
      "(mu) id: sd(f[c]:g[u])",
      "(mu) id: sd(f[a]:g[v])",
      "(mu) id: sd(f[b]:g[v])",
      "(mu) id: sd(f[c]:g[v])"
    )
  )
  expect_identical(
    sd_labels(~ 1 + diag(0 + f * g | id)),
    c(
      "(mu) id: sd(f[b])",
      "(mu) id: sd(f[c])",
      "(mu) id: sd(g[v])",
      "(mu) id: sd(f[b]:g[v])",
      "(mu) id: sd(f[c]:g[v])"
    )
  )
})


test_that("prior_factor_levels() sets the factor metadata of formula factor terms", {

  level_names <- c("low", "mid", "high")
  data <- data.frame(f = factor(rep(level_names, 2), levels = level_names))
  designs <- list(
    treatment   = stats::contr.treatment(3),
    independent = contr.independent(level_names),
    meandif     = contr.meandif(level_names),
    orthonormal = contr.orthonormal(level_names)
  )
  for(contrast in names(designs)){
    distribution <- if(contrast %in% c("meandif", "orthonormal")) "mnormal" else "normal"
    prior <- prior_factor(distribution, list(0, 1), contrast = contrast)

    named <- prior_factor_levels(prior, level_names)
    expect_identical(attr(named, "levels"), 3L, info = contrast)
    expect_identical(attr(named, "level_names"), level_names, info = contrast)
    expect_equal(attr(named, "factor_design"), unname(designs[[contrast]]),
                 ignore_attr = TRUE, info = contrast)
    expect_identical(attr(named, "factor_cell_names"), level_names, info = contrast)

    # bound to its parameter, the prior has the design of the formula term
    bound <- .complete_factor_metadata(named, "mu_f")
    formula_prior <- JAGS_formula(
      ~ f, "mu", data,
      list(intercept = prior("normal", list(0, 1)), f = prior)
    )$prior_list$mu_f
    expect_equal(attr(bound, "factor_design"), attr(formula_prior, "factor_design"),
                 ignore_attr = TRUE, info = contrast)
    expect_identical(attr(bound, "factor_terms"), "mu_f", info = contrast)
    expect_identical(unname(attr(bound, "factor_contrasts")),
                     unname(attr(formula_prior, "factor_contrasts")), info = contrast)
    expect_identical(names(attr(bound, "factor_contrasts")), "mu_f", info = contrast)

    # levels given by their number are named 1..K
    counted <- prior_factor_levels(prior, 3)
    expect_identical(attr(counted, "level_names"), c("1", "2", "3"), info = contrast)
    expect_identical(attr(counted, "factor_design"), attr(named, "factor_design"), info = contrast)
  }

  # the components of a spike-and-slab factor prior carry the same metadata
  spike <- prior_factor_levels(
    prior_spike_and_slab(prior_factor("mnormal", list(0, 1), contrast = "meandif")),
    level_names
  )
  for(component in c(list(spike), as.list(spike))){
    expect_equal(attr(component, "factor_design"), unname(designs$meandif), ignore_attr = TRUE)
    expect_identical(attr(component, "level_names"), level_names)
  }
  expect_identical(
    attr(.complete_factor_metadata(spike, "p1")[[1L]], "factor_terms"),
    "p1"
  )

  # an ordered prior is bound to the fitted parameter's nodes
  ordered <- prior_factor_levels(prior_ordered(prior("normal", list(0, 1))), level_names)
  expect_no_warning(density(ordered, n_points = 11))
  bound <- .complete_factor_metadata(ordered, "mu_o")
  expect_identical(attr(bound, "ordered_metadata")$parameter_name, "mu_o")
  expect_equal(attr(bound, "coefficient_dim"), 2)

  expect_error(prior_factor_levels(prior("normal", list(0, 1)), 3),
               "'prior' must be a factor prior distribution.", fixed = TRUE)
  expect_error(prior_factor_levels(prior_factor("normal", list(0, 1), contrast = "treatment"), 1),
               "A factor prior with 'contr.treatment' contrasts requires at least two levels.",
               fixed = TRUE)
  expect_identical(
    attr(prior_factor_levels(prior_factor("normal", list(0, 1), contrast = "independent"), 1),
         "factor_design"),
    matrix(1)
  )
  expect_error(prior_factor_levels(prior_factor("normal", list(0, 1), contrast = "treatment"),
                                   c("a", "a")),
               "The 'levels' argument must contain unique, nonempty level names.", fixed = TRUE)
  expect_error(prior_factor_levels(prior_factor("normal", list(0, 1), contrast = "treatment"), NA),
               "The 'levels' argument cannot contain NA/NaN values.", fixed = TRUE)
})

test_that("factor priors without their factor levels stop instead of being counted", {

  skip_if_not_installed("runjags")
  message <- paste0(
    "The factor prior of 'p1' has no complete factor-level metadata. Set its ",
    "levels with 'prior_factor_levels()', or specify the factor in a formula."
  )
  counted <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  attr(counted, "levels") <- 3
  for(prior in list(counted, prior_spike_and_slab(counted))){
    expect_error(
      JAGS_fit("model{}", prior_list = list(p1 = prior), chains = 1,
               adapt = 50, burnin = 50, sample = 100, seed = 1),
      message, fixed = TRUE
    )
  }

  # coordinate names come from the factor design, never from the number of
  # coordinates (a treatment interaction with as many coordinates as cells was
  # named by the full cell grid, otherwise by the grid beyond the references)
  treatment <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(treatment, "levels") <- 3
  attr(treatment, "level_names") <- c("a", "b", "c")
  expect_error(.bt_label_prior_column_names("p1", treatment), message, fixed = TRUE)
  expect_identical(
    .bt_label_prior_column_names("p1", prior_factor_levels(treatment, c("a", "b", "c"))),
    c("p1[b]", "p1[c]")
  )
})

test_that("random-effect SD factor priors carry the design of their term", {

  data <- data.frame(
    g = factor(rep(c("a", "b", "c"), 4)),
    h = factor(rep(c("u", "v"), each = 6)),
    x = seq(-1, 1, length.out = 12),
    id = factor(rep(1:4, each = 3))
  )
  sd_priors <- function(formula){
    prior_list <- JAGS_formula(
      formula, "mu", data, list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(
        sd = prior("normal", list(0, 1), list(0, Inf))
      ))
    )$prior_list
    prior_list <- .complete_factor_metadata_prior_list(prior_list)
    prior_list[vapply(prior_list, is.prior.factor, logical(1))]
  }
  coordinate_names <- function(prior_list){
    lapply(names(prior_list), function(parameter){
      .bt_label_prior_column_names(parameter, prior_list[[parameter]])
    })
  }

  # treatment slopes code the levels after the reference level
  treatment <- sd_priors(~ 1 + (1 + g || id))
  expect_equal(attr(treatment$mu__xREx__id_g, "factor_design"),
               unname(stats::contr.treatment(3)), ignore_attr = TRUE)
  expect_identical(coordinate_names(treatment),
                   list(c("mu__xREx__id_g[b]", "mu__xREx__id_g[c]")))

  # factor-by-continuous slopes code every level (the completion from the
  # coordinate count failed with "invalid 'times' value" for this term)
  by_level <- sd_priors(~ 1 + (1 + g:x || id))
  expect_equal(attr(by_level$mu__xREx__id_g__xXx__x, "factor_design"), diag(3))
  expect_identical(coordinate_names(by_level), list(c(
    "mu__xREx__id_g[a]__xXx__x", "mu__xREx__id_g[b]__xXx__x", "mu__xREx__id_g[c]__xXx__x"
  )))

  # an interaction without its main effects codes every cell, one with them
  # the cells beyond both reference levels
  cells <- sd_priors(~ 1 + (0 + g:h || id))
  expect_equal(attr(cells$mu__xREx__id_g__xXx__h, "factor_design"), diag(6))
  hierarchical <- sd_priors(~ 1 + diag(0 + g * h | id))
  expect_equal(
    attr(hierarchical$mu__xREx__id_g__xXx__h, "factor_design"),
    kronecker(stats::contr.treatment(2), stats::contr.treatment(3)),
    ignore_attr = TRUE
  )
  expect_identical(coordinate_names(hierarchical)[[3L]],
                   c("mu__xREx__id_g[b]__xXx__h[v]", "mu__xREx__id_g[c]__xXx__h[v]"))
  for(prior in c(treatment, by_level, cells, hierarchical)){
    expect_true(.bt_factor_metadata_complete(prior))
  }
})

test_that(".format_factor_level_parameter_names rejects interaction metadata mismatch", {

  expect_error(
    BayesTools:::.format_factor_level_parameter_names(
      "mu_a__xXx__year__xXx__b",
      list(a = c("a1", "a2"), missing = c("m1", "m2", "m3")),
      n_parameters = 4
    ),
    "factor metadata do not match"
  )
})

test_that(".generate_prior_sample_matrix errors on unsupported prior RNGs", {

  unsupported_prior <- list(distribution = "unsupported")
  class(unsupported_prior) <- c("prior", "prior.unsupported")

  expect_error(
    BayesTools:::.generate_prior_sample_matrix(
      list(theta = unsupported_prior),
      n_samples = 4
    ),
    "Could not generate samples for prior 'theta'",
    fixed = TRUE
  )
})

test_that("random SD unscaling ignores obsolete term-name metadata", {

  posterior <- matrix(
    1,
    nrow = 2,
    ncol = 2,
    dimnames = list(NULL, c("mu__xREx__study_x[1]", "mu__xREx__study_x[2]"))
  )
  formula_scale <- list(mu_x = list(mean = 0, sd = 2))
  attr(formula_scale, "random_effect_terms") <- c("mu__xREx__study_x" = "x")

  expect_identical(
    BayesTools:::.apply_random_sd_unscale(
      posterior = posterior,
      random_sd_cols = colnames(posterior),
      formula_scale = formula_scale,
      prefix = "mu"
    ),
    posterior
  )
})


test_that(".transform_factor_contrasts transforms orthonormal to differences", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create mock samples with orthonormal contrasts
  model_samples <- matrix(rnorm(300), ncol = 3)
  colnames(model_samples) <- c("group[1]", "group[2]", "group[3]")

  # Create a factor prior with levels attribute (as would be set by JAGS_formula)
  prior_obj <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  attr(prior_obj, "levels") <- 4  # 4 levels total (orthonormal has K-1 parameters for K levels)
  attr(prior_obj, "level_names") <- c("A", "B", "C", "D")  # Should be a vector, not a list
  attr(prior_obj, "factor_terms") <- ".factor"
  attr(prior_obj, "factor_contrasts") <- c(.factor = "contr.orthonormal")
  
  prior_list <- list(group = prior_obj)

  expect_silent(
    result <- BayesTools:::.transform_factor_contrasts(
      model_samples, prior_list, transform_factors = TRUE
    )
  )

  expect_message(
    BayesTools:::.transform_factor_contrasts(
      model_samples,
      prior_list,
      transform_factors = TRUE,
      transformations = list(group = list(fun = exp, arg = list()))
    ),
    "transformation was applied"
  )

  # Should have 4 columns after transformation (one per level)
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), paste0("group[dif: ", c("A", "B", "C", "D"), "]"))
})


test_that(".transform_factor_contrasts handles multi-factor transformed interactions", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3"))
  )

  for(interaction_contrast in c("orthonormal", "meandif")){
    formula_result <- JAGS_formula(
      formula = ~ a * b,
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
        b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
        "a:b"     = prior_factor("mnormal", list(0, 1), contrast = interaction_contrast)
      )
    )
    interaction_prior <- formula_result$prior_list$mu_a__xXx__b

    model_samples <- matrix(seq_len(20), nrow = 10, ncol = 2)
    colnames(model_samples) <- paste0("mu_a__xXx__b[", 1:2, "]")

    transformed <- suppressMessages(BayesTools:::.transform_factor_contrasts(
      model_samples,
      list(mu_a__xXx__b = interaction_prior),
      transform_factors = TRUE
    ))
    expected <- model_samples %*% t(attr(interaction_prior, "factor_design"))

    expect_equal(unname(transformed), unname(expected))
    expected_names <- paste0(
      "mu_a[dif: ",
      rep(c("a1", "a2"), times = 3),
      "]__xXx__b[dif: ",
      rep(c("b1", "b2", "b3"), each = 2),
      "]"
    )
    expect_equal(
      colnames(transformed),
      expected_names
    )
    expect_equal(anyDuplicated(colnames(transformed)), 0L)
  }
})

test_that(".transform_factor_contrasts handles one-coefficient interactions", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2"), levels = c("b1", "b2"))
  )
  formula_result <- JAGS_formula(
    formula = ~ a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b"     = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_a__xXx__b

  expect_equal(BayesTools:::.JAGS_prior_factor_names("mu_a__xXx__b", interaction_prior), "mu_a__xXx__b")

  model_samples <- matrix(seq_len(10), nrow = 10, ncol = 1)
  colnames(model_samples) <- "mu_a__xXx__b"

  transformed <- suppressMessages(BayesTools:::.transform_factor_contrasts(
    model_samples,
    list(mu_a__xXx__b = interaction_prior),
    transform_factors = TRUE
  ))
  expected <- model_samples %*% t(attr(interaction_prior, "factor_design"))

  expect_equal(unname(transformed), unname(expected))
  expect_equal(ncol(transformed), 4L)
  expect_equal(
    colnames(transformed),
    paste0(
      "mu_a[dif: ",
      rep(c("a1", "a2"), times = 2),
      "]__xXx__b[dif: ",
      rep(c("b1", "b2"), each = 2),
      "]"
    )
  )
})

test_that(".transform_factor_contrasts ignores unrelated transformations", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2"), levels = c("b1", "b2"))
  )
  formula_result <- JAGS_formula(
    formula = ~ a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b"     = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_a__xXx__b

  model_samples <- matrix(seq_len(10), nrow = 5, ncol = 2)
  colnames(model_samples) <- c("mu_intercept", "mu_a__xXx__b")

  transformed <- suppressMessages(BayesTools:::.transform_factor_contrasts(
    model_samples,
    list(
      mu_intercept = formula_result$prior_list$mu_intercept,
      mu_a__xXx__b = interaction_prior
    ),
    transform_factors = TRUE,
    transformations = list(mu_intercept = list(fun = exp))
  ))

  expect_equal(unname(transformed[, -1, drop = FALSE]), unname(model_samples[, 2, drop = FALSE] %*% t(attr(interaction_prior, "factor_design"))))
})

test_that(".transform_factor_contrasts reconstructs multi-factor designs from metadata", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3"))
  )
  formula_result <- JAGS_formula(
    formula = ~ a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b"     = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_a__xXx__b
  expected_design <- attr(interaction_prior, "factor_design")

  model_samples <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(model_samples) <- paste0("mu_a__xXx__b[", 1:2, "]")

  metadata_only_prior <- interaction_prior
  attr(metadata_only_prior, "factor_design") <- NULL
  attr(metadata_only_prior, "factor_cell_names") <- NULL
  transformed <- suppressMessages(BayesTools:::.transform_factor_contrasts(
    model_samples,
    list(mu_a__xXx__b = metadata_only_prior),
    transform_factors = TRUE
  ))
  expect_equal(unname(transformed), unname(model_samples %*% t(expected_design)))

  unnamed_contrast_prior <- metadata_only_prior
  attr(unnamed_contrast_prior, "factor_contrasts") <- unname(attr(unnamed_contrast_prior, "factor_contrasts"))
  transformed_unnamed <- suppressMessages(BayesTools:::.transform_factor_contrasts(
    model_samples,
    list(mu_a__xXx__b = unnamed_contrast_prior),
    transform_factors = TRUE
  ))
  expect_equal(unname(transformed_unnamed), unname(model_samples %*% t(expected_design)))

  inferred_contrast_prior <- metadata_only_prior
  attr(inferred_contrast_prior, "factor_contrasts") <- NULL
  model_samples_full <- cbind(
    mu_a = seq_len(nrow(model_samples)),
    `mu_b[1]` = seq_len(nrow(model_samples)) + 10,
    `mu_b[2]` = seq_len(nrow(model_samples)) + 20,
    model_samples
  )
  expect_error(
    suppressMessages(BayesTools:::.transform_factor_contrasts(
      model_samples_full,
      list(
        mu_a = formula_result$prior_list$mu_a,
        mu_b = formula_result$prior_list$mu_b,
        mu_a__xXx__b = inferred_contrast_prior
      ),
      transform_factors = TRUE
    )),
    "Factor contrast metadata is missing"
  )
})

test_that(".transform_factor_contrasts validates multi-factor metadata", {

  prior_obj <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  attr(prior_obj, "levels") <- 3
  attr(prior_obj, "level_names") <- list(a = c("a1", "a2"), b = c("b1", "b2"))
  attr(prior_obj, "interaction") <- TRUE
  attr(prior_obj, "factor_terms") <- c("a", "b")
  attr(prior_obj, "factor_contrasts") <- c(a = "contr.orthonormal")

  model_samples <- matrix(seq_len(10), nrow = 5, ncol = 2)
  colnames(model_samples) <- paste0("mu_a__xXx__b[", 1:2, "]")

  expect_error(
    BayesTools:::.factor_term_design_from_metadata(prior_obj),
    "incomplete"
  )

  missing_contrasts_prior <- prior_obj
  attr(missing_contrasts_prior, "factor_contrasts") <- NULL
  expect_error(
    BayesTools:::.factor_term_design_from_metadata(missing_contrasts_prior),
    "missing"
  )

  missing_levels_prior <- prior_obj
  attr(missing_levels_prior, "level_names") <- NULL
  attr(missing_levels_prior, "levels") <- NULL
  expect_error(
    BayesTools:::.factor_term_design_from_metadata(missing_levels_prior),
    "level names"
  )

  mismatch_prior <- prior_obj
  attr(mismatch_prior, "factor_contrasts") <- c(a = "contr.orthonormal", b = "contr.orthonormal")
  attr(mismatch_prior, "factor_design") <- diag(3)
  expect_error(
    BayesTools:::.transform_factor_contrasts(
      model_samples,
      list(mu_a__xXx__b = mismatch_prior),
      transform_factors = TRUE
    ),
    "has 3 coefficient columns"
  )
})

test_that("plain factor priors are canonicalized before mixed posterior transformation", {

  factor_prior <- prior_factor_levels(
    prior_factor("mnormal", list(0, 1), contrast = "meandif"),
    c("low", "mid", "high")
  )

  canonical_prior <- BayesTools:::.complete_factor_metadata(factor_prior, "p1")
  expect_equal(attr(canonical_prior, "factor_terms"), "p1")
  expect_equal(attr(canonical_prior, "factor_contrasts"), c(p1 = "contr.meandif"))
  expect_equal(
    attr(canonical_prior, "factor_design"),
    contr.meandif(c("low", "mid", "high"))
  )
  expect_equal(attr(canonical_prior, "factor_cell_names"), c("low", "mid", "high"))

  posterior <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(posterior) <- paste0("p1[", 1:2, "]")
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- list(p1 = factor_prior)
  fit <- attach_test_parameter_map(fit)

  mixed <- as_mixed_posteriors(fit, parameters = "p1")
  expect_equal(attr(mixed$p1, "factor_terms"), "p1")
  expect_equal(attr(mixed$p1, "factor_contrasts"), c(p1 = "contr.meandif"))

  transformed <- transform_factor_samples(mixed)$p1
  expect_equal(
    as.vector(transformed),
    as.vector(posterior %*% t(contr.meandif(c("low", "mid", "high"))))
  )
  expect_equal(colnames(transformed), paste0("p1[dif: ", c("low", "mid", "high"), "]"))
})

test_that("factor sample transformations still reject incomplete mixed posterior metadata", {

  incomplete_samples <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(incomplete_samples) <- paste0("p1[", 1:2, "]")
  attr(incomplete_samples, "levels") <- 2
  attr(incomplete_samples, "level_names") <- c("low", "mid", "high")
  attr(incomplete_samples, "meandif") <- TRUE
  class(incomplete_samples) <- c(
    "mixed_posteriors",
    "mixed_posteriors.factor",
    "mixed_posteriors.vector",
    "matrix"
  )

  expect_error(
    transform_factor_samples(list(p1 = incomplete_samples)),
    "Factor contrast metadata is missing"
  )
})

test_that("as_mixed_posteriors propagates multi-factor contrast metadata", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3"))
  )
  formula_result <- JAGS_formula(
    formula = ~ a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b"     = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_a__xXx__b

  posterior <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(posterior) <- paste0("mu_a__xXx__b[", 1:2, "]")
  fit <- coda::mcmc(complete_test_posterior(posterior, formula_result$prior_list))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  fit <- attach_test_parameter_map(fit)

  mixed <- as_mixed_posteriors(fit, parameters = "mu_a__xXx__b")

  expect_equal(attr(mixed$mu_a__xXx__b, "factor_terms"), attr(interaction_prior, "factor_terms"))
  expect_equal(attr(mixed$mu_a__xXx__b, "factor_contrasts"), attr(interaction_prior, "factor_contrasts"))
  expect_equal(attr(mixed$mu_a__xXx__b, "factor_design"), attr(interaction_prior, "factor_design"))
  expect_equal(attr(mixed$mu_a__xXx__b, "factor_cell_names"), attr(interaction_prior, "factor_cell_names"))
  mixed$mu_a__xXx__b <- .bt_meta_set(mixed$mu_a__xXx__b, "support", stats::setNames(
    rep(
      list(BayesTools:::.posterior_support_new(c(-1, 1), source = "test")),
      ncol(mixed$mu_a__xXx__b)
    ),
    colnames(mixed$mu_a__xXx__b)
  ))

  transformed <- transform_factor_samples(mixed)$mu_a__xXx__b
  expect_equal(
    as.vector(transformed),
    as.vector(posterior %*% t(attr(interaction_prior, "factor_design")))
  )
  expect_null(.bt_meta_get(transformed, "support"))
})


test_that(".filter_parameters removes spike at 0 priors", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    delta = prior("point", list(0)),  # spike at 0
    tau = prior("normal", list(1, 1))
  )

  # With remove_spike_0 = TRUE
  result <- BayesTools:::.filter_parameters(prior_list, remove_spike_0 = TRUE)
  expect_true("delta" %in% result)
  expect_false("mu" %in% result)
  expect_false("tau" %in% result)

  # With remove_spike_0 = FALSE
  result <- BayesTools:::.filter_parameters(prior_list, remove_spike_0 = FALSE)
  expect_equal(length(result), 0)
})

test_that(".filter_parameters keeps point priors with expression locations", {

  prior_list <- list(
    a = prior("normal", list(0, 1)),
    b = prior("point", list(location = expression(a))),
    c = prior("point", list(0))
  )

  # the derived point b is not a structural spike at zero (no coercion error)
  result <- BayesTools:::.filter_parameters(prior_list, remove_spike_0 = TRUE)
  expect_identical(result, "c")

  # mixed posteriors treat the derived point as having unknown support
  set.seed(10)
  a <- rnorm(50)
  fit <- coda::mcmc(cbind(a = a, b = a))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- prior_list[c("a", "b")]
  fit <- attach_test_parameter_map(fit)
  samples <- as_mixed_posteriors(fit, c("a", "b"))
  expect_null(.bt_meta_get(samples$b, "support"))
  expect_equal(as.numeric(marginal_posterior(samples, "b")), a)
})


test_that(".filter_parameters removes character specified parameters", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    delta = prior("normal", list(0, 1)),
    tau = prior("normal", list(1, 1))
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = c("mu", "tau"), remove_spike_0 = FALSE)
  expect_true("mu" %in% result)
  expect_true("tau" %in% result)
  expect_false("delta" %in% result)
})


test_that(".filter_parameters removes non-formula parameters when TRUE", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create priors with formula attributes
  prior_formula <- prior("normal", list(0, 1))
  attr(prior_formula, "parameter") <- "y"

  prior_list <- list(
    intercept = prior_formula,
    sigma = prior("normal", list(1, 1))  # no formula attribute
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = TRUE, remove_spike_0 = FALSE)
  expect_true("sigma" %in% result)
  expect_false("intercept" %in% result)
})


test_that(".filter_parameters removes formula-specific parameters", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create priors with different formula attributes
  prior_y <- prior("normal", list(0, 1))
  attr(prior_y, "parameter") <- "y"

  prior_x <- prior("normal", list(0, 1))
  attr(prior_x, "parameter") <- "x"

  prior_list <- list(
    intercept_y = prior_y,
    slope_y = prior_y,
    intercept_x = prior_x,
    sigma = prior("normal", list(1, 1))  # no formula attribute
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_formulas = "y", remove_spike_0 = FALSE)
  expect_true("intercept_y" %in% result)
  expect_true("slope_y" %in% result)
  expect_false("intercept_x" %in% result)
  expect_false("sigma" %in% result)
})


test_that(".filter_parameters keeps only specified parameters", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    delta = prior("normal", list(0, 1)),
    tau = prior("normal", list(1, 1))
  )

  result <- BayesTools:::.filter_parameters(prior_list, keep_parameters = "mu", remove_spike_0 = FALSE)
  expect_false("mu" %in% result)
  expect_true("delta" %in% result)
  expect_true("tau" %in% result)
})


test_that(".filter_parameters keeps only specified formulas", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create priors with different formula attributes
  prior_y <- prior("normal", list(0, 1))
  attr(prior_y, "parameter") <- "y"

  prior_x <- prior("normal", list(0, 1))
  attr(prior_x, "parameter") <- "x"

  prior_list <- list(
    intercept_y = prior_y,
    slope_y = prior_y,
    intercept_x = prior_x,
    sigma = prior("normal", list(1, 1))  # no formula attribute
  )

  result <- BayesTools:::.filter_parameters(prior_list, keep_formulas = "y", remove_spike_0 = FALSE)
  expect_false("intercept_y" %in% result)
  expect_false("slope_y" %in% result)
  expect_true("intercept_x" %in% result)
  expect_true("sigma" %in% result)
})


test_that(".filter_parameters combines keep_parameters and keep_formulas", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create priors with different formula attributes
  prior_y <- prior("normal", list(0, 1))
  attr(prior_y, "parameter") <- "y"

  prior_x <- prior("normal", list(0, 1))
  attr(prior_x, "parameter") <- "x"

  prior_list <- list(
    intercept_y = prior_y,
    slope_y = prior_y,
    intercept_x = prior_x,
    sigma = prior("normal", list(1, 1))  # no formula attribute
  )

  # Keep formula "y" and parameter "sigma"
  result <- BayesTools:::.filter_parameters(prior_list, keep_parameters = "sigma", keep_formulas = "y", remove_spike_0 = FALSE)
  expect_false("intercept_y" %in% result)
  expect_false("slope_y" %in% result)
  expect_false("sigma" %in% result)
  expect_true("intercept_x" %in% result)
})


test_that(".filter_parameters removes bias-related parameters when bias is removed", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with PET component (simulating bias)
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, 1), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = bias_prior
  )

  # When bias is removed, PET should also be removed
  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = "bias", remove_spike_0 = FALSE)
  expect_true("bias" %in% result)
  expect_true("PET" %in% result)
  expect_false("mu" %in% result)
})


test_that(".filter_parameters removes bias-related parameters when bias contains PEESE", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with PEESE component
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PEESE("normal", list(0, 1), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = bias_prior
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = "bias", remove_spike_0 = FALSE)
  expect_true("bias" %in% result)
  expect_true("PEESE" %in% result)
  expect_false("mu" %in% result)
})


test_that(".filter_parameters removes bias-related parameters when bias contains weightfunction", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with weightfunction component
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = bias_prior
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = "bias", remove_spike_0 = FALSE)
  expect_true("bias" %in% result)
  expect_true("omega" %in% result)
  expect_false("mu" %in% result)
})


test_that(".filter_parameters keeps bias-related parameters when bias is kept", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with PET component
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, 1), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(1, 1)),
    bias = bias_prior
  )

  # When only bias is kept, mu and tau should be removed, but PET should be kept
  result <- BayesTools:::.filter_parameters(prior_list, keep_parameters = "bias", remove_spike_0 = FALSE)
  expect_false("bias" %in% result)
  expect_false("PET" %in% result)
  expect_true("mu" %in% result)
  expect_true("tau" %in% result)
})


test_that(".filter_parameters handles non-mixture bias priors", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a single PET prior named bias (not a mixture)
  bias_prior <- prior_PET("normal", list(0, 1))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = bias_prior
  )

  result <- BayesTools:::.filter_parameters(prior_list, remove_parameters = "bias", remove_spike_0 = FALSE)
  expect_true("bias" %in% result)
  expect_true("PET" %in% result)
  expect_false("mu" %in% result)
})


test_that(".filter_parameters removes bias-related parameters when bias is not in keep list", {
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Create a mixture prior with PET component
  bias_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, 1), prior_weights = 1)
  ))
  
  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(1, 1)),
    bias = bias_prior
  )

  # When only mu is kept, bias should be removed along with PET
  result <- BayesTools:::.filter_parameters(prior_list, keep_parameters = "mu", remove_spike_0 = FALSE)
  expect_false("mu" %in% result)
  expect_true("bias" %in% result)
  expect_true("PET" %in% result)
  expect_true("tau" %in% result)
})


