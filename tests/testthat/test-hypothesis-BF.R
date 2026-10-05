skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Hypothesis Bayes Factors
# ============================================================================ #

.hypothesis_marginal_posterior_for_test <- function(samples, prior_density){

  class(samples) <- c("marginal_posterior.simple", "marginal_posterior", class(samples))
  samples <- .bt_meta_set(samples, "prior_density", prior_density)
  samples <- .bt_meta_set(samples, "atoms", posterior_atom_attribute())

  samples
}


.posterior_density_for_test <- function(x, y, method = "iwmde",
                                        density_method = "precomputed", ...){
  posterior_density_attribute(x = x, y = y, method = method,
                              density_method = density_method, ...)
}


.posterior_ordinate_for_test <- function(value, ordinate, method = "qCMDE",
                                         density_method = "precomputed", ...){
  posterior_ordinate_attribute(value = value, ordinate = ordinate,
                               method = method,
                               density_method = density_method, ...)
}


.hypothesis_prior_density_for_test <- function(){

  BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("normal", list(mean = 0, sd = 1))),
    weights    = c(theta = 1),
    n_grid     = 1024
  )
}


.hypothesis_log_odds_var_for_test <- function(left, right){

  n       <- length(left)
  left    <- as.numeric(left)
  right   <- as.numeric(right)
  p_left  <- mean(left)
  p_right <- mean(right)

  stats::var(left) / (n * p_left^2) +
    stats::var(right) / (n * p_right^2) -
    2 * stats::cov(left, right) / (n * p_left * p_right)
}


.hypothesis_log_prob_var_for_test <- function(values){

  n      <- length(values)
  values <- as.numeric(values)
  p      <- mean(values)

  stats::var(values) / (n * p^2)
}


test_that("hypothesis_BF requires a finite seed", {

  arguments <- list(
    posterior = c(-2, -1, 1, 2),
    prior = c(-2, -1, 1, 2),
    hypothesis = "theta > 0",
    parameter = "theta"
  )

  for(seed in c(Inf, -Inf)){
    expect_error(
      do.call(hypothesis_BF, c(arguments, list(seed = seed))),
      "'seed' must be finite.",
      fixed = TRUE
    )
  }
  expect_s3_class(
    do.call(hypothesis_BF, c(arguments, list(seed = 1.5))),
    "BayesTools_hypothesis_BF"
  )
})

test_that("transformed factor levels preserve joint affine prior provenance", {

  formula <- JAGS_formula(
    ~ fac, "mu", data.frame(fac = factor(c("A", "B"))),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      fac = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  fit <- coda::mcmc(cbind(mu_intercept = seq(-1, 1, length.out = 201),
                          mu_fac = seq(-2, 2, length.out = 201)))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula$prior_list
  fit <- attach_test_parameter_map(fit)
  mixed <- as_mixed_posteriors(fit, "mu_fac")
  raw <- marginal_posterior(mixed, "mu_fac", use_formula = FALSE,
                             prior_samples = TRUE)
  affine <- marginal_posterior(
    mixed, "mu_fac", use_formula = FALSE, prior_samples = TRUE,
    transformation = "lin", transformation_arguments = list(a = 3, b = 2)
  )
  target <- hypothesis_linear_target(affine, "mu_fac[B] - mu_fac[A] = 0", "mu_fac")
  ordinate <- prior_density_ordinate(.bt_meta_get(target$posterior, "prior_density"), 0)
  expect_equal(ordinate$log_density, stats::dnorm(0, sd = 2, log = TRUE), tolerance = 1e-12)
  expect_identical(.bt_meta_get(affine$B, "linear_offset"), 3)
  expect_equal(.bt_meta_get(affine$B, "linear_weights")[["mu_fac"]], 2)

  raw_region <- hypothesis_BF(raw, hypothesis = "mu_fac[B] > 0", seed = 82, columns = "all")
  affine_region <- hypothesis_BF(affine, hypothesis = "mu_fac[B] > 3", seed = 82, columns = "all")
  expect_identical(affine_region$prior, raw_region$prior)
  expect_identical(attr(affine_region, "raw_BF"), attr(raw_region, "raw_BF"))

  nonlinear <- marginal_posterior(
    mixed, "mu_fac", use_formula = FALSE, prior_samples = TRUE,
    transformation = "exp"
  )
  message <- "Joint prior information is unavailable for nonlinear transformed level hypotheses. Use untransformed levels or a direct scalar hypothesis."
  expect_error(hypothesis_linear_target(nonlinear, "mu_fac[B] - mu_fac[A] = 0", "mu_fac"),
               message, fixed = TRUE)
  expect_error(hypothesis_BF(nonlinear, hypothesis = "mu_fac[B] > mu_fac[A]", seed = 82),
               message, fixed = TRUE)
  expect_s3_class(hypothesis_BF(nonlinear, hypothesis = "mu_fac[B] = 1", density_method = "normal"),
                  "BayesTools_hypothesis_BF")
})

test_that("explicit hypothesis seeds preserve the caller's RNG state", {

  withr::local_seed(714)
  original_seed <- .Random.seed
  arguments <- list(
    posterior = seq(-2, 2, length.out = 201),
    prior = prior("normal", list(0, 1)),
    hypothesis = "theta > 0",
    parameter = "theta",
    seed = 19
  )
  first <- do.call(hypothesis_BF, arguments)
  expect_identical(.Random.seed, original_seed)
  expect_identical(do.call(hypothesis_BF, arguments), first)
  expect_identical(.Random.seed, original_seed)
  expect_error(
    hypothesis_BF(environment(), hypothesis = "theta > 0", seed = 1),
    "Unsupported posterior input"
  )
  expect_identical(.Random.seed, original_seed)

  rm(".Random.seed", envir = .GlobalEnv)
  do.call(hypothesis_BF, arguments)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("prior density construction failures remain visible in hypothesis tests", {

  testthat::local_mocked_bindings(
    .prior_linear_combination_density = function(...) {
      stop("Prior density construction failed.", call. = FALSE)
    }
  )
  expect_error(
    hypothesis_BF(seq(-2, 2, length.out = 201),
                  prior = prior("normal", list(0, 1)),
                  hypothesis = "theta = 0", parameter = "theta"),
    "Prior density construction failed.",
    fixed = TRUE
  )
})


test_that("hypothesis_BF computes point-null Savage-Dickey from numeric draws", {

  set.seed(1)
  prior     <- stats::rnorm(12000, mean = 0, sd = 1)
  posterior <- stats::rnorm(12000, mean = 0.4, sd = 1.2)

  expect_warning(out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0",
    parameter  = "theta"
  ), class = "BayesTools_inexact_ordinate")

  expected <- stats::dnorm(0, mean = 0, sd = 1) /
    stats::dnorm(0, mean = 0.4, sd = 1.2)

  expect_equal(attr(out, "raw_BF"), expected, tolerance = 0.08)
  expect_s3_class(attr(out, "hypothesis_ast"), "BayesTools_hypothesis_ast")
  expect_null(attr(out, "parsed", exact = TRUE))

  expect_warning(equivalent <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta != 0",
    parameter  = "theta"
  ), class = "BayesTools_inexact_ordinate")

  expect_equal(attr(equivalent, "raw_BF"), attr(out, "raw_BF"),
               tolerance = 1e-12)

  expect_warning(explicit_inverse <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0 vs theta != 0",
    parameter  = "theta"
  ), class = "BayesTools_inexact_ordinate")

  expect_equal(attr(explicit_inverse, "raw_BF"), 1 / attr(out, "raw_BF"),
               tolerance = 1e-12)
})


test_that("hypothesis_BF warns for exact point masses in raw point-null draws", {

  continuous_draws <- seq(-3, 3, length.out = 400)
  posterior_spike <- c(rep(0, 20), seq(-3, 3, length.out = 380))
  prior_spike <- c(rep(0, 20), seq(-3, 3, length.out = 380))

  expect_warning(expect_warning(
    hypothesis_BF(
      posterior  = posterior_spike,
      prior      = continuous_draws,
      hypothesis = "theta = 0",
      parameter  = "theta"
    ),
    "posterior draws exactly match"
  ), class = "BayesTools_inexact_ordinate")

  expect_warning(expect_warning(
    hypothesis_BF(
      posterior  = continuous_draws,
      prior      = prior_spike,
      hypothesis = "theta = 0",
      parameter  = "theta"
    ),
    "prior draws exactly match"
  ), class = "BayesTools_inexact_ordinate")
})


test_that("hypothesis_BF KDEs retain Gaussian tails beyond sample grids", {

  posterior <- seq(-1, 1, length.out = 101)
  prior     <- seq(-1.2, 1.2, length.out = 101)
  null      <- 3

  posterior_bw <- stats::density(posterior)[["bw"]]
  prior_bw     <- stats::density(prior)[["bw"]]
  expected_posterior <- mean(stats::dnorm(
    null,
    mean = posterior,
    sd   = posterior_bw
  ))
  expected_prior <- mean(stats::dnorm(
    null,
    mean = prior,
    sd   = prior_bw
  ))

  expect_warning(
    posterior_height <- BayesTools:::.hypothesis_sample_density_height(
      posterior, null, "posterior"
    ),
    "posterior samples do not span"
  )
  expect_warning(
    prior_height <- BayesTools:::.hypothesis_sample_density_height(
      prior, null, "prior"
    ),
    "prior samples do not span"
  )
  expect_equal(posterior_height, expected_posterior, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(prior_height, expected_prior, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_true(is.list(attr(posterior_height, "kde_extrapolation")))
  expect_true(is.list(attr(prior_height, "kde_extrapolation")))
  expect_gt(posterior_height, 0)
  expect_gt(prior_height, 0)

  out <- suppressWarnings(hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 3",
    parameter  = "theta",
    columns    = "all"
  ))
  expect_equal(out[["posterior"]], expected_posterior, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(out[["prior"]], expected_prior, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(
    attr(out, "raw_BF"),
    expected_prior / expected_posterior,
    tolerance = 1e-12
  )
  expect_true(is.finite(attr(out, "raw_BF")))
})


test_that("hypothesis_BF marginal KDEs retain tails beyond evaluation grids", {

  prior_density <- .hypothesis_prior_density_for_test()
  posterior_samples <- seq(-1, 1, length.out = 101)
  posterior <- .hypothesis_marginal_posterior_for_test(
    posterior_samples,
    prior_density
  )
  null <- 2

  posterior_bw <- stats::density(posterior_samples)[["bw"]]
  expected_height <- mean(stats::dnorm(
    null,
    mean = posterior_samples,
    sd   = posterior_bw
  ))
  expect_gt(null, max(stats::density(posterior_samples)[["x"]]))

  expect_warning(
    height <- BayesTools:::.Savage_Dickey_BF.kd(posterior, null),
    "Gaussian kernel tails"
  )
  expect_equal(height, expected_height, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_gt(height, 0)

  supported_samples <- seq(.01, 1, length.out = 101)
  supported_bw <- stats::density(supported_samples)[["bw"]]
  expected_supported_height <-
    mean(stats::dnorm(null, mean = supported_samples, sd = supported_bw)) +
    mean(stats::dnorm(null, mean = -supported_samples, sd = supported_bw))
  expect_warning(
    supported_height <- BayesTools:::.Savage_Dickey_BF.kd(
      supported_samples,
      null,
      support = c(0, Inf)
    ),
    "Gaussian kernel tails"
  )
  expect_equal(
    as.numeric(supported_height),
    expected_supported_height,
    tolerance = 1e-12
  )
  expect_gt(supported_height, 0)
  expect_true(attr(supported_height, "boundary_reflection"))

  out <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = "theta = 2",
    parameter  = "theta",
    columns    = "all"
  )
  expect_equal(out[["posterior"]], expected_height, tolerance = 1e-12)
  expect_true(is.finite(attr(out, "raw_BF")))
  expect_match(attr(out, "warnings"), "Posterior samples do not span")
})


test_that("hypothesis_BF rejects prior point mass metadata for point nulls", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      theta = prior_spike_and_slab(
        prior("normal", list(mean = 0, sd = 1)),
        prior_inclusion = prior("point", list(location = 0.5))
      )
    ),
    weights = c(theta = 1),
    n_grid = 1024
  )

  expect_error(
    BayesTools:::.hypothesis_prior_density_height(prior_density, 0),
    "point mass in the prior"
  )
})


test_that("hypothesis_BF routes normal point-null density requests", {

  set.seed(11)
  prior_density <- .hypothesis_prior_density_for_test()
  posterior <- .hypothesis_marginal_posterior_for_test(
    stats::rnorm(8000, mean = 0.35, sd = 1.15),
    prior_density
  )

  out <- hypothesis_BF(
    posterior      = posterior,
    hypothesis     = "theta = 0",
    parameter      = "theta",
    density_method = "normal",
    columns        = "all"
  )
  expected <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = 0,
    normal_approximation = TRUE,
    silent               = TRUE
  )

  expect_equal(attr(out, "raw_BF"), as.numeric(expected), tolerance = 1e-12)
  expect_equal(out[["method"]], "Savage-Dickey (normal)")
})


test_that("hypothesis_BF uses boundary-reflected KDE for marginal point nulls", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("beta", list(alpha = 1, beta = 1))),
    weights    = c(theta = 1),
    n_grid     = 4096
  )
  posterior <- .hypothesis_marginal_posterior_for_test(
    seq(.001, .999, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "support", posterior_support_attribute(c(0, 1)))

  expected <- Savage_Dickey_BF(
    posterior,
    null_hypothesis      = 0,
    normal_approximation = FALSE,
    silent               = TRUE
  )
  out <- hypothesis_BF(
    posterior      = posterior,
    hypothesis     = "theta = 0",
    parameter      = "theta",
    density_method = "KDE",
    columns        = "all"
  )

  expect_true(attr(expected, "posterior_density_boundary_reflection"))
  expect_equal(attr(out, "raw_BF"), as.numeric(expected), tolerance = 1e-12)
  expect_equal(out[["method"]], "Savage-Dickey")
})


test_that("hypothesis_BF applies normal point-null routing to list levels", {

  set.seed(12)
  prior_density <- .hypothesis_prior_density_for_test()
  level_a <- .hypothesis_marginal_posterior_for_test(
    stats::rnorm(5000, mean = 0.2, sd = 1.1),
    prior_density
  )
  level_b <- .hypothesis_marginal_posterior_for_test(
    stats::rnorm(5000, mean = -0.2, sd = 0.9),
    prior_density
  )
  posterior <- list(a = level_a, b = level_b)
  class(posterior) <- c("marginal_posterior.factor", "marginal_posterior", "list")

  explicit <- hypothesis_BF(
    posterior      = posterior,
    hypothesis     = "theta[ a ] = 0",
    parameter      = "theta",
    density_method = "normal",
    columns        = "all"
  )
  expected <- Savage_Dickey_BF(
    level_a,
    null_hypothesis      = 0,
    normal_approximation = TRUE,
    silent               = TRUE
  )

  expect_equal(attr(explicit, "raw_BF"), as.numeric(expected), tolerance = 1e-12)
  expect_equal(explicit[["method"]], "Savage-Dickey (normal)")

  expanded <- hypothesis_BF(
    posterior      = posterior,
    hypothesis     = "theta = 0",
    parameter      = "theta",
    density_method = "normal",
    columns        = "all"
  )

  expect_equal(nrow(expanded), 2L)
  expect_equal(expanded[["method"]], rep("Savage-Dickey (normal)", 2))
})


test_that("hypothesis_BF normal compound point expressions do not use KDE", {

  set.seed(13)
  prior     <- stats::rnorm(7000, mean = 0, sd = 1)
  posterior <- stats::rnorm(7000, mean = 0.3, sd = 1.2)

  expect_warning(out <- hypothesis_BF(
    posterior      = posterior,
    prior          = prior,
    hypothesis     = "theta + 0 = 0",
    parameter      = "theta",
    density_method = "normal",
    columns        = "all"
  ), class = "BayesTools_inexact_ordinate")

  expected_prior <- stats::dnorm(0, mean = mean(prior), sd = stats::sd(prior))
  expected_posterior <- stats::dnorm(
    0,
    mean = mean(posterior),
    sd   = stats::sd(posterior)
  )

  expect_equal(attr(out, "raw_BF"), expected_prior / expected_posterior)
  expect_equal(out[["prior"]], expected_prior)
  expect_equal(out[["posterior"]], expected_posterior)
  expect_equal(out[["method"]], "Savage-Dickey (normal)")
  expect_false(identical(out[["method"]], "kernel Savage-Dickey"))
})


test_that("hypothesis parser helpers expose point and level references", {

  parsed <- hypothesis_parse_point_reference(c(
    "theta[ a ] = 0",
    "theta + 0 != 1",
    "theta == 0 VS theta > -1"
  ))

  expect_equal(parsed[["symbol"]][1], "theta[a]")
  expect_equal(parsed[["parameter"]][1], "theta")
  expect_equal(parsed[["level"]][1], "a")
  expect_true(parsed[["direct"]][1])
  expect_false(parsed[["direct"]][2])
  expect_equal(parsed[["operator"]][3], "=")

  expect_error(
    hypothesis_parse_point_reference("theta + 0 = 0", allow_compound = FALSE),
    "not a direct"
  )
  expect_false(
    hypothesis_parse_point_reference("theta = other")[["direct"]]
  )

  refs <- hypothesis_parse_level_reference(c("theta[ a ]", "theta"))
  expect_equal(refs[["symbol"]][1], "theta[a]")
  expect_equal(refs[["parameter"]][1], "theta")
  expect_equal(refs[["level"]][1], "a")
  expect_true(refs[["direct"]][1])
  expect_false(refs[["direct"]][2])
  expect_equal(hypothesis_normalize_level_references("theta[ a ] = 0"),
               "`theta[a]` = 0")

  parsed_vs_level <- hypothesis_parse_point_reference("theta[a vs b] = 0")
  expect_equal(parsed_vs_level[["symbol"]], "theta[a vs b]")
  expect_equal(parsed_vs_level[["level"]], "a vs b")
})


test_that("hypothesis_BF returns compact BayesTools table by default", {

  set.seed(101)
  prior     <- stats::rnorm(4000)
  posterior <- stats::rnorm(4000, mean = 0.2)

  expect_warning(out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0",
    parameter  = "theta"
  ), class = "BayesTools_inexact_ordinate")

  expect_s3_class(out, "BayesTools_table")
  expect_s3_class(out, "BayesTools_hypothesis_BF")
  expect_equal(colnames(out), c("Alternative", "Null", "BF", "BF_error"))
  expect_equal(attr(out, "type"),
               c("hypothesis_label", "hypothesis_label", "BF", "BF_error"))
  expect_null(attr(out, "footnotes"))
  expect_equal(attr(out[["BF_error"]], "name"), "error%(BF)")
  expect_equal(out[["Alternative"]], "theta != 0")
  expect_equal(out[["Null"]], "theta = 0")

  printed <- utils::capture.output(print(out))
  expect_true(any(grepl("Alternative:", printed, fixed = TRUE)))
  expect_true(any(grepl("Null:", printed, fixed = TRUE)))
  expect_true(any(grepl("error%(BF)", printed, fixed = TRUE)))

  log_out <- update(out, logBF = TRUE)
  expect_equal(as.numeric(log_out[["BF"]]), log(attr(out, "raw_BF")),
               tolerance = 1e-12)
  expect_equal(attr(log_out[["BF_error"]], "name"), "error%(BF)")
  expect_true(attr(log_out, "logBF"))
  expect_false(attr(log_out, "BF01"))

  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      prior      = prior,
      hypothesis = "theta = 0",
      parameter  = "theta",
      columns    = c("all", "typo")
    ),
    "cannot be combined"
  )
  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      prior      = prior,
      hypothesis = "theta = 0",
      parameter  = "theta",
      columns    = "prior_probability"
    ),
    "Unknown 'columns' value"
  )

  # each column has one spelling: its column name
  expect_warning(selected <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0",
    parameter  = "theta",
    columns    = c("method", "prior")
  ), class = "BayesTools_inexact_ordinate")
  expect_equal(
    colnames(selected),
    c("Alternative", "Null", "BF", "BF_error", "method", "prior")
  )
  for(spelling in c("Method", "computation_method", "Prior", "error%(BF)",
                    "alternative", "null")){
    expect_error(
      hypothesis_BF(
        posterior  = posterior,
        prior      = prior,
        hypothesis = "theta = 0",
        parameter  = "theta",
        columns    = spelling
      ),
      paste0(
        "Unknown 'columns' value: '", spelling, "'. Columns are selected by ",
        "their names: 'Alternative', 'Null', 'BF', 'BF_error', 'prior', ",
        "'posterior', 'method', or 'default' or 'all'."
      ),
      fixed = TRUE
    )
  }

  expect_warning(detailed <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0",
    parameter  = "theta",
    columns    = "all"
  ), class = "BayesTools_inexact_ordinate")

  expect_equal(
    colnames(detailed),
    c("Alternative", "Null", "BF", "BF_error", "prior", "posterior", "method")
  )
  expect_match(attr(detailed, "footnotes"), "diagnostic values", fixed = TRUE)
  expect_match(attr(detailed, "footnotes"), "not always probabilities",
               fixed = TRUE)

  printed_detailed <- utils::capture.output(print(detailed))
  expect_true(any(grepl("diagnostic values", printed_detailed, fixed = TRUE)))

  expect_warning(multi <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = c("theta = 0", "theta > 0 vs theta < 0"),
    parameter  = "theta"
  ), class = "BayesTools_inexact_ordinate")
  expect_equal(attr(multi[2, , drop = FALSE], "raw_BF"),
               attr(multi, "raw_BF")[2])

  expect_null(attr(multi[c(1, 1), , drop = FALSE], "raw_BF"))

  mismatched <- multi[1, , drop = FALSE]
  rownames(mismatched) <- "missing"
  attr(mismatched, "raw_BF") <- attr(multi, "raw_BF")
  mismatched <- BayesTools:::.subset_table_hypothesis_attributes(multi, mismatched)
  expect_null(attr(mismatched, "raw_BF"))
})


test_that("hypothesis_BF row names identify repeated statements on one quantity", {

  set.seed(1)
  posterior <- stats::rnorm(4000, mean = 0.3, sd = 0.1)
  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior("normal", list(mean = 0, sd = 1)),
    hypothesis = c("theta > 0", "theta = 0.3", "theta < -1"),
    parameter  = "theta",
    seed       = 1
  )
  expect_identical(rownames(out), c("theta (1)", "theta (2)", "theta (3)"))
  expect_identical(
    out[["Alternative"]],
    c("theta > 0", "theta != 0.3", "theta < -1")
  )
  # Warnings and row subsets are keyed by the same row names.
  expect_identical(names(attr(out, "warnings")), "theta (3)")
  expect_equal(attr(out["theta (2)", , drop = FALSE], "raw_BF"),
               attr(out, "raw_BF")[2])
  expect_null(attr(out[c(1, 1), , drop = FALSE], "raw_BF"))

  single <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior("normal", list(mean = 0, sd = 1)),
    hypothesis = "theta > 0",
    parameter  = "theta"
  )
  expect_identical(rownames(single), "theta")

  draws <- data.frame(mu = posterior, phi = stats::rnorm(4000, 0.1, 0.1))
  prior_draws <- data.frame(mu = stats::rnorm(4000), phi = stats::rnorm(4000))
  table <- hypothesis_BF(
    posterior  = draws,
    prior      = prior_draws,
    hypothesis = c("mu > phi", "mu - phi > 0.5")
  )
  expect_identical(rownames(table), c("draws (1)", "draws (2)"))
})


test_that("hypothesis_BF accepts scalar BayesTools prior objects", {

  set.seed(2)
  posterior <- stats::rnorm(8000, mean = 0.3, sd = 1.1)

  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior("normal", list(mean = 0, sd = 1)),
    hypothesis = "theta = 0",
    parameter  = "theta",
    columns    = "all",
    seed       = 3
  )

  expect_equal(out[["prior"]], stats::dnorm(0), tolerance = 1e-3)
  expect_equal(
    attr(out, "raw_BF"),
    out[["prior"]] / out[["posterior"]],
    tolerance = 1e-12
  )
})


test_that("hypothesis_BF accepts mixture and spike-and-slab prior objects", {

  set.seed(7)
  posterior <- stats::rnorm(20000, mean = 0.3, sd = 0.2)
  mixture <- prior_mixture(list(
    prior("normal", list(mean = 0, sd = 1), prior_weights = 1),
    prior("normal", list(mean = 0, sd = 0.3), prior_weights = 3)
  ), is_null = c(FALSE, FALSE))
  spike_and_slab <- prior_spike_and_slab(
    prior("normal", list(mean = 0, sd = 1)),
    prior_inclusion = prior("spike", list(0.3))
  )

  # Reference: analytic mixture densities; 1e-4 is the grid-refinement
  # criterion of the deterministic prior density.
  point <- hypothesis_BF(posterior, mixture, hypothesis = "theta = 0",
                         parameter = "theta", columns = "all", seed = 1)
  expect_equal(as.numeric(point[["prior"]]),
               .25 * stats::dnorm(0) + .75 * stats::dnorm(0, sd = 0.3),
               tolerance = 1e-4)
  slab <- hypothesis_BF(posterior, spike_and_slab, hypothesis = "theta = 0.3",
                        parameter = "theta", columns = "all", seed = 1)
  expect_equal(as.numeric(slab[["prior"]]), 0.3 * stats::dnorm(0.3),
               tolerance = 1e-4)
  expect_error(
    hypothesis_BF(posterior, spike_and_slab, hypothesis = "theta = 0",
                  parameter = "theta", seed = 1),
    "point mass in the prior"
  )

  # Region masses of non-simple prior objects come from 20000 prior draws;
  # compare the prior odds with the analytic value within 5 Monte Carlo SEs.
  region_cases <- list(
    list(prior = mixture,
         mass  = .25 * stats::pnorm(0.2, lower.tail = FALSE) +
           .75 * stats::pnorm(0.2, sd = 0.3, lower.tail = FALSE)),
    list(prior = spike_and_slab,
         mass  = 0.3 * stats::pnorm(0.2, lower.tail = FALSE))
  )
  for(case in region_cases){
    region <- hypothesis_BF(posterior, case$prior, hypothesis = "theta > 0.2",
                            parameter = "theta", columns = "all", seed = 1)
    odds    <- case$mass / (1 - case$mass)
    odds_se <- sqrt(case$mass * (1 - case$mass) / 20000) / (1 - case$mass)^2
    expect_lt(abs(region[["prior"]] - odds), 5 * odds_se)
    expect_true(is.finite(attr(region, "raw_BF")))
  }
})


test_that("hypothesis_BF uses exact scalar prior density beyond its grid", {

  theta_prior <- prior("normal", list(mean = 0, sd = 1))
  prior_grid <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = theta_prior),
    weights    = c(theta = 1)
  )
  expect_gt(4, max(prior_grid[["density"]][["x"]]))

  posterior <- seq(3.5, 4.5, length.out = 401)
  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = theta_prior,
    hypothesis = "theta = 4",
    parameter  = "theta",
    seed       = 14,
    columns    = "all"
  )

  expect_equal(out[["prior"]], stats::dnorm(4), tolerance = 1e-15)
  expect_gt(out[["posterior"]], 0)
  expect_true(is.finite(attr(out, "raw_BF")))
})


test_that("hypothesis_BF uses analytic scalar prior probabilities for simple regions", {

  posterior <- c(
    seq(-1, 0.2, length.out = 30),
    seq(0.3, 2, length.out = 70)
  )

  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior("normal", list(mean = 0, sd = 1)),
    hypothesis = "theta > 0.25",
    parameter  = "theta",
    columns    = "all",
    seed       = 1
  )
  out_seed2 <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior("normal", list(mean = 0, sd = 1)),
    hypothesis = "theta > 0.25",
    parameter  = "theta",
    columns    = "all",
    seed       = 2
  )

  prior_left      <- stats::pnorm(0.25, lower.tail = FALSE)
  prior_right     <- stats::pnorm(0.25)
  posterior_left  <- mean(posterior > 0.25)
  posterior_right <- mean(posterior <= 0.25)
  expected        <- (posterior_left / posterior_right) /
    (prior_left / prior_right)

  expect_equal(out[["prior"]], prior_left / prior_right, tolerance = 1e-12)
  expect_equal(attr(out, "raw_BF"), expected, tolerance = 1e-12)
  expect_equal(attr(out_seed2, "raw_BF"), expected, tolerance = 1e-12)
})


test_that("hypothesis_BF uses prior-posterior odds for region hypotheses", {

  prior <- data.frame(theta = c(
    seq(-0.75, -0.25, length.out = 25),
    seq(0.05, 0.95, length.out = 75)
  ))
  posterior <- data.frame(theta = c(
    seq(-0.75, -0.25, length.out = 20),
    seq(0.05, 0.95, length.out = 80)
  ))

  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta > 0.25",
    columns    = "all"
  )

  prior_left      <- mean(prior[["theta"]] > 0.25)
  prior_right     <- mean(prior[["theta"]] <= 0.25)
  posterior_left  <- mean(posterior[["theta"]] > 0.25)
  posterior_right <- mean(posterior[["theta"]] <= 0.25)
  expected        <- (posterior_left / posterior_right) /
    (prior_left / prior_right)

  expect_equal(attr(out, "raw_BF"), expected, tolerance = 1e-12)
  expect_equal(out[["method"]], "prior-posterior odds")

  expected_error <- 100 * sqrt(
    .hypothesis_log_odds_var_for_test(posterior[["theta"]] > 0.25,
                                      posterior[["theta"]] <= 0.25) +
      .hypothesis_log_odds_var_for_test(prior[["theta"]] > 0.25,
                                        prior[["theta"]] <= 0.25)
  )

  expect_equal(as.numeric(out[["BF_error"]]), expected_error, tolerance = 1e-12)
})


test_that("hypothesis_BF region error omits variance of exact prior masses", {

  set.seed(1)
  posterior   <- stats::rnorm(20000, mean = 0.3, sd = 0.1)
  prior_draws <- stats::rnorm(20000)
  theta_prior <- prior("normal", list(mean = 0, sd = 1))

  # The prior masses of 'theta > .2' and its complement come from the prior
  # object's distribution function, so only posterior indicators vary.
  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = theta_prior,
    hypothesis = "theta > 0.2",
    parameter  = "theta",
    seed       = 1
  )
  expect_equal(
    as.numeric(out[["BF_error"]]),
    100 * sqrt(.hypothesis_log_odds_var_for_test(posterior > 0.2,
                                                 posterior <= 0.2)),
    tolerance = 1e-12
  )

  quantity <- BayesTools:::.hypothesis_quantity_from_draws(
    posterior    = data.frame(theta = posterior),
    prior        = data.frame(theta = prior_draws),
    label        = "theta",
    parameter    = "theta",
    prior_object = theta_prior
  )
  statement <- hypothesis_parse("theta > 0.2 vs abs(theta) < 1")$statements[[1L]]
  # Only the compound side's prior mass is estimated from prior draws.
  expect_equal(
    BayesTools:::.hypothesis_region_odds_BF_error_percent(
      quantity, statement$left, statement$right
    ),
    100 * sqrt(
      .hypothesis_log_odds_var_for_test(posterior > 0.2, abs(posterior) < 1) +
        .hypothesis_log_prob_var_for_test(abs(prior_draws) < 1)
    ),
    tolerance = 1e-12
  )
  expect_identical(
    BayesTools:::.hypothesis_region_log_mass_mc_var(
      quantity, statement$left, prior = TRUE
    ),
    0
  )
  expect_equal(
    BayesTools:::.hypothesis_region_log_mass_mc_var(
      quantity, statement$right, prior = TRUE
    ),
    .hypothesis_log_prob_var_for_test(abs(prior_draws) < 1),
    tolerance = 1e-12
  )
})


test_that("hypothesis_BF respects inclusive and exclusive boundaries", {

  prior     <- data.frame(theta = c(-1, 0, 1, 2))
  posterior <- data.frame(theta = c(0, 0, 1, 1))

  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta >= 0 vs theta > 0"
  )

  expected <- (mean(posterior[["theta"]] >= 0) / mean(posterior[["theta"]] > 0)) /
    (mean(prior[["theta"]] >= 0) / mean(prior[["theta"]] > 0))

  expect_equal(attr(out, "raw_BF"), expected, tolerance = 1e-12)
})


test_that("hypothesis_BF estimates region BF error for overlapping regions", {

  prior     <- data.frame(theta = seq(-2, 2, length.out = 101))
  posterior <- data.frame(theta = seq(-1.5, 1.5, length.out = 101))

  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta > -0.5 vs theta < 0.5"
  )

  expected_error <- 100 * sqrt(
    .hypothesis_log_odds_var_for_test(posterior[["theta"]] > -0.5,
                                      posterior[["theta"]] < 0.5) +
      .hypothesis_log_odds_var_for_test(prior[["theta"]] > -0.5,
                                        prior[["theta"]] < 0.5)
  )

  expect_true(any(posterior[["theta"]] > -0.5 & posterior[["theta"]] < 0.5))
  expect_equal(as.numeric(out[["BF_error"]]), expected_error, tolerance = 1e-12)
})


test_that("hypothesis_BF supports interval complements", {

  prior     <- data.frame(theta = seq(-1, 1, length.out = 101))
  posterior <- data.frame(theta = seq(-0.8, 0.8, length.out = 101))

  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta > -0.5 & theta < 0.5"
  )

  prior_left      <- mean(prior[["theta"]] > -0.5 & prior[["theta"]] < 0.5)
  prior_right     <- 1 - prior_left
  posterior_left  <- mean(posterior[["theta"]] > -0.5 & posterior[["theta"]] < 0.5)
  posterior_right <- 1 - posterior_left
  expected        <- (posterior_left / posterior_right) /
    (prior_left / prior_right)

  expect_equal(attr(out, "raw_BF"), expected, tolerance = 1e-12)
})


test_that("hypothesis_BF composes point-vs-region tests by transitivity", {

  set.seed(4)
  prior     <- stats::rnorm(10000)
  posterior <- stats::rnorm(10000, mean = 0.25, sd = 1.1)

  expect_warning(point <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0",
    parameter  = "theta"
  ), class = "BayesTools_inexact_ordinate")
  expect_warning(transitive <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0 vs theta > 0",
    parameter  = "theta"
  ), class = "BayesTools_inexact_ordinate")

  region_unrestricted <- mean(posterior > 0) / mean(prior > 0)
  expect_equal(
    attr(transitive, "raw_BF"),
    (1 / attr(point, "raw_BF")) / region_unrestricted,
    tolerance = 1e-12
  )
  expect_true(is.na(transitive[["BF_error"]]))
})


test_that("hypothesis_BF reports boundary-valued posterior region evidence", {

  prior     <- data.frame(theta = c(-2, -1, 0, 1))
  posterior <- data.frame(theta = c(-3, -2, -1, 0))

  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta > 0 vs theta <= 0"
  )

  expect_equal(attr(out, "raw_BF"), 0)
  expect_true(is.na(out[["BF_error"]]))
  expect_match(attr(out, "warnings"), "Posterior region mass is zero")
})


test_that("hypothesis_BF reports undefined odds when both posterior regions are empty", {

  prior     <- data.frame(theta = c(-2, -1.5, 1.5, 2))
  posterior <- data.frame(theta = c(-0.5, 0, 0.5))

  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta > 1 vs theta < -1",
    columns    = "all"
  )

  expect_true(is.na(attr(out, "raw_BF")))
  expect_false(is.nan(attr(out, "raw_BF")))
  expect_true(is.na(out[["posterior"]]))
  expect_true(is.na(out[["BF_error"]]))
  expect_match(attr(out, "warnings"), "Both posterior region masses are zero")
})


test_that("hypothesis_BF rejects zero or one prior region mass", {

  prior     <- data.frame(theta = c(1, 2, 3))
  posterior <- data.frame(theta = c(1, 2, 3))

  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      prior      = prior,
      hypothesis = "theta > 0"
    ),
    "Prior region mass"
  )
})


test_that("explicit region comparisons accept an encompassing prior region", {

  set.seed(1)
  half_normal <- prior("normal", list(mean = 0, sd = 1), list(0, Inf))
  posterior   <- abs(stats::rnorm(20000, mean = 0.3, sd = 0.1))

  # Prior masses are analytic: P(theta > 0) = 1, P(theta > .5) = 2 pnorm(-.5).
  expected <- mean(posterior > 0.5) / (2 * stats::pnorm(-0.5))
  nested <- hypothesis_BF(
    posterior  = posterior,
    prior      = half_normal,
    hypothesis = "theta > 0.5 vs theta > 0",
    seed       = 1
  )
  expect_equal(attr(nested, "raw_BF"), expected, tolerance = 1e-12)

  point <- hypothesis_BF(
    posterior  = posterior,
    prior      = half_normal,
    hypothesis = "theta = 0.2",
    seed       = 1
  )
  transitive <- hypothesis_BF(
    posterior  = posterior,
    prior      = half_normal,
    hypothesis = "theta = 0.2 vs theta > 0",
    seed       = 1
  )
  expect_equal(attr(transitive, "raw_BF"), 1 / attr(point, "raw_BF"),
               tolerance = 1e-12)

  # Implicit statements keep the complement guard; zero mass stays invalid.
  expect_error(
    hypothesis_BF(posterior, half_normal, hypothesis = "theta > 0"),
    "complement has zero prior mass"
  )
  expect_error(
    hypothesis_BF(posterior, half_normal,
                  hypothesis = "theta > 0.5 vs theta < 0"),
    "is zero or non-finite"
  )

  # The same comparisons on a deterministic grid density (L16 makes the
  # encompassing grid mass exactly one).
  marginal <- .hypothesis_marginal_posterior_for_test(
    posterior,
    BayesTools:::.prior_linear_combination_density(
      prior_list = list(theta = half_normal),
      weights    = c(theta = 1)
    )
  )
  grid_nested <- hypothesis_BF(
    posterior  = marginal,
    hypothesis = "theta > 0.5 vs theta > 0",
    parameter  = "theta"
  )
  expect_equal(attr(grid_nested, "raw_BF"), expected, tolerance = 1e-4)
  expect_error(
    hypothesis_BF(marginal, hypothesis = "theta > 0", parameter = "theta"),
    "complement has zero prior mass"
  )
})


test_that("hypothesis_BF evaluates transformed and non-syntactic quantities", {

  posterior <- data.frame(
    `Level A` = c(1, 3, 0, 4),
    `Level B` = c(0, 1, 1, 3),
    check.names = FALSE
  )
  prior <- data.frame(
    `Level A` = c(0, 2, 3, 0),
    `Level B` = c(1, 1, 1, 2),
    check.names = FALSE
  )

  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "`Level A` - `Level B` > 0 vs `Level A` - `Level B` < 0"
  )

  posterior_diff <- posterior[["Level A"]] - posterior[["Level B"]]
  prior_diff     <- prior[["Level A"]] - prior[["Level B"]]
  expected       <- (mean(posterior_diff > 0) / mean(posterior_diff < 0)) /
    (mean(prior_diff > 0) / mean(prior_diff < 0))

  expect_equal(attr(out, "raw_BF"), expected, tolerance = 1e-12)
})


test_that("hypothesis_BF rejects duplicate quantity names before evaluation", {

  posterior <- data.frame(
    first  = c(-1, 0, 1),
    second = c(10, 11, 12),
    check.names = FALSE
  )
  names(posterior) <- c("theta", "theta")

  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      hypothesis = "theta > 0"
    ),
    "posterior.*duplicate quantity names.*theta"
  )

  posterior_matrix <- as.matrix(posterior)
  expect_error(
    hypothesis_BF(
      posterior  = posterior_matrix,
      hypothesis = "theta > 0"
    ),
    "posterior.*duplicate quantity names.*theta"
  )

  prior <- data.frame(
    first  = c(-1, 0, 1),
    second = c(10, 11, 12),
    check.names = FALSE
  )
  names(prior) <- c("theta", "theta")

  expect_error(
    hypothesis_BF(
      posterior  = data.frame(theta = c(-1, 0, 1)),
      prior      = prior,
      hypothesis = "theta > 0"
    ),
    "prior.*duplicate quantity names.*theta"
  )

  marginal_with_prior <- .hypothesis_marginal_posterior_for_test(
    seq(-2, 2, length.out = 201),
    .hypothesis_prior_density_for_test()
  )
  expect_s3_class(
    hypothesis_BF(
      posterior  = marginal_with_prior,
      prior      = prior,
      hypothesis = "theta = 0",
      parameter  = "theta"
    ),
    "BayesTools_hypothesis_BF"
  )

  marginal <- list(
    alternate = structure(
      c(-1, 0, 1),
      class = c("marginal_posterior.simple", "numeric")
    ),
    alternate = structure(
      c(1, 2, 3),
      class = c("marginal_posterior.simple", "numeric")
    )
  )
  class(marginal) <- c(
    "list", "marginal_posterior.factor", "marginal_posterior"
  )
  attr(marginal, "parameter") <- "mu"

  expect_error(
    hypothesis_BF(
      posterior  = marginal,
      hypothesis = "mu[alternate] > 0"
    ),
    "posterior.*duplicate quantity names.*alternate"
  )
})


test_that("hypothesis_BF rejects non-finite posterior draws before computation", {

  for(nonfinite in c(Inf, -Inf, NaN)){
    expect_error(
      hypothesis_BF(
        posterior  = c(-1, nonfinite, 1),
        hypothesis = "theta > 0",
        parameter  = "theta"
      ),
      "Posterior draws must contain only finite values"
    )
  }

  expect_error(
    hypothesis_BF(
      posterior = data.frame(
        theta    = c(-1, 0, 1),
        nuisance = c(0, Inf, 0)
      ),
      hypothesis = "theta > 0"
    ),
    "Posterior draws must contain only finite values"
  )

  marginal <- structure(
    c(-1, -Inf, 1),
    class = c(
      "marginal_posterior.simple", "marginal_posterior", "numeric"
    )
  )
  expect_error(
    hypothesis_BF(
      posterior  = marginal,
      hypothesis = "theta > 0",
      parameter  = "theta"
    ),
    "Posterior draws must contain only finite values"
  )
})


test_that("hypothesis_BF evaluates explicit marginal posterior level comparisons", {

  context <- BayesTools:::.prior_density_context(
    prior_list   = list(
      alt  = prior("normal", list(mean = 0, sd = 1)),
      rand = prior("normal", list(mean = 0, sd = 1))
    ),
    column_names = c("alt", "rand"),
    n_grid       = 128
  )
  posterior <- list(
    alternate = .bt_meta_update(
      structure(c(rep(1, 80), rep(-1, 20)), class = c("marginal_posterior.simple", "numeric")),
      linear_weights = matrix(c(1, 0), nrow = 1,
                              dimnames = list(NULL, c("alt", "rand")))
    ),
    random = .bt_meta_update(
      structure(rep(0, 100), class = c("marginal_posterior.simple", "numeric")),
      linear_weights = c(alt = 0, rand = 1)
    )
  )
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter")             <- "mu_alloc"
  posterior <- .bt_meta_set(posterior, "prior_context", context)

  out <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = "mu_alloc[alternate] > mu_alloc[random]",
    columns    = "all",
    seed       = 11
  )

  expect_equal(out[["posterior"]], 4, tolerance = 1e-12)
  expect_equal(attr(out, "raw_BF"), out[["posterior"]] / out[["prior"]],
               tolerance = 1e-12)
  expect_true(out[["prior"]] > .90 && out[["prior"]] < 1.10)
  expect_true(is.finite(out[["BF_error"]]))
  expect_equal(out[["method"]], "prior-posterior odds")
})


test_that("hypothesis_BF references level names that contain brackets", {

  context <- BayesTools:::.prior_density_context(
    prior_list   = list(
      a = prior("normal", list(mean = 0, sd = 1)),
      b = prior("normal", list(mean = 0, sd = 1))
    ),
    column_names = c("a", "b"),
    n_grid       = 128
  )
  make_posterior <- function(level_names){
    set.seed(2)
    posterior <- list(
      .bt_meta_update(
        structure(stats::rnorm(4000, 0.5, 0.2), class = c("marginal_posterior.simple", "numeric")),
        linear_weights = c(a = 1, b = 0),
        atoms = posterior_atom_attribute()
      ),
      .bt_meta_update(
        structure(stats::rnorm(4000, 0, 0.2), class = c("marginal_posterior.simple", "numeric")),
        linear_weights = c(a = 0, b = 1),
        atoms = posterior_atom_attribute()
      )
    )
    names(posterior) <- level_names
    class(posterior) <- c("list", "marginal_posterior.factor",
                          "marginal_posterior")
    attr(posterior, "parameter")             <- "mu"
    posterior <- .bt_meta_set(posterior, "prior_context", context)
    posterior
  }

  reference <- hypothesis_BF(
    make_posterior(c("A", "B")),
    hypothesis = "mu[A] > mu[B]",
    seed       = 1,
    columns    = "all"
  )
  # cut() levels; backticks or the parameter catalog's escaped component
  # form ("(0,1%5D", optionally quoted) reference them.
  posterior <- make_posterior(c("(0,1]", "(1,2]"))
  for(hypothesis in c(
    "`mu[(0,1]]` > `mu[(1,2]]`",
    "mu[\"(0,1%5D\"] > mu[\"(1,2%5D\"]",
    "mu[(0,1%5D] > `mu[(1,2]]`"
  )){
    out <- hypothesis_BF(posterior, hypothesis = hypothesis, seed = 1,
                         columns = "all")
    expect_equal(attr(out, "raw_BF"), attr(reference, "raw_BF"),
                 tolerance = 1e-12, info = hypothesis)
    expect_equal(out[["method"]], "prior-posterior odds", info = hypothesis)
  }
  expect_identical(
    hypothesis_parse_level_reference("`mu[(0,1]]`")[c("parameter", "level")],
    data.frame(parameter = "mu", level = "(0,1]", stringsAsFactors = FALSE)
  )
  expect_error(
    hypothesis_BF(posterior, hypothesis = "`mu[(0,2]]` > `mu[(1,2]]`"),
    "unknown level '(0,2]'",
    fixed = TRUE,
    class = "BayesTools_parameter_not_found"
  )

  contrast <- hypothesis_linear_target(
    posterior, "`mu[(0,1]]` - `mu[(1,2]]` = 0.1", "mu"
  )
  reference_contrast <- hypothesis_linear_target(
    make_posterior(c("A", "B")), "mu[A] - mu[B] = 0.1", "mu"
  )
  expect_identical(contrast$weights, reference_contrast$weights)
  expect_equal(as.numeric(contrast$posterior),
               as.numeric(reference_contrast$posterior))

  # an unknown level of a linear target is an unresolved reference, as in
  # hypothesis_BF()
  condition <- tryCatch(
    hypothesis_linear_target(
      posterior, "`mu[(0,2]]` - `mu[(1,2]]` = 0.1", "mu"
    ),
    error = identity
  )
  expect_identical(
    class(condition),
    c("BayesTools_parameter_not_found",
      "BayesTools_parameter_resolution_error", "error", "condition")
  )
  expect_identical(
    conditionMessage(condition),
    "Hypothesis references unknown level '(0,2]' for parameter 'mu'."
  )
  expect_identical(condition$alias, "mu[(0,2]]")
  expect_identical(condition$available, c("mu[(0,1]]", "mu[(1,2]]"))
})


test_that("hypothesis_BF rejects conditional level comparisons with different conditionals", {

  context <- BayesTools:::.prior_density_context(
    prior_list   = list(
      alt  = prior("normal", list(mean = 0, sd = 1)),
      rand = prior("normal", list(mean = 0, sd = 1))
    ),
    column_names = c("alt", "rand"),
    n_grid       = 128
  )
  posterior <- list(
    alternate = .bt_meta_update(
      structure(c(rep(1, 80), rep(-1, 20)), class = c("marginal_posterior.simple", "numeric")),
      linear_weights = c(alt = 1, rand = 0),
      condition = list(effective_conditional = "mu_intercept")
    ),
    random = .bt_meta_update(
      structure(rep(0, 100), class = c("marginal_posterior.simple", "numeric")),
      linear_weights = c(alt = 0, rand = 1),
      condition = list(effective_conditional = c("mu_intercept", "mu_alloc"))
    )
  )
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter")             <- "mu_alloc"
  posterior <- .bt_meta_set(posterior, "prior_context", context)

  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      hypothesis = "mu_alloc[alternate] > mu_alloc[random]"
    ),
    "different conditional posterior subsets"
  )
})


test_that("hypothesis_BF treats same conditional labels with AND and OR as different events", {

  prior_list <- list(
    theta = prior_spike_and_slab(
      prior("normal", list(mean = 0, sd = 1)),
      prior_inclusion = prior("point", list(location = 0.5))
    ),
    phi = prior_spike_and_slab(
      prior("normal", list(mean = 0, sd = 1)),
      prior_inclusion = prior("point", list(location = 0.5))
    )
  )
  and_event <- BayesTools:::.condition_event(
    prior_list       = prior_list,
    conditional      = c("theta", "phi"),
    conditional_rule = "AND"
  )
  or_event <- BayesTools:::.condition_event(
    prior_list       = prior_list,
    conditional      = c("theta", "phi"),
    conditional_rule = "OR"
  )
  posterior <- list(
    and = .bt_meta_update(
      structure(c(rep(1, 80), rep(-1, 20)), class = c("marginal_posterior.simple", "numeric")),
      linear_weights = c(theta = 1, phi = 0),
      condition = list(effective_conditional = c("theta", "phi"), effective_conditional_rule = "AND", condition_key = and_event[["condition_key"]], resolved_condition_event = and_event)
    ),
    or = .bt_meta_update(
      structure(rep(0, 100), class = c("marginal_posterior.simple", "numeric")),
      linear_weights = c(theta = 0, phi = 1),
      condition = list(effective_conditional = c("theta", "phi"), effective_conditional_rule = "OR", condition_key = or_event[["condition_key"]], resolved_condition_event = or_event)
    )
  )
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter") <- "mu_alloc"

  expect_false(identical(
    .bt_meta_condition(posterior[["and"]], "resolved_condition_event"),
    .bt_meta_condition(posterior[["or"]], "resolved_condition_event")
  ))
  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      hypothesis = "mu_alloc[and] > mu_alloc[or]"
    ),
    "different conditional posterior subsets"
  )
})


test_that("hypothesis_BF uses child precomputed density for explicit level point null", {

  prior_density <- .hypothesis_prior_density_for_test()
  alternate <- .hypothesis_marginal_posterior_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  alternate <- .bt_meta_set(alternate, "posterior_ordinate", .posterior_ordinate_for_test(
    value       = 0,
    ordinate    = 0.50,
    method      = "IWMDE",
    diagnostics = list(relative_mcse = 0.03)
  ))
  posterior <- list(alternate = alternate)
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter") <- "mu_alloc"

  out <- hypothesis_BF(
    posterior      = posterior,
    hypothesis     = "mu_alloc[alternate] = 0",
    columns        = "all",
    density_method = "precomputed"
  )

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) / 0.50

  expect_equal(attr(out, "raw_BF"), expected, tolerance = 1e-12)
  expect_equal(out[["posterior"]], 0.50, tolerance = 1e-12)
  expect_equal(as.numeric(out[["BF_error"]]), 3, tolerance = 1e-12)
  expect_equal(out[["method"]], "Savage-Dickey (precomputed)")
})


test_that("hypothesis_linear_target compiles an exact atom-free linear target", {

  context <- BayesTools:::.prior_density_context(
    prior_list   = list(
      alt  = prior("normal", list(mean = 0, sd = 1)),
      rand = prior("normal", list(mean = 0, sd = 1))
    ),
    column_names = c("alt", "rand"),
    n_grid       = 128
  )
  posterior <- list(
    alternate = .bt_meta_update(
      structure(c(rep(1, 80), rep(-1, 20)), class = c("marginal_posterior.simple", "numeric")),
      linear_weights = c(alt = 1, rand = 0),
      atoms = posterior_atom_attribute()
    ),
    random = .bt_meta_update(
      structure(rep(0, 100), class = c("marginal_posterior.simple", "numeric")),
      linear_weights = c(alt = 0, rand = 1),
      atoms = posterior_atom_attribute()
    )
  )
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter")             <- "mu_alloc"
  posterior <- .bt_meta_set(posterior, "prior_context", context)

  # the two-level contrast (the former hypothesis_level_contrast()) keeps its
  # numbers: weights (1, -1), the N(0, 2) prior ordinate at 0
  target <- hypothesis_linear_target(
    posterior  = posterior,
    hypothesis = paste(
      "mu_alloc[alternate] > mu_alloc[random] vs",
      "mu_alloc[alternate] = mu_alloc[random]"
    ),
    parameter  = "mu_alloc"
  )

  expect_equal(as.numeric(target$posterior),
               as.numeric(posterior$alternate - posterior$random))
  expect_identical(target$weights, c(alt = 1, rand = -1))
  expect_identical(target$parameter, ".BayesTools_linear_target")
  expect_identical(
    hypothesis_render(target$hypothesis),
    paste(
      ".BayesTools_linear_target > 0 vs",
      ".BayesTools_linear_target = 0"
    )
  )
  expect_identical(posterior_metadata(target$posterior, "linear_weights"), c(alt = 1, rand = -1))
  expect_identical(.bt_meta_get(target$posterior, "linear_offset"), 0)
  expect_identical(posterior_metadata(target$posterior, "prior_context"), context)
  expect_true(posterior_atoms_free(target$posterior))
  ordinate <- prior_density_ordinate(
    .bt_meta_get(target$posterior, "prior_density"),
    0
  )
  expect_true(ordinate$exact)
  expect_equal(ordinate$log_density, stats::dnorm(0, sd = sqrt(2), log = TRUE),
               tolerance = 1e-8)
  out <- hypothesis_BF(
    posterior  = target$posterior,
    hypothesis = target$hypothesis,
    parameter  = target$parameter,
    columns    = "all"
  )
  expect_true(is.finite(attr(out, "raw_BF")))
  expect_identical(out$method, "transitive Savage-Dickey")

  precise <- hypothesis_linear_target(
    posterior = posterior,
    hypothesis = paste(
      "mu_alloc[alternate] - mu_alloc[random] > 123456789.123 vs",
      "mu_alloc[alternate] - mu_alloc[random] = 0.123456789"
    ),
    parameter = "mu_alloc"
  )
  precise_statement <- precise$hypothesis$statements[[1L]]
  expect_identical(precise_statement$left$expression$right$value, 123456789.123)
  expect_identical(precise_statement$right$value, 0.123456789)

  # scaled levels: 2 alt - rand ~ N(0, 5)
  scaled <- hypothesis_linear_target(
    posterior  = posterior,
    hypothesis = paste(
      "2 * mu_alloc[alternate] - mu_alloc[random] > 0 vs",
      "2 * mu_alloc[alternate] - mu_alloc[random] = 0.5"
    ),
    parameter  = "mu_alloc"
  )
  expect_identical(scaled$weights, c(alt = 2, rand = -1))
  expect_equal(as.numeric(scaled$posterior),
               2 * as.numeric(posterior$alternate) - as.numeric(posterior$random))
  expect_identical(scaled$hypothesis$statements[[1L]]$right$value, 0.5)
  expect_equal(
    prior_density_ordinate(.bt_meta_get(scaled$posterior, "prior_density"), .5)$log_density,
    stats::dnorm(.5, sd = sqrt(5), log = TRUE),
    tolerance = 1e-8
  )

  # the average of two levels with a constant moved to the value:
  # (alt + rand) / 2 + 1 = 1.25  <=>  (alt + rand) / 2 = 0.25, N(0, 1/2)
  average <- hypothesis_linear_target(
    posterior  = posterior,
    hypothesis = "(mu_alloc[alternate] + mu_alloc[random]) / 2 + 1 = 1.25",
    parameter  = "mu_alloc"
  )
  expect_identical(average$weights, c(alt = .5, rand = .5))
  expect_identical(average$hypothesis$statements[[1L]]$left$value, 0.25)
  expect_equal(
    prior_density_ordinate(.bt_meta_get(average$posterior, "prior_density"), .25)$log_density,
    stats::dnorm(.25, sd = sqrt(.5), log = TRUE),
    tolerance = 1e-8
  )

  expect_error(
    hypothesis_linear_target(
      posterior  = posterior,
      hypothesis = paste(
        "mu_alloc[alternate] > mu_alloc[random] vs",
        "mu_alloc[alternate] + mu_alloc[random] = 0"
      ),
      parameter  = "mu_alloc"
    ),
    "Hypothesis statements must all use the same linear combination of levels.",
    fixed = TRUE
  )
  expect_error(
    hypothesis_linear_target(
      posterior  = posterior,
      hypothesis = "mu_alloc[alternate] * mu_alloc[random] = 0",
      parameter  = "mu_alloc"
    ),
    "A linear target must be a linear combination of levels of 'mu_alloc' and numbers.",
    fixed = TRUE
  )
  expect_error(
    hypothesis_linear_target(
      posterior  = posterior,
      hypothesis = "mu_alloc[alternate] - mu[random] = 0",
      parameter  = "mu_alloc"
    ),
    "A linear target may reference levels of only 'mu_alloc'.",
    fixed = TRUE
  )

  posterior_with_baseline <- posterior
  posterior_with_baseline$baseline <- .bt_meta_update(
    structure(rep(0, 100), class = c("marginal_posterior.simple", "numeric")),
    linear_weights = c(alt = 0, rand = 0),
    atoms = posterior_atom_attribute(
      data.frame(x = 0, mass = 1)
    )
  )
  baseline_target <- hypothesis_linear_target(
    posterior  = posterior_with_baseline,
    hypothesis = paste(
      "mu_alloc[random] > mu_alloc[baseline] vs",
      "mu_alloc[random] = mu_alloc[baseline]"
    ),
    parameter  = "mu_alloc"
  )
  expect_equal(
    as.vector(baseline_target$posterior),
    as.vector(posterior$random)
  )

  posterior_without_declaration <- posterior
  posterior_without_declaration$alternate <- .bt_meta_set(posterior_without_declaration$alternate, "atoms", NULL)
  expect_error(
    hypothesis_linear_target(
      posterior  = posterior_without_declaration,
      hypothesis = paste(
        "mu_alloc[alternate] > mu_alloc[random] vs",
        "mu_alloc[alternate] = mu_alloc[random]"
      ),
      parameter  = "mu_alloc"
    ),
    "posterior-atom declarations"
  )
})

test_that("prior_ordinate_status reports the point-hypothesis exactness rule per value", {

  columns <- c("value", "eligible", "condition", "reason", "continuous_behavior")
  point_mass <- paste0(
    "There is a point mass in the prior at the exact null hypothesis value. ",
    "The Savage-Dickey density ratio is invalid."
  )
  spike <- prior_spike_and_slab(
    prior("normal", list(0, 1)),
    prior_inclusion = prior("spike", list(.5))
  )
  spike_density <- BayesTools:::.prior_linear_combination_density(
    list(x = spike), c(x = 1)
  )
  status <- prior_ordinate_status(spike_density, c(0, .5, -1))
  expect_s3_class(status, "data.frame")
  expect_identical(names(status), columns)
  expect_identical(status$value, c(0, .5, -1))
  expect_identical(status$eligible, c(FALSE, TRUE, TRUE))
  expect_identical(status$condition,
                   c("BayesTools_point_mass_at_null", NA, NA))
  expect_identical(status$reason, c(point_mass, NA, NA))
  # the slab is regular at the atom
  expect_identical(status$continuous_behavior, rep("regular", 3L))

  gamma_half <- prior("gamma", list(.5, 1))
  status <- prior_ordinate_status(gamma_half, c(0, -1, 1),
                                  labels = c("s = 0", "s = -1", "s = 1"))
  expect_identical(status$eligible, c(FALSE, FALSE, TRUE))
  expect_identical(status$condition, c(
    "BayesTools_infinite_ordinate", "BayesTools_zero_ordinate", NA
  ))
  expect_identical(status$reason, c(
    "Prior density at point hypothesis 's = 0' is infinite, so the Savage-Dickey density ratio is undefined.",
    "Prior density at point hypothesis 's = -1' is zero, so the Savage-Dickey density ratio is undefined.",
    NA
  ))
  expect_identical(status$continuous_behavior, c("infinite", "zero", "regular"))
  # default labels are the values
  expect_match(prior_ordinate_status(gamma_half, 0)$reason,
               "point hypothesis '0' is infinite", fixed = TRUE)

  # the stopping rule of hypothesis_BF() and Savage_Dickey_BF() is the first
  # ineligible row, with the same class and message
  t3 <- prior("t", list(0, 1, 3))
  three_t <- BayesTools:::.prior_linear_combination_density(
    list(a = t3, b = t3, c = t3), c(a = 1, b = 1, c = 1)
  )
  nonnegative_spike <- prior_spike_and_slab(
    prior("normal", list(0, 1), list(0, Inf)),
    prior_inclusion = prior("spike", list(.5))
  )
  undefined <- BayesTools:::.prior_linear_combination_density(
    list(x = nonnegative_spike), c(x = 1), output_transformation = "exp_lin",
    output_transformation_arguments = list(a = 0, b = 2)
  )
  cases <- list(
    list(spike_density, 0, "BayesTools_point_mass_at_null"),
    list(gamma_half, 0, "BayesTools_infinite_ordinate"),
    list(gamma_half, -1, "BayesTools_zero_ordinate"),
    list(three_t, 0, "BayesTools_inexact_ordinate"),
    list(undefined, 0, "BayesTools_point_mass_at_null")
  )
  for(case in cases){
    status <- prior_ordinate_status(case[[1L]], case[[2L]], labels = "theta = v")
    expect_identical(status$condition, case[[3L]])
    expect_false(status$eligible)
    condition <- tryCatch(
      BayesTools:::.hypothesis_check_prior_ordinate(case[[1L]], case[[2L]], "theta = v"),
      error = function(e) e
    )
    expect_identical(class(condition)[[1L]], status$condition)
    expect_s3_class(condition, "BayesTools_hypothesis_ordinate")
    expect_identical(conditionMessage(condition), status$reason)
  }
  expect_identical(prior_ordinate_status(three_t, 0)$continuous_behavior, "unknown")
  expect_true(prior_ordinate_status(undefined, 1)$eligible)
  expect_equal(exp(prior_density_ordinate(undefined, 1)$log_density), stats::dnorm(1) / 2, tolerance = 1e-12)
  # an eligible value returns its ordinate: the slab density times its weight
  expect_equal(
    BayesTools:::.hypothesis_check_prior_ordinate(spike_density, .5, "x = 0.5")$log_density,
    stats::dnorm(.5, log = TRUE) + log(.5),
    tolerance = 1e-12
  )

  expect_error(prior_ordinate_status(spike_density, Inf),
               "The 'values' argument must contain only finite values.", fixed = TRUE)
  expect_error(prior_ordinate_status(spike_density, numeric()),
               "The 'values' argument must contain at least one value.", fixed = TRUE)
  expect_error(prior_ordinate_status(spike_density, NULL),
               "The 'values' argument cannot be NULL.", fixed = TRUE)
  expect_error(prior_ordinate_status(spike_density, c(0, 1), labels = "a"),
               "The 'labels' argument must have length '2'.", fixed = TRUE)
  expect_error(prior_ordinate_status(stats::rnorm(10), 0),
               "The 'prior_density' argument must be a BayesTools prior or prior_linear_density object.",
               fixed = TRUE)
})

test_that("point-hypothesis refusals of regular ordinates without a value give the ordinate's reason", {

  # a gamma(3, 0.7) density at the subnormal value 1e-320 has no value (the
  # full-precision rule of prior_density_ordinate()); the refusal states that
  # reason instead of a missing structural value
  gamma_prior <- prior("gamma", list(3, .7))
  status <- prior_ordinate_status(gamma_prior, 1e-320, labels = "s = v")
  expect_false(status$eligible)
  expect_identical(status$condition, "BayesTools_inexact_ordinate")
  expect_identical(status$reason, paste0(
    "Prior density at point hypothesis 's = v' is unavailable: its regular prior ordinate has no ",
    "value (The value at which the density is evaluated is not representable at full precision ",
    "in ordinary floating-point arithmetic). Test a region hypothesis instead."
  ))
  condition <- tryCatch(
    BayesTools:::.hypothesis_check_prior_ordinate(gamma_prior, 1e-320, "s = v"),
    error = function(e) e
  )
  expect_s3_class(condition, "BayesTools_inexact_ordinate")
  expect_identical(conditionMessage(condition), status$reason)
})

test_that("linear target refusals are classed with their reason", {

  normal_context <- BayesTools:::.prior_density_context(
    prior_list   = list(
      alt  = prior("normal", list(mean = 0, sd = 1)),
      rand = prior("normal", list(mean = 0, sd = 1))
    ),
    column_names = c("alt", "rand"),
    n_grid       = 128
  )
  spike_context <- BayesTools:::.prior_density_context(
    prior_list   = list(
      alt  = prior_spike_and_slab(
        prior("normal", list(mean = 0, sd = 1)),
        prior_inclusion = prior("point", list(location = 0.5))
      ),
      rand = prior("point", list(location = 0))
    ),
    column_names = c("alt", "rand"),
    n_grid       = 128
  )
  level <- function(values, weights, context = NULL){
    out <- .bt_meta_update(
      structure(values, class = c("marginal_posterior.simple", "numeric")),
      linear_weights = weights,
      atoms = posterior_atom_attribute()
    )
    .bt_meta_set(out, "prior_context", context)
  }
  factor_posterior <- function(levels, context = NULL){
    class(levels) <- c("list", "marginal_posterior.factor", "marginal_posterior")
    attr(levels, "parameter") <- "mu"
    .bt_meta_set(levels, "prior_context", context)
  }
  levels <- list(
    a = level(c(rep(1, 80), rep(-1, 20)), c(alt = 1, rand = 0)),
    b = level(rep(0, 100), c(alt = 0, rand = 1))
  )
  expect_reason <- function(posterior, reason, message){
    condition <- tryCatch(
      hypothesis_linear_target(posterior, "mu[a] - mu[b] = 0", "mu"),
      error = function(e) e
    )
    expect_s3_class(condition, "BayesTools_linear_target_unavailable")
    expect_s3_class(condition, "BayesTools_hypothesis_target")
    expect_identical(condition$reason, reason)
    expect_identical(conditionMessage(condition), message)
  }

  expect_reason(
    factor_posterior(levels),
    "prior_context",
    "A valid joint prior context is required for a linear target."
  )
  partial <- levels
  partial$a <- .bt_meta_set(partial$a, "prior_context", normal_context)
  expect_reason(
    factor_posterior(partial, normal_context),
    "prior_context",
    "Level comparisons require prior contexts for all conditional levels."
  )
  undeclared <- levels
  undeclared$a <- .bt_meta_set(undeclared$a, "atoms", NULL)
  expect_reason(
    factor_posterior(undeclared, normal_context),
    "atom_declarations",
    "Complete structural posterior-atom declarations are required for a linear target."
  )
  expect_reason(
    factor_posterior(levels, spike_context),
    "posterior_atoms",
    "The linear target prior is not structurally atom-free."
  )
  # level weights on a column the joint context does not contain
  foreign <- levels
  foreign$b <- level(rep(0, 100), c(alt = 0, rand = 0, other = 1))
  expect_reason(
    factor_posterior(foreign, normal_context),
    "prior_context",
    "Linear prior weights reference columns not available in the joint prior context: other."
  )
})

test_that("hypothesis_BF uses parent precomputed metadata for level point nulls", {

  prior_density <- .hypothesis_prior_density_for_test()
  alternate <- .hypothesis_marginal_posterior_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  random <- .hypothesis_marginal_posterior_for_test(
    seq(-2, 2, length.out = 301),
    prior_density
  )
  alternate <- .bt_meta_set(alternate, "posterior_ordinate", .posterior_ordinate_for_test(
    value    = 1,
    ordinate = 100,
    method   = "wrong-null"
  ))
  posterior <- list(alternate = alternate, random = random)
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter") <- "mu_alloc"
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", list(
    .posterior_ordinate_for_test(
      parameter   = "mu_alloc[alternate]",
      value       = 0,
      ordinate    = 0.50,
      method      = "IWMDE",
      diagnostics = list(relative_mcse = 0.03)
    ),
    .posterior_ordinate_for_test(
      parameter   = "mu_alloc[random]",
      value       = 0,
      ordinate    = 0.25,
      method      = "IWMDE",
      diagnostics = list(relative_mcse = 0.04)
    )
  ))

  explicit <- hypothesis_BF(
    posterior      = posterior,
    hypothesis     = "mu_alloc[alternate] = 0",
    columns        = "all",
    density_method = "precomputed"
  )
  expanded <- hypothesis_BF(
    posterior      = posterior,
    hypothesis     = "mu_alloc = 0",
    parameter      = "mu_alloc",
    columns        = "all",
    density_method = "precomputed"
  )

  expected_alternate <- BayesTools:::.prior_linear_density_height(prior_density, 0) / 0.50
  expected_random <- BayesTools:::.prior_linear_density_height(prior_density, 0) / 0.25

  expect_equal(attr(explicit, "raw_BF"), expected_alternate, tolerance = 1e-12)
  expect_equal(explicit[["posterior"]], 0.50, tolerance = 1e-12)
  expect_equal(as.numeric(explicit[["BF_error"]]), 3, tolerance = 1e-12)
  expect_equal(expanded[["posterior"]], c(0.50, 0.25), tolerance = 1e-12)
  expect_equal(attr(expanded, "raw_BF"), c(expected_alternate, expected_random),
               tolerance = 1e-12)
  expect_equal(expanded[["method"]], rep("Savage-Dickey (precomputed)", 2))
})

test_that("hypothesis_BF uses parent posterior_densities for explicit level point null", {

  prior_density <- .hypothesis_prior_density_for_test()
  alternate <- .hypothesis_marginal_posterior_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  random <- .hypothesis_marginal_posterior_for_test(
    seq(-2, 2, length.out = 301),
    prior_density
  )
  posterior <- list(alternate = alternate, random = random)
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter") <- "mu_alloc"
  posterior <- .bt_meta_set(posterior, "posterior_densities", list(list(
    .posterior_density_for_test(
      parameter = "mu_alloc[alternate]",
      x         = seq(-1, 1, length.out = 101),
      y         = rep(0.50, 101),
      method    = "qCMDE"
    ),
    .posterior_density_for_test(
      parameter = "mu_alloc[random]",
      x         = seq(-1, 1, length.out = 101),
      y         = rep(0.25, 101),
      method    = "qCMDE"
    )
  )))

  out <- hypothesis_BF(
    posterior      = posterior,
    hypothesis     = "mu_alloc[alternate] = 0",
    columns        = "all",
    density_method = "precomputed"
  )

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) / 0.50

  expect_equal(attr(out, "raw_BF"), expected, tolerance = 1e-12)
  expect_equal(out[["posterior"]], 0.50, tolerance = 1e-12)
  expect_equal(out[["method"]], "Savage-Dickey (precomputed)")
})


test_that("hypothesis_BF infers marginal_inference parameter from bracket syntax", {

  context <- BayesTools:::.prior_density_context(
    prior_list   = list(
      alt  = prior("normal", list(mean = 0, sd = 1)),
      rand = prior("normal", list(mean = 0, sd = 1))
    ),
    column_names = c("alt", "rand"),
    n_grid       = 128
  )
  posterior <- list(
    alternate = .bt_meta_update(
      structure(c(rep(1, 75), rep(-1, 25)), class = c("marginal_posterior.simple", "numeric")),
      linear_weights = c(alt = 1, rand = 0)
    ),
    random = .bt_meta_update(
      structure(rep(0, 100), class = c("marginal_posterior.simple", "numeric")),
      linear_weights = c(alt = 0, rand = 1)
    )
  )
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter")             <- "mu_alloc"
  posterior <- .bt_meta_set(posterior, "prior_context", context)
  inference <- list(
    averaged    = list(mu_alloc = posterior),
    conditional = list(mu_alloc = posterior),
    inference   = list()
  )
  class(inference) <- c("list", "marginal_inference")

  out <- hypothesis_BF(
    posterior  = inference,
    hypothesis = "mu_alloc[alternate] > mu_alloc[random]",
    columns    = "all",
    seed       = 12
  )

  expect_equal(out[["posterior"]], 3, tolerance = 1e-12)
  expect_equal(out[["method"]], "prior-posterior odds")

  # a 'parameter' the inference does not contain is an unresolved reference
  condition <- tryCatch(
    hypothesis_BF(
      posterior  = inference,
      hypothesis = "mu_alloc[alternate] > mu_alloc[random]",
      parameter  = "mu_other"
    ),
    error = identity
  )
  expect_identical(
    class(condition),
    c("BayesTools_parameter_not_found",
      "BayesTools_parameter_resolution_error", "error", "condition")
  )
  expect_identical(
    conditionMessage(condition),
    "Parameter 'mu_other' is not available in 'posterior'."
  )
  expect_identical(condition$alias, "mu_other")
  expect_identical(condition$available, "mu_alloc")
})


test_that("hypothesis_BF samples level priors from mixture and conditional contexts", {

  alt_1  <- BayesTools:::.set_prior_model_weight(prior("normal", list(-1, 1)), 1)
  alt_2  <- BayesTools:::.set_prior_model_weight(prior("normal", list(1, 1)), 1)
  rand_1 <- BayesTools:::.set_prior_model_weight(prior("normal", list(0, 1)), 1)
  rand_2 <- BayesTools:::.set_prior_model_weight(prior("normal", list(0, 1)), 1)
  mixture_context <- BayesTools:::.prior_density_model_mixture_context(
    prior_list   = list(alt = list(alt_1, alt_2), rand = list(rand_1, rand_2)),
    column_names = c("alt", "rand"),
    n_grid       = 128
  )
  conditional_context <- BayesTools:::.prior_density_build_context(
    prior_list        = list(
      alt  = prior("normal", list(mean = 0, sd = 1)),
      rand = prior("normal", list(mean = 0, sd = 1))
    ),
    column_names      = c("alt", "rand"),
    conditional       = "alt",
    conditional_rule  = "OR",
    n_grid            = 128
  )

  for(context in list(mixture_context, conditional_context)){
    posterior <- list(
      alternate = .bt_meta_update(
        structure(c(rep(1, 75), rep(-1, 25)), class = c("marginal_posterior.simple", "numeric")),
        linear_weights = c(alt = 1, rand = 0)
      ),
      random = .bt_meta_update(
        structure(rep(0, 100), class = c("marginal_posterior.simple", "numeric")),
        linear_weights = c(alt = 0, rand = 1)
      )
    )
    class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
    attr(posterior, "parameter")             <- "mu_alloc"
    posterior <- .bt_meta_set(posterior, "prior_context", context)

    out <- hypothesis_BF(
      posterior  = posterior,
      hypothesis = "mu_alloc[alternate] > mu_alloc[random]",
      columns    = "all",
      seed       = 13
    )

    expect_true(is.finite(attr(out, "raw_BF")))
    expect_equal(out[["method"]], "prior-posterior odds")
  }
})


test_that("hypothesis_BF rejects missing nonzero level weight columns", {

  context <- BayesTools:::.prior_density_context(
    prior_list   = list(alt = prior("normal", list(mean = 0, sd = 1))),
    column_names = "alt",
    n_grid       = 128
  )
  posterior <- list(
    alternate = .bt_meta_update(
      structure(c(rep(1, 75), rep(-1, 25)), class = c("marginal_posterior.simple", "numeric")),
      linear_weights = c(alt = 1)
    ),
    random = .bt_meta_update(
      structure(rep(0, 100), class = c("marginal_posterior.simple", "numeric")),
      linear_weights = c(rand = 1)
    )
  )
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter")             <- "mu_alloc"
  posterior <- .bt_meta_set(posterior, "prior_context", context)

  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      hypothesis = "mu_alloc[alternate] > mu_alloc[random]",
      seed       = 14
    ),
    "columns not available"
  )

  unsupported_context <- BayesTools:::.prior_density_context(
    prior_list   = list(alt = prior("normal", list(mean = 0, sd = 1))),
    column_names = c("alt", "rand"),
    n_grid       = 128
  )
  posterior <- .bt_meta_set(posterior, "prior_context", unsupported_context)

  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      hypothesis = "mu_alloc[alternate] > mu_alloc[random]",
      seed       = 15
    ),
    "No prior distribution"
  )

  posterior[["random"]] <- .bt_meta_set(posterior[["random"]], "linear_weights", c(alt = 0, rand = 0))
  out <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = "mu_alloc[alternate] > mu_alloc[random]",
    seed       = 16
  )
  expect_true(is.finite(attr(out, "raw_BF")))

  posterior[["random"]] <- .bt_meta_set(posterior[["random"]], "linear_weights", c(alt = 1, rand = 0))
  testthat::local_mocked_bindings(
    .generate_transformed_prior_samples = function(..., n_samples) {
      matrix(numeric(), nrow = n_samples, ncol = 0L)
    }
  )
  expect_error(
    hypothesis_BF(posterior, hypothesis = "mu_alloc[alternate] > mu_alloc[random]",
                  seed = 16),
    "Linear prior weights reference columns not available in the joint prior context: alt.",
    fixed = TRUE
  )
})


test_that("hypothesis_BF rejects incompatible point-vs-region expressions", {

  set.seed(5)
  prior     <- stats::rnorm(1000)
  posterior <- stats::rnorm(1000)

  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      prior      = prior,
      hypothesis = "theta = 1 vs theta^2 > 0 & theta^2 < 4",
      parameter  = "theta"
    ),
    "same scalar expression"
  )
})


test_that("hypothesis_BF rejects unsafe or ambiguous expressions", {

  prior     <- data.frame(theta = c(-1, 0, 1))
  posterior <- data.frame(theta = c(-1, 0, 1))

  expect_equal(
    BayesTools:::.hypothesis_eval_expression("qlogis(plogis(theta))", posterior),
    posterior[["theta"]],
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.hypothesis_eval_condition("plogis(theta) >= 0.5", posterior),
    c(FALSE, TRUE, TRUE)
  )

  expect_error(
    hypothesis_BF(posterior, prior, "theta <- 1"),
    "Use '=' or '=='"
  )
  expect_error(
    hypothesis_BF(posterior, prior, "theta<-1"),
    "Use '=' or '=='"
  )
  valid_negative <- hypothesis_BF(
    posterior  = data.frame(theta = c(-2, 0, 2)),
    prior      = data.frame(theta = c(-2, 0, 2)),
    hypothesis = "theta < -1"
  )
  expect_true(is.finite(attr(valid_negative, "raw_BF")))

  expect_error(
    hypothesis_BF(posterior, prior, "theta == 0 & theta > -1"),
    "Equality constraints"
  )
  expect_error(
    hypothesis_BF(posterior, prior, "system(theta) > 0"),
    "Unsupported hypothesis expression"
  )
  expect_error(
    hypothesis_BF(posterior, prior, "missing > 0"),
    "unknown quantity",
    class = "BayesTools_parameter_not_found"
  )
  # unknown quantities are unresolved references, as in the catalog
  for(hypothesis in c("missing > 0", "missing = 0", "theta + missing = 0")){
    condition <- tryCatch(
      hypothesis_BF(posterior, prior, hypothesis),
      error = identity
    )
    expect_identical(
      class(condition),
      c("BayesTools_parameter_not_found",
        "BayesTools_parameter_resolution_error", "error", "condition"),
      info = hypothesis
    )
    expect_identical(
      conditionMessage(condition),
      "Hypothesis expression references unknown quantity 'missing'.",
      info = hypothesis
    )
    expect_identical(condition$alias, "missing", info = hypothesis)
    expect_identical(condition$available, names(posterior), info = hypothesis)
  }
})


test_that("hypothesis_BF rejects unknown marginal posterior levels", {

  posterior <- list(alternate = structure(
    rnorm(10),
    class = c("marginal_posterior.simple", "numeric")
  ))
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter") <- "mu_alloc"

  condition <- tryCatch(
    hypothesis_BF(
      posterior  = posterior,
      hypothesis = "mu_alloc[alternate] > mu_alloc[random]"
    ),
    error = identity
  )
  expect_identical(
    class(condition),
    c("BayesTools_parameter_not_found",
      "BayesTools_parameter_resolution_error", "error", "condition")
  )
  expect_identical(
    conditionMessage(condition),
    "Hypothesis references unknown level 'random' for parameter 'mu_alloc'."
  )
  expect_identical(condition$alias, "mu_alloc[random]")
  expect_identical(condition$available, "mu_alloc[alternate]")
})


test_that("hypothesis_BF rejects degenerate point-null KDE inputs", {

  expect_error(
    hypothesis_BF(
      posterior  = rep(0, 10),
      prior      = stats::rnorm(100),
      hypothesis = "theta = 0",
      parameter  = "theta"
    ),
    "degenerate samples"
  )
})


test_that("hypothesis_BF rejects precomputed density requests on raw draws", {

  set.seed(6)
  posterior <- stats::rnorm(100)
  prior <- stats::rnorm(100)

  expect_error(
    hypothesis_BF(
      posterior      = posterior,
      prior          = prior,
      hypothesis     = "theta = 0",
      parameter      = "theta",
      density_method = "precomputed"
    ),
    "requires valid posterior density metadata"
  )

  expect_error(
    hypothesis_BF(
      posterior      = data.frame(theta = posterior),
      prior          = data.frame(theta = prior),
      hypothesis     = "theta + 0 = 0",
      density_method = "precomputed"
    ),
    "raw draws or compound expressions"
  )
})


test_that("hypothesis_BF reuses stored IWMDE ordinates and BF error", {

  prior_density <- .hypothesis_prior_density_for_test()
  posterior <- .hypothesis_marginal_posterior_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", .posterior_ordinate_for_test(
    value       = c(0, 0.5),
    ordinate    = c(0.50, 0.25),
    method      = "IWMDE",
    diagnostics = list(relative_mcse = c(0.03, 0.07))
  ))

  out <- hypothesis_BF(
    posterior      = posterior,
    hypothesis     = "theta = 0.5",
    parameter      = "theta",
    columns        = "all",
    density_method = "precomputed"
  )

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0.5) / 0.25

  expect_equal(attr(out, "raw_BF"), expected, tolerance = 1e-12)
  expect_equal(out[["posterior"]], 0.25, tolerance = 1e-12)
  expect_equal(as.numeric(out[["BF_error"]]), 7, tolerance = 1e-12)
  expect_equal(out[["method"]], "Savage-Dickey (precomputed)")
})


test_that("hypothesis_BF propagates point and region error for point-vs-region tests", {

  prior_density <- .hypothesis_prior_density_for_test()
  samples       <- seq(-3, 3, length.out = 301)
  posterior     <- .hypothesis_marginal_posterior_for_test(samples, prior_density)
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", .posterior_ordinate_for_test(
    value       = 0,
    ordinate    = 0.50,
    method      = "IWMDE",
    diagnostics = list(relative_mcse = 0.03)
  ))

  out <- hypothesis_BF(
    posterior      = posterior,
    hypothesis     = "theta = 0 vs theta > 0",
    parameter      = "theta",
    columns        = "all",
    density_method = "precomputed"
  )

  expected_region_error <- 100 *
    sqrt(.hypothesis_log_prob_var_for_test(samples > 0))
  expected_error <- 100 * sqrt(
    (3 / 100)^2 + (expected_region_error / 100)^2
  )

  expect_equal(as.numeric(out[["BF_error"]]), expected_error, tolerance = 1e-12)
  expect_equal(out[["method"]], "transitive Savage-Dickey")
})


test_that("hypothesis_BF accepts marginal posterior subclasses without base class", {

  prior_density <- .hypothesis_prior_density_for_test()
  posterior <- .hypothesis_marginal_posterior_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  class(posterior) <- setdiff(class(posterior), "marginal_posterior")
  posterior <- .bt_meta_set(posterior, "posterior_ordinate", .posterior_ordinate_for_test(
    value       = 0,
    ordinate    = 0.50,
    method      = "IWMDE",
    diagnostics = list(relative_mcse = 0.02)
  ))

  out <- hypothesis_BF(
    posterior      = posterior,
    hypothesis     = "theta = 0",
    parameter      = "theta",
    density_method = "precomputed"
  )

  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0) / 0.50

  expect_equal(attr(out, "raw_BF"), expected, tolerance = 1e-12)
  expect_equal(as.numeric(out[["BF_error"]]), 2, tolerance = 1e-12)
})


test_that("hypothesis_BF reuses stored qCMDE density and BF error", {

  prior_density <- .hypothesis_prior_density_for_test()
  posterior <- .hypothesis_marginal_posterior_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  stored_x <- seq(-4, 4, length.out = 401)
  stored_y <- stats::dnorm(stored_x, mean = 0.25, sd = 1.1)
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x           = stored_x,
    y           = stored_y,
    method      = "qCMDE",
    diagnostics = list(
      bf_value          = 0.25,
      bf_relative_mcse = 0.04
    )
  ))

  out <- hypothesis_BF(
    posterior      = posterior,
    hypothesis     = "theta = 0.25",
    parameter      = "theta",
    columns        = "all",
    density_method = "precomputed"
  )

  posterior_height <- stats::approx(stored_x, stored_y, xout = 0.25)[["y"]]
  expected <- BayesTools:::.prior_linear_density_height(prior_density, 0.25) /
    posterior_height

  expect_equal(attr(out, "raw_BF"), expected, tolerance = 1e-12)
  expect_equal(out[["posterior"]], posterior_height, tolerance = 1e-12)
  expect_equal(as.numeric(out[["BF_error"]]), 4, tolerance = 1e-12)
})

test_that("precomputed posterior densities carry no point masses", {

  prior_density <- .hypothesis_prior_density_for_test()
  posterior <- .hypothesis_marginal_posterior_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  expect_error(
    .posterior_density_for_test(
      x            = seq(-1, 1, length.out = 101),
      y            = rep(1, 101),
      method       = "point-mass",
      point_masses = data.frame(x = 0, mass = .2)
    ),
    "Posterior densities do not carry 'point_masses'",
    fixed = TRUE
  )

  expect_error(
    .bt_meta_set(posterior, "posterior_density", list(
      x            = seq(-1, 1, length.out = 101),
      y            = rep(1, 101),
      method       = "raw-list",
      point_masses = list(x = 0, mass = .2)
    )),
    "Posterior density metadata must be created with 'posterior_density_attribute()'.",
    fixed = TRUE
  )
})


test_that("hypothesis_BF rejects precomputed point density missing the null", {

  prior_density <- .hypothesis_prior_density_for_test()
  posterior <- .hypothesis_marginal_posterior_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = seq(2, 3, length.out = 101),
    y      = rep(1, 101),
    method = "qCMDE"
  ))

  expect_error(
    hypothesis_BF(
      posterior      = posterior,
      hypothesis     = "theta = 0",
      parameter      = "theta",
      columns        = "all",
      density_method = "precomputed"
    ),
    "Stored posterior density does not span",
    fixed = TRUE
  )
})


test_that("hypothesis_BF rejects zero prior density at point null", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("beta", list(alpha = 2, beta = 2))),
    weights    = c(theta = 1),
    n_grid     = 1024
  )
  posterior <- .hypothesis_marginal_posterior_for_test(
    seq(.001, .999, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "support", BayesTools:::.posterior_support_new(c(0, 1), source = "test"))

  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      hypothesis = "theta = 0",
      parameter  = "theta"
    ),
    "Prior density at the null hypothesis value is zero or non-finite|Prior density at point hypothesis"
  )
})


test_that("hypothesis_BF precomputed point null rejects density grid missing support boundary", {

  prior_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(theta = prior("beta", list(alpha = 1, beta = 1))),
    weights    = c(theta = 1),
    n_grid     = 1024
  )
  posterior <- .hypothesis_marginal_posterior_for_test(
    seq(.001, .999, length.out = 301),
    prior_density
  )
  posterior <- .bt_meta_set(posterior, "support", BayesTools:::.posterior_support_new(c(0, 1), source = "test"))
  posterior <- .bt_meta_set(posterior, "posterior_density", .posterior_density_for_test(
    x      = seq(.25, .75, length.out = 101),
    y      = rep(1, 101),
    method = "qCMDE"
  ))

  expect_error(
    hypothesis_BF(
      posterior      = posterior,
      hypothesis     = "theta = 0",
      parameter      = "theta",
      columns        = "all",
      density_method = "precomputed"
    ),
    "Stored posterior density does not span",
    fixed = TRUE
  )
})


test_that("point hypotheses need an exact regular prior ordinate on every route", {

  # Non-regular or inexact prior ordinates stop with classed conditions
  # (each also of class BayesTools_hypothesis_ordinate).
  set.seed(21)
  draws <- stats::rnorm(4000, .3, .5)
  expect_ordinate_class <- function(expr, class){
    condition <- tryCatch(expr, error = function(e) e)
    expect_s3_class(condition, class)
    expect_s3_class(condition, "BayesTools_hypothesis_ordinate")
  }

  # prior objects of numeric draws: the atom of a spike-and-slab prior, the
  # infinite gamma(1/2) density at 0 and the zero half-normal density below 0
  spike <- prior_spike_and_slab(prior("normal", list(0, 1)),
                                prior_inclusion = prior("spike", list(.5)))
  expect_ordinate_class(hypothesis_BF(draws, spike, "theta = 0"),
                        "BayesTools_point_mass_at_null")
  expect_ordinate_class(hypothesis_BF(abs(draws), prior("gamma", list(.5, 1)), "theta = 0"),
                        "BayesTools_infinite_ordinate")
  expect_ordinate_class(hypothesis_BF(draws, prior("normal", list(0, 1), list(0, Inf)), "theta = -0.5"),
                        "BayesTools_zero_ordinate")

  # prior densities of marginal posteriors: three t terms have no structural
  # ordinate (a grid) and a density grid without provenance none at all
  three_t <- BayesTools:::.prior_linear_combination_density(
    list(a = prior("t", list(0, 1, 3)), b = prior("t", list(0, 1, 3)), c = prior("t", list(0, 1, 3))),
    c(a = 1, b = 1, c = 1)
  )
  expect_ordinate_class(
    hypothesis_BF(.hypothesis_marginal_posterior_for_test(draws, three_t), hypothesis = "theta = 0",
                  parameter = "theta"),
    "BayesTools_inexact_ordinate"
  )
  grid <- structure(list(density = list(x = c(-3, 0, 3), y = c(0, 1 / 3, 0), mass = 1),
                         points = data.frame(x = numeric(), p = numeric())),
                    class = c("prior_linear_density", "prior_density"))
  expect_ordinate_class(
    hypothesis_BF(.hypothesis_marginal_posterior_for_test(draws, grid), hypothesis = "theta = 0",
                  parameter = "theta"),
    "BayesTools_inexact_ordinate"
  )
  # A positive power retains the zero atom and has an exact continuous ordinate at 1.
  nonnegative_spike <- prior_spike_and_slab(prior("normal", list(0, 1), list(0, Inf)),
                                            prior_inclusion = prior("spike", list(.5)))
  undefined <- BayesTools:::.prior_linear_combination_density(
    list(x = nonnegative_spike), c(x = 1), output_transformation = "exp_lin",
    output_transformation_arguments = list(a = 0, b = 2)
  )
  transformed_BF <- hypothesis_BF(
    .hypothesis_marginal_posterior_for_test(abs(draws), undefined), hypothesis = "theta = 1",
    parameter = "theta", columns = "all"
  )
  expect_true(is.finite(transformed_BF$BF) && transformed_BF$BF > 0)
  expect_equal(as.numeric(transformed_BF$prior), stats::dnorm(1) / 2, tolerance = 1e-12)
  expect_ordinate_class(
    Savage_Dickey_BF(.hypothesis_marginal_posterior_for_test(draws, BayesTools:::.prior_linear_combination_density(
      list(x = spike), c(x = 1))), null_hypothesis = 0, silent = TRUE),
    "BayesTools_point_mass_at_null"
  )

  # an affine expression of a prior object has its exact density (2 theta ~
  # N(0, 2) at 1, i.e. dnorm(1 / 2) / 2); nonlinear expressions are refused
  affine <- hypothesis_BF(draws, prior("normal", list(0, 1)), "2 * theta = 1", columns = "all")
  expect_identical(affine$method, "kernel Savage-Dickey")
  expect_equal(as.numeric(affine$prior), stats::dnorm(.5) / 2, tolerance = 1e-14)
  expect_ordinate_class(hypothesis_BF(draws, prior("normal", list(0, 1)), "exp(theta) = 1"),
                        "BayesTools_inexact_ordinate")

  # user-supplied prior draws keep the kernel estimate of the prior ordinate,
  # with one classed warning per call (draw-only inputs have no structural
  # prior density)
  prior_draws <- stats::rnorm(4000)
  raw_warning <- tryCatch(
    hypothesis_BF(draws, prior_draws, c("theta = 0", "theta = 0.5"), columns = "all"),
    warning = function(w) w
  )
  expect_s3_class(raw_warning, "BayesTools_inexact_ordinate")
  expect_s3_class(raw_warning, "BayesTools_hypothesis_ordinate")
  expect_match(conditionMessage(raw_warning), "'theta = 0', 'theta = 0.5'", fixed = TRUE)
  expect_warning(raw <- hypothesis_BF(draws, prior_draws, "theta = 0", columns = "all"),
                 class = "BayesTools_inexact_ordinate")
  expect_identical(raw$method, "kernel Savage-Dickey")
  expect_equal(
    as.numeric(raw$prior),
    as.numeric(BayesTools:::.hypothesis_sample_density_height(prior_draws, 0, "prior")),
    tolerance = 1e-14
  )
})


test_that("linear level expressions use the joint prior density for normal ordinates", {

  # treatment levels B and C have iid N(0, 1) coefficients, so B - C ~ N(0,
  # sqrt(2)) and B - A = B (A is the structural reference level 0); the
  # normal method approximates only the posterior ordinate
  formula <- JAGS_formula(
    ~ fac, "mu", data.frame(fac = factor(c("A", "B", "C"))),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      fac = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  set.seed(3)
  fit <- coda::mcmc(cbind(mu_intercept = stats::rnorm(2000, .1, .3),
                          "mu_fac[1]" = stats::rnorm(2000, .4, .3),
                          "mu_fac[2]" = stats::rnorm(2000, -.2, .3)))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula$prior_list
  fit <- attach_test_parameter_map(fit)
  posterior <- marginal_posterior(as_mixed_posteriors(fit, "mu_fac"), "mu_fac",
                                  use_formula = FALSE, prior_samples = TRUE)
  draws <- fit[, "mu_fac[1]"] - fit[, "mu_fac[2]"]
  for(method in c("KDE", "normal")){
    contrast <- hypothesis_BF(posterior, hypothesis = "mu_fac[B] - mu_fac[C] = 0",
                              density_method = method, columns = "all")
    expect_equal(as.numeric(contrast$prior), stats::dnorm(0, 0, sqrt(2)), tolerance = 1e-14)
    reference <- hypothesis_BF(posterior, hypothesis = "mu_fac[B] - mu_fac[A] = 0.1",
                               density_method = method, columns = "all")
    expect_equal(as.numeric(reference$prior), stats::dnorm(.1), tolerance = 1e-14)
  }
  expect_equal(as.numeric(contrast$posterior),
               stats::dnorm(0, mean(draws), stats::sd(draws)), tolerance = 1e-12)
  expect_identical(contrast$method, "Savage-Dickey (normal)")

  # precomputed ordinates do not exist for expression draws
  expect_error(
    hypothesis_BF(posterior, hypothesis = "mu_fac[B] - mu_fac[C] = 0",
                  density_method = "precomputed"),
    "raw draws or compound expressions"
  )
})


test_that("rejected prior-ordinate quadratures stop with the inexact class", {

  # N(0, 1e-3) * Cauchy(0, 1) is a pure scale mixture whose quadrature away
  # from its offset is rejected by its diagnostics ('probably divergent')
  priors <- list(beta = prior("normal", list(0, 1e-3)), sigma = prior("cauchy", list(0, 1)))
  attr(priors$beta, "multiply_by") <- "sigma"
  density <- .prior_linear_combination_density(priors, c(beta = 1))
  ordinate <- prior_density_ordinate(density, .3)
  expect_identical(ordinate$behavior, "regular")
  expect_false(ordinate$exact)
  expect_true(is.na(ordinate$log_density))
  expect_false(ordinate$provenance$integration$converged)
  posterior <- .hypothesis_marginal_posterior_for_test(stats::rnorm(1000, .01, .05), density)
  condition <- tryCatch(
    hypothesis_BF(posterior, hypothesis = "theta = 0.3", parameter = "theta"),
    error = function(e) e
  )
  expect_s3_class(condition, "BayesTools_inexact_ordinate")
  expect_s3_class(condition, "BayesTools_hypothesis_ordinate")
  expect_match(conditionMessage(condition), "rejected by its diagnostics", fixed = TRUE)
})

test_that("N28 symbolic and centered point-versus-region targets agree", {

  posterior <- data.frame(theta = seq(-1, 2, length.out = 201),
                           phi = seq(-2, 1, length.out = 201)^3)
  prior <- data.frame(theta = seq(-2, 2, length.out = 201),
                       phi = rev(seq(-3, 3, length.out = 201)))
  run <- function(statement){
    suppressWarnings(hypothesis_BF(posterior, prior, statement, columns = "all"))
  }
  symbolic <- tryCatch(run("theta = phi vs theta > phi"), error = identity)
  expect_s3_class(symbolic, "BayesTools_hypothesis_BF")
  centered <- run("theta - phi = 0 vs theta - phi > 0")
  if(!inherits(symbolic, "error")){
    expect_equal(as.numeric(symbolic$BF), as.numeric(centered$BF), tolerance = 1e-14)
    expect_identical(symbolic$method, centered$method)
    expect_equal(symbolic$BF_error, centered$BF_error, tolerance = 1e-14)
  }
  expect_error(run("theta = phi vs phi < theta"))
  expect_error(run("theta = phi vs theta > 0 & phi > 0"))
})

test_that("N33 hypothesis exports retain precision, warning rows and selected BF bounds", {

  result <- hypothesis_BF(seq(-1, 1, length.out = 201), prior("normal", list(0, 1)),
    c("theta > 0", "theta < 0"), parameter = "theta", columns = "all")
  expected_names <- rownames(result)
  attr(result, "warnings") <- setNames(c("first", "second", "global", "other"),
                                       c(expected_names[1], expected_names[1], "", "other"))
  raw <- c(0.98039215686274506, 2.001)
  attr(raw, "bound_operator") <- c(">", "<")
  result$BF <- format_BF(raw)
  result[["component"]] <- c("visible", "visible")
  for(logBF in c(FALSE, TRUE)) for(BF01 in c(FALSE, TRUE)){
    table <- result
    table$BF <- format_BF(raw, logBF = logBF, BF01 = BF01)
    exported <- as.data.frame(table)
    expect_identical(class(exported), "data.frame")
    expect_identical(exported, data.frame(table))
    expect_identical(names(exported)[1:3], c("component", "parameter", "warning"))
    expect_identical(exported$parameter[1:2], expected_names)
    expect_identical(exported$warning, c("first second", NA_character_, "global", "other"))
    name <- paste0(if(logBF) "log" else "", if(BF01) "BF01" else "BF10")
    expect_identical(exported[[name]][1:2], as.numeric(table$BF))
    expect_identical(exported$BF_bound_operator[1:2], attr(table$BF, "bound_operator"))
    selected <- table[2:1, , drop = FALSE]
    expect_identical(attr(selected$BF, "bound_operator"), rev(attr(table$BF, "bound_operator")))
    expect_identical(as.data.frame(selected)$BF_bound_operator[1:2], rev(attr(table$BF, "bound_operator")))
    updated <- update(selected, logBF = logBF, BF01 = BF01)
    expect_identical(as.data.frame(updated)$BF_bound_operator[1:2], attr(updated$BF, "bound_operator"))
    expect_false(any(grepl("hypothesis/warnings|component.1|BF_bound_operator", capture.output(print(table)))))
  }
  attr(result$BF, "bound_operator") <- c(NA_character_, NA_character_)
  expect_false("BF_bound_operator" %in% names(as.data.frame(result)))
  renamed <- result
  rownames(renamed) <- c("consumer first", "consumer second")
  attr(renamed, "rownames") <- FALSE
  expect_identical(as.data.frame(renamed)$parameter[1:2], rownames(renamed))
  removed <- result[, setdiff(names(result), "BF"), drop = FALSE]
  expect_false("BF10" %in% names(as.data.frame(removed)))
  collision <- result
  collision$BF <- format_BF(raw)
  collision$BF10 <- c(17, 18)
  collision$BF_bound_operator <- c("user left", "user right")
  collision <- collision[, c("BF10", setdiff(names(collision), "BF10")), drop = FALSE]
  exported_collision <- as.data.frame(collision)
  expect_identical(exported_collision$BF10[1:2], as.numeric(collision$BF))
  expect_identical(exported_collision$BF10.1[1:2], c(17, 18))
  expect_identical(exported_collision$BF_bound_operator[1:2], c(">", "<"))
  expect_identical(exported_collision$BF_bound_operator.1[1:2], c("user left", "user right"))
  expect_error(as.data.frame(result, row.names = "one"))
})

test_that("N34 a genuine compiled combination prints its persisted quantity label", {

  formula <- JAGS_formula(~ fac, "mu", data.frame(fac = factor(c("A", "B", "C"))),
    prior_list = list(intercept = prior("normal", list(0, 1)),
      fac = prior_factor("normal", list(0, 1), contrast = "treatment")))
  fit <- coda::mcmc(cbind(mu_intercept = seq(-1, 1, length.out = 201),
    "mu_fac[1]" = seq(-1, 2, length.out = 201), "mu_fac[2]" = seq(-2, 1, length.out = 201)^3))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula$prior_list
  fit <- attach_test_parameter_map(fit)
  posterior <- marginal_posterior(as_mixed_posteriors(fit, "mu_fac"), "mu_fac",
                                  use_formula = FALSE, prior_samples = TRUE)
  target <- hypothesis_linear_target(posterior, "mu_fac[B] - mu_fac[C] = 0", "mu_fac")
  expected <- parameter_labels(.bt_meta_get(target$posterior, "quantities")$label_parts, style = "table")
  result <- hypothesis_BF(target$posterior, hypothesis = target$hypothesis, parameter = target$parameter)
  expect_identical(rownames(result), unname(expected))
  expect_false(any(grepl(".BayesTools_linear_target", rownames(result), fixed = TRUE)))
})

test_that("R116 N05 direct ordinate matching filters original parameter and condition leaves", {

  prior_density <- .prior_linear_combination_density(list(theta = prior("normal", list(0, 1))), c(theta = 1), n_grid = 128)
  posterior <- .hypothesis_marginal_posterior_for_test(seq(-2, 2, length.out = 101), prior_density)
  attr(posterior, "parameter") <- "theta"
  leaf <- function(value, parameter = "theta", conditional = NULL, height = .5){
    posterior_ordinate_attribute(value, height, "qCMDE", "precomputed", diagnostics = list(relative_mcse = .1), parameter = parameter, conditional = conditional)
  }
  original <- posterior_ordinate_append(leaf(0), leaf(1, "phi", height = 100))
  posterior_metadata(posterior, "posterior_ordinate") <- original
  status <- .posterior_ordinate_direct_status(posterior, null_hypothesis = 1)
  expect_identical(status[c("present", "relevant", "valid")], list(present = TRUE, relevant = TRUE, valid = FALSE))
  expect_null(status$value)
  expect_null(.posterior_ordinate_direct_attribute(posterior, null_hypothesis = 1))
  expect_error(Savage_Dickey_BF(posterior, 1, silent = TRUE, density_method = "precomputed"),
    "Precomputed posterior ordinate metadata is present but invalid for the requested null hypothesis.", fixed = TRUE)
  selected <- .posterior_ordinate_direct_attribute(posterior, null_hypothesis = 0)
  expect_identical(selected, original$ordinates[[1L]])
  expect_identical(.posterior_ordinate_direct_attribute(posterior), selected)
  expect_equal(as.numeric(Savage_Dickey_BF(posterior, 0, silent = TRUE, density_method = "precomputed")), stats::dnorm(0) / .5, tolerance = 1e-12)
  compatible <- posterior_ordinate_append(leaf(0), leaf(1, height = .8))
  posterior_metadata(posterior, "posterior_ordinate") <- compatible
  matched <- .posterior_ordinate_direct_attribute(posterior, null_hypothesis = 1)
  expect_s3_class(matched, "BayesTools_posterior_ordinates")
  expect_length(matched$ordinates, 2L)
  expect_equal(.posterior_ordinate_from_attribute(matched, 0)$y, .5)
  expect_equal(.posterior_ordinate_from_attribute(matched, 1)$y, .8)
  expect_equal(.posterior_ordinate_from_attribute(matched, 1)$diagnostics$relative_mcse, .1)
  conditioned <- .bt_meta_update(posterior, condition = list(conditional = "theta", conditional_rule = "AND"))
  posterior_metadata(conditioned, "posterior_ordinate") <- posterior_ordinate_append(leaf(0, conditional = "theta"), leaf(1, conditional = "phi"))
  expect_null(.posterior_ordinate_direct_attribute(conditioned, null_hypothesis = 1))
  expect_true(.posterior_ordinate_direct_status(conditioned, null_hypothesis = 1)$relevant)
  expect_equal(.posterior_ordinate_direct_status(conditioned, null_hypothesis = 0)$value$y, .5)
  posterior_metadata(posterior, "posterior_ordinate") <- leaf(0, "phi")
  expect_identical(.posterior_ordinate_direct_status(posterior, null_hypothesis = 0)[c("present", "relevant", "valid")], list(present = TRUE, relevant = FALSE, valid = FALSE))
  posterior_metadata(posterior, "posterior_ordinate") <- NULL
  expect_identical(.posterior_ordinate_direct_status(posterior, null_hypothesis = 0)[c("present", "relevant", "valid")], list(present = FALSE, relevant = FALSE, valid = FALSE))
})

test_that("R116 N05 source collection retains only compatible raw ordinate leaves", {

  theta <- posterior_ordinate_attribute(0, .5, "qCMDE", "precomputed", parameter = "theta")
  phi <- posterior_ordinate_attribute(1, 100, "qCMDE", "precomputed", parameter = "phi")
  aggregate <- posterior_ordinate_append(theta, phi)
  aggregate$parameter <- "theta"
  aggregate$conditional <- "irrelevant-parent"
  matched <- .posterior_ordinate_from_sources(list(aggregate), "theta")
  expect_identical(matched, theta)
  expect_null(.posterior_ordinate_from_sources(list(aggregate), "theta", null_hypothesis = 1))
  expect_identical(.posterior_ordinate_from_sources(list(aggregate), "theta", null_hypothesis = 0), theta)
  unlabeled <- posterior_ordinate_attribute(0, .6, "qCMDE", "precomputed")
  expect_identical(.posterior_ordinate_from_sources(list(list(theta = unlabeled)), "theta"), unlabeled)
  expect_null(.posterior_ordinate_from_sources(list(list(theta = phi)), "theta"))
  expect_identical(.posterior_ordinate_from_sources(list(list(theta = aggregate)), "theta"), theta)
  invalid <- aggregate
  invalid$ordinates[[2L]]$ordinate <- -1
  expect_error(.posterior_ordinate_from_sources(list(invalid), "theta"), "Posterior ordinate metadata is invalid", fixed = TRUE)
  expect_error(posterior_ordinate_append(theta, theta), "cannot contain duplicate values", fixed = TRUE)
})
