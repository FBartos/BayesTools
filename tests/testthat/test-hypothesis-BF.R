skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Hypothesis Bayes Factors
# ============================================================================ #

.hypothesis_marginal_posterior_for_test <- function(samples, prior_density){

  class(samples) <- c("marginal_posterior.simple", "marginal_posterior", class(samples))
  attr(samples, "prior_density") <- prior_density

  samples
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


test_that("hypothesis_BF computes point-null Savage-Dickey from numeric draws", {

  set.seed(1)
  prior     <- stats::rnorm(12000, mean = 0, sd = 1)
  posterior <- stats::rnorm(12000, mean = 0.4, sd = 1.2)

  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0",
    parameter  = "theta"
  )

  expected <- stats::dnorm(0, mean = 0, sd = 1) /
    stats::dnorm(0, mean = 0.4, sd = 1.2)

  expect_equal(attr(out, "raw_BF"), expected, tolerance = 0.08)

  equivalent <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta != 0",
    parameter  = "theta"
  )

  expect_equal(attr(equivalent, "raw_BF"), attr(out, "raw_BF"),
               tolerance = 1e-12)

  explicit_inverse <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0 vs theta != 0",
    parameter  = "theta"
  )

  expect_equal(attr(explicit_inverse, "raw_BF"), 1 / attr(out, "raw_BF"),
               tolerance = 1e-12)
})


test_that("hypothesis_BF returns compact BayesTools table by default", {

  set.seed(101)
  prior     <- stats::rnorm(4000)
  posterior <- stats::rnorm(4000, mean = 0.2)

  out <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0",
    parameter  = "theta"
  )

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

  detailed <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0",
    parameter  = "theta",
    columns    = "all"
  )

  expect_equal(
    colnames(detailed),
    c("Alternative", "Null", "BF", "BF_error", "prior", "posterior", "method")
  )
  expect_match(attr(detailed, "footnotes"), "diagnostic values", fixed = TRUE)
  expect_match(attr(detailed, "footnotes"), "not always probabilities",
               fixed = TRUE)

  printed_detailed <- utils::capture.output(print(detailed))
  expect_true(any(grepl("diagnostic values", printed_detailed, fixed = TRUE)))

  multi <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = c("theta = 0", "theta > 0 vs theta < 0"),
    parameter  = "theta"
  )
  expect_equal(attr(multi[2, , drop = FALSE], "raw_BF"),
               attr(multi, "raw_BF")[2])
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

  point <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0",
    parameter  = "theta"
  )
  transitive <- hypothesis_BF(
    posterior  = posterior,
    prior      = prior,
    hypothesis = "theta = 0 vs theta > 0",
    parameter  = "theta"
  )

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
    alternate = structure(
      c(rep(1, 80), rep(-1, 20)),
      class          = c("marginal_posterior.simple", "numeric"),
      linear_weights = matrix(c(1, 0), nrow = 1,
                              dimnames = list(NULL, c("alt", "rand")))
    ),
    random = structure(
      rep(0, 100),
      class          = c("marginal_posterior.simple", "numeric"),
      linear_weights = c(alt = 0, rand = 1)
    )
  )
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter")             <- "mu_alloc"
  attr(posterior, "prior_density_context") <- context

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
    alternate = structure(
      c(rep(1, 80), rep(-1, 20)),
      class                 = c("marginal_posterior.simple", "numeric"),
      linear_weights        = c(alt = 1, rand = 0),
      effective_conditional = "mu_intercept"
    ),
    random = structure(
      rep(0, 100),
      class                 = c("marginal_posterior.simple", "numeric"),
      linear_weights        = c(alt = 0, rand = 1),
      effective_conditional = c("mu_intercept", "mu_alloc")
    )
  )
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter")             <- "mu_alloc"
  attr(posterior, "prior_density_context") <- context

  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      hypothesis = "mu_alloc[alternate] > mu_alloc[random]"
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
  attr(alternate, "posterior_ordinate") <- list(
    value       = 0,
    ordinate    = 0.50,
    method      = "IWMDE",
    diagnostics = list(relative_mcse = 0.03)
  )
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
    alternate = structure(
      c(rep(1, 75), rep(-1, 25)),
      class          = c("marginal_posterior.simple", "numeric"),
      linear_weights = c(alt = 1, rand = 0)
    ),
    random = structure(
      rep(0, 100),
      class          = c("marginal_posterior.simple", "numeric"),
      linear_weights = c(alt = 0, rand = 1)
    )
  )
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter")             <- "mu_alloc"
  attr(posterior, "prior_density_context") <- context
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
      alternate = structure(
        c(rep(1, 75), rep(-1, 25)),
        class          = c("marginal_posterior.simple", "numeric"),
        linear_weights = c(alt = 1, rand = 0)
      ),
      random = structure(
        rep(0, 100),
        class          = c("marginal_posterior.simple", "numeric"),
        linear_weights = c(alt = 0, rand = 1)
      )
    )
    class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
    attr(posterior, "parameter")             <- "mu_alloc"
    attr(posterior, "prior_density_context") <- context

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
    alternate = structure(
      c(rep(1, 75), rep(-1, 25)),
      class          = c("marginal_posterior.simple", "numeric"),
      linear_weights = c(alt = 1)
    ),
    random = structure(
      rep(0, 100),
      class          = c("marginal_posterior.simple", "numeric"),
      linear_weights = c(rand = 1)
    )
  )
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter")             <- "mu_alloc"
  attr(posterior, "prior_density_context") <- context

  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      hypothesis = "mu_alloc[alternate] > mu_alloc[random]",
      seed       = 14
    ),
    "columns not available"
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
    "unknown quantity"
  )
})


test_that("hypothesis_BF rejects unknown marginal posterior levels", {

  posterior <- list(alternate = structure(
    rnorm(10),
    class = c("marginal_posterior.simple", "numeric")
  ))
  class(posterior) <- c("list", "marginal_posterior.factor", "marginal_posterior")
  attr(posterior, "parameter") <- "mu_alloc"

  expect_error(
    hypothesis_BF(
      posterior  = posterior,
      hypothesis = "mu_alloc[alternate] > mu_alloc[random]"
    ),
    "unknown level"
  )
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


test_that("hypothesis_BF reuses stored IWMDE ordinates and BF error", {

  prior_density <- .hypothesis_prior_density_for_test()
  posterior <- .hypothesis_marginal_posterior_for_test(
    seq(-3, 3, length.out = 301),
    prior_density
  )
  attr(posterior, "posterior_ordinate") <- list(
    value       = c(0, 0.5),
    ordinate    = c(0.50, 0.25),
    method      = "IWMDE",
    diagnostics = list(relative_mcse = c(0.03, 0.07))
  )

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
  attr(posterior, "posterior_ordinate") <- list(
    value       = 0,
    ordinate    = 0.50,
    method      = "IWMDE",
    diagnostics = list(relative_mcse = 0.03)
  )

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
  attr(posterior, "posterior_ordinate") <- list(
    value       = 0,
    ordinate    = 0.50,
    method      = "IWMDE",
    diagnostics = list(relative_mcse = 0.02)
  )

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
  attr(posterior, "posterior_density") <- list(
    x           = stored_x,
    y           = stored_y,
    method      = "qCMDE",
    diagnostics = list(
      bf_value          = 0.25,
      bf_relative_mcse = 0.04
    )
  )

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
