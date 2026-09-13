skip_if_not_test_profile(c("unit", "fixture"))

# TEST FILE: Weightfunction prior redesign
# ============================================================================ #
#
# PURPOSE:
#   Contract tests for the unified prior_weightfunction() API, independent
#   weight priors, one/two-sided mapping, JAGS code generation, and bridge
#   sampling parameter extraction.
#
# TAGS: @priors, @jags, @weightfunctions
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

test_that("selection models declare three fixed source choices without changing weight priors", {

  default <- selection_model()
  expect_identical(unclass(default), list(
    estimate_random_effects = "integrate", other_random_effects = "condition",
    known_sampling_variance = "integrate",
    weight_rule = "product", group = NULL
  ))
  expect_identical(selection_model_spec(prior_weightfunction()), default)

  weights <- list(
    wf_cumulative(c(2, 3)), wf_fixed(c(1, 0)),
    wf_independent(prior("gamma", list(2, 1))),
    wf_independent(prior("normal", list(0, 1)), scale = "log_omega")
  )
  for(weight in weights){
    reference <- prior_weightfunction(steps = .05, weights = weight, prior_weights = 3)
    for(estimate in c("integrate", "condition")){
      for(other in c("condition", "integrate")){
        for(sampling in c("condition", "integrate")){
          for(rule in c("product", "best")){
            model <- selection_model(estimate, other, sampling, rule, group = paper_id)
            candidate <- prior_weightfunction(steps = .05, weights = weight,
                                               prior_weights = 3, model = model)
            expect_identical(selection_model_spec(candidate), model)
            expect_identical(candidate[names(candidate) != "model"],
                             reference[names(reference) != "model"])
          }
        }
      }
    }
    set.seed(31)
    reference_draws <- rng(reference, 20)
    set.seed(31)
    expect_identical(rng(candidate, 20), reference_draws)
    expect_identical(JAGS_to_monitor(list(omega = candidate)),
                     JAGS_to_monitor(list(omega = reference)))
    expect_identical(JAGS_get_inits(list(omega = candidate), chains = 1, seed = 31),
                     JAGS_get_inits(list(omega = reference), chains = 1, seed = 31))
    expect_identical(.JAGS_bridgesampling_posterior_info.weightfunction(candidate),
                     .JAGS_bridgesampling_posterior_info.weightfunction(reference))
  }
})

test_that("selection group references survive deferred and wrapper capture", {

  paper_id <- seq_len(3)
  expect_identical(selection_model(group = paper_id)$group, "paper_id")
  expect_identical(selection_model(group = `paper id`)$group, "paper id")
  expect_identical(selection_model(group = "paper id")$group, "paper id")
  expect_null(selection_model(group = NULL)$group)

  wrapper <- function(group = paper_id) selection_model(group = group)
  nested_wrapper <- function(column) wrapper(group = column)
  explicit_wrapper <- function(group) selection_model(group = {{group}})
  forced_wrapper <- function(group) {
    force(group)
    selection_model(group = group)
  }
  expect_identical(wrapper()$group, "paper_id")
  expect_identical(wrapper(paper_id)$group, "paper_id")
  expect_identical(nested_wrapper(`paper id`)$group, "paper id")
  expect_identical(explicit_wrapper(paper_id)$group, "paper_id")
  expect_identical(forced_wrapper("paper_id")$group, "paper_id")
  expect_error(forced_wrapper(1:3), "'group' must be a data-column name", fixed = TRUE)

  column_name <- "paper id"
  stored <- do.call(selection_model, list(group = column_name))
  expect_identical(stored, selection_model(group = `paper id`))
  expect_identical(unserialize(serialize(stored, NULL)), stored)
  expect_identical(attributes(stored), list(
    names = c("estimate_random_effects", "other_random_effects",
              "known_sampling_variance", "weight_rule", "group"),
    class = "selection_model"
  ))
  expect_true(all(vapply(unclass(stored), is.character, logical(1))))

  evaluated <- FALSE
  expect_error(selection_model(group = { evaluated <- TRUE; "paper_id" }),
               "'group' must be a data-column name", fixed = TRUE)
  expect_false(evaluated)
  expect_error(selection_model(group = ""), "'group' must be a data-column name", fixed = TRUE)
  expect_error(selection_model(group = NA_character_), "'group' must be a data-column name", fixed = TRUE)
  expect_error(do.call(selection_model, list(group = c("paper_id", "study"))),
               "'group' must be a data-column name", fixed = TRUE)
  circular <- function(group = other, other = group) selection_model(group = group)
  expect_error(circular(), "'group' contains a circular wrapper argument reference. Forward a data-column name.", fixed = TRUE)
})

test_that("selection model validation rejects malformed and obsolete specifications", {

  expect_error(selection_model(estimate_random_effects = "exact"), "estimate_random_effects", fixed = TRUE)
  expect_error(selection_model(other_random_effects = "approximate"), "other_random_effects", fixed = TRUE)
  expect_error(selection_model(known_sampling_variance = "approximate"), "known_sampling_variance", fixed = TRUE)
  expect_error(selection_model(known_sampling_variance = .5), "known_sampling_variance", fixed = TRUE)
  expect_error(selection_model(random_effects = "condition"), "unused argument", fixed = TRUE)
  expect_error(selection_model(known_covariance = "condition"), "unused argument", fixed = TRUE)
  expect_error(selection_model(weight_rule = "maximum_weight"), "weight_rule", fixed = TRUE)
  expect_error(selection_model(estimate_random_effects = c("condition", "integrate")), "length", fixed = TRUE)
  expect_error(selection_model(known_sampling_variance = NA_character_), "cannot contain NA/NaN", fixed = TRUE)
  expect_identical(withVisible(check_selection_model(selection_model())),
                   list(value = selection_model(), visible = FALSE))
  expect_error(check_selection_model(NULL, name = "selection"),
               "'selection' must be a specification from 'selection_model()'.", fixed = TRUE)
  expect_error(prior_weightfunction(model = NULL),
               "'model' must be a specification from 'selection_model()'.", fixed = TRUE)
  malformed <- prior_weightfunction()
  malformed$model$weight_rule <- "maximum_weight"
  expect_error(selection_model_spec(malformed), "weight_rule", fixed = TRUE)
  expect_error(prior_mixture(list(prior_none(), malformed)), "weight_rule", fixed = TRUE)
  expect_error(prior_bias(selection = malformed), "weight_rule", fixed = TRUE)
  missing_model <- prior_weightfunction()
  missing_model$model <- NULL
  expect_error(selection_model_spec(missing_model),
               "'prior$model' must be a specification from 'selection_model()'.", fixed = TRUE)
})

test_that("bias composition preserves each child selection model and its odds", {

  first <- prior_weightfunction(steps = .05, prior_weights = 2,
    model = selection_model(group = paper_id))
  second <- prior_weightfunction(steps = .05, weights = wf_fixed(c(1, 1.5)),
    prior_weights = 3,
    model = selection_model("condition", "integrate", "integrate", "best", group = study_id))
  combined <- prior_bias(selection = second, phacking = prior_phacking(), prior_weights = 5)
  branches <- list(prior_none(prior_weights = 7), first, combined)
  before <- serialize(branches, NULL)
  mixture <- prior_mixture(branches, is_null = c(TRUE, FALSE, FALSE),
                           components = c("null", "product", "best"))
  expected <- list(NULL, first$model, second$model)
  expect_identical(lapply(branches, selection_model_spec), expected)
  expect_identical(lapply(mixture, selection_model_spec), expected)
  expect_identical(lapply(unserialize(serialize(mixture, NULL)), selection_model_spec), expected)
  expect_identical(attr(mixture, "prior_weights"), c(7, 2, 5))
  expect_identical(attr(mixture, "components"), c("null", "product", "best"))
  expect_identical(vapply(mixture, attr, character(1), which = "component"),
                   c("null", "product", "best"))
  expect_identical(lapply(lapply(branches, .selection_branch_info), function(branch){
    selection_model_spec(branch$selection)
  }), expected)
  expect_error(selection_model_spec(mixture),
    "A selection-model specification is unavailable for a mixture as a whole. Inspect each prior branch with 'selection_model_spec()'.", fixed = TRUE)
  expect_error(prior_mixture(list(first, prior("normal", list(0, 1)))),
               "Publication-bias prior mixtures", fixed = TRUE)
  expect_identical(serialize(branches, NULL), before)
  expect_null(selection_model_spec(prior_bias(phacking = prior_phacking())))
  expect_null(selection_model_spec(prior_PET("normal", list(0, 1))))
})

test_that("selection model printing distinguishes source choices from weight-prior labels", {

  model <- selection_model("integrate", "condition", "integrate", "best", group = `paper id`)
  expect_identical(utils::capture.output(print(model)), c(
    "Selection model:",
    "  Estimate random effects: integrate (average effects before normalization).",
    "  Other random effects: condition (retain unknown effects during normalization).",
    "  Known sampling error: integrate (average full error vector before normalization).",
    "  Weight rule: best (Weight of the best p-value).",
    "  Group: 'paper id' (unresolved data column).",
    "  Sources are resolved when model data are bound."
  ))
  expect_identical(utils::capture.output(print(selection_model(
    known_sampling_variance = "condition"
  ))), c(
    "Selection model:",
    "  Estimate random effects: integrate (average effects before normalization).",
    "  Other random effects: condition (retain unknown effects during normalization).",
    "  Known sampling error: condition (retain full unknown error vector during normalization).",
    "  Weight rule: product (Product of estimate weights).",
    "  Group: automatic (resolved when data are bound).",
    "  Sources are resolved when model data are bound."
  ))
  candidate <- prior_weightfunction(steps = .05, model = model)
  reference <- prior_weightfunction(steps = .05)
  expect_identical(print(candidate, silent = TRUE),
                   "omega[one-sided: .05] ~ CumDirichlet(1, 1)")
  expect_identical(print(candidate, plot = TRUE), print(reference, plot = TRUE))
  expect_identical(utils::capture.output(print(candidate, inline = TRUE)),
                   "omega[one-sided: .05] ~ CumDirichlet(1, 1)")
  combined <- prior_bias(selection = candidate)
  expect_true(any(grepl("Weight of the best p-value", utils::capture.output(print(combined)), fixed = TRUE)))
  mixture <- prior_mixture(list(reference, combined))
  printed <- utils::capture.output(print(mixture))
  expect_true(all(c("Prior branch 1", "Prior branch 2") %in% printed))
  expect_true(any(grepl("Product of estimate weights", printed, fixed = TRUE)))
  expect_true(any(grepl("Weight of the best p-value", printed, fixed = TRUE)))
})

test_that("prior_weightfunction stores canonical geometry and weight priors", {

  wf <- prior_weightfunction(
    side = "one-sided",
    steps = c(.025, .05),
    weights = wf_cumulative(c(1, 2, 3)),
    prior_weights = 2
  )

  expect_true(is.prior.weightfunction(wf))
  expect_equal(wf$distribution, "weightfunction")
  expect_equal(wf$side, "one-sided")
  expect_equal(wf$steps, c(.025, .05))
  expect_equal(wf$bins$lower, c(0, .025, .05))
  expect_equal(wf$bins$upper, c(.025, .05, 1))
  expect_equal(wf$bins$reference, c(TRUE, FALSE, FALSE))
  expect_equal(wf$weights$type, "cumulative")
  expect_equal(wf$weights$alpha, c(1, 2, 3))
  expect_equal(wf$prior_weights, 2)

  samples <- rng(wf, 1000)
  expect_true(all(samples[,1] == 1))
})

test_that("legacy monotone weightfunction helpers use the canonical reference-first orientation", {

  alpha <- c(2, 4)
  q <- .5

  expect_equal(
    mdone.sided(q, alpha = alpha),
    matrix(c(0, stats::dbeta(q, 4, 2)), nrow = 1)
  )
  expect_equal(
    mpone.sided(q, alpha = alpha),
    matrix(c(0, stats::pbeta(q, 4, 2)), nrow = 1)
  )
  expect_equal(
    mqone.sided(q, alpha = alpha),
    matrix(c(1, stats::qbeta(q, 4, 2)), nrow = 1)
  )

  set.seed(11)
  helper_samples <- rone.sided(5000, alpha = alpha)
  expect_true(all(helper_samples[,1] == 1))
  expect_equal(mean(helper_samples[,2]), 4 / 6, tolerance = .02)

  prior <- prior_weightfunction("one-sided", c(.05), wf_cumulative(alpha))
  set.seed(11)
  prior_samples <- rng(prior, 5000)
  expect_equal(unname(helper_samples), unname(prior_samples), tolerance = 1e-12)
})

test_that("weightfunction constructors validate independent scales", {

  expect_silent(wf_independent(prior("beta", list(1, 1))))
  expect_silent(wf_independent(prior("gamma", list(2, 1))))
  expect_error(
    wf_independent(prior("normal", list(0, 1))),
    "non-negative support"
  )

  expect_silent(wf_independent(
    prior("normal", list(0, 1)),
    scale = "log_omega"
  ))
  expect_error(wf_cumulative(c(1, Inf)), "finite")
  expect_error(wf_fixed(c(1, Inf)), "finite")

  expect_error(
    prior_weightfunction("one-sided", c(.05), wf_fixed(c(.9, .5))),
    "reference-bin"
  )
  expect_silent(prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, 0))))

  fixed_above_one <- prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, 1.5)))
  expect_equal(range(fixed_above_one), c(0, 1.5))
  expect_equal(unname(rng(fixed_above_one, 2)), matrix(c(1, 1, 1.5, 1.5), nrow = 2))
  expect_equal(rone.sided_fixed(2, omega = c(1, 1.5)), matrix(c(1, 1, 1.5, 1.5), nrow = 2))
})

test_that("weightfunctions_mapping expands two-sided priors onto one-sided cuts", {

  one_sided <- prior_weightfunction(
    "one-sided",
    c(.025, .05),
    wf_cumulative(c(1, 1, 1))
  )
  two_sided <- prior_weightfunction(
    "two-sided",
    c(.05),
    wf_fixed(c(1, .5))
  )

  expect_equal(
    weightfunctions_mapping(list(one_sided, two_sided), cuts_only = TRUE, one_sided = TRUE),
    c(0, .025, .05, .975, 1)
  )
  expect_equal(
    weightfunctions_mapping(list(one_sided, two_sided), one_sided = TRUE),
    list(c(1L, 2L, 3L, 3L), c(1L, 2L, 2L, 1L))
  )
})

test_that("JAGS generation uses component-local omega for bias mixtures", {

  fixed <- prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .5)))
  independent <- prior_weightfunction(
    "one-sided",
    c(.025, .05),
    wf_independent(prior("beta", list(2, 3)))
  )
  bias <- prior_mixture(list(prior_none(), fixed, independent))

  syntax <- JAGS_add_priors("model{}", list(bias = bias))

  expect_match(syntax, "omega_component_1\\[1\\] <- 1")
  expect_match(syntax, "omega_local_component_2\\[2\\] <- 0.5")
  expect_match(syntax, "omega_local_component_3\\[2\\] ~ dbeta\\(2,3\\)")
  expect_match(syntax, "omega\\[1\\] <- omega_component_1\\[1\\] \\* equals\\(bias_indicator, 1\\)")
  expect_false(grepl("eta2omega", syntax, fixed = TRUE))
})

test_that("binary cumulative weights use their exact beta marginal", {

  cumulative <- prior_weightfunction(
    "one-sided",
    .05,
    wf_cumulative(c(2, 4))
  )
  syntax <- JAGS_add_priors("model{}", list(omega = cumulative))
  expect_match(syntax, "omega_ratio ~ dbeta\\(4, 2\\)")
  expect_match(syntax, "omega\\[2\\] <- omega_ratio")
  expect_false(grepl("eta\\[|std_eta", syntax))

  inits <- JAGS_get_inits(list(omega = cumulative), chains = 2, seed = 1)
  expect_true(all(vapply(inits, function(x){
    is.finite(x$omega_ratio) && x$omega_ratio > 0 && x$omega_ratio < 1
  }, logical(1))))
  expect_equal(JAGS_to_monitor(list(omega = cumulative)), "omega")

  posterior_info <- .JAGS_bridgesampling_posterior_info.weightfunction(cumulative)
  expect_equal(as.vector(posterior_info), "omega[2]")
  expect_equal(attr(posterior_info, "lb"), c("omega[2]" = 0))
  expect_equal(attr(posterior_info, "ub"), c("omega[2]" = 1))

  samples <- c("omega[2]" = .6)
  expect_equal(
    JAGS_marglik_priors(samples, list(omega = cumulative)),
    stats::dbeta(.6, 4, 2, log = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    JAGS_marglik_parameters(samples, list(omega = cumulative))$omega,
    c(1, .6)
  )
  expect_equal(
    JAGS_marglik_priors(c("omega[2]" = 1.1), list(omega = cumulative)),
    -Inf
  )
  expect_error(
    JAGS_marglik_parameters(
      c("omega[2]" = 1.1),
      list(omega = cumulative)
    ),
    "out-of-support binary cumulative weightfunction coordinate"
  )
  expect_error(
    JAGS_marglik_priors(numeric(), list(omega = cumulative)),
    "does not contain the monitored binary cumulative weightfunction parameter"
  )

  eta <- c(.8, 1.2)
  total <- sum(eta)
  omega <- eta[2] / total
  gamma_density_with_jacobian <-
    sum(stats::dgamma(eta, shape = c(2, 4), rate = 1, log = TRUE)) +
    log(total)
  factorized_density <-
    stats::dbeta(omega, 4, 2, log = TRUE) +
    stats::dgamma(total, 6, rate = 1, log = TRUE)
  expect_equal(
    gamma_density_with_jacobian,
    factorized_density,
    tolerance = 1e-12
  )
})

test_that("JAGS bridge helpers use natural latent weight parameters", {

  cumulative <- prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)))
  independent <- prior_weightfunction(
    "one-sided",
    c(.025, .05),
    wf_independent(prior("beta", list(2, 3)))
  )
  independent_gamma <- prior_weightfunction(
    "one-sided",
    c(.05),
    wf_independent(prior("gamma", list(shape = 9, rate = 3)))
  )
  log_independent <- prior_weightfunction(
    "one-sided",
    c(.05),
    wf_independent(
      prior("normal", list(0, 1)),
      "log_omega"
    )
  )
  two_sided <- prior_weightfunction(
    "two-sided",
    c(.05, .10),
    wf_cumulative(c(1, 2, 3))
  )

  expect_equal(as.vector(.JAGS_bridgesampling_posterior_info.weightfunction(cumulative)), paste0("eta[", 1:3, "]"))
  expect_equal(as.vector(.JAGS_bridgesampling_posterior_info.weightfunction(independent)), paste0("omega[", 2:3, "]"))
  expect_equal(as.vector(.JAGS_bridgesampling_posterior_info.weightfunction(independent_gamma)), "omega[2]")
  expect_equal(attr(.JAGS_bridgesampling_posterior_info.weightfunction(independent_gamma), "ub"), c("omega[2]" = Inf))
  expect_equal(as.vector(.JAGS_bridgesampling_posterior_info.weightfunction(log_independent)), "log_omega[2]")
  expect_equal(attr(.JAGS_bridgesampling_posterior_info.weightfunction(log_independent), "ub"), c("log_omega[2]" = Inf))

  samples <- c("eta[1]" = 1, "eta[2]" = 2, "eta[3]" = 3, "omega[2]" = 1.2, "omega[3]" = .1, "log_omega[2]" = .5)
  expect_equal(JAGS_marglik_parameters(samples, list(omega = cumulative))$omega, c(1, 5/6, 1/2))
  expect_equal(JAGS_marglik_parameters(samples, list(omega = independent))$omega, c(1, 1.2, .1))
  expect_equal(JAGS_marglik_parameters(samples, list(omega = independent_gamma))$omega, c(1, 1.2))
  expect_equal(
    JAGS_marglik_priors(samples, list(omega = independent_gamma)),
    mlpdf(independent_gamma$weights$prior, 1.2),
    tolerance = 1e-12
  )
  expect_equal(JAGS_marglik_parameters(samples, list(omega = log_independent))$omega, c(1, exp(.5)))
  expect_equal(
    JAGS_marglik_priors(samples, list(omega = log_independent)),
    mlpdf(log_independent$weights$prior, .5),
    tolerance = 1e-12
  )
  edge_samples <- samples
  edge_samples[["eta[1]"]] <- 0
  expect_equal(JAGS_marglik_priors(edge_samples, list(omega = cumulative)), -Inf)
  expect_error(
    JAGS_marglik_parameters(edge_samples, list(omega = cumulative)),
    "out-of-support positive auxiliary coordinate"
  )

  samples_two_sided <- c("eta[1]" = 1, "eta[2]" = 2, "eta[3]" = 3)
  expect_equal(
    JAGS_marglik_parameters(samples_two_sided, list(omega = two_sided))$omega,
    c(1, 5/6, 1/2, 5/6, 1),
    tolerance = 1e-12
  )
})

test_that("JAGS bridge posterior rejects monitored deterministic owned aliases", {

  cumulative <- prior_weightfunction("one-sided", c(.025), wf_cumulative(c(1, 2)))
  cumulative_posterior <- matrix(
    c(1, 2, 1, 1 / 3),
    nrow = 1,
    dimnames = list(NULL, c("eta[1]", "eta[2]", "omega[1]", "std_eta[1]"))
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      cumulative_posterior,
      prior_list = list(omega = cumulative),
      add_parameters = "omega[1]",
      add_bounds = list(lb = c("omega[1]" = 0), ub = c("omega[1]" = Inf))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      cumulative_posterior,
      prior_list = list(omega = cumulative),
      add_parameters = "std_eta[1]",
      add_bounds = list(lb = c("std_eta[1]" = 0), ub = c("std_eta[1]" = 1))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )
  two_sided <- prior_weightfunction("two-sided", c(.05, .10), wf_cumulative(c(1, 2, 3)))
  expect_error(
    JAGS_bridgesampling_posterior(
      cumulative_posterior,
      prior_list = list(omega = two_sided),
      add_parameters = "omega[5]",
      add_bounds = list(lb = c("omega[5]" = 0), ub = c("omega[5]" = Inf))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )

  phacking <- prior_phacking(form = "linear")
  phacking_posterior <- matrix(
    c(.2, 1, 1.5),
    nrow = 1,
    dimnames = list(NULL, c("alpha", "phack_kind", "phack_z_source[1]"))
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      phacking_posterior,
      prior_list = list(phacking = phacking),
      add_parameters = "phack_kind",
      add_bounds = list(lb = c(phack_kind = 0), ub = c(phack_kind = 2))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      phacking_posterior,
      prior_list = list(phacking = phacking),
      add_parameters = "omega[1]",
      add_bounds = list(lb = c("omega[1]" = 0), ub = c("omega[1]" = Inf))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      phacking_posterior,
      prior_list = list(phacking = phacking),
      add_parameters = "phack_z_source[1]",
      add_bounds = list(lb = c("phack_z_source[1]" = -Inf), ub = c("phack_z_source[1]" = Inf))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )

  theta <- prior_factor("invgamma", list(2, 1), contrast = "independent")
  attr(theta, "levels") <- 2
  theta_posterior <- matrix(
    c(1, 0.5),
    nrow = 1,
    dimnames = list(NULL, c("theta[1]", "theta[2]"))
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      theta_posterior,
      prior_list = list(theta = theta),
      add_parameters = "theta[1]",
      add_bounds = list(lb = c("theta[1]" = 0), ub = c("theta[1]" = Inf))
    ),
    "BayesTools-owned",
    fixed = TRUE
  )
})

test_that("heterogeneous bias mixtures map cumulative, omega, log-omega, fixed, and null weightfunctions", {

  bias <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)), prior_weights = 1),
    prior_weightfunction("one-sided", c(.05, .10), wf_independent(prior("gamma", list(shape = 9, rate = 3))), prior_weights = 1),
    prior_weightfunction("one-sided", c(.025), wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega"), prior_weights = 1),
    prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .4)), prior_weights = 1)
  ))

  cuts <- weightfunctions_mapping(bias[sapply(bias, is.prior.weightfunction)], cuts_only = TRUE, one_sided = TRUE)
  mapping <- weightfunctions_mapping(bias[sapply(bias, is.prior.weightfunction)], one_sided = TRUE)

  expect_equal(cuts, c(0, .025, .05, .10, .975, 1))
  expect_equal(mapping[[1]], c(1L, 2L, 3L, 3L, 3L))
  expect_equal(mapping[[2]], c(1L, 1L, 2L, 3L, 3L))
  expect_equal(mapping[[3]], c(1L, 2L, 2L, 2L, 2L))
  expect_equal(mapping[[4]], c(1L, 2L, 2L, 2L, 1L))

  syntax <- JAGS_add_priors("model{}", list(bias = bias))

  expect_match(syntax, "eta_component_2\\[1\\] ~ dgamma\\(1, 1\\)")
  expect_match(syntax, "omega_local_component_3\\[2\\] ~ dgamma\\(9,3\\)")
  expect_match(syntax, "log_omega_component_4\\[2\\] ~ dnorm")
  expect_match(syntax, "omega_local_component_4\\[2\\] <- exp\\(log_omega_component_4\\[2\\]\\)")
  expect_match(syntax, "omega_local_component_5\\[2\\] <- 0.4")
  expect_match(syntax, "omega_component_5\\[5\\] <- omega_local_component_5\\[1\\]")
  expect_match(syntax, "omega\\[3\\] <- omega_component_1\\[3\\] \\* equals\\(bias_indicator, 1\\)")
  expect_false(grepl("eta2omega", syntax, fixed = TRUE))

  set.seed(13)
  prior_samples <- rng(bias, 600)
  components <- attr(prior_samples, "components")

  expect_equal(colnames(prior_samples), paste0("omega[", 1:5, "]"))
  expect_true(all(prior_samples[components == 1, ] == 1))
  expect_true(all(prior_samples[components == 5, "omega[2]"] == .4))
  expect_true(all(prior_samples[components == 5, "omega[5]"] == 1))
  expect_gt(mean(prior_samples[components == 3, "omega[3]"] > 1), .90)
  expect_gt(mean(prior_samples[components == 4, "omega[2]"] > 1), .95)
})

test_that("JAGS syntax and fitting allow independent omega weights above one", {

  skip_if_not_test_profile("fixture")
  skip_if_not_installed("rjags")
  skip_if_missing_fits(c("fit_wf_independent_gamma", "fit_wf_independent_log"))

  omega_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("gamma", list(shape = 9, rate = 3)))
  )
  log_prior <- prior_weightfunction(
    "one-sided", c(.05),
    wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega")
  )

  omega_syntax <- JAGS_add_priors("model{}", list(omega = omega_prior))
  log_syntax   <- JAGS_add_priors("model{}", list(omega = log_prior))

  expect_match(omega_syntax, "omega\\[2\\] ~ dgamma\\(9,3\\)")
  expect_match(log_syntax, "log_omega\\[2\\] ~ dnorm")
  expect_match(log_syntax, "omega\\[2\\] <- exp\\(log_omega\\[2\\]\\)")

  omega_fit <- readRDS(file.path(temp_fits_dir, "fit_wf_independent_gamma.RDS"))
  log_fit   <- readRDS(file.path(temp_fits_dir, "fit_wf_independent_log.RDS"))

  omega_samples <- as.matrix(.fit_to_posterior(omega_fit))
  log_samples   <- as.matrix(.fit_to_posterior(log_fit))

  expect_true("omega[2]" %in% colnames(omega_samples))
  expect_true("omega[2]" %in% colnames(log_samples))
  expect_true("log_omega[2]" %in% colnames(log_samples))
  expect_gt(mean(omega_samples[, "omega[2]"] > 1), .90)
  expect_gt(mean(log_samples[, "omega[2]"] > 1), .90)
  expect_equal(
    unname(log_samples[, "omega[2]"]),
    unname(exp(log_samples[, "log_omega[2]"])),
    tolerance = 1e-8
  )
})

test_that("JAGS fits heterogeneous bias mixtures with omega and log-omega weights above one", {

  skip_if_not_test_profile("fixture")
  skip_if_not_installed("rjags")
  skip_if_missing_fits("fit_bias_heterogeneous_wf")

  fit <- readRDS(file.path(temp_fits_dir, "fit_bias_heterogeneous_wf.RDS"))

  posterior <- as.matrix(.fit_to_posterior(fit))
  expect_true(all(c("bias_indicator", paste0("omega[", 1:5, "]")) %in% colnames(posterior)))
  expect_true(all(1:5 %in% posterior[, "bias_indicator"]))

  indicator <- posterior[, "bias_indicator"]
  expect_true(all(posterior[indicator == 1, paste0("omega[", 1:5, "]")] == 1))
  expect_true(all(posterior[indicator == 2, paste0("omega[", 1:5, "]")] <= 1))
  expect_gt(mean(posterior[indicator == 3, "omega[3]"] > 1), .90)
  expect_gt(mean(posterior[indicator == 4, "omega[2]"] > 1), .95)
  expect_true(all(abs(posterior[indicator == 5, "omega[2]"] - .4) < 1e-8))
  expect_true(all(abs(posterior[indicator == 5, "omega[5]"] - 1) < 1e-8))

  mixed <- as_mixed_posteriors(fit, parameters = "bias", conditional = "omega")
  expect_equal(colnames(mixed$bias), c("omega[0,0.025]", "omega[0.025,0.05]", "omega[0.05,0.1]", "omega[0.1,0.975]", "omega[0.975,1]"))
  expect_true(all(attr(mixed$bias, "models_ind") %in% 2:5))
  expect_gt(mean(mixed$bias[attr(mixed$bias, "models_ind") == 3, "omega[0.05,0.1]"] > 1), .90)
  expect_gt(mean(mixed$bias[attr(mixed$bias, "models_ind") == 4, "omega[0.025,0.05]"] > 1), .95)
})

test_that("JAGS fits full bias mixtures with PET, PEESE, and heterogeneous weightfunctions", {

  skip_if_not_test_profile("fixture")
  skip_if_not_installed("rjags")
  skip_if_missing_fits("fit_bias_petpeese_hetero_wf")

  bias <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, .4), prior_weights = 1),
    prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)), prior_weights = 1),
    prior_weightfunction("one-sided", c(.05, .10), wf_independent(prior("gamma", list(shape = 9, rate = 3))), prior_weights = 1),
    prior_PEESE("gamma", list(shape = 3, rate = 2), prior_weights = 1),
    prior_weightfunction("one-sided", c(.025), wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega"), prior_weights = 1),
    prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .4)), prior_weights = 1)
  ))

  omega_names <- c("omega[0,0.025]", "omega[0.025,0.05]", "omega[0.05,0.1]", "omega[0.1,0.975]", "omega[0.975,1]")

  syntax <- JAGS_add_priors("model{}", list(bias = bias))
  expect_match(syntax, "PET <- PET_1 \\* equals\\(bias_indicator, 2\\)")
  expect_match(syntax, "eta_component_3\\[1\\] ~ dgamma\\(1, 1\\)")
  expect_match(syntax, "omega_local_component_4\\[2\\] ~ dgamma\\(9,3\\)")
  expect_match(syntax, "PEESE <- PEESE_1 \\* equals\\(bias_indicator, 5\\)")
  expect_match(syntax, "log_omega_component_6\\[2\\] ~ dnorm")
  expect_match(syntax, "omega_local_component_7\\[2\\] <- 0.4")
  expect_match(syntax, "omega_component_2\\[3\\] <- 1")
  expect_match(syntax, "omega_component_5\\[3\\] <- 1")
  expect_match(syntax, "omega\\[3\\] <- omega_component_1\\[3\\] \\* equals\\(bias_indicator, 1\\) \\+ omega_component_2\\[3\\] \\* equals\\(bias_indicator, 2\\)")

  fit <- readRDS(file.path(temp_fits_dir, "fit_bias_petpeese_hetero_wf.RDS"))

  posterior <- as.matrix(.fit_to_posterior(fit))
  indicator <- posterior[, "bias_indicator"]
  expect_true(all(c("bias_indicator", paste0("omega[", 1:5, "]"), "PET", "PEESE") %in% colnames(posterior)))
  expect_true(all(1:7 %in% indicator))

  expect_true(all(posterior[indicator %in% c(1, 2, 5), paste0("omega[", 1:5, "]")] == 1))
  expect_true(all(posterior[indicator == 3, paste0("omega[", 1:5, "]")] <= 1))
  expect_gt(mean(posterior[indicator == 4, "omega[3]"] > 1), .90)
  expect_gt(mean(posterior[indicator == 6, "omega[2]"] > 1), .95)
  expect_true(all(abs(posterior[indicator == 7, "omega[2]"] - .4) < 1e-8))
  expect_true(all(abs(posterior[indicator == 7, "omega[5]"] - 1) < 1e-8))

  expect_true(all(abs(posterior[indicator != 2, "PET"]) < 1e-8))
  expect_true(any(posterior[indicator == 2, "PET"] > 0))
  expect_true(all(abs(posterior[indicator != 5, "PEESE"]) < 1e-8))
  expect_true(any(posterior[indicator == 5, "PEESE"] > 0))

  mixed_all <- as_mixed_posteriors(fit, parameters = "bias")
  expect_equal(colnames(mixed_all$bias), c(omega_names, "PET", "PEESE"))
  expect_true(all(1:7 %in% attr(mixed_all$bias, "models_ind")))

  mixed_omega <- as_mixed_posteriors(fit, parameters = "bias", conditional = "omega")
  expect_equal(colnames(mixed_omega$bias), omega_names)
  expect_true(all(attr(mixed_omega$bias, "models_ind") %in% c(3, 4, 6, 7)))
  expect_gt(mean(mixed_omega$bias[attr(mixed_omega$bias, "models_ind") == 4, "omega[0.05,0.1]"] > 1), .90)
  expect_gt(mean(mixed_omega$bias[attr(mixed_omega$bias, "models_ind") == 6, "omega[0.025,0.05]"] > 1), .95)

  mixed_pet <- as_mixed_posteriors(fit, parameters = "bias", conditional = "PET")
  expect_equal(colnames(mixed_pet$bias), "PET")
  expect_true(all(attr(mixed_pet$bias, "models_ind") == 2))
  expect_true(all(mixed_pet$bias[, "PET"] > 0))

  mixed_peese <- as_mixed_posteriors(fit, parameters = "bias", conditional = "PEESE")
  expect_equal(colnames(mixed_peese$bias), "PEESE")
  expect_true(all(attr(mixed_peese$bias, "models_ind") == 5))
  expect_true(all(mixed_peese$bias[, "PEESE"] > 0))

  mixed_petpeese <- as_mixed_posteriors(fit, parameters = "bias", conditional = "PETPEESE")
  expect_equal(colnames(mixed_petpeese$bias), c("PET", "PEESE"))
  expect_true(all(attr(mixed_petpeese$bias, "models_ind") %in% c(2, 5)))
  expect_true(all(mixed_petpeese$bias[attr(mixed_petpeese$bias, "models_ind") == 2, "PEESE"] == 0))
  expect_true(all(mixed_petpeese$bias[attr(mixed_petpeese$bias, "models_ind") == 5, "PET"] == 0))

  table_samples <- suppressWarnings(runjags_estimates_table(
    fit,
    conditional        = TRUE,
    return_samples     = TRUE,
    remove_diagnostics = TRUE
  ))
  expect_true(all(c(omega_names, "PET", "PEESE") %in% colnames(table_samples)))
  expect_true(any(is.na(table_samples[, "PET"])))
  expect_true(any(table_samples[, "PET"] > 0, na.rm = TRUE))
  expect_true(any(is.na(table_samples[, "PEESE"])))
  expect_true(any(table_samples[, "PEESE"] > 0, na.rm = TRUE))
  expect_true(any(is.na(table_samples[, "omega[0.05,0.1]"])))
  expect_true(any(table_samples[, "omega[0.05,0.1]"] > 1, na.rm = TRUE))
})

test_that("omega diagnostics reject bias mixtures without weightfunctions", {

  fit <- structure(list(), class = c("runjags", "BayesTools_fit"))
  attr(fit, "prior_list") <- list(
    bias = prior_mixture(list(
      prior_PET("normal", list(0, 1), prior_weights = 1),
      prior_PEESE("normal", list(0, 1), prior_weights = 1)
    ))
  )

  expect_error(
    JAGS_diagnostics_density(fit, parameter = "omega"),
    "at least one weightfunction component"
  )
})
