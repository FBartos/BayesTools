skip_if_not_test_profile("unit")

# A fitted object with the metadata JAGS_fit() stores: the prior list, the
# monitored 'add_parameters' and model syntax, and the parameter map whose
# convergence roles the check reads.
.mock_convergence_fit <- function(chain_1, chain_2 = NULL, prior_list = list(),
                                  add_parameters = NULL, model_syntax = NULL,
                                  data_names = NULL){

  chains <- list(coda::mcmc(chain_1))
  if(!is.null(chain_2)){
    chains <- c(chains, list(coda::mcmc(chain_2)))
  }
  fit <- list(
    mcmc = do.call(coda::mcmc.list, chains),
    summary.pars = list(mutate = NULL)
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list
  attr(fit, "add_parameters") <- add_parameters
  attr(fit, "model_syntax") <- model_syntax
  fit <- BayesTools:::.bt_attach_parameter_map(fit, data_names = data_names)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  BayesTools:::.bt_attach_fit_contract(fit)
}

.convergence_states <- function(result){

  diagnostics <- attr(result, "diagnostics")
  stats::setNames(diagnostics$state, diagnostics$parameter)
}

.convergence_roles <- function(fit){

  coordinates <- parameter_coordinates(fit)
  stats::setNames(coordinates$convergence_role, coordinates$coordinate_name)
}

.convergence_weight_draws <- function(n){
  pmin(pmax(stats::rbeta(n, 4, 2), 1e-3), 1 - 1e-3)
}


test_that("JAGS_check_convergence preserves internal indicator and inclusion substrings", {

  parameters <- c(
    "theta_indicator_effect",
    "theta_indicator_effect[1]",
    "theta_inclusion_effect",
    "theta_inclusion_effect[1]"
  )
  for(parameter in parameters){
    chain_1 <- matrix(
      rep(-5, 100),
      ncol = 1,
      dimnames = list(NULL, parameter)
    )
    chain_2 <- matrix(
      rep(5, 100),
      ncol = 1,
      dimnames = list(NULL, parameter)
    )
    prior_list <- list(prior("normal", list(0, 1)))
    names(prior_list) <- sub("\\[.*$", "", parameter)

    convergence <- JAGS_check_convergence(
      .mock_convergence_fit(chain_1, chain_2, prior_list),
      max_Rhat = 1.05,
      min_ESS = NULL,
      max_error = NULL,
      max_SD_error = NULL
    )

    expect_false(convergence, info = parameter)
    expect_match(attr(convergence, "errors"), "R-hat", info = parameter)
  }
})

test_that("convergence roles come from declarations, not name suffixes", {

  set.seed(1)
  n <- 200
  # A spike-and-slab prior declares its indicator and inclusion probability.
  priors <- list(mu = prior_spike_and_slab(
    prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  ))
  chain <- function(state){
    cbind(
      mu_indicator = state,
      mu_inclusion = stats::runif(n),
      mu = stats::rnorm(n),
      mu_variable = stats::rnorm(n)
    )
  }
  fit <- .mock_convergence_fit(chain(0), chain(1), priors)
  expect_identical(
    .convergence_roles(fit),
    c(mu_indicator = "indicator", mu_inclusion = "auxiliary",
      mu = "sampled", mu_variable = "sampled")
  )
  check <- function(fit, ...){
    JAGS_check_convergence(
      fit, max_Rhat = 1.05, min_ESS = NULL, max_error = NULL,
      max_SD_error = NULL, ...
    )
  }
  expect_true(check(fit, check_indicators = FALSE))
  indicators <- check(fit, check_indicators = TRUE)
  expect_false(indicators)
  expect_match(attr(indicators, "errors"), "R-hat .* for 'mu_indicator'")

  # Continuous priors whose names end in '_indicator' or '_inclusion' are
  # ordinary sampled parameters: they are checked, never rejected as
  # indicators and never skipped as inclusion probabilities.
  named <- list(
    dose_indicator = prior("normal", list(0, 1)),
    p_inclusion = prior("beta", list(1, 1))
  )
  draws <- function(stuck){
    cbind(
      dose_indicator = stats::rnorm(n),
      p_inclusion = if(stuck) rep(0.3, n) else stats::rbeta(n, 1, 1)
    )
  }
  named_fit <- .mock_convergence_fit(draws(FALSE), draws(TRUE), named)
  expect_identical(
    .convergence_roles(named_fit),
    c(dose_indicator = "sampled", p_inclusion = "sampled")
  )
  named_check <- JAGS_check_convergence(
    named_fit, max_Rhat = 1.2, min_ESS = 50, max_error = NULL,
    max_SD_error = NULL
  )
  expect_false(named_check)
  expect_identical(
    .convergence_states(named_check),
    c(dose_indicator = "assessable", p_inclusion = "assessable")
  )
  expect_match(attr(named_check, "errors"), "for 'p_inclusion'", fixed = TRUE)
})

test_that("indicator diagnostics are invariant to categorical labels", {

  set.seed(12)
  chain_1_states <- sample(1:3, 2000, replace = TRUE, prob = c(.55, .3, .15))
  chain_2_states <- sample(1:3, 2000, replace = TRUE, prob = c(.55, .3, .15))
  product_prior <- list(theta = prior_mixture(list(
    prior("normal", list(0, 1)),
    prior("normal", list(1, 1)),
    prior("normal", list(2, 1))
  )))

  check_states <- function(chain_1, chain_2){
    fit <- .mock_convergence_fit(
      cbind(theta = stats::rnorm(2000), theta_indicator = chain_1),
      cbind(theta = stats::rnorm(2000), theta_indicator = chain_2),
      product_prior
    )
    JAGS_check_convergence(
      fit,
      max_Rhat = 1.2,
      min_ESS = 1,
      max_error = 1,
      max_SD_error = 1,
      check_indicators = TRUE,
      monitor = "theta_indicator"
    )
  }

  original <- check_states(chain_1_states, chain_2_states)
  expect_true(original)
  original_diagnostics <- attr(original, "diagnostics")
  original_diagnostics <- original_diagnostics[
    grepl("theta_indicator", original_diagnostics$parameter, fixed = TRUE),
  ]
  expect_equal(
    original_diagnostics$parameter,
    paste0("theta_indicator (state ", 1:3, ")")
  )
  expect_true(all(original_diagnostics$state == "assessable"))

  labels <- c(30, -4, 8)
  relabeled <- check_states(labels[chain_1_states], labels[chain_2_states])
  expect_true(relabeled)
  relabeled_diagnostics <- attr(relabeled, "diagnostics")
  relabeled_diagnostics <- relabeled_diagnostics[
    grepl("theta_indicator", relabeled_diagnostics$parameter, fixed = TRUE),
  ]
  original_state_order <- match(labels, sort(labels))
  diagnostic_columns <- c("Rhat", "ESS", "MCMC_error", "MCMC_SD_error")
  expect_equal(
    original_diagnostics[, diagnostic_columns],
    relabeled_diagnostics[original_state_order, diagnostic_columns],
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
})

test_that("binary indicator diagnostics are scale and label invariant", {

  set.seed(14)
  states_1 <- stats::rbinom(2000, 1, .4)
  states_2 <- stats::rbinom(2000, 1, .4)
  product_prior <- list(theta = prior_spike_and_slab(
    prior("normal", list(0, 1))
  ))
  check_states <- function(chain_1, chain_2){
    fit <- .mock_convergence_fit(
      cbind(theta = stats::rnorm(2000), theta_indicator = chain_1),
      cbind(theta = stats::rnorm(2000), theta_indicator = chain_2),
      product_prior
    )
    JAGS_check_convergence(
      fit,
      max_Rhat = 2,
      min_ESS = 1,
      max_error = 1,
      max_SD_error = 1,
      check_indicators = TRUE,
      monitor = "theta_indicator"
    )
  }

  original <- check_states(states_1, states_2)
  relabeled <- check_states(ifelse(states_1 == 1, -10, 40),
                            ifelse(states_2 == 1, -10, 40))
  diagnostic_columns <- c("Rhat", "ESS", "MCMC_error", "MCMC_SD_error")
  expect_true(original)
  expect_true(relabeled)
  expect_equal(
    attr(original, "diagnostics")[, diagnostic_columns],
    attr(relabeled, "diagnostics")[, diagnostic_columns],
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
})

test_that("only prior-supported one-state indicators are structural", {

  set.seed(13)
  chain_1 <- cbind(theta = stats::rnorm(100), theta_indicator = 1)
  chain_2 <- cbind(theta = stats::rnorm(100), theta_indicator = 1)
  structural_prior <- list(theta = prior_mixture(list(
    prior("normal", list(0, 1))
  )))
  structural <- JAGS_check_convergence(
    .mock_convergence_fit(chain_1, chain_2, structural_prior),
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL,
    check_indicators = TRUE,
    monitor = "theta_indicator"
  )
  expect_true(structural)
  expect_equal(
    attr(structural, "diagnostics")$state,
    c("not_requested", "structural_constant")
  )

  sampled_prior <- list(theta = prior_mixture(list(
    prior("normal", list(0, 1)),
    prior("normal", list(1, 1))
  )))
  sampled <- JAGS_check_convergence(
    .mock_convergence_fit(chain_1, chain_2, sampled_prior),
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL,
    check_indicators = TRUE,
    monitor = "theta_indicator"
  )
  expect_false(sampled)
  expect_equal(attr(sampled, "diagnostics")$state[[2L]], "not_assessable")

  # A spike-and-slab prior whose inclusion is a point at 0 or 1 fixes its
  # indicator and inclusion probability.
  point_inclusion <- list(mu = prior_spike_and_slab(
    prior("normal", list(0, 1)),
    prior_inclusion = prior("spike", list(1))
  ))
  draws <- function(){
    cbind(mu_indicator = 1, mu_inclusion = 1, mu = stats::rnorm(100),
          mu_variable = stats::rnorm(100))
  }
  fixed <- .mock_convergence_fit(draws(), draws(), point_inclusion)
  expect_identical(
    unname(.convergence_roles(fixed)[c("mu_indicator", "mu_inclusion")]),
    c("structural", "structural")
  )
  # At inclusion 1 the effective coordinate is the slab and stays sampled.
  expect_identical(unname(.convergence_roles(fixed)["mu"]), "sampled")

  # At inclusion 0 the effective coordinate is the spike, a constant zero:
  # structural like the indicator and the inclusion probability, while the
  # slab coordinate keeps its own role.
  zero_inclusion <- list(theta = prior_spike_and_slab(
    prior("normal", list(0, 1)),
    prior_inclusion = prior("point", list(0))
  ))
  zero_draws <- function(){
    cbind(theta = 0, theta_indicator = 0, theta_inclusion = 0,
          theta_variable = stats::rnorm(200))
  }
  zero <- .mock_convergence_fit(zero_draws(), zero_draws(), zero_inclusion)
  expect_identical(
    .convergence_roles(zero),
    c(theta = "structural", theta_indicator = "structural",
      theta_inclusion = "structural", theta_variable = "sampled")
  )
  zero_result <- JAGS_check_convergence(
    zero, max_Rhat = NULL, min_ESS = 1, max_error = NULL, max_SD_error = NULL
  )
  expect_true(zero_result)
  expect_identical(.convergence_states(zero_result)[["theta"]], "structural_constant")
  expect_identical(.convergence_states(zero_result)[["theta_variable"]], "assessable")
})

test_that("sampled constants are not assessable but point priors are structural", {

  constant_chain <- matrix(
    rep(1, 20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  )
  sampled <- JAGS_check_convergence(
    .mock_convergence_fit(
      constant_chain, constant_chain,
      list(mu = prior("normal", list(0, 1)))
    ),
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL
  )
  expect_false(sampled)
  expect_match(attr(sampled, "errors"), "R-hat.*not assessable")
  sampled_diagnostics <- attr(sampled, "diagnostics")
  expect_s3_class(
    sampled_diagnostics,
    "BayesTools_convergence_diagnostics"
  )
  expect_equal(sampled_diagnostics$state, "not_assessable")

  structural <- JAGS_check_convergence(
    .mock_convergence_fit(
      constant_chain, constant_chain,
      list(mu = prior("point", list(1)))
    ),
    max_Rhat = 1.05,
    min_ESS = 10,
    max_error = 0.01,
    max_SD_error = 0.05
  )
  expect_true(structural)
  expect_equal(
    attr(structural, "diagnostics")$state,
    "structural_constant"
  )
})

test_that("a point prior with an expression location is derived", {

  set.seed(15)
  draws <- function(){
    a <- stats::rnorm(200)
    cbind(a = a, b = a)
  }
  priors <- list(
    a = prior("normal", list(0, 1)),
    b = prior("point", list(location = expression(a)))
  )
  fit <- .mock_convergence_fit(draws(), draws(), priors)
  expect_identical(.convergence_roles(fit), c(a = "sampled", b = "derived"))
  # The parameter map classifies it as a derived coordinate as well.
  coordinates <- parameter_coordinates(fit)
  expect_identical(coordinates$monitor_status[coordinates$coordinate_name == "b"], "sampled")
  expect_true(is.na(coordinates$fixed_value[coordinates$coordinate_name == "b"]))

  result <- JAGS_check_convergence(
    fit, max_Rhat = 1.2, min_ESS = 1, max_error = NULL, max_SD_error = NULL
  )
  expect_true(result)
  expect_identical(
    .convergence_states(result),
    c(a = "assessable", b = "not_requested")
  )
  requested <- JAGS_check_convergence(
    fit, max_Rhat = 1.2, min_ESS = 1, max_error = NULL, max_SD_error = NULL,
    monitor = "b"
  )
  expect_identical(
    .convergence_states(requested),
    c(a = "not_requested", b = "assessable")
  )
})

test_that("a supplied prior list must be the fit's own", {

  set.seed(47)
  chain <- function() cbind(mu = stats::rnorm(100), tau = abs(stats::rnorm(100)))
  priors <- list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(0, 1), list(0, Inf))
  )
  fit <- .mock_convergence_fit(chain(), chain(), priors)
  check <- function(prior_list){
    JAGS_check_convergence(
      fit, prior_list = prior_list, max_Rhat = 1.2, min_ESS = 1,
      max_error = NULL, max_SD_error = NULL
    )
  }

  expected <- check(NULL)
  expect_identical(check(attr(fit, "prior_list")), expected)
  expect_identical(check(priors), expected)

  mismatch <- paste0(
    "'prior_list' differs from the prior list stored with 'fit'. ",
    "The convergence classification comes from the fit's coordinate roles; ",
    "omit 'prior_list' or supply 'attr(fit, \"prior_list\")'."
  )
  # A different prior, or a list that would have classified 'tau' as a
  # structural constant, is not the fitted prior list.
  expect_error(check(priors["mu"]), mismatch, fixed = TRUE)
  expect_error(
    check(list(mu = priors$mu, tau = prior("point", list(1)))),
    mismatch,
    fixed = TRUE
  )
  expect_error(check(list()), mismatch, fixed = TRUE)
})

test_that("explicit convergence monitors distinguish omitted parameters", {

  set.seed(42)
  chain_1 <- cbind(mu = rnorm(100), tau = rep(1, 100))
  chain_2 <- cbind(mu = rnorm(100), tau = rep(1, 100))
  priors <- list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(0, 1))
  )
  fit <- .mock_convergence_fit(chain_1, chain_2, priors)

  selected <- JAGS_check_convergence(
    fit,
    max_Rhat = NULL,
    min_ESS = 1,
    max_error = NULL,
    max_SD_error = NULL,
    monitor = "mu"
  )
  expect_true(selected)
  diagnostics <- attr(selected, "diagnostics")
  expect_equal(
    setNames(diagnostics$state, diagnostics$parameter),
    c(mu = "assessable", tau = "not_requested")
  )

  all_parameters <- JAGS_check_convergence(
    fit,
    max_Rhat = NULL,
    min_ESS = 1,
    max_error = NULL,
    max_SD_error = NULL
  )
  expect_false(all_parameters)
  expect_equal(
    attr(all_parameters, "diagnostics")$state,
    c("assessable", "not_assessable")
  )
  expect_error(
    JAGS_check_convergence(
      fit,
      monitor = "misspelled"
    ),
    "requested convergence monitor 'misspelled' is not available"
  )
})

test_that("targeted convergence checks retain requested product-space indicators", {

  set.seed(45)
  chain_1 <- cbind(
    mu = stats::rnorm(200),
    mu_indicator = 0,
    mu_inclusion = 0
  )
  chain_2 <- cbind(
    mu = stats::rnorm(200),
    mu_indicator = 1,
    mu_inclusion = 1
  )
  priors <- list(mu = prior_spike_and_slab(
    prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  ))
  fit <- .mock_convergence_fit(chain_1, chain_2, priors)

  targeted <- JAGS_check_convergence(
    fit,
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL,
    check_indicators = TRUE,
    monitor = "mu"
  )

  # The indicator is stuck at 0 in one chain and at 1 in the other: its pooled
  # draws vary, so R-hat is the value Inf that fails the criterion.
  expect_false(targeted)
  diagnostics <- attr(targeted, "diagnostics")
  expect_equal(
    setNames(diagnostics$state, diagnostics$parameter),
    c(
      mu = "assessable",
      mu_indicator = "assessable",
      mu_inclusion = "not_requested"
    )
  )
  expect_match(attr(targeted, "errors"), "R-hat Inf for 'mu_indicator'",
               fixed = TRUE, all = FALSE)

  explicit_inclusion <- JAGS_check_convergence(
    fit,
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL,
    check_indicators = FALSE,
    monitor = c("mu", "mu_inclusion")
  )
  expect_false(explicit_inclusion)
  expect_equal(
    attr(explicit_inclusion, "diagnostics")$state,
    c("assessable", "not_requested", "assessable")
  )
  expect_match(attr(explicit_inclusion, "errors"), "R-hat Inf for 'mu_inclusion'",
               fixed = TRUE, all = FALSE)
})

test_that("a base convergence monitor selects every indexed element", {

  chain_1 <- cbind(
    mu = 1:20,
    "beta[1]" = 1:20,
    "beta[2]" = 20:1
  )
  chain_2 <- cbind(
    mu = 2:21,
    "beta[1]" = 2:21,
    "beta[2]" = 21:2
  )
  result <- JAGS_check_convergence(
    .mock_convergence_fit(chain_1, chain_2, list(
      mu = prior("normal", list(0, 1)),
      beta = prior("mnormal", list(mean = 0, sd = 1, K = 2))
    )),
    max_Rhat = NULL,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL,
    monitor = "beta"
  )

  expect_true(result)
  diagnostics <- attr(result, "diagnostics")
  expect_equal(
    diagnostics$parameter[diagnostics$state == "assessable"],
    c("beta[1]", "beta[2]")
  )
  expect_equal(
    diagnostics$state[diagnostics$parameter == "mu"],
    "not_requested"
  )
})

test_that("empty convergence selections do not claim convergence", {

  set.seed(43)
  chain_1 <- cbind(mu = rnorm(20))
  chain_2 <- cbind(mu = rnorm(20))
  result <- JAGS_check_convergence(
    .mock_convergence_fit(chain_1, chain_2, list(mu = prior("normal", list(0, 1)))),
    monitor = character()
  )

  expect_identical(as.vector(result), logical())
  expect_s3_class(
    attr(result, "diagnostics"),
    "BayesTools_convergence_diagnostics"
  )
  expect_equal(nrow(attr(result, "diagnostics")), 0L)
})

test_that("empty available convergence parameters are vacuously TRUE", {

  set.seed(44)
  chain_1 <- cbind(mu = rnorm(20))
  chain_2 <- cbind(mu = rnorm(20))
  # Drop every sample column via add_parameters removal with no structural priors.
  result <- JAGS_check_convergence(
    .mock_convergence_fit(chain_1, chain_2),
    add_parameters = "mu",
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL
  )

  expect_true(result)
  expect_equal(nrow(attr(result, "diagnostics")), 0L)
})

test_that("not-assessable diagnostics require an explicit opt-in to ignore", {

  constant_chain <- matrix(
    rep(1, 20),
    ncol = 1,
    dimnames = list(NULL, "mu")
  )
  result <- JAGS_check_convergence(
    .mock_convergence_fit(
      constant_chain, constant_chain,
      list(mu = prior("normal", list(0, 1)))
    ),
    max_Rhat = 1.05,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL,
    allow_not_assessable = TRUE
  )

  expect_true(result)
  expect_equal(
    attr(result, "diagnostics")$state,
    "not_assessable"
  )
  expect_null(attr(result, "errors"))
})

test_that("chains stuck at different values fail the criteria and are not ignored", {

  # Each chain is constant but at a different value: the pooled draws vary, so
  # R-hat is Inf and ESS is 0. These are values that fail the criteria, not
  # undefined results that allow_not_assessable may ignore.
  stuck_low <- matrix(rep(-5, 200), ncol = 1, dimnames = list(NULL, "mu"))
  stuck_high <- matrix(rep(5, 200), ncol = 1, dimnames = list(NULL, "mu"))
  stuck_fit <- .mock_convergence_fit(
    stuck_low, stuck_high, list(mu = prior("normal", list(0, 1)))
  )
  for(allow in c(TRUE, FALSE)){
    result <- JAGS_check_convergence(stuck_fit, allow_not_assessable = allow)
    expect_false(result)
    errors <- attr(result, "errors")
    expect_match(errors, "R-hat Inf for 'mu' is larger than the set target (1.05).",
                 fixed = TRUE, all = FALSE)
    expect_match(errors, "ESS 0 for 'mu' is lower than the set target (500).",
                 fixed = TRUE, all = FALSE)
    # the time-series SE of constant chains is 0, which measures no error: the
    # MCMC error criteria are not assessable (F322) and only these are ignored
    # on request; R-hat and ESS stay values
    expect_identical(
      any(grepl("MCMC error, MCMC SD error diagnostics are not assessable for 'mu'.", errors, fixed = TRUE)),
      !allow
    )
    diagnostics <- attr(result, "diagnostics")
    expect_identical(diagnostics$state, "not_assessable")
    expect_identical(diagnostics$Rhat, Inf)
    expect_identical(diagnostics$ESS, 0)
    expect_true(is.na(diagnostics$MCMC_error))
    expect_true(is.na(diagnostics$MCMC_SD_error))
  }

  # one criterion at a time
  r_hat_only <- JAGS_check_convergence(
    stuck_fit, max_Rhat = 1.05, min_ESS = NULL, max_error = NULL,
    max_SD_error = NULL, allow_not_assessable = TRUE
  )
  expect_false(r_hat_only)
  expect_match(attr(r_hat_only, "errors"), "R-hat Inf for 'mu'", fixed = TRUE)
  ess_only <- JAGS_check_convergence(
    stuck_fit, max_Rhat = NULL, min_ESS = 500, max_error = NULL,
    max_SD_error = NULL, allow_not_assessable = TRUE
  )
  expect_false(ess_only)
  expect_match(attr(ess_only, "errors"), "ESS 0 for 'mu'", fixed = TRUE)

  # a model indicator stuck in different components per chain
  indicator_priors <- list(mu = prior_mixture(list(
    prior("normal", list(0, 1)),
    prior("normal", list(1, 1))
  )))
  indicator_fit <- .mock_convergence_fit(
    cbind(mu_indicator = rep(1, 200)), cbind(mu_indicator = rep(2, 200)),
    indicator_priors
  )
  stuck_indicator <- JAGS_check_convergence(
    indicator_fit, check_indicators = TRUE, allow_not_assessable = TRUE
  )
  expect_false(stuck_indicator)
  expect_match(attr(stuck_indicator, "errors"), "R-hat Inf for 'mu_indicator",
               fixed = TRUE, all = FALSE)

  # a column that is constant when pooled stays not assessable, unchanged
  constant <- matrix(rep(1, 200), ncol = 1, dimnames = list(NULL, "mu"))
  constant_fit <- .mock_convergence_fit(
    constant, constant, list(mu = prior("normal", list(0, 1)))
  )
  ignored <- JAGS_check_convergence(constant_fit, allow_not_assessable = TRUE)
  expect_true(ignored)
  expect_identical(attr(ignored, "diagnostics")$state, "not_assessable")
  expect_true(is.na(attr(ignored, "diagnostics")$Rhat))
  expect_true(is.na(attr(ignored, "diagnostics")$ESS))
  refused <- JAGS_check_convergence(constant_fit)
  expect_false(refused)
  expect_match(attr(refused, "errors"), "not assessable for 'mu'", fixed = TRUE, all = FALSE)
})

test_that("chains stuck at different values do not pass the MCMC error criteria (F322)", {

  # An error-only check used to accept them: the time-series SE of constant
  # chains is exactly 0 although the pooled draws vary.
  stuck_fit <- .mock_convergence_fit(
    matrix(rep(-5, 200), ncol = 1, dimnames = list(NULL, "mu")),
    matrix(rep(5, 200), ncol = 1, dimnames = list(NULL, "mu")),
    list(mu = prior("normal", list(0, 1)))
  )
  criteria <- list(
    error    = list(max_error = 0.01, max_SD_error = NULL),
    SD_error = list(max_error = NULL, max_SD_error = 0.05),
    both     = list(max_error = 0.01, max_SD_error = 0.05)
  )
  unavailable <- c(
    error    = "MCMC error diagnostic is not assessable for 'mu'.",
    SD_error = "MCMC SD error diagnostic is not assessable for 'mu'.",
    both     = "MCMC error, MCMC SD error diagnostics are not assessable for 'mu'."
  )
  for(name in names(criteria)){
    check <- function(...){
      do.call(JAGS_check_convergence, c(
        list(stuck_fit, max_Rhat = NULL, min_ESS = NULL), criteria[[name]], list(...)
      ))
    }
    result <- check()
    expect_false(result, info = name)
    expect_identical(attr(result, "errors"), unavailable[[name]], info = name)
    diagnostics <- attr(result, "diagnostics")
    expect_identical(diagnostics$state, "not_assessable", info = name)
    expect_true(is.na(diagnostics$MCMC_error) && is.na(diagnostics$MCMC_SD_error), info = name)

    # a request to ignore the undefined results ignores them
    ignored <- check(allow_not_assessable = TRUE)
    expect_true(ignored, info = name)
    expect_identical(attr(ignored, "diagnostics")$state, "not_assessable", info = name)
  }

  # chains that vary have a positive time-series SE and are assessed as before
  set.seed(7)
  mixed_fit <- .mock_convergence_fit(
    matrix(stats::rnorm(400), ncol = 1, dimnames = list(NULL, "mu")),
    matrix(stats::rnorm(400), ncol = 1, dimnames = list(NULL, "mu")),
    list(mu = prior("normal", list(0, 1)))
  )
  result <- JAGS_check_convergence(mixed_fit, max_Rhat = NULL, min_ESS = NULL, max_error = 1, max_SD_error = 1)
  expect_true(result)
  diagnostics <- attr(result, "diagnostics")
  expect_identical(diagnostics$state, "assessable")
  expect_gt(diagnostics$MCMC_error, 0)
  expect_gt(diagnostics$MCMC_SD_error, 0)
})

test_that("weightfunction reference omega bins are structural constants", {

  set.seed(21)
  n <- 200
  cumulative <- prior_weightfunction("one-sided", .05, wf_cumulative(c(2, 4)))
  chain_1 <- cbind(
    mu = stats::rnorm(n),
    "omega[1]" = 1,
    "omega[2]" = pmin(pmax(stats::rbeta(n, 4, 2), 1e-3), 1 - 1e-3)
  )
  chain_2 <- cbind(
    mu = stats::rnorm(n),
    "omega[1]" = 1,
    "omega[2]" = pmin(pmax(stats::rbeta(n, 4, 2), 1e-3), 1 - 1e-3)
  )
  priors <- list(
    mu = prior("normal", list(0, 1)),
    omega = cumulative
  )
  fit <- .mock_convergence_fit(chain_1, chain_2, priors)
  expect_identical(
    .convergence_roles(fit),
    c(mu = "sampled", "omega[1]" = "structural", "omega[2]" = "sampled")
  )

  conv <- JAGS_check_convergence(
    fit,
    max_Rhat = 1.2,
    min_ESS = 1,
    max_error = 1,
    max_SD_error = 1
  )
  expect_true(conv)
  diagnostics <- attr(conv, "diagnostics")
  expect_equal(
    setNames(diagnostics$state, diagnostics$parameter),
    c(
      mu = "assessable",
      "omega[0,0.05]" = "structural_constant",
      "omega[0.05,1]" = "assessable"
    )
  )

  monitored <- JAGS_check_convergence(
    fit,
    max_Rhat = 1.2,
    min_ESS = 1,
    max_error = 1,
    max_SD_error = 1,
    monitor = "omega"
  )
  expect_true(monitored)
  expect_equal(
    attr(monitored, "diagnostics")$state[attr(monitored, "diagnostics")$parameter == "omega[0,0.05]"],
    "structural_constant"
  )

  mix <- prior_mixture(list(prior_none(), cumulative))
  mixture <- JAGS_check_convergence(
    .mock_convergence_fit(
      cbind(mu = chain_1[, "mu"], "omega[1]" = 1, "omega[2]" = chain_1[, "omega[2]"]),
      cbind(mu = chain_2[, "mu"], "omega[1]" = 1, "omega[2]" = chain_2[, "omega[2]"]),
      list(mu = priors$mu, bias = mix)
    ),
    max_Rhat = 1.2,
    min_ESS = 1,
    max_error = 1,
    max_SD_error = 1
  )
  expect_true(mixture)
  mixture_diagnostics <- attr(mixture, "diagnostics")
  expect_equal(
    mixture_diagnostics$state[mixture_diagnostics$parameter == "omega[1]"],
    "structural_constant"
  )
  expect_equal(
    mixture_diagnostics$state[mixture_diagnostics$parameter == "omega[2]"],
    "assessable"
  )
})

test_that("fixed weightfunction omega bins are structural constants", {

  set.seed(22)
  n <- 80
  fixed <- prior_weightfunction("one-sided", .05, wf_fixed(c(1, .5)))
  chain_1 <- cbind(
    mu = stats::rnorm(n),
    "omega[1]" = 1,
    "omega[2]" = .5
  )
  chain_2 <- cbind(
    mu = stats::rnorm(n),
    "omega[1]" = 1,
    "omega[2]" = .5
  )
  conv <- JAGS_check_convergence(
    .mock_convergence_fit(chain_1, chain_2, list(
      mu = prior("normal", list(0, 1)),
      omega = fixed
    )),
    max_Rhat = NULL,
    min_ESS = NULL,
    max_error = NULL,
    max_SD_error = NULL
  )
  expect_true(conv)
  expect_equal(
    setNames(attr(conv, "diagnostics")$state, attr(conv, "diagnostics")$parameter),
    c(
      mu = "assessable",
      "omega[0,0.05]" = "structural_constant",
      "omega[0.05,1]" = "structural_constant"
    )
  )
})

test_that("composed bias priors keep every constant omega bin structural", {

  # Column names and constants are those of prior-only JAGS fits of the same
  # priors: a composed bias prior expands its selection onto the one-sided
  # grid and is summarized under renamed bins, a bias mixture keeps raw bins.
  set.seed(31)
  n <- 200
  mu_prior <- prior("normal", list(0, 1))
  check <- function(bias, make, ...){
    JAGS_check_convergence(
      .mock_convergence_fit(make(), make(), list(mu = mu_prior, bias = bias)),
      max_Rhat = 1.2,
      min_ESS = 1,
      max_error = 1,
      max_SD_error = 1,
      ...
    )
  }

  two_sided <- check(
    prior_bias(selection = prior_weightfunction(
      "two-sided", .05, wf_cumulative(c(1, 1))
    )),
    function() cbind(
      mu = stats::rnorm(n),
      "omega[1]" = 1,
      "omega[2]" = .convergence_weight_draws(n),
      "omega[3]" = 1
    )
  )
  expect_true(two_sided)
  expect_equal(
    .convergence_states(two_sided),
    c(
      mu = "assessable",
      "omega[0,0.025]" = "structural_constant",
      "omega[0.025,0.975]" = "assessable",
      "omega[0.975,1]" = "structural_constant"
    )
  )

  fixed <- check(
    prior_bias(selection = prior_weightfunction(
      "one-sided", c(.025, .05), wf_fixed(c(1, .7, .4))
    )),
    function() cbind(
      mu = stats::rnorm(n),
      "omega[1]" = 1,
      "omega[2]" = .7,
      "omega[3]" = .4
    )
  )
  expect_true(fixed)
  expect_equal(
    .convergence_states(fixed),
    c(
      mu = "assessable",
      "omega[0,0.025]" = "structural_constant",
      "omega[0.025,0.05]" = "structural_constant",
      "omega[0.05,1]" = "structural_constant"
    )
  )

  mixture_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("two-sided", .05, wf_cumulative(c(1, 1)), prior_weights = 1),
    prior_weightfunction("two-sided", c(.05, .10), wf_cumulative(c(1, 1, 1)), prior_weights = 1)
  ), is_null = c(TRUE, FALSE, FALSE))
  mixture <- check(
    mixture_prior,
    function() cbind(
      mu = stats::rnorm(n),
      bias_indicator = rep(1:3, length.out = n),
      "omega[1]" = 1,
      "omega[2]" = .convergence_weight_draws(n),
      "omega[3]" = .convergence_weight_draws(n),
      "omega[4]" = .convergence_weight_draws(n),
      "omega[5]" = 1
    )
  )
  expect_true(mixture)
  expect_equal(
    .convergence_states(mixture)[paste0("omega[", 1:5, "]")],
    stats::setNames(
      c("structural_constant", "assessable", "assessable",
        "assessable", "structural_constant"),
      paste0("omega[", 1:5, "]")
    )
  )

  # A shared constant across fixed and reference bins is structural; a bin
  # that differs between branches is sampled through the indicator.
  shared_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one-sided", .025, wf_fixed(c(1, 1)), prior_weights = 1),
    prior_weightfunction("two-sided", .05, wf_cumulative(c(1, 1)), prior_weights = 1)
  ), is_null = c(TRUE, FALSE, FALSE))
  expect_identical(
    BayesTools:::.bt_convergence_role_constant_bins(shared_prior),
    c("omega[1]", "omega[3]")
  )
  differing_prior <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one-sided", .025, wf_fixed(c(1, .5)), prior_weights = 1)
  ), is_null = c(TRUE, FALSE))
  expect_identical(
    BayesTools:::.bt_convergence_role_constant_bins(differing_prior),
    "omega[1]"
  )

  # Without a step selection, the single omega coordinate is constant.
  phacking <- check(
    prior_bias(phacking = prior_phacking()),
    function() cbind(
      mu = stats::rnorm(n),
      omega = 1,
      pi_null = .convergence_weight_draws(n)
    ),
    monitor = c("mu", "omega", "pi_null")
  )
  expect_true(phacking)
  expect_equal(
    .convergence_states(phacking),
    c(mu = "assessable", omega = "structural_constant", pi_null = "assessable")
  )
})

test_that("declared p-hacking kinds are structural convergence constants", {

  # Columns and constants of prior-only JAGS fits: the reported pi_null is
  # sampled, alpha is dropped as unreported, and omega and phack_kind are
  # fixed by the declared p-hacking form.
  set.seed(34)
  n <- 200
  mu_prior <- prior("normal", list(0, 1))
  check <- function(bias, make){
    JAGS_check_convergence(
      .mock_convergence_fit(make(), make(), list(mu = mu_prior, bias = bias)),
      max_Rhat = 1.2,
      min_ESS = 1,
      max_error = 1,
      max_SD_error = 1
    )
  }
  single_draws <- function(){
    alpha <- stats::rbeta(n, 1, 1)
    cbind(
      mu = stats::rnorm(n),
      omega = 1,
      alpha = alpha,
      phack_kind = 1,
      pi_null = 0.1 * alpha
    )
  }

  for(bias in list(prior_phacking(), prior_bias(phacking = prior_phacking()))){
    single <- check(bias, single_draws)
    expect_true(single)
    expect_null(attr(single, "errors"))
    expect_equal(
      .convergence_states(single),
      c(
        mu = "assessable",
        omega = "structural_constant",
        phack_kind = "structural_constant",
        pi_null = "assessable"
      )
    )
  }

  # Mixture branches that declare the same form share the constant kind.
  same_kind <- prior_mixture(list(
    prior_phacking(source = .25, prior_weights = 1),
    prior_phacking(source = .5, prior_weights = 1)
  ))
  same <- check(same_kind, function(){
    cbind(
      mu = stats::rnorm(n),
      bias_indicator = rep(1:2, length.out = n),
      omega = 1,
      alpha = stats::rbeta(n, 1, 1),
      phack_kind = 1,
      pi_null = stats::rbeta(n, 1, 9)
    )
  })
  expect_true(same)
  expect_identical(
    unname(.convergence_states(same)["phack_kind"]),
    "structural_constant"
  )

  # Kinds that vary across branches (including branches without p-hacking,
  # which declare 0) follow the mixture indicator and stay assessable.
  varying_forms <- prior_mixture(list(
    prior_phacking(form = "linear", prior_weights = 1),
    prior_phacking(form = "quadratic", prior_weights = 1)
  ))
  with_none <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_phacking(prior_weights = 1)
  ), is_null = c(TRUE, FALSE))
  expect_false(BayesTools:::.bt_convergence_role_constant_phacking_kind(varying_forms))
  expect_false(BayesTools:::.bt_convergence_role_constant_phacking_kind(with_none))
  varying <- check(varying_forms, function(){
    indicator <- sample(1:2, n, replace = TRUE)
    cbind(
      mu = stats::rnorm(n),
      bias_indicator = indicator,
      omega = 1,
      alpha = stats::rbeta(n, 1, 1),
      phack_kind = indicator,
      pi_null = stats::rbeta(n, 1, 9)
    )
  })
  expect_true(varying)
  expect_identical(unname(.convergence_states(varying)["phack_kind"]), "assessable")
})

test_that("ordered priors with point totals keep their constants structural", {

  set.seed(32)
  n <- 200
  data <- data.frame(f = ordered(
    rep(c("low", "mid", "high"), each = 2),
    levels = c("low", "mid", "high")
  ))
  ordered_prior_list <- function(f_prior){
    JAGS_formula(
      ~ 1 + f,
      parameter = "mu",
      data = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        f = f_prior
      )
    )$prior_list
  }
  eta_names <- paste0("prior_par_eta_mu_f_ordered_alloc_f_1[", 1:2, "]")
  check <- function(prior_list, coefficients, total, allocation = TRUE){
    make <- function(){
      chain <- cbind(
        mu_intercept = stats::rnorm(n),
        "mu_f[1]" = coefficients[[1L]](),
        "mu_f[2]" = coefficients[[2L]](),
        mu_f_ordered_total = total
      )
      if(allocation){
        chain <- cbind(
          chain,
          matrix(stats::rgamma(2L * n, 1), ncol = 2L, dimnames = list(NULL, eta_names))
        )
      }
      chain
    }
    JAGS_check_convergence(
      .mock_convergence_fit(make(), make(), prior_list),
      max_Rhat = 1.2,
      min_ESS = 1,
      max_error = 1,
      max_SD_error = 1
    )
  }
  constant <- function(value) function() rep(value, n)
  sampled  <- function() stats::rnorm(n)

  zero_total <- check(
    ordered_prior_list(prior_ordered(prior("point", list(0)))),
    list(constant(0), constant(0)),
    total = 0
  )
  expect_true(zero_total)
  expect_equal(
    .convergence_states(zero_total)[c("mu_f[1]", "mu_f[2]", "mu_f_ordered_total", eta_names)],
    stats::setNames(
      c(rep("structural_constant", 3L), rep("assessable", 2L)),
      c("mu_f[1]", "mu_f[2]", "mu_f_ordered_total", eta_names)
    )
  )

  nonzero_total <- check(
    ordered_prior_list(prior_ordered(prior("point", list(0.5)))),
    list(sampled, sampled),
    total = 0.5
  )
  expect_true(nonzero_total)
  expect_equal(
    .convergence_states(nonzero_total)[c("mu_f[1]", "mu_f[2]", "mu_f_ordered_total")],
    c("mu_f[1]" = "assessable", "mu_f[2]" = "assessable",
      mu_f_ordered_total = "structural_constant")
  )

  fixed_split <- check(
    ordered_prior_list(prior_ordered(prior("point", list(0.5)), allocation = c(0.3, 0.7))),
    list(constant(0.15), constant(0.35)),
    total = 0.5,
    allocation = FALSE
  )
  expect_true(fixed_split)
  expect_equal(
    unname(.convergence_states(fixed_split)[c("mu_f[1]", "mu_f[2]", "mu_f_ordered_total")]),
    rep("structural_constant", 3L)
  )

  # A sampled total keeps all its coordinates sampled.
  sampled_total <- ordered_prior_list(prior_ordered(prior("normal", list(0, 1))))
  sampled_roles <- .convergence_roles(.mock_convergence_fit(
    cbind(mu_intercept = sampled(), "mu_f[1]" = sampled(), "mu_f[2]" = sampled(),
          mu_f_ordered_total = sampled()),
    cbind(mu_intercept = sampled(), "mu_f[1]" = sampled(), "mu_f[2]" = sampled(),
          mu_f_ordered_total = sampled()),
    sampled_total
  ))
  expect_true(all(sampled_roles == "sampled"))
})

test_that("explicit convergence monitors reach additional monitored parameters", {

  set.seed(33)
  make <- function() cbind(mu = stats::rnorm(200), theta = stats::rnorm(200))
  fit <- .mock_convergence_fit(make(), make(), list(mu = prior("normal", list(0, 1))))

  # Parameters listed in 'add_parameters' are excluded from the default selection.
  default <- JAGS_check_convergence(
    fit, add_parameters = "theta",
    max_Rhat = 1.2, min_ESS = 1, max_error = 1, max_SD_error = 1
  )
  expect_true(default)
  expect_equal(.convergence_states(default), c(mu = "assessable"))

  requested <- JAGS_check_convergence(
    fit, add_parameters = "theta", monitor = "theta",
    max_Rhat = 1.2, min_ESS = 1, max_error = 1, max_SD_error = 1
  )
  expect_true(requested)
  expect_equal(
    .convergence_states(requested),
    c(mu = "not_requested", theta = "assessable")
  )
  expect_error(
    JAGS_check_convergence(
      fit, add_parameters = "theta", monitor = "theta[2]"
    ),
    "The requested convergence monitor 'theta[2]' is not available in the fitted model.",
    fixed = TRUE
  )

  expect_silent(BayesTools:::.bt_convergence_validate_monitor_names(
    c("theta", "mu", "bias_indicator (state 2)", "omega[0,0.05]"),
    c("mu", "theta", "bias_indicator", "omega")
  ))
  expect_error(
    BayesTools:::.bt_convergence_validate_monitor_names("thetaa", c("mu", "theta")),
    "The requested convergence monitor 'thetaa' is not monitored by the model.",
    fixed = TRUE
  )
})

test_that("unavailable monitors carry the parameter-not-found class and add_parameters refuses NA (F212)", {

  before_sampling <- tryCatch(
    BayesTools:::.bt_convergence_validate_monitor_names("thetaa", c("mu", "theta")),
    error = function(e) e
  )
  expect_s3_class(before_sampling, "BayesTools_parameter_not_found")
  expect_identical(before_sampling$alias, "thetaa")
  expect_identical(before_sampling$available, c("mu", "theta"))

  after_sampling <- tryCatch(
    BayesTools:::.bt_convergence_resolve_monitor("theta[3]", c("theta[1]", "theta[2]", "mu")),
    error = function(e) e
  )
  expect_s3_class(after_sampling, "BayesTools_parameter_not_found")
  expect_identical(after_sampling$alias, "theta[3]")
  expect_identical(after_sampling$available, c("mu", "theta[1]", "theta[2]"))

  # a missing parameter name is an argument error, not an internal failure
  expect_error(
    JAGS_check_convergence(structure(list(), class = "runjags"), add_parameters = NA),
    "add_parameters",
    fixed = TRUE
  )
  expect_error(
    JAGS_check_convergence(structure(list(), class = "runjags"), add_parameters = c("mu", NA)),
    "add_parameters",
    fixed = TRUE
  )
})

test_that("monitored deterministic nodes are structural only without stochastic ancestors", {

  set.seed(36)
  make <- function(d = .731){
    draws <- cbind(
      mu = stats::rnorm(200), theta = 2, eta = .5,
      L11 = -1, L21 = 2, d = d
    )
    colnames(draws)[4:5] <- c("L[1,1]", "L[2,1]")
    draws
  }
  priors <- list(mu = prior("normal", list(0, 1)))
  syntax <- paste(
    "model{",
    "  mu ~ dnorm(0, 1)",
    "  theta <- 2 * one",
    "  eta ~ dnorm(0, 1)",
    "  for(i in 1:2){ L[i, 1] = pow(-1, i) * i }",
    "  z ~ dnorm(0, 1)",
    "  d <- z",
    "}",
    sep = "\n"
  )
  mock <- function(chain_1, chain_2){
    .mock_convergence_fit(
      chain_1, chain_2, priors,
      add_parameters = c("theta", "eta", "L", "d"),
      model_syntax = syntax,
      data_names = "one"
    )
  }
  check <- function(fit, ...){
    JAGS_check_convergence(
      fit, max_Rhat = 1.2, min_ESS = 1, max_error = NULL,
      max_SD_error = NULL, ...
    )
  }

  # 'theta' reads only data and 'L' only its loop index: both are constants
  # of the model. 'eta' is stochastic, and 'd' copies the unmonitored
  # stochastic 'z': identical draws in every chain are a stuck sampler, not a
  # constant.
  fit <- mock(make(), make())
  expect_identical(
    .convergence_roles(fit),
    c(mu = "sampled", theta = "structural", eta = "sampled",
      "L[1,1]" = "structural", "L[2,1]" = "structural", d = "sampled")
  )
  informed <- check(fit)
  expect_false(informed)
  expect_equal(
    .convergence_states(informed),
    c(mu = "assessable", theta = "structural_constant", eta = "not_assessable",
      "L[1,1]" = "structural_constant", "L[2,1]" = "structural_constant",
      d = "not_assessable")
  )
  expect_match(attr(informed, "errors"), "not assessable for 'eta'", fixed = TRUE, all = FALSE)
  expect_match(attr(informed, "errors"), "not assessable for 'd'", fixed = TRUE, all = FALSE)
  expect_true(check(fit, allow_not_assessable = TRUE))

  # Without the data names, an undeclared name is not known to be constant.
  unknown <- .mock_convergence_fit(
    make(), make(), priors,
    add_parameters = c("theta", "eta", "L", "d"),
    model_syntax = syntax
  )
  expect_identical(unname(.convergence_roles(unknown)["theta"]), "sampled")
})

test_that("monitored fully observed data are structural constants", {

  set.seed(50)
  syntax <- "model{\n  for(i in 1:N){ y[i] ~ dnorm(mu + b * x[i], 1) }\n}"
  priors <- list(
    mu = prior("normal", list(0, 1)),
    b = prior("normal", list(0, 1))
  )
  x <- stats::rnorm(10)
  x_names <- paste0("x[", 1:10, "]")
  make <- function(){
    draws <- cbind(
      mu = stats::rnorm(200), b = stats::rnorm(200), N = 10,
      matrix(rep(x, each = 200), ncol = 10L)
    )
    colnames(draws)[4:13] <- x_names
    draws
  }
  check <- function(fit){
    JAGS_check_convergence(
      fit, max_Rhat = 1.2, min_ESS = 1, max_error = NULL, max_SD_error = NULL
    )
  }

  # 'N' and 'x' are data that the syntax reads but never defines.
  fit <- .mock_convergence_fit(
    make(), make(), priors, add_parameters = c("N", "x"),
    model_syntax = syntax, data_names = c("y", "x", "N")
  )
  expect_true(all(.convergence_roles(fit)[c("N", x_names)] == "structural"))
  result <- check(fit)
  expect_true(result)
  expect_true(all(.convergence_states(result)[c("N", x_names)] == "structural_constant"))

  # Without the data names, such names are not known to be constant.
  unknown <- .mock_convergence_fit(
    make(), make(), priors, add_parameters = c("N", "x"),
    model_syntax = syntax
  )
  expect_identical(unname(.convergence_roles(unknown)["N"]), "sampled")
  expect_false(check(unknown))
})

test_that("the JAGS syntax graph reads whole statements", {

  graph <- BayesTools:::.bt_jags_syntax_graph(paste(
    "data{ one <- 1; K <- length(v) }",
    "model{",
    "  for(i in 1:N){ y[i] ~ dnorm(m[i], 1); m[i] <- a + b * x[i] }",
    "  u[1:2] ~ ddirch(alpha) # d ~ dnorm(0, 1)",
    "  c.d ~ dgamma(1, 1) T(0.1,)",
    "  e = 2 * c.d",
    "  logit(p) <- q",
    "  q ~ dnorm(0, 1)",
    "  f <- h +",
    "    g",
    "  s[1,2] <- s0; s[2,1] <- 0",
    "}",
    sep = "\n"
  ))
  expect_setequal(graph$stochastic, c("y", "u", "c.d", "q"))
  expect_setequal(graph$data, c("one", "K"))
  expect_setequal(names(graph$deterministic), c("m", "e", "p", "f", "s"))
  expect_setequal(graph$deterministic$m, c("a", "b", "x"))
  expect_identical(graph$deterministic$e, "c.d")
  expect_identical(graph$deterministic$p, "q")
  expect_setequal(graph$deterministic$f, c("h", "g"))
  expect_identical(graph$deterministic$s, "s0")

  constants <- c("h", "g", "s0", "x")
  expect_identical(
    BayesTools:::.bt_jags_stochastic_ancestry(
      c("f", "s", "e", "p", "m", "one", "y"),
      graph = graph,
      constants = constants
    ),
    c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE)
  )

  # The two-character operators are single tokens: a comparison is part of
  # the right-hand side of its '<-' statement, never a relation of its own.
  expect_identical(
    BayesTools:::.bt_jags_syntax_tokens("a <- b <= c >= d == e != f && g || h < -1"),
    c("a", "<-", "b", "<=", "c", ">=", "d", "==", "e", "!=", "f", "&&", "g",
      "||", "h", "<", "-", "1")
  )
  comparison_graph <- function(statement, extra = ""){
    BayesTools:::.bt_jags_syntax_graph(paste0(
      "model{\n  theta ~ dnorm(0, 1)\n  ", statement, "\n", extra, "}"
    ))
  }
  for(statement in c(
    "flag <- theta == 0.25", "flag <- 0.25 == theta", "flag <- theta != 0",
    "flag <- 0 != theta", "flag <- theta >= 0", "flag <- 0 >= theta",
    "flag <- theta <= 0", "flag <- 0 <= theta",
    "flag <- (theta >= 0) && (theta <= 1)", "flag <- (theta > 0) || (theta == -1)"
  )){
    graph <- comparison_graph(statement)
    expect_identical(graph$stochastic, "theta")
    expect_identical(names(graph$deterministic), "flag")
    expect_setequal(graph$deterministic$flag, "theta")
  }
  # A lone '=' on the left of a comparison stays a relation, and an indexed
  # data operand of a comparison does not become a node of its own.
  graph <- comparison_graph("flag <- x[1] == theta", "  xm <- 2 * x[1]\n  r = theta\n")
  expect_setequal(names(graph$deterministic), c("flag", "xm", "r"))
  expect_setequal(graph$deterministic$flag, c("x", "theta"))
  expect_identical(graph$deterministic$xm, "x")
  expect_identical(graph$deterministic$r, "theta")

  # Convergence roles: a monitored comparison of a sampled node varies with
  # it, and data compared with data stays structural.
  roles <- function(statement, data_names = character(), extra = ""){
    syntax <- paste0(
      "model{\n  theta ~ dnorm(0, 1)\n  ", statement, "\n", extra, "}"
    )
    monitored <- c("flag", if(nzchar(extra)) "xm")
    coordinates <- BayesTools:::.bt_build_parameter_coordinates(
      columns = c("theta", monitored),
      prior_list = list(theta = prior("normal", list(0, 1))),
      add_parameters = monitored,
      model_syntax = syntax,
      data_names = data_names
    )
    stats::setNames(coordinates$convergence_role, coordinates$coordinate_name)
  }
  for(statement in c(
    "flag <- theta == 0.25", "flag <- 0 <= theta", "flag <- theta != 0",
    "flag <- 0 >= theta"
  )){
    expect_identical(roles(statement), c(theta = "sampled", flag = "sampled"))
  }
  expect_identical(
    roles("flag <- x[1] == 2", data_names = "x"),
    c(theta = "sampled", flag = "structural")
  )
  expect_identical(
    roles("flag <- x[1] == theta", data_names = "x", extra = "  xm <- 2 * x[1]\n"),
    c(theta = "sampled", flag = "sampled", xm = "structural")
  )
})

test_that("formula monitors are sampled when stochastic and derived otherwise", {

  set.seed(34)
  data <- data.frame(
    g = factor(rep(paste0("g", 1:4), each = 4)),
    x = stats::rnorm(16)
  )
  us_formula <- function(parameterization){
    JAGS_formula(
      ~ 1 + x + random(1 + x | g, name = "g", covariance = "us"),
      "mu",
      data,
      list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(g = random_block(
        sd = prior("normal", list(0, 1), list(0, Inf)),
        parameterization = parameterization,
        monitor = random_monitor(coefficients = TRUE, lkj_primitives = TRUE)
      ))
    )
  }
  stem <- "mu__xREx__g_xRE_"
  # 'theta' is a user-defined deterministic node and 'mu' the formula output;
  # both depend on sampled nodes and stay checked.
  roles <- function(formula){
    columns <- c(
      "mu_intercept", "mu_x", "mu__xREx__g_intercept", "mu__xREx__g_x",
      "theta", "mu[1]",
      paste0(stem, "CORx_L[", c("1,1", "2,1", "1,2", "2,2"), "]"),
      paste0(stem, "CORx_R[", c("1,1", "2,1", "1,2", "2,2"), "]"),
      paste0(stem, "CORx_lkj_u[1]"), paste0(stem, "CORx_lkj_cpc[1]"),
      paste0(stem, "Zx[1,", 1:2, "]"), paste0(stem, "COEFx[1,", 1:2, "]")
    )
    coordinates <- BayesTools:::.bt_build_parameter_coordinates(
      columns = columns,
      prior_list = formula$prior_list,
      formula_design = list(mu = formula$formula_design),
      add_parameters = c("theta", "mu", formula$add_parameters),
      model_syntax = paste0(
        "model{\n", formula$formula_syntax,
        "theta <- 2 * mu_intercept # not ~ stochastic\n}"
      ),
      data_names = names(formula$data)
    )
    stats::setNames(coordinates$convergence_role, coordinates$coordinate_name)
  }
  expected <- function(latent, coefficients){
    c(
      mu_intercept = "sampled", mu_x = "sampled",
      mu__xREx__g_intercept = "sampled", mu__xREx__g_x = "sampled",
      theta = "sampled", "mu[1]" = "sampled",
      stats::setNames(
        c("structural", "derived", "structural", "derived",
          "structural", "derived", "derived", "structural"),
        paste0(stem, c(
          paste0("CORx_L[", c("1,1", "2,1", "1,2", "2,2"), "]"),
          paste0("CORx_R[", c("1,1", "2,1", "1,2", "2,2"), "]")
        ))
      ),
      stats::setNames(c("sampled", "derived"), paste0(stem, c("CORx_lkj_u[1]", "CORx_lkj_cpc[1]"))),
      stats::setNames(rep(latent, 2L), paste0(stem, "Zx[1,", 1:2, "]")),
      stats::setNames(rep(coefficients, 2L), paste0(stem, "COEFx[1,", 1:2, "]"))
    )
  }

  # Non-centered: latent effects and LKJ primitives are sampled; correlation
  # matrices, their Cholesky factor, CPCs, and coefficients are derived, with
  # the unit diagonals and zero upper Cholesky entry structural.
  expect_identical(roles(us_formula("noncentered")), expected("sampled", "derived"))

  # Centered: coefficients are sampled and the latent effects derived.
  expect_identical(roles(us_formula("centered")), expected("derived", "sampled"))
})

test_that("random-effect SD spike-and-slab roles follow the declared inclusion", {

  set.seed(49)
  n <- 200
  data <- data.frame(id = factor(rep(paste0("g", 1:4), each = 2)))
  stem <- "mu__xREx__id_intercept"
  mock <- function(inclusion, indicator){
    formula <- JAGS_formula(
      ~ 1 + (1 || id), "mu", data,
      list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = prior_spike_and_slab(
        prior("normal", list(0, 1), list(0, Inf)),
        prior_inclusion = inclusion
      )))
    )
    chain <- function(){
      sd <- abs(stats::rnorm(n))
      states <- indicator(n)
      draws <- cbind(
        mu_intercept = stats::rnorm(n),
        states,
        stats::runif(n),
        sd * states,
        sd,
        matrix(stats::rnorm(4L * n), ncol = 4L)
      )
      colnames(draws) <- c(
        "mu_intercept",
        paste0(stem, c("_indicator", "_inclusion", "", "_variable")),
        paste0("mu__xREx__id_xRE_Zx[", 1:4, ",1]")
      )
      draws
    }
    fit <- list(
      mcmc = coda::mcmc.list(coda::mcmc(chain()), coda::mcmc(chain())),
      summary.pars = list(mutate = NULL)
    )
    class(fit) <- c("runjags", "BayesTools_fit")
    attr(fit, "prior_list") <- formula$prior_list
    attr(fit, "formula_design") <- list(mu = formula$formula_design)
    attr(fit, "add_parameters") <- formula$add_parameters
    attr(fit, "model_syntax") <- paste0("model{\n", formula$formula_syntax, "\n}")
    fit <- BayesTools:::.bt_attach_parameter_map(fit, data_names = names(formula$data))
    fit <- BayesTools:::.bt_attach_draw_geometry(fit)
    BayesTools:::.bt_attach_fit_contract(fit)
  }
  roles <- function(fit){
    .convergence_roles(fit)[paste0(stem, c("_indicator", "_inclusion", "_variable"))]
  }
  check <- function(fit){
    JAGS_check_convergence(
      fit, max_Rhat = 1.2, min_ESS = 1, max_error = NULL, max_SD_error = NULL,
      check_indicators = TRUE
    )
  }

  # An inclusion probability declared as a point at 1 fixes the indicator in
  # every draw: a structural constant, as for the fixed-effect spike and slab.
  point <- mock(prior("spike", list(1)), function(n) rep(1, n))
  expect_identical(
    unname(roles(point)),
    c("structural", "structural", "sampled")
  )
  point_check <- check(point)
  expect_true(point_check)
  expect_identical(
    unname(.convergence_states(point_check)[paste0(stem, "_indicator")]),
    "structural_constant"
  )

  # A random inclusion probability keeps the indicator checked and the
  # probability auxiliary; a stuck indicator is not assessable.
  beta <- mock(prior("beta", list(1, 1)), function(n) rep(1, n))
  expect_identical(
    unname(roles(beta)),
    c("indicator", "auxiliary", "sampled")
  )
  beta_check <- check(beta)
  expect_false(beta_check)
  expect_identical(
    unname(.convergence_states(beta_check)[paste0(stem, "_indicator")]),
    "not_assessable"
  )
})

test_that("one chain warns and retains its other convergence criteria", {

  set.seed(44)
  chain <- cbind(mu = rnorm(2000))
  fit <- .mock_convergence_fit(chain, prior_list = list(mu = prior("normal", list(0, 1))))
  expect_warning(result <- JAGS_check_convergence(
    fit,
    max_Rhat = 1.05,
    min_ESS = 100,
    max_error = NULL,
    max_SD_error = NULL
  ), "Only one chain was run. R-hat cannot be computed; checking the remaining enabled convergence criteria.", fixed = TRUE)

  expect_true(result)
  expect_null(attr(result, "errors"))
  expect_true(is.na(attr(result, "diagnostics")$Rhat))
  expect_equal(attr(result, "diagnostics")$state, "assessable")
  expect_warning(failed <- JAGS_check_convergence(
    fit,
    min_ESS = 1e6,
    max_error = NULL,
    max_SD_error = NULL
  ), "Only one chain was run", fixed = TRUE)
  expect_false(failed)
  expect_match(attr(failed, "errors"), "ESS")
})

test_that("fail_fast stops computing after a failed parameter", {

  set.seed(46)
  chain <- cbind(mu = rnorm(100), tau = rnorm(100))
  fit <- .mock_convergence_fit(chain, chain, list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(0, 1))
  ))
  calls <- 0L
  original <- BayesTools:::.bt_convergence_parameter_diagnostics
  testthat::local_mocked_bindings(
    .bt_convergence_parameter_diagnostics = function(...){
      calls <<- calls + 1L
      original(...)
    },
    .package = "BayesTools"
  )
  result <- JAGS_check_convergence(
    fit, min_ESS = 1e6, max_Rhat = NULL, max_error = NULL,
    max_SD_error = NULL, fail_fast = TRUE
  )
  expect_false(result)
  expect_identical(calls, 1L)
  expect_equal(attr(result, "diagnostics")$state,
               c("assessable", "not_checked"))
  expect_length(attr(result, "errors"), 1L)
  complete <- JAGS_check_convergence(
    fit, min_ESS = 1e6, max_Rhat = NULL, max_error = NULL,
    max_SD_error = NULL, fail_fast = FALSE
  )
  expect_false(complete)
  expect_identical(calls, 3L)
  expect_equal(attr(complete, "diagnostics")$state,
               c("assessable", "assessable"))
})

test_that("computed failures are reported when another diagnostic of the row is not assessable (F277)", {

  diagnostics <- data.frame(
    parameter = c("a", "b"),
    state = c("not_assessable", "assessable"),
    Rhat = c(Inf, 1.01),
    ESS = c(NA_real_, 900),
    MCMC_error = c(NA_real_, 0.001),
    MCMC_SD_error = c(NA_real_, 0.01),
    stringsAsFactors = FALSE
  )
  failures <- function(allow_not_assessable){
    BayesTools:::.bt_convergence_failures(
      diagnostics, max_Rhat = 1.05, min_ESS = 500, max_error = 0.01,
      max_SD_error = 0.05, allow_not_assessable = allow_not_assessable
    )
  }

  # the computed R-hat of 'a' fails whatever allow_not_assessable says
  expect_identical(failures(TRUE), "R-hat Inf for 'a' is larger than the set target (1.05).")
  strict <- failures(FALSE)
  expect_length(strict, 2L)
  expect_match(strict[[1L]], "ESS, MCMC error, MCMC SD error diagnostics are not assessable for 'a'.", fixed = TRUE)
  expect_identical(strict[[2L]], failures(TRUE))

  # a row whose computed diagnostics pass reports nothing when undefined ones are allowed
  diagnostics$Rhat[1L] <- 1.01
  expect_identical(failures(TRUE), character())
})
