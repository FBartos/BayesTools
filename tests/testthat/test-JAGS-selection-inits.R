skip_if_not_test_profile("unit")

test_that("cumulative selection initializations stay strictly positive", {

  expect_error(wf_cumulative(c(1e-300, 1e-300)),
               "The 'alpha' must be equal or higher than 0.01.", fixed = TRUE)
  cumulative <- prior_weightfunction(
    side = "one-sided",
    steps = .05,
    weights = wf_cumulative(c(1, 1))
  )

  direct_inits <- JAGS_get_inits(
    prior_list = list(bias = cumulative),
    chains = 1,
    seed = 1
  )
  for(chain_inits in direct_inits){
    expect_true(is.finite(chain_inits$omega_ratio))
    expect_true(chain_inits$omega_ratio > 0)
    expect_true(chain_inits$omega_ratio < 1)
  }

  mixture <- prior_mixture(list(
    prior_none(),
    cumulative
  ))
  mixture_inits <- JAGS_get_inits(
    prior_list = list(bias = mixture),
    chains = 1,
    seed = 1
  )
  for(chain_inits in mixture_inits){
    expect_true(is.finite(chain_inits$omega_ratio_component_2))
    expect_true(chain_inits$omega_ratio_component_2 > 0)
    expect_true(chain_inits$omega_ratio_component_2 < 1)
  }

  set.seed(1)
  component_inits <- BayesTools:::.JAGS_init.weightfunction(
    cumulative,
    component_id = 2
  )
  expect_true(is.finite(component_inits$omega_ratio_component_2))
  expect_true(component_inits$omega_ratio_component_2 > 0)
  expect_true(component_inits$omega_ratio_component_2 < 1)
})

test_that("independent weight functions emit JAGS inits for free bins", {

  omega <- prior_weightfunction(
    side = "one-sided",
    steps = c(.025, .05),
    weights = wf_independent(prior("gamma", list(2, 1)))
  )
  log_omega <- prior_weightfunction(
    side = "one-sided",
    steps = .05,
    weights = wf_independent(prior("normal", list(0, 1)), scale = "log_omega")
  )

  set.seed(1)
  omega_inits <- JAGS_get_inits(list(bias = omega), chains = 1, seed = 1)[[1]]
  expect_equal(length(omega_inits$omega), 3L)
  expect_true(is.na(omega_inits$omega[[1L]]))
  expect_true(all(is.finite(omega_inits$omega[-1L])))
  expect_true(all(omega_inits$omega[-1L] > 0))

  log_inits <- JAGS_get_inits(list(bias = log_omega), chains = 1, seed = 1)[[1]]
  expect_equal(length(log_inits$log_omega), 2L)
  expect_true(is.na(log_inits$log_omega[[1L]]))
  expect_true(is.finite(log_inits$log_omega[[2L]]))

  # In a mixture the component is initialized on its local node:
  # `omega_component_2` is deterministic and expanded onto the global cuts.
  mixture <- prior_mixture(list(prior_none(), omega))
  mixture_inits <- JAGS_get_inits(list(bias = mixture), chains = 1, seed = 1)[[1]]
  expect_null(mixture_inits$omega_component_2)
  expect_equal(length(mixture_inits$omega_local_component_2), 3L)
  expect_true(is.na(mixture_inits$omega_local_component_2[[1L]]))
  expect_true(all(is.finite(mixture_inits$omega_local_component_2[-1L])))
})

test_that("Gamma initialization fallback preserves representable median proportions", {

  set.seed(1)
  expect_warning(
    initialization <- BayesTools:::.JAGS_positive_gamma_initialization(
      rep(.001, 4),
      "test latent Gamma variables"
    ),
    "deterministic, order-one rescaling"
  )
  expect_identical(initialization, rep(1, 4))

  set.seed(1)
  expect_warning(
    expect_error(
      BayesTools:::.JAGS_positive_gamma_initialization(
        c(1e-323, 2e-323),
        "test latent Gamma variables"
      ),
      "medians are not representable"
    ),
    "deterministic, order-one rescaling"
  )
})

test_that("a two-bin cumulative initialization falls back to the median of the weight (F222)", {

  weights  <- function(alpha) prior_weightfunction("one-sided", .05, wf_cumulative(alpha))
  fallback <- "deterministic, order-one rescaling"
  inits    <- function(alpha, seed){
    JAGS_get_inits(list(omega = weights(alpha)), chains = 1, seed = seed)[[1L]]$omega_ratio
  }

  # the seeds at which the drawn Gamma pair rounds the ratio to 0 or 1 refused
  # their initial value (147 of 400 at the constructor minimum): the interior
  # median of Beta(a, a) with the existing warning instead
  for(seed in c(1, 2, 3, 5, 8, 11, 12, 14)){
    expect_warning(value <- inits(c(.01, .01), seed), fallback, info = seed)
    expect_equal(value, .5, tolerance = 1e-12, info = seed)
  }
  # a seed with an interior ratio keeps the drawn value and warns of nothing
  expect_no_warning(value <- inits(c(.01, .01), 4))
  expect_gt(abs(value - .5), 1e-6)

  # every seed returns an interior value and warns exactly where it fell back
  for(alpha in list(c(.01, .01), c(.02, .02), c(.05, .05), c(.1, .1))){
    fallbacks <- 0L
    values <- vapply(1:200, function(seed){
      withCallingHandlers(
        inits(alpha, seed),
        warning = function(w){
          expect_match(conditionMessage(w), fallback, fixed = TRUE)
          fallbacks <<- fallbacks + 1L
          invokeRestart("muffleWarning")
        }
      )
    }, numeric(1))
    expect_true(all(values > 0 & values < 1), info = paste(alpha, collapse = ", "))
    expect_identical(fallbacks, sum(abs(values - .5) < 1e-12), info = paste(alpha, collapse = ", "))
  }

  # the component of a mixture takes the same fallback
  mixture <- prior_mixture(list(prior_none(), weights(c(.01, .01))))
  for(seed in 1:20){
    value <- suppressWarnings(
      JAGS_get_inits(list(bias = mixture), chains = 1, seed = seed)[[1L]]$omega_ratio_component_2
    )
    expect_true(value > 0 && value < 1)
  }

  # unequal shapes start at the median of Beta(alpha[2], alpha[1]) ...
  expect_warning(value <- inits(c(.01, .02), 1), fallback)
  expect_equal(value, stats::qbeta(.5, .02, .01), tolerance = 1e-12)
  # ... and only a median that is not an interior value keeps the refusal
  for(alpha in list(c(.01, 1), c(.01, 100))){
    expect_warning(
      expect_error(inits(alpha, 1), "is not representable as a finite interior cumulative weight", fixed = TRUE),
      fallback
    )
  }
})

test_that("equal Gamma shapes retain the all-underflow median fallback", {

  set.seed(1)
  expect_warning(
    initialization <- BayesTools:::.JAGS_positive_gamma_initialization(
      rep(1e-300, 2), "test latent Gamma variables"
    ),
    "deterministic, order-one rescaling"
  )
  expect_identical(initialization, rep(1, 2))
})


# Initial values and model syntax are generated by separate code paths, and a
# name that drifts between them is not a soft failure: JAGS aborts with
# "Dimension mismatch in values supplied for <node>" whenever the node it does
# find has a different length. These helpers assert the general invariant --
# every initialized node is stochastic in the syntax, and lengths agree --
# rather than any particular set of names.
.selection_init_stochastic_nodes <- function(syntax){

  lines <- trimws(unlist(strsplit(paste(syntax, collapse = "\n"), "\n")))
  lines <- lines[grepl("~", lines, fixed = TRUE)]
  targets <- trimws(sub("~.*$", "", lines))
  targets <- targets[nzchar(targets)]

  indexed <- grepl("[", targets, fixed = TRUE)
  out <- list()
  for(target in targets[!indexed]){
    out[[target]] <- max(c(out[[target]], 1L))
  }
  for(target in targets[indexed]){
    parts <- strsplit(target, "[", fixed = TRUE)[[1]]
    name  <- parts[[1]]
    index <- suppressWarnings(as.integer(
      strsplit(parts[[2]], "]", fixed = TRUE)[[1]][[1]]
    ))
    if(!is.na(index)){
      out[[name]] <- max(c(out[[name]], index))
    }
  }

  out
}

expect_selection_inits_match_syntax <- function(prior_list){

  syntax <- JAGS_add_priors("model{\n}", prior_list)
  inits  <- JAGS_get_inits(prior_list = prior_list, chains = 1, seed = 1)[[1]]
  nodes  <- .selection_init_stochastic_nodes(syntax)

  for(name in setdiff(names(inits), c(".RNG.seed", ".RNG.name"))){
    expect_true(
      name %in% names(nodes),
      info = paste0("initialized node '", name, "' is not stochastic in the model syntax")
    )
    if(name %in% names(nodes)){
      expect_gte(nodes[[name]], length(inits[[name]]))
    }
  }

  invisible(inits)
}


test_that("selection initial values name stochastic nodes of the matching length", {

  independent_omega <- prior_weightfunction(
    "one-sided", c(.05, .10), wf_independent(prior("gamma", list(shape = 9, rate = 3)))
  )
  independent_log <- prior_weightfunction(
    "one-sided", c(.025),
    wf_independent(prior("normal", list(mean = log(1.5), sd = .15)), "log_omega")
  )

  # standalone priors: no cut expansion
  expect_selection_inits_match_syntax(list(bias = independent_omega))
  expect_selection_inits_match_syntax(list(bias = independent_log))

  # mixtures expand every component onto the union of all cut points, so the
  # public `omega_component_*` arrays are longer than the local ones
  mixture <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one-sided", c(.025, .05), wf_cumulative(c(1, 2, 3)), prior_weights = 1),
    independent_omega,
    independent_log,
    prior_weightfunction("two-sided", c(.05), wf_fixed(c(1, .4)), prior_weights = 1)
  ))
  inits <- expect_selection_inits_match_syntax(list(bias = mixture))

  # the initialized omega node is the local one, never the expanded target
  expect_true("omega_local_component_3" %in% names(inits))
  expect_false("omega_component_3" %in% names(inits))
  expect_true("log_omega_component_4" %in% names(inits))
})
