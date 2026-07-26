skip_if_not_test_profile("unit")

test_that("factor prior samplers use one RNG stream across coefficients", {

  n_samples <- 40
  null_prior <- prior("spike", list(0), prior_weights = 1)

  for(contrast in c("treatment", "independent")){
    factor_prior <- prior_factor(
      "normal",
      list(0, 0.2),
      contrast = contrast,
      prior_weights = 3
    )
    attr(factor_prior, "levels") <- 3
    n_coefficients <- .get_prior_factor_levels(factor_prior)
    priors <- list(null_prior, factor_prior)

    set.seed(611)
    expected_counts <- .prior_mixture_sample_counts(c(.25, .75), n_samples)
    expected_mixture <- lapply(seq_len(n_coefficients), function(i){
      .mix_priors.simple(
        priors,
        paste0("fac[", i, "]"),
        seed = NULL,
        n_samples = n_samples,
        sample_counts = expected_counts
      )
    })
    expected_mixture <- do.call(cbind, lapply(expected_mixture, as.numeric))
    expected_mixture_state <- .Random.seed

    mixture <- .mix_priors.factor(
      priors,
      "fac",
      seed = 611,
      n_samples = n_samples
    )
    mixture_state <- .Random.seed

    expect_equal(
      matrix(as.numeric(mixture), nrow = n_samples),
      expected_mixture
    )
    expect_identical(mixture_state, expected_mixture_state)
    expect_false(identical(
      as.numeric(mixture[, 1]),
      as.numeric(mixture[, 2])
    ))

    set.seed(612)
    expected_single <- replicate(
      n_coefficients,
      rng(factor_prior, n_samples, transform_factor_samples = FALSE)
    )
    expected_single_state <- .Random.seed

    single <- .as_mixed_priors.factor(
      factor_prior,
      "fac",
      seed = 612,
      n_samples = n_samples
    )
    single_state <- .Random.seed

    expect_equal(
      matrix(as.numeric(single), nrow = n_samples),
      expected_single
    )
    expect_identical(single_state, expected_single_state)
    expect_false(identical(
      as.numeric(single[, 1]),
      as.numeric(single[, 2])
    ))
  }
})

test_that("prior sampler orchestrators continue the caller RNG stream", {

  n_samples <- 30
  theta_prior <- prior("normal", list(0, 1))
  eta_prior <- prior("normal", list(2, 0.5))

  set.seed(621)
  expected_theta <- .as_mixed_priors.simple(
    theta_prior,
    "theta",
    seed = NULL,
    n_samples = n_samples
  )
  expected_eta <- .as_mixed_priors.simple(
    eta_prior,
    "eta",
    seed = NULL,
    n_samples = n_samples
  )
  expected_single_state <- .Random.seed

  set.seed(621)
  single <- .as_mixed_priors(
    list(theta = theta_prior, eta = eta_prior),
    seed = NULL,
    n_samples = n_samples
  )
  single_state <- .Random.seed

  expect_equal(as.numeric(single$theta), as.numeric(expected_theta))
  expect_equal(as.numeric(single$eta), as.numeric(expected_eta))
  expect_identical(single_state, expected_single_state)

  set.seed(622)
  expected_theta <- .mix_priors.simple(
    list(theta_prior),
    "theta",
    seed = NULL,
    n_samples = n_samples
  )
  expected_eta <- .mix_priors.simple(
    list(eta_prior),
    "eta",
    seed = NULL,
    n_samples = n_samples
  )
  expected_mixture_state <- .Random.seed

  set.seed(622)
  mixture <- .mix_priors(
    list(theta = list(theta_prior), eta = list(eta_prior)),
    seed = NULL,
    n_samples = n_samples
  )
  mixture_state <- .Random.seed

  expect_equal(as.numeric(mixture$theta), as.numeric(expected_theta))
  expect_equal(as.numeric(mixture$eta), as.numeric(expected_eta))
  expect_identical(mixture_state, expected_mixture_state)
})

test_that("seeded spike-and-slab priors do not restart the slab RNG", {

  n_samples <- 80
  variable_prior <- prior("normal", list(0, 1))
  inclusion_prior <- prior("point", list(0.5))
  spike_and_slab <- prior_spike_and_slab(
    variable_prior,
    prior_inclusion = inclusion_prior
  )

  set.seed(631)
  expected_inclusion <- stats::rbinom(
    n_samples,
    size = 1,
    prob = rng(inclusion_prior, n_samples)
  )
  expected_variable <- rng(
    variable_prior,
    n_samples,
    transform_factor_samples = FALSE
  )
  expected_state <- .Random.seed

  samples <- .as_mixed_priors.spike_and_slab(
    spike_and_slab,
    "theta",
    seed = 631,
    n_samples = n_samples
  )
  samples_state <- .Random.seed

  expect_identical(attr(samples, "models_ind"), expected_inclusion)
  expect_equal(
    as.numeric(samples),
    as.numeric(expected_variable * expected_inclusion)
  )
  expect_identical(samples_state, expected_state)
})
