skip_if_not_test_profile("unit")

test_that("R133 raw affine dependence refuses bounded bridge coordinates before sampling", {
  x2 <- seq(1, 4, length.out = 200)
  dependent <- cbind(x1 = 2 * x2 + 1, x2 = x2)
  sampler_calls <- 0L
  testthat::local_mocked_bindings(bridge_sampler = function(...){
    sampler_calls <<- sampler_calls + 1L
    structure(list(logml = -1, niter = 1L, method = "normal", mcse_logml = .1), class = "bridge")
  }, .package = "bridgesampling")
  for(bounds in list(c(0, 10), c(0, Inf), c(-Inf, 10), c(-Inf, Inf))){
    posterior <- dependent
    attr(posterior, "lb") <- c(x1 = bounds[1], x2 = bounds[1])
    attr(posterior, "ub") <- c(x1 = bounds[2], x2 = bounds[2])
    expect_error(.bt_JAGS_bridge_check_coordinate_rank(posterior), "rank-deficient (rank 1 of 2", fixed = TRUE)
    expect_error(JAGS_bridgesampling(coda::as.mcmc(dependent),
      log_posterior = function(parameters, data) 0, data = list(),
      prior_list = list(x1 = prior("normal", list(0, 1), truncation = list(bounds[1], bounds[2]))),
      add_parameters = "x2", add_bounds = list(lb = c(x2 = bounds[1]), ub = c(x2 = bounds[2]))),
      "rank-deficient (rank 1 of 2", fixed = TRUE)
  }
  expect_identical(sampler_calls, 0L)
  full_rank <- dependent
  full_rank[, "x2"] <- dependent[c(seq(2, 200, 2), seq(1, 199, 2)), "x2"]
  expect_s3_class(JAGS_bridgesampling(coda::as.mcmc(full_rank),
    log_posterior = function(parameters, data) 0, data = list(),
    prior_list = list(x1 = prior("normal", list(0, 1), truncation = list(0, 10))),
    add_parameters = "x2", add_bounds = list(lb = c(x2 = 0), ub = c(x2 = 10))), "BayesTools_marglik")
  expect_identical(sampler_calls, 1L)
})

test_that("R133 bridge rank checks preserve transformed-only and ordinary controls", {
  t <- seq(1, 2, length.out = 200)
  transformed_dependence <- cbind(x1 = exp(t), x2 = exp(2 * t))
  attr(transformed_dependence, "lb") <- c(x1 = 0, x2 = 0)
  attr(transformed_dependence, "ub") <- c(x1 = Inf, x2 = Inf)
  expect_error(.bt_JAGS_bridge_check_coordinate_rank(transformed_dependence), "rank-deficient (rank 1 of 2", fixed = TRUE)
  full_rank <- cbind(x1 = t, x2 = sin(t))
  expect_true(.bt_JAGS_bridge_check_coordinate_rank(full_rank))
  expect_true(.bt_JAGS_bridge_check_coordinate_rank(cbind(x1 = t, constant = 1)))
  expect_error(.bt_JAGS_bridge_check_coordinate_rank(full_rank[1:2, ]), "2 posterior draws span only rank 1 of 2", fixed = TRUE)
  expect_error(.bt_JAGS_bridge_check_coordinate_rank(cbind(x1 = t, x2 = 2 * t)), "rank-deficient (rank 1 of 2", fixed = TRUE)
})

test_that("R143 independent bridge coordinates retain rank at large magnitudes", {

  independent <- cbind(x1 = c(1, 2, 3, 4), x2 = c(1, 4, 2, 3))
  # These centered columns are independent; their first two rows have
  # determinant -3. This moderate-scale reference also uses the existing QR.
  expect_identical(qr(scale(independent), tol = 1e-7)$rank, 2L)
  sampler_calls <- 0L
  sampler_arguments <- NULL
  testthat::local_mocked_bindings(bridge_sampler = function(...){
    sampler_calls <<- sampler_calls + 1L
    sampler_arguments <<- list(...)
    structure(list(logml = -1, niter = 1L, method = "normal", mcse_logml = .1), class = "bridge")
  }, .package = "bridgesampling")

  for(coordinate_scale in c(1, 1e160)){
    coordinates <- independent * coordinate_scale
    original <- coordinates
    expect_true(all(is.finite(coordinates)))
    expect_error(.bt_JAGS_bridge_check_matrix_rank(coordinates), NA)
    expect_identical(coordinates, original)
    for(lower in c(0, -Inf)){
      posterior <- coordinates
      attr(posterior, "lb") <- c(x1 = lower, x2 = lower)
      attr(posterior, "ub") <- c(x1 = Inf, x2 = Inf)
      original_posterior <- posterior
      expect_error(.bt_JAGS_bridge_check_coordinate_rank(posterior), NA)
      expect_identical(posterior, original_posterior)
      sampler_calls <- 0L
      sampler_arguments <- NULL
      expect_error(JAGS_bridgesampling(coda::as.mcmc(coordinates),
        log_posterior = function(parameters, data) 0, data = list(),
        prior_list = list(x1 = prior("normal", list(0, coordinate_scale), truncation = list(lower, Inf))),
        add_parameters = "x2", add_bounds = list(lb = c(x2 = lower), ub = c(x2 = Inf))), NA)
      expect_identical(sampler_calls, 1L)
      if(!is.null(sampler_arguments)){
        expect_identical(sampler_arguments$samples, coda::as.mcmc(original_posterior))
        expect_identical(sampler_arguments$lb, attr(original_posterior, "lb"))
        expect_identical(sampler_arguments$ub, attr(original_posterior, "ub"))
      }
      expect_identical(coordinates, original)
    }
  }

  # Numerical conditioning preserves the accepted QR tolerance and the
  # exclusion of genuinely constant or nonfinite columns.
  for(epsilon in c(1e-8, 1e-5)){
    near_affine <- cbind(x1 = independent[, "x1"],
      x2 = independent[, "x1"] + epsilon * independent[, "x2"])
    expected_rank <- if(epsilon == 1e-8) 1L else 2L
    expect_identical(qr(scale(near_affine), tol = 1e-7)$rank, expected_rank)
    for(coordinate_scale in c(1, 1e160)){
      if(expected_rank == 1L){
        condition <- tryCatch(.bt_JAGS_bridge_check_matrix_rank(near_affine * coordinate_scale),
          error = identity)
        expect_s3_class(condition, "BayesTools_marglik_rank_deficiency")
        if(inherits(condition, "error")){
          expect_match(conditionMessage(condition), "rank-deficient (rank 1 of 2", fixed = TRUE)
          expect_null(conditionCall(condition))
        }
      }else{
        expect_error(.bt_JAGS_bridge_check_matrix_rank(near_affine * coordinate_scale), NA)
      }
    }
  }
  expect_true(.bt_JAGS_bridge_check_matrix_rank(cbind(x1 = independent[, "x1"],
    zero = 0, constant = 7, nonfinite = c(1, 2, Inf, 4))))
})

test_that("R143 tiny affine bridge coordinates refuse before bound transformations", {

  x2 <- seq(1, 4, length.out = 200)
  moderate <- cbind(x1 = 2 * x2 + 1, x2 = x2)
  expect_identical(qr(scale(moderate), tol = 1e-7)$rank, 1L)
  # Log transformation removes the affine relation, so only the transformed
  # check would wrongly accept this supplied deterministic coordinate.
  expect_identical(qr(scale(log(moderate)), tol = 1e-7)$rank, 2L)
  sampler_calls <- 0L
  testthat::local_mocked_bindings(bridge_sampler = function(...){
    sampler_calls <<- sampler_calls + 1L
    structure(list(logml = -1, niter = 1L, method = "normal", mcse_logml = .1), class = "bridge")
  }, .package = "bridgesampling")

  for(coordinate_scale in c(1, 1e-200)){
    scaled_x2 <- coordinate_scale * x2
    coordinates <- cbind(x1 = 2 * scaled_x2 + coordinate_scale, x2 = scaled_x2)
    original <- coordinates
    expect_true(all(is.finite(coordinates)))
    expect_identical(qr(scale(log(coordinates)), tol = 1e-7)$rank, 2L)
    expect_error(.bt_JAGS_bridge_check_matrix_rank(coordinates),
      "rank-deficient (rank 1 of 2", fixed = TRUE)
    expect_identical(coordinates, original)
    for(lower in c(0, -Inf)){
      posterior <- coordinates
      attr(posterior, "lb") <- c(x1 = lower, x2 = lower)
      attr(posterior, "ub") <- c(x1 = Inf, x2 = Inf)
      original_posterior <- posterior
      expect_error(.bt_JAGS_bridge_check_coordinate_rank(posterior),
        "rank-deficient (rank 1 of 2", fixed = TRUE)
      expect_identical(posterior, original_posterior)
      sampler_calls <- 0L
      expect_error(JAGS_bridgesampling(coda::as.mcmc(coordinates),
        log_posterior = function(parameters, data) 0, data = list(),
        prior_list = list(x1 = prior("normal", list(0, coordinate_scale), truncation = list(lower, Inf))),
        add_parameters = "x2", add_bounds = list(lb = c(x2 = lower), ub = c(x2 = Inf))),
        "rank-deficient (rank 1 of 2", fixed = TRUE)
      expect_identical(sampler_calls, 0L)
      expect_identical(coordinates, original)
    }
  }
})

test_that("bridge rank refusals expose diagnostic classes and exact context", {

  sampler_calls <- 0L
  testthat::local_mocked_bindings(bridge_sampler = function(...){
    sampler_calls <<- sampler_calls + 1L
    structure(list(logml = -1, niter = 1L, method = "normal"), class = "bridge")
  }, .package = "bridgesampling")
  public_rank <- function(coordinates){
    tryCatch(JAGS_bridgesampling(
      coda::as.mcmc(coordinates),
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = setNames(rep(list(prior("normal", list(0, 1))),
        ncol(coordinates)), colnames(coordinates))
    ), error = identity)
  }
  t <- seq(-1, 1, length.out = 20L)
  cases <- list(
    list(coordinates = cbind(x1 = c(-1, 1), x2 = c(1, -1)),
      reason = "insufficient_draws", n_draws = 2L,
      message = paste0("The bridge-sampling target was rejected by diagnostics: ",
        "2 posterior draws span only rank 1 of 2 varying bridge coordinates. ",
        "Increase the number of posterior draws, for example with a larger ",
        "'sample' in JAGS_fit() or with JAGS_extend().")),
    list(coordinates = cbind(x1 = t, x2 = 2 * t),
      reason = "dependent_coordinates", n_draws = 20L,
      message = paste0("The bridge-sampling target was rejected by diagnostics: ",
        "the posterior draws of the bridge coordinates are rank-deficient ",
        "(rank 1 of 2 varying coordinates; linearly dependent coordinate(s): 'x2'). ",
        "Only free stochastic nodes whose prior density 'log_posterior' adds ",
        "can be supplied through 'add_parameters'; reconstruct deterministic ",
        "nodes, such as row-shaped SD sources computed from other parameters, ",
        "with parameter_source(values = ..., inputs = ...) or inside 'log_posterior'."))
  )
  for(case in cases){
    condition <- public_rank(case$coordinates)
    expect_s3_class(condition, "BayesTools_marglik_rank_deficiency")
    expect_s3_class(condition, "BayesTools_marglik_input")
    expect_identical(condition$reason, case$reason)
    expect_identical(condition$rank, 1L)
    expect_identical(condition$n_coordinates, 2L)
    expect_identical(condition$n_draws, case$n_draws)
    expect_identical(condition$dependent_coordinates, "x2")
    expect_identical(conditionMessage(condition), case$message)
    expect_null(conditionCall(condition))
  }
  expect_identical(sampler_calls, 0L)

  # The fields use the filtered diagnostic matrix, not the supplied width.
  coordinates <- cbind(x1 = t, x2 = 2 * t, x3 = 3 * t, x4 = 4 * t,
    x5 = 5 * t, constant = 7, nonfinite = c(Inf, rep(1, 19)))
  condition <- tryCatch(BayesTools:::.bt_JAGS_bridge_check_matrix_rank(coordinates),
    error = identity)
  expect_identical(condition$n_coordinates, 5L)
  expect_identical(condition$dependent_coordinates, paste0("x", 2:5))
  expect_match(conditionMessage(condition), "'x2', 'x3', 'x4', ...", fixed = TRUE)
  colnames(coordinates) <- NULL
  condition <- tryCatch(BayesTools:::.bt_JAGS_bridge_check_matrix_rank(coordinates),
    error = identity)
  expect_identical(condition$dependent_coordinates, character())

  # Reasons describe the observed sample diagnostic, not an exclusive cause.
  a <- c(0, 1, 0, 1)
  b <- c(0, 0, 1, 1)
  mixed <- cbind(a = a, b = b, c = a + b, d = 2 * a + b, e = a + 2 * b)
  condition <- public_rank(mixed)
  expect_s3_class(condition, "BayesTools_marglik_rank_deficiency")
  expect_identical(condition$reason, "insufficient_draws")
  expect_identical(condition$rank, 2L)
  expect_identical(condition$n_coordinates, 5L)
  expect_identical(condition$n_draws, 4L)
  expect_identical(condition$dependent_coordinates, c("c", "d", "e"))
  duplicated <- cbind(x1 = rep(c(0, 1), each = 3L),
    x2 = rep(c(0, 2), each = 3L))
  condition <- public_rank(duplicated)
  expect_s3_class(condition, "BayesTools_marglik_rank_deficiency")
  expect_identical(condition$reason, "dependent_coordinates")
  expect_identical(condition$rank, 1L)
  expect_identical(condition$n_coordinates, 2L)
  expect_identical(condition$n_draws, 6L)
  expect_identical(condition$dependent_coordinates, "x2")
  expect_identical(sampler_calls, 0L)
})
