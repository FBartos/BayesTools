skip_if_not_test_profile("cran")

# The only tests CRAN runs (tests/testthat.R): the package and its native code
# load, every registered method resolves, and the JAGS module compiles and
# samples. Full behavior is covered by the other profiles.

test_that("native routines and S3 methods are registered", {

  dll <- getLoadedDLLs()[["BayesTools"]]
  expect_false(is.null(dll))
  registered <- getDLLRegisteredRoutines(dll)[[".Call"]]
  expect_setequal(.BayesTools_native_symbols(), names(registered))

  # columns: generic, class, method, and the package of a foreign generic
  methods  <- getNamespaceInfo("BayesTools", "S3methods")
  resolved <- vapply(seq_len(nrow(methods)), function(i){
    owner <- if(is.na(methods[i, 4L])) "BayesTools" else methods[i, 4L]
    !is.null(getS3method(methods[i, 1L], methods[i, 2L], optional = TRUE,
                         envir = asNamespace(owner)))
  }, logical(1))
  expect_identical(methods[!resolved, 3L], character())

  exports <- getNamespaceExports("BayesTools")
  missing <- exports[!vapply(exports, exists, logical(1),
                             envir = asNamespace("BayesTools"), inherits = FALSE)]
  expect_identical(missing, character())
})

test_that("the JAGS module compiles, samples, and bridges a minimal model", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  skip_if_not_installed("bridgesampling")
  expect_true(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = TRUE)))

  # the priors and the LKJ block use every distribution and function of the module
  prior_list <- list(
    theta = prior("moment", list(mode = .5)),
    phi   = prior("invmoment", list(mode = .5, df = 6)),
    sigma = prior("invgamma", list(2, 1))
  )
  lkj <- JAGS_lkj_corr_cholesky("Omega", K = 3, eta = 1, include_correlation = TRUE)
  ordered <- "Omega_L_ordered[1:9] <- bt_lkj_cholesky_ordered(Omega_lkj_u[1:3], 3, position[1:3])"

  fit <- suppressWarnings(JAGS_fit(
    model_syntax   = paste0("model{\n", lkj$syntax, ordered, "\n}"),
    data           = list(position = c(3, 1, 2)),
    prior_list     = prior_list,
    add_parameters = c(lkj$monitor, "Omega_lkj_u", "Omega_L_ordered"),
    chains = 1, adapt = 100, burnin = 100, sample = 500,
    silent = TRUE, seed = 1
  ))
  expect_s3_class(fit, "BayesTools_fit")
  expect_true("BayesTools" %in% attr(fit, "jags_modules"))
  samples <- as.matrix(fit$mcmc)
  expect_true(all(c("theta", "phi", "sigma") %in% colnames(samples)))
  expect_true(all(is.finite(samples)))
  expect_true(all(samples[, "sigma"] > 0))
  # the ordered factor of the module is the R kernel (JAGS flattens row by row)
  u <- samples[, paste0("Omega_lkj_u[", 1:3, "]"), drop = FALSE]
  expected <- .bt_lkj_cholesky_cpc_u_to_L_ordered(u, K = 3, position = c(3, 1, 2))
  expect_identical(
    unname(samples[, paste0("Omega_L_ordered[", 1:9, "]"), drop = FALSE]),
    matrix(aperm(expected, c(1L, 3L, 2L)), nrow = nrow(u))
  )

  marglik <- JAGS_bridgesampling(
    fit           = fit,
    log_posterior = function(parameters, data) 0,
    data          = list(),
    prior_list    = prior_list,
    seed          = 1
  )
  expect_s3_class(marglik, "BayesTools_marglik")
  expect_true(is.finite(marglik$logml))
})
