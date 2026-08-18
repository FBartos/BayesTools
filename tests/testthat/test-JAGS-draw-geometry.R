skip_if_not_test_profile("unit")

test_that("draw geometry records exact chain-major MCMC timing", {

  chains <- coda::mcmc.list(
    coda::mcmc(matrix(1:6, nrow = 3, dimnames = list(NULL, c("theta", "private"))),
               start = 11, thin = 2),
    coda::mcmc(matrix(7:12, nrow = 3, dimnames = list(NULL, c("theta", "private"))),
               start = 11, thin = 2)
  )
  geometry <- .bt_draw_geometry_from_chains(chains)

  expect_s3_class(geometry, "BayesTools_draw_geometry")
  expect_identical(geometry$schema_version, 1L)
  expect_identical(geometry$total_draws, 6L)
  expect_identical(geometry$chain_order, 1:2)
  expect_identical(geometry$chains$iterations, c(3L, 3L))
  expect_identical(geometry$chains$start, c(11L, 11L))
  expect_identical(geometry$chains$end, c(15L, 15L))
  expect_identical(geometry$chains$thin, c(2L, 2L))
  expect_identical(geometry$chains$draw_start, c(1L, 4L))
  expect_identical(geometry$chains$draw_end, c(3L, 6L))

  malformed <- geometry
  malformed$chains$end[1L] <- 16L
  expect_error(.bt_validate_draw_geometry(malformed), "Refit the model")
})

test_that("coordinate materialization adds fixed values and removes internals", {

  chains <- coda::mcmc.list(
    coda::mcmc(
      matrix(c(1, 0, 2, 0, 3, 0), ncol = 2, byrow = TRUE,
             dimnames = list(NULL, c("theta", .bt_backend_anchor_name))),
      start = 5,
      thin = 3
    ),
    coda::mcmc(
      matrix(c(4, 0, 5, 0, 6, 0), ncol = 2, byrow = TRUE,
             dimnames = list(NULL, c("theta", .bt_backend_anchor_name))),
      start = 5,
      thin = 3
    )
  )
  fit <- chains
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- list(
    theta = prior("normal", list(0, 1)),
    fixed = prior("point", list(-4))
  )
  attr(fit, "backend_anchor") <- .bt_backend_anchor_name
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = colnames(chains[[1L]]),
    prior_list = attr(fit, "prior_list"),
    backend_anchor = .bt_backend_anchor_name
  )
  attr(fit, "draw_geometry") <- .bt_draw_geometry_from_chains(chains)
  fit <- .bt_attach_fit_contract(fit)

  materialized <- JAGS_materialize_draws(fit)
  expect_s3_class(materialized, "mcmc.list")
  expect_identical(colnames(materialized[[1L]]), c("theta", "fixed"))
  expect_identical(as.numeric(materialized[[1L]][, "fixed"]), rep(-4, 3))
  expect_identical(attr(materialized[[1L]], "mcpar"), c(5, 11, 3))
  expect_false(.bt_backend_anchor_name %in% colnames(materialized[[1L]]))

  fixed_only <- JAGS_materialize_draws(fit, "fixed")
  expect_identical(dim(fixed_only[[1L]]), c(3L, 1L))
  expect_identical(as.numeric(fixed_only[[2L]][, 1L]), rep(-4, 3))
})

test_that("replacement draws refresh fitted draw geometry", {

  original <- coda::mcmc.list(coda::mcmc(
    matrix(1:6, ncol = 1L, dimnames = list(NULL, "theta"))
  ))
  fit <- structure(
    list(mcmc = original),
    class = c("runjags", "BayesTools_fit")
  )
  attr(fit, "parameter_map") <- .bt_build_parameter_map(columns = "theta")
  attr(fit, "draw_geometry") <- .bt_draw_geometry_from_chains(original)
  fit <- .bt_attach_fit_contract(fit)
  replacement <- coda::mcmc.list(coda::mcmc(
    matrix(1:4, ncol = 1L, dimnames = list(NULL, "theta"))
  ))

  replaced <- JAGS_with_draws(fit, replacement)

  expect_equal(nrow(replaced[["mcmc"]][[1L]]), 4L)
  expect_equal(JAGS_draw_geometry(replaced)$total_draws, 4L)
  expect_equal(nrow(fit[["mcmc"]][[1L]]), 6L)
})

test_that("zero-public materialization retains chain and iteration dimensions", {

  chains <- coda::mcmc.list(coda::mcmc(
    matrix(0, nrow = 4, ncol = 1,
           dimnames = list(NULL, .bt_backend_anchor_name)),
    start = 2,
    thin = 2
  ))
  fit <- chains
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "backend_anchor") <- .bt_backend_anchor_name
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = .bt_backend_anchor_name,
    backend_anchor = .bt_backend_anchor_name
  )
  attr(fit, "draw_geometry") <- .bt_draw_geometry_from_chains(chains)
  fit <- .bt_attach_fit_contract(fit)

  materialized <- JAGS_materialize_draws(fit)
  expect_identical(length(materialized), 1L)
  expect_identical(dim(materialized[[1L]]), c(4L, 0L))
  expect_identical(attr(materialized[[1L]], "mcpar"), c(2, 8, 2))
})

test_that("private backend anchor is inserted only for an empty monitor set", {

  anchored <- .bt_add_backend_anchor(
    "model { theta <- 0 }",
    data = list(),
    prior_list = list(),
    add_parameters = character(),
    monitor = ""
  )
  expect_identical(anchored$backend_anchor, .bt_backend_anchor_name)
  expect_identical(anchored$monitor, .bt_backend_anchor_name)
  expect_match(anchored$model_syntax, .bt_backend_anchor_name, fixed = TRUE)

  ordinary <- .bt_add_backend_anchor(
    "model { theta <- 0 }",
    data = list(),
    prior_list = list(),
    add_parameters = "theta",
    monitor = "theta"
  )
  expect_null(ordinary$backend_anchor)
  expect_identical(ordinary$model_syntax, "model { theta <- 0 }")

  expect_error(
    .bt_add_backend_anchor(
      paste0("model { ", .bt_backend_anchor_name, " <- 1 }"),
      data = list(), prior_list = list(), add_parameters = character(), monitor = ""
    ),
    "reserved BayesTools backend node"
  )
})
