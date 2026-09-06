skip_if_not_test_profile("unit")

test_that("random group dependencies are structural and row aligned", {

  terms <- list(
    list(block_name = "study", group_map = c(1L, 1L, 2L, 2L)),
    list(block_name = "estimate", group_map = 1:4)
  )
  block_diagonal <- kronecker(diag(2L), matrix(TRUE, 2L, 2L)) != 0
  expect_identical(random_effects_dependency_matrix(terms, 4L), block_diagonal)
  expect_identical(random_effects_dependency_matrix(terms, 4L, "estimate"), diag(TRUE, 4L))
  expect_identical(random_effects_dependency_matrix(terms, 4L, character()), diag(TRUE, 4L))
  terms[[1L]]$group_covariance <- list(type = "known", kernel = matrix(c(1, .3, .3, 1), 2L))
  expect_identical(random_effects_dependency_matrix(terms, 4L), matrix(TRUE, 4L, 4L))
  terms[[1L]]$group_covariance$kernel <- diag(2L)
  expect_identical(random_effects_dependency_matrix(terms, 4L), block_diagonal)
  expect_error(
    random_effects_dependency_matrix(terms, 4L, "missing"),
    "Requested random-effect dependency blocks must be unique existing block names.",
    fixed = TRUE
  )
  terms[[1L]]$group_map[[1L]] <- 1.5
  expect_error(
    random_effects_dependency_matrix(terms, 4L),
    "Random-effect grouping metadata are invalid for dependency construction.",
    fixed = TRUE
  )
})

test_that("factor covariance uses the shared group and row-scale geometry", {

  X <- cbind(1, c(-1, 0, 1))
  Z <- cbind(X * c(1, 1, 0), X * c(0, 0, 1))
  G <- matrix(c(.04, .02, .02, .10), 2L)
  L <- matrix(c(.2, .1, 0, .3), 2L)
  for(type in c("group", "row_group", "known_group")){
    kernel <- if(type == "known_group") matrix(c(1, .4, .4, 1), 2L) else diag(2L)
    scale <- if(type == "row_group") c(1, 2, .5) else rep(1, 3L)
    plan <- list(type = type, model_matrix = X, group_map = c(1L, 1L, 2L),
                 coefficient_structure = "dense", group_covariance = kernel)
    state <- list(coefficient_factor = L, row_scale = scale)
    factors <- list(representation = "factor_state", factor_plans = list(plan),
                    factor_states = list(state), row_blocks = list(1:3))
    expected <- (Z * scale) %*% kronecker(kernel, G) %*% t(Z * scale)
    observed <- random_effects_marginal_factor_vcov(factors)
    expect_identical(dim(observed), c(1L, 3L, 3L))
    expect_equal(observed[1L, , ], unname(expected), tolerance = 1e-14)
    factors$factor_states[[1L]]$coefficient_factor[,] <- 0
    expect_identical(random_effects_marginal_factor_vcov(factors), array(0, c(1L, 3L, 3L)))
  }
})

.rank_one_factor_fixture <- function(){

  group_matrix <- matrix(
    c(1, 0,
      0, 1,
      1, 0,
      0, 1),
    nrow = 4L,
    byrow = TRUE
  )
  plans <- list(
    study = list(
      type                  = "group",
      model_matrix          = group_matrix,
      group_map             = c(1L, 1L, 2L, 2L),
      coefficient_structure = "dense"
    ),
    estimate = list(
      type                  = "group",
      model_matrix          = matrix(1, nrow = 4L, ncol = 1L),
      group_map             = seq_len(4L),
      coefficient_structure = "diagonal"
    )
  )
  states <- lapply(c(1, 1.5), function(scale){
    list(
      study = list(coefficient_factor = matrix(
        c(0.8, 0, -0.25, 0.5) * scale,
        nrow = 2L,
        byrow = TRUE
      )),
      estimate = list(coefficient_factor = matrix(0.2 * scale, 1L, 1L))
    )
  })
  factors <- list(
    factor_plans  = plans,
    factor_states = states,
    row_blocks    = list(1:2, 3:4),
    metadata = list(
      n_draws         = length(states),
      n_rows          = 4L,
      included_blocks = names(plans)
    )
  )
  class(factors) <- c(
    "BayesTools_random_effects_marginal_factor_states",
    "list"
  )
  factors
}


test_that("certified factor reduction reproduces compiled covariance", {

  factors <- .rank_one_factor_fixture()
  factorized <- random_effects_marginal_diagonal_factor(factors)
  dense <- random_effects_marginal_factor_vcov(factors)
  expect_identical(dim(dense), c(2L, 4L, 4L))
  expect_equal(dense[2L, , ], 1.5^2 * dense[1L, , ], tolerance = 1e-14)

  expect_s3_class(
    factorized,
    "BayesTools_random_effects_marginal_diagonal_factor"
  )
  expect_identical(factorized$ranks, c(1L, 1L))
  expect_identical(factorized$row_blocks, factors$row_blocks)
  for(draw in seq_along(factors$factor_states)){
    direct <- matrix(0, nrow = 4L, ncol = 4L)
    for(block in seq_along(factors$factor_plans)){
      plan  <- factors$factor_plans[[block]]
      state <- factors$factor_states[[draw]][[block]]
      basis <- plan$model_matrix %*% state$coefficient_factor
      direct <- direct + tcrossprod(basis) *
        outer(plan$group_map, plan$group_map, "==")
    }
    for(rows in factors$row_blocks){
      block <- which(vapply(
        factors$row_blocks,
        identical,
        logical(1),
        rows
      ))
      factor_loading <- factorized$loadings[[block]][
        draw,
        ,
        ,
        drop = FALSE
      ]
      dim(factor_loading) <- c(length(rows), factorized$ranks[[block]])
      factor_reconstructed <- diag(factorized$diagonal[draw, rows]) +
        tcrossprod(factor_loading)
      expect_equal(
        factor_reconstructed,
        direct[rows, rows, drop = FALSE],
        tolerance = 1e-14
      )
    }
  }
  expect_lt(factorized$loadings[[1L]][1L, 2L, 1L], 0)

  bridge_state <- list(
    representation = "factor_state",
    factor_plans   = factors$factor_plans,
    factor_states  = factors$factor_states[[1L]],
    row_blocks     = factors$row_blocks
  )
  bridge_factorized <- random_effects_marginal_diagonal_factor(bridge_state)
  expect_equal(
    bridge_factorized$diagonal,
    factorized$diagonal[1L, , drop = FALSE]
  )
  expect_equal(
    bridge_factorized$loadings[[1L]],
    factorized$loadings[[1L]][1L, , , drop = FALSE]
  )
})


test_that("certified factor reduction preserves higher structural rank", {

  factors <- .rank_one_factor_fixture()
  factors$factor_plans <- list(higher_rank = list(
    type                  = "group",
    model_matrix          = diag(3L),
    group_map             = rep(1L, 3L),
    coefficient_structure = "dense"
  ))
  factors$factor_states <- list(list(higher_rank = list(
    coefficient_factor = diag(3L)
  )))
  factors$row_blocks <- list(1:3)
  factors$metadata$n_draws <- 1L
  factors$metadata$n_rows <- 3L
  factors$metadata$included_blocks <- "higher_rank"

  factorized <- random_effects_marginal_diagonal_factor(factors)
  expect_identical(factorized$ranks, 2L)
  loading <- factorized$loadings[[1L]][1L, , , drop = FALSE]
  dim(loading) <- c(3L, 2L)
  expect_equal(
    diag(factorized$diagonal[1L, ]) + tcrossprod(loading),
    diag(3L),
    tolerance = 1e-14
  )

  factors$factor_plans[[1L]]$type <- character()
  expect_error(
    random_effects_marginal_diagonal_factor(factors),
    "unsupported factor structure",
    class = "BayesTools_random_effects_marginal_factor_unavailable"
  )
})


test_that("batched factor reduction preserves row scales and zero draws", {

  factors <- .rank_one_factor_fixture()
  factors$factor_plans$study$type <- "row_group"
  factors$factor_states[[1L]]$study$row_scale <- c(1, 0, 2, .5)
  factors$factor_states[[2L]]$study$row_scale <- c(.5, 2, 0, 1)
  reduced <- random_effects_marginal_diagonal_factor(factors)
  for(draw in 1:2){
    state <- factors$factor_states[[draw]]
    basis <- (factors$factor_plans$study$model_matrix %*%
      state$study$coefficient_factor) * state$study$row_scale
    expected <- tcrossprod(basis) * outer(c(1L, 1L, 2L, 2L),
                                         c(1L, 1L, 2L, 2L), "==") +
      diag(as.numeric(state$estimate$coefficient_factor)^2, 4L)
    for(block in 1:2){
      rows <- factors$row_blocks[[block]]
      loading <- matrix(reduced$loadings[[block]][draw, , ], ncol = 1L)
      expect_equal(diag(reduced$diagonal[draw, rows]) + tcrossprod(loading),
                   expected[rows, rows], tolerance = 1e-14)
    }
  }

  # One row and one coefficient must keep the draw axis, including a zero SD.
  factors$factor_plans <- list(single = list(
    type = "row_group", model_matrix = matrix(1, 1L, 1L),
    group_map = 1L, coefficient_structure = "diagonal"
  ))
  factors$factor_states <- lapply(c(0, 2, 3), function(sd){
    list(single = list(coefficient_factor = matrix(sd, 1L, 1L), row_scale = 2))
  })
  factors$row_blocks <- list(1L)
  factors$metadata <- list(n_draws = 3L, n_rows = 1L, included_blocks = "single")
  reduced <- random_effects_marginal_diagonal_factor(factors)
  expect_identical(reduced$diagonal, matrix(c(0, 16, 36), ncol = 1L))
  expect_identical(dim(reduced$loadings[[1L]]), c(3L, 1L, 0L))
})


test_that("bridge factor reduction caches only invariant contract plans", {

  factors <- .rank_one_factor_fixture()
  contract_id <- new.env(parent = emptyenv())
  bridge_state <- list(
    representation = "factor_state",
    contract_id    = contract_id,
    factor_plans   = factors$factor_plans,
    factor_states  = factors$factor_states[[1L]],
    row_blocks     = factors$row_blocks
  )
  cache <- new.env(parent = emptyenv())
  plan_calls <- 0L
  original_plan <- .bt_random_effect_marginal_diagonal_factor_plan
  testthat::local_mocked_bindings(
    .bt_random_effect_marginal_diagonal_factor_plan = function(...) {
      plan_calls <<- plan_calls + 1L
      original_plan(...)
    },
    .package = "BayesTools"
  )

  first <- random_effects_marginal_diagonal_factor(
    bridge_state,
    cache = cache
  )
  bridge_state$factor_states <- lapply(
    bridge_state$factor_states,
    function(state) {
      state$coefficient_factor <- 2 * state$coefficient_factor
      state
    }
  )
  second <- random_effects_marginal_diagonal_factor(
    bridge_state,
    cache = cache
  )

  expect_identical(plan_calls, 1L)
  expect_equal(second$loadings[[1L]], 2 * first$loadings[[1L]])
  expect_equal(second$diagonal, 4 * first$diagonal)

  bridge_state$contract_id <- new.env(parent = emptyenv())
  expect_error(
    random_effects_marginal_diagonal_factor(bridge_state, cache = cache),
    "contract changed",
    fixed = TRUE
  )
  expect_error(
    random_effects_marginal_diagonal_factor(bridge_state, cache = list()),
    "'cache' must be NULL or an environment.",
    fixed = TRUE
  )
})
