skip_if_not_test_profile("unit")

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
