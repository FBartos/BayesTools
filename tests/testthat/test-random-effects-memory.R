skip_if_not_test_profile("unit")

.random_memory_option <- "BayesTools.random_effects_memory_limit_bytes"

.random_memory_restore_option <- function(value){
  do.call(options, stats::setNames(list(value), .random_memory_option))
}

.random_memory_formula <- function(){

  data <- data.frame(
    f = factor(rep(letters[1:6], each = 2L), levels = letters[1:6]),
    id = factor(rep(c("g1", "g2", "g3"), length.out = 12L))
  )
  result <- JAGS_formula(
    formula = ~ 1 + random(1 + f | id, name = "id", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior(
          "normal",
          list(0, 1),
          truncation = list(lower = 0, upper = Inf)
        )
      )
    )
  )

  list(data = data, result = result)
}

test_that("random-effect memory estimates expose payload and peak components", {

  estimate <- BayesTools:::.bt_random_effect_design_memory_estimate(
    n_rows = 100L,
    n_columns = 20L,
    n_groups = 5L,
    monitor_policy = random_monitor(
      latent = TRUE,
      coefficients = FALSE
    ),
    structure = "us",
    compile_mode = "sampled"
  )

  expect_equal(estimate$payload_bytes, 8 * 100 * 20)
  expect_gt(estimate$peak_bytes, estimate$payload_bytes)
  expect_named(
    estimate$components,
    c("design_working", "group_state", "covariance_state", "row_maps")
  )
  expect_match(
    BayesTools:::.bt_random_effect_format_bytes(16 * 1024^3),
    "16.0 GiB",
    fixed = TRUE
  )

  output <- BayesTools:::.bt_random_effect_output_memory_estimate(
    operation = "test covariance",
    n_rows = 100L,
    n_draws = 20L,
    covariance = TRUE
  )
  expect_equal(output$payload_bytes, 8 * 20 * 100^2)
  expect_equal(output$peak_bytes, 3 * output$payload_bytes)
})

test_that("random-effect memory option is validated and can disable the guard", {

  old_limit <- getOption(.random_memory_option)
  on.exit(.random_memory_restore_option(old_limit), add = TRUE)

  .random_memory_restore_option(NULL)
  expect_equal(
    BayesTools:::.bt_random_effect_memory_limit_bytes(),
    16 * 1024^3
  )

  for(value in list(0, -1, NA_real_, c(1, 2), "16 GiB")){
    .random_memory_restore_option(value)
    expect_error(
      BayesTools:::.bt_random_effect_memory_limit_bytes(),
      "must be a positive numeric scalar or Inf",
      fixed = TRUE
    )
  }

  .random_memory_restore_option(Inf)
  estimate <- BayesTools:::.bt_random_effect_output_memory_estimate(
    operation = "test output",
    n_rows = 1e8,
    n_draws = 1e8
  )
  expect_silent(
    BayesTools:::.bt_random_effect_check_memory(
      estimate,
      alternative = "Reduce the test size."
    )
  )
})

test_that("formula design memory guard stops before the full dense allocation", {

  old_limit <- getOption(.random_memory_option)
  on.exit(.random_memory_restore_option(old_limit), add = TRUE)
  original_design <- BayesTools:::.bt_random_effect_design_matrix
  full_design_calls <- 0L
  testthat::local_mocked_bindings(
    .bt_random_effect_design_matrix = function(formula, data, ...){
      if(nrow(data) > 1L){
        full_design_calls <<- full_design_calls + 1L
      }
      original_design(formula = formula, data = data, ...)
    },
    .package = "BayesTools"
  )

  .random_memory_restore_option(512)
  expect_error(
    .random_memory_formula(),
    regexp = paste0(
      "Estimated peak memory for random-effect design construction for block ",
      "'id' exceeds option '", .random_memory_option, "'. ",
      "Dimensions: 12 rows x 6 columns; 3 groups; primary payload: ",
      "576 bytes"
    ),
    fixed = TRUE
  )
  expect_equal(full_design_calls, 0L)

  .random_memory_restore_option(Inf)
  expect_s3_class(.random_memory_formula()$result$formula_design,
                  "BayesTools_formula_design")
  expect_equal(full_design_calls, 1L)
})

test_that("output guard reports the conservative peak and safer alternative", {

  old_limit <- getOption(.random_memory_option)
  on.exit(.random_memory_restore_option(old_limit), add = TRUE)
  .random_memory_restore_option(Inf)
  specification <- .random_memory_formula()
  random_term <- specification$result$formula_design$random_effects[[1L]]
  sd_names <- unique(random_term$sd_parameter_names)
  posterior <- matrix(
    1,
    nrow = 1L,
    ncol = length(sd_names),
    dimnames = list(NULL, sd_names)
  )

  .random_memory_restore_option(1000)
  expect_error(
    random_effects_marginal_vcov(
      fit = specification$result$formula_design,
      posterior_samples = posterior,
      prior_list = specification$result$prior_list
    ),
    "Use diagonal_only = TRUE",
    fixed = TRUE
  )
  diagonal <- random_effects_marginal_vcov(
    fit = specification$result$formula_design,
    posterior_samples = posterior,
    prior_list = specification$result$prior_list,
    diagonal_only = TRUE
  )
  expect_identical(dim(diagonal$samples), c(1L, nrow(specification$data)))

  .random_memory_restore_option(1e6)
  estimate <- BayesTools:::.bt_random_effect_output_memory_estimate(
    operation = "random-effect marginal covariance",
    n_rows = 100L,
    n_draws = 20L,
    covariance = TRUE
  )

  expect_error(
    BayesTools:::.bt_random_effect_check_memory(
      estimate = estimate,
      block_name = "id",
      alternative = "Use diagonal_only = TRUE."
    ),
    regexp = paste0(
      "primary payload: 1,600,000 bytes",
      ".*conservative peak estimate: 4,800,000 bytes",
      ".*configured ceiling: 1,000,000 bytes",
      ".*Use diagonal_only = TRUE"
    )
  )
})
