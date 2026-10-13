skip_if_not_test_profile("unit")

test_that("CAR admission includes factors before either ordinary or row SD reconstruction", {

  old <- getOption("BayesTools.random_effects_memory_limit_bytes")
  on.exit(do.call(options, list(BayesTools.random_effects_memory_limit_bytes = old)), add = TRUE)
  options(BayesTools.random_effects_memory_limit_bytes = Inf)
  data <- data.frame(time = c(0, .5, 2), id = rep("a", 3L))
  compile <- function(row_source){
    block <- if(row_source){
      random_block(sd_source = random_sd_source("tau", shape = "row"),
        cor = prior("normal", list(0, .5)))
    }else{
      random_block(sd = prior("normal", list(0, 1), list(0, Inf)),
        cor = prior("normal", list(0, .5)))
    }
    JAGS_formula(~ 1 + car(time | id), "mu", data,
      list(intercept = prior("normal", list(0, 1))), prior_random = prior_random(id = block))
  }
  fixtures <- lapply(c(FALSE, TRUE), compile)
  original <- BayesTools:::.bt_random_effect_cholesky_draws
  factor_calls <- 0L
  testthat::local_mocked_bindings(
    .bt_random_effect_cholesky_draws = function(...){
      factor_calls <<- factor_calls + 1L
      original(...)
    }, .package = "BayesTools")
  for(i in seq_along(fixtures)){
    fixture <- fixtures[[i]]
    term <- fixture$formula_design$random_effects[[1L]]
    sd <- if(i == 1L) stats::setNames(2, unique(term$sd_parameter_names)) else
      stats::setNames(c(2, 3, 4), paste0("tau[", 1:3, "]"))
    values <- c(sd, stats::setNames(atanh(.5), term$correlation$sample_name))
    posterior <- matrix(rep(values, each = 2L), 2L, dimnames = list(NULL, names(values)))
    call <- function(diagonal){
      random_effects_marginal_vcov(fixture$formula_design, posterior_samples = posterior,
        prior_list = fixture$prior_list, diagonal_only = diagonal)
    }
    options(BayesTools.random_effects_memory_limit_bytes = 200)
    calls_before <- factor_calls
    expect_error(call(TRUE), "3 coefficient/time columns", fixed = TRUE)
    expect_error(call(FALSE), "CAR factor storage is still required", fixed = TRUE)
    expect_identical(factor_calls, calls_before)
    options(BayesTools.random_effects_memory_limit_bytes = Inf)
    correlation <- .5^abs(outer(data$time, data$time, "-"))
    scales <- if(i == 1L) rep(2, 3L) else c(2, 3, 4)
    expected <- correlation * tcrossprod(scales)
    expect_equal(unname(call(FALSE)$samples[1L, , ]), expected, tolerance = 1e-14)
    expect_equal(unname(call(TRUE)$samples[1L, ]), diag(expected), tolerance = 1e-14)
  }
  for(diagonal in c(FALSE, TRUE)){
    estimate <- BayesTools:::.bt_random_effect_output_memory_estimate("test", 3L, 2L,
      covariance = TRUE, diagonal_only = diagonal, car_columns = 3L)
    payload <- 8 * 2 * if(diagonal) 3 else 9
    expect_identical(estimate$components,
      c(output_working = 3 * payload, car_factor_working = 2 * 8 * 2 * 9,
        car_contraction_working = 4 * 8 * 3 * 3))
    expect_equal(estimate$peak_bytes, sum(estimate$components))
  }
  expect_identical(BayesTools:::.bt_random_effect_output_memory_estimate("test", 3L, 2L)$components,
    c(output_working = 3 * 8 * 3 * 2))
  expect_identical(BayesTools:::.bt_random_effect_output_memory_estimate("test", 1e200, 1e200,
    car_columns = 1e200)$peak_bytes, Inf)
})

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
  cell_design_calls <- 0L
  testthat::local_mocked_bindings(
    .bt_random_effect_design_matrix = function(formula, data, ...){
      # the dense design of all 12 data rows; the SD factor prior of term 'f'
      # evaluates its design once at the 6 level cells
      if(nrow(data) == 12L){
        full_design_calls <<- full_design_calls + 1L
      }else if(nrow(data) > 1L){
        expect_identical(nrow(data), 6L)
        cell_design_calls <<- cell_design_calls + 1L
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
  expect_equal(cell_design_calls, 0L)

  .random_memory_restore_option(Inf)
  expect_s3_class(.random_memory_formula()$result$formula_design,
                  "BayesTools_formula_design")
  expect_equal(full_design_calls, 1L)
  expect_equal(cell_design_calls, 1L)
})

test_that("factor-level design memory guard stops before the dense cell design", {

  old_limit <- getOption(.random_memory_option)
  on.exit(.random_memory_restore_option(old_limit), add = TRUE)
  # 100 rows of a 20 x 20 factor interaction: the block has 400 columns, so
  # the data design is 100 x 400 (0.3 MiB, 1.3 MiB peak with its working
  # copies) and the design at every level cell 400 x 400 (1.2 MiB, 4.9 MiB)
  set.seed(1)
  data <- data.frame(
    g  = factor(sample(sprintf("g%02d", 1:20), 100L, replace = TRUE), levels = sprintf("g%02d", 1:20)),
    h  = factor(sample(sprintf("h%02d", 1:20), 100L, replace = TRUE), levels = sprintf("h%02d", 1:20)),
    id = factor(rep(1:5, length.out = 100L))
  )
  compile <- function(){
    JAGS_formula(
      ~ 1 + (0 + g:h || id), "mu", data,
      list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(
        sd = prior("normal", list(0, 1), list(0, Inf))
      ))
    )
  }
  original_design <- BayesTools:::.bt_random_effect_design_matrix
  cell_design_calls <- 0L
  testthat::local_mocked_bindings(
    .bt_random_effect_design_matrix = function(formula, data, ...){
      if(nrow(data) == 400L){
        cell_design_calls <<- cell_design_calls + 1L
      }
      original_design(formula = formula, data = data, ...)
    },
    .package = "BayesTools"
  )

  .random_memory_restore_option(2 * 1024^2)
  expect_error(
    compile(),
    paste0(
      "Estimated peak memory for random-effect factor-level design construction ",
      "for block 'id' exceeds option '", .random_memory_option, "'. ",
      "Dimensions: 400 level cells x 400 columns; primary payload: ",
      "1,280,000 bytes (1.22 MiB); conservative peak estimate: ",
      "5,126,400 bytes (4.89 MiB); configured ceiling: ",
      "2,097,152 bytes (2.00 MiB). This guard cannot promise that smaller ",
      "allocations will succeed. Reduce the number of levels of the factors ",
      "in the random-effect term or the number of random-effect columns, or ",
      "simplify the random-effect design."
    ),
    fixed = TRUE
  )
  expect_identical(cell_design_calls, 0L)

  .random_memory_restore_option(Inf)
  result <- compile()
  expect_identical(cell_design_calls, 1L)
  sd_prior <- result$prior_list$mu__xREx__id_g__xXx__h
  expect_identical(dim(attr(sd_prior, "factor_design")), c(400L, 400L))
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

test_that("alternate marginal covariance builders pass the shared memory ceiling", {

  old_limit <- getOption(.random_memory_option)
  on.exit(.random_memory_restore_option(old_limit), add = TRUE)
  .random_memory_restore_option(Inf)
  data <- data.frame(id = factor(seq_len(20L)))
  result <- JAGS_formula(
    formula = ~ 1 + (1 | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("normal", list(0, 1), list(0, Inf)))
    ),
    random_effects_compile = random_effects_compile(marginalized = "id")
  )
  random_term <- result$formula_design$random_effects[[1L]]
  sd_name <- unique(random_term$sd_parameter_names)
  posterior <- matrix(2, nrow = 1L, ncol = 1L, dimnames = list(NULL, sd_name))
  factors <- random_effects_marginal_factor_states(
    result$formula_design,
    posterior_samples = posterior,
    prior_list = result$prior_list,
    row_blocks = as.list(seq_len(20L))
  )
  primary <- function(){
    random_effects_marginal_vcov(
      result$formula_design,
      posterior_samples = posterior,
      prior_list = result$prior_list
    )
  }
  dense <- function() random_effects_marginal_factor_vcov(factors)
  variance <- function() random_effects_marginal_variance_factors(result$formula_design)

  reference_dense <- dense()
  reference_variance <- variance()
  expect_identical(dim(reference_dense), c(1L, 20L, 20L))
  expect_equal(unname(reference_dense[1L, , ]), diag(4, 20L), tolerance = 1e-12)
  expect_equal(unname(reference_variance$blocks$id$row_multiplier), rep(1, 20L))

  # one draw of 20 x 20 doubles: the shared estimate is three payloads
  estimate <- BayesTools:::.bt_random_effect_output_memory_estimate(
    "test", n_rows = 20L, n_draws = 1L, covariance = TRUE
  )
  expect_equal(estimate$peak_bytes, 3 * 8 * 20^2)

  .random_memory_restore_option(100)
  expect_error(primary(), "random-effect marginal covariance for block 'id' exceeds option", fixed = TRUE)
  expect_error(
    dense(),
    paste0("Estimated peak memory for random_effects_marginal_factor_vcov() ",
           "exceeds option 'BayesTools.random_effects_memory_limit_bytes'. ",
           "Dimensions: 1 draws x 20 rows x 20 rows; primary payload: 3,200 bytes"),
    fixed = TRUE
  )
  expect_error(
    variance(),
    paste0("Estimated peak memory for random_effects_marginal_variance_factors() ",
           "for block 'id' exceeds option ",
           "'BayesTools.random_effects_memory_limit_bytes'."),
    fixed = TRUE
  )
  expect_error(dense(), "Raise the option", fixed = TRUE)
  expect_error(variance(), "Raise the option", fixed = TRUE)

  # the alternate builders sit on the threshold of the shared estimate
  .random_memory_restore_option(estimate$peak_bytes - 1)
  expect_error(dense(), "exceeds option", fixed = TRUE)
  expect_error(variance(), "exceeds option", fixed = TRUE)
  .random_memory_restore_option(estimate$peak_bytes)
  expect_identical(dense(), reference_dense)
  expect_identical(variance()$blocks, reference_variance$blocks)
  .random_memory_restore_option(1e6)
  expect_no_error(primary())
})

# Runs `call` with the memory guard recording the peak estimate it admitted and,
# where R can profile memory, the largest vector the call allocated and the
# total of its allocations of at least `threshold` bytes. Allocations are
# profiled rather than read from gc()'s max-used column, which is the heap
# size at the time of a collection: it follows the collector's trigger and the
# garbage of small allocations, not the arrays the call keeps alive.
.random_memory_measure <- function(call, threshold = 1000){

  admitted <- numeric()
  original <- BayesTools:::.bt_random_effect_check_memory
  profile <- isTRUE(unname(capabilities("profmem")))
  log <- tempfile(fileext = ".log")
  on.exit(unlink(log), add = TRUE)
  value <- testthat::with_mocked_bindings(
    {
      if(profile){
        utils::Rprofmem(log, threshold = threshold)
      }
      tryCatch(call(), finally = if(profile) utils::Rprofmem(NULL))
    },
    .bt_random_effect_check_memory = function(estimate, ...){
      admitted <<- c(admitted, estimate$peak_bytes)
      original(estimate, ...)
    },
    .package = "BayesTools"
  )
  sizes <- if(profile){
    suppressWarnings(as.numeric(sub(" .*", "", readLines(log))))
  }else{
    NA_real_
  }
  largest <- if(profile) max(0, sizes, na.rm = TRUE) else NA_real_
  allocated <- if(profile) sum(sizes, na.rm = TRUE) else NA_real_

  list(
    value = value, admitted = admitted, largest = largest, allocated = allocated,
    profile = profile
  )
}

.random_memory_covariance <- function(result, posterior, diagonal_only){

  random_effects_marginal_vcov(
    result$formula_design,
    posterior_samples = posterior,
    prior_list = result$prior_list,
    diagonal_only = diagonal_only
  )
}

.random_memory_known_group <- function(n_rows){

  data <- data.frame(id = factor(rep(c("a", "b"), length.out = n_rows)))
  kernel <- diag(2L)
  dimnames(kernel) <- list(c("a", "b"), c("a", "b"))
  JAGS_formula(
    random_effects_formula(
      ~ 1 | id,
      group_covariance = random_group_covariance(kernel, scale = "none")
    ),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("normal", list(0, 1), list(0, Inf)))
    )
  )
}

test_that("known-group diagonal output is linear in the rows and covered by its estimate", {

  old_limit <- getOption(.random_memory_option)
  on.exit(.random_memory_restore_option(old_limit), add = TRUE)
  .random_memory_restore_option(Inf)
  n_rows <- 1000L
  result <- .random_memory_known_group(n_rows)
  random_term <- result$formula_design$random_effects[[1L]]
  posterior <- matrix(
    2, nrow = 1L, ncol = 1L,
    dimnames = list(NULL, unique(random_term$sd_parameter_names))
  )

  # admitted under a ceiling far below the rows x rows matrices
  .random_memory_restore_option(50000)
  measured <- .random_memory_measure(
    function() .random_memory_covariance(result, posterior, diagonal_only = TRUE)
  )
  expect_identical(dim(measured$value$samples), c(1L, n_rows))
  expect_equal(unname(measured$value$samples[1L, ]), rep(4, n_rows))
  expect_length(measured$admitted, 1L)
  expect_lt(measured$admitted, 50000)
  if(measured$profile){
    expect_lte(measured$largest, measured$admitted)
  }

  # the full output keeps the rows x rows matrices and budgets their working copies
  .random_memory_restore_option(Inf)
  small_result <- .random_memory_known_group(200L)
  full <- .random_memory_measure(
    function() .random_memory_covariance(small_result, posterior, diagonal_only = FALSE)
  )
  expect_identical(dim(full$value$samples), c(1L, 200L, 200L))
  # the kernel expansion, mask and covariance copies come on top of 3 outputs
  expect_gte(full$admitted, (3 * 8 + 44) * 200^2)
  if(full$profile){
    expect_lte(full$largest, full$admitted)
  }
  .random_memory_restore_option(1500000)
  expect_error(
    .random_memory_covariance(small_result, posterior, diagonal_only = FALSE),
    "exceeds option 'BayesTools.random_effects_memory_limit_bytes'",
    fixed = TRUE
  )
})

test_that("the output estimate bounds the large arrays a dense factor path allocates", {

  old_limit <- getOption(.random_memory_option)
  on.exit(.random_memory_restore_option(old_limit), add = TRUE)
  .random_memory_restore_option(Inf)

  # A US block of 30 coefficient columns and 2,000 draws on 2 rows: the
  # Cholesky draws and the correlation draws (draws x 30 x 30) are the only
  # arrays of their size. Reading the Cholesky draws, validating them and
  # forming the correlation draws allocate no further full copy.
  n_columns <- 30L
  n_draws <- 2000L
  data <- data.frame(
    id = factor(c("a", "a")),
    f = factor(c("level1", "level2"), levels = paste0("level", seq_len(n_columns)))
  )
  result <- JAGS_formula(
    ~ 1 + us(0 + f | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("normal", list(0, 1), list(0, Inf)),
        contrasts = c(f = "independent")
      )
    )
  )
  term <- result$formula_design$random_effects[[1L]]
  values <- c(
    stats::setNames(rep(1, n_columns), unique(term$sd_parameter_names)),
    stats::setNames(
      as.vector(diag(n_columns)),
      as.vector(.bt_random_effect_cholesky_names(term, n_columns))
    )
  )
  posterior <- matrix(rep(values, each = n_draws), n_draws, dimnames = list(NULL, names(values)))
  array_bytes <- 8 * n_draws * n_columns^2
  measured <- .random_memory_measure(
    function() .random_memory_covariance(result, posterior, diagonal_only = FALSE),
    threshold = array_bytes / 4
  )
  expect_identical(dim(measured$value$samples), c(n_draws, 2L, 2L))
  expect_equal(unname(measured$value$samples[n_draws, , ]), diag(2L), tolerance = 1e-14)
  expect_length(measured$admitted, 1L)
  expect_gte(measured$admitted, 2 * array_bytes)
  if(measured$profile){
    # every allocation of a quarter of a draws x 30 x 30 array or more, together,
    # stays within what the call was admitted with
    expect_lte(measured$allocated, measured$admitted)
  }

  # A marginalized random intercept with 2,000 rows: the rows x rows
  # covariance is the only array of its size; its kernel and product, and the
  # copies of its diagonal check, are not formed.
  n_rows <- 2000L
  variance_result <- JAGS_formula(
    ~ 1 + (1 | id),
    parameter = "mu",
    data = data.frame(id = factor(seq_len(n_rows))),
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("normal", list(0, 1), list(0, Inf)))
    ),
    random_effects_compile = random_effects_compile(marginalized = "id")
  )
  measured <- .random_memory_measure(
    function() random_effects_marginal_variance_factors(variance_result$formula_design),
    threshold = 8 * n_rows^2 / 4
  )
  expect_length(measured$admitted, 1L)
  expect_true(measured$value$blocks$id$row_covariance_diagonal)
  expect_identical(measured$value$blocks$id$max_off_diagonal, 0)
  expect_equal(unname(measured$value$blocks$id$row_multiplier), rep(1, n_rows))
  if(measured$profile){
    expect_lte(measured$allocated, measured$admitted)
  }
})

test_that("marginal variance factors build the row covariance of groups and kernels column by column", {

  old_limit <- getOption(.random_memory_option)
  on.exit(.random_memory_restore_option(old_limit), add = TRUE)
  .random_memory_restore_option(Inf)
  base_covariance <- BayesTools:::.bt_random_effect_marginal_variance_base_covariance
  diagonal_status <- BayesTools:::.bt_random_effect_row_covariance_diagonal_status

  # groups of repeated rows: the same-group indicator times the products of
  # the row multipliers
  multiplier <- c(1, -2, 3, 4, -5)
  group_map <- c(1L, 1L, 2L, 3L, 3L)
  reference <- outer(group_map, group_map, "==") * tcrossprod(multiplier)
  data <- data.frame(id = factor(c("a", "a", "b", "c", "c")), x = multiplier)
  repeated <- JAGS_formula(
    ~ 1 + (0 + x | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("normal", list(0, 1), list(0, Inf)))
    ),
    random_effects_compile = random_effects_compile(marginalized = "id")
  )
  factors <- random_effects_marginal_variance_factors(
    repeated$formula_design, require_diagonal = FALSE
  )
  block <- repeated$formula_design$random_effects[[1L]]
  covariance <- base_covariance(block, matrix(multiplier, ncol = 1L), group_map)
  expect_identical(unname(covariance), unname(reference))
  expect_false(factors$blocks$id$row_covariance_diagonal)
  expect_identical(factors$blocks$id$max_off_diagonal, 20)
  expect_equal(unname(factors$blocks$id$row_multiplier), multiplier^2)

  # a known group kernel replaces the same-group indicator
  kernel <- matrix(
    c(4, 1, .5, 1, 9, 2, .5, 2, 16), 3L,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  known <- JAGS_formula(
    random_effects_formula(
      ~ 0 + x | id, group_covariance = random_group_covariance(kernel, scale = "none")
    ),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("normal", list(0, 1), list(0, Inf)))
    ),
    random_effects_compile = random_effects_compile(marginalized = "id")
  )
  known_block <- known$formula_design$random_effects[[1L]]
  expect_identical(
    unname(base_covariance(known_block, matrix(multiplier, ncol = 1L), group_map)),
    unname(kernel[group_map, group_map] * tcrossprod(multiplier))
  )

  # the diagonal status reads the off-diagonal entries without a copy of the matrix
  status <- diagonal_status(reference)
  expect_false(status$is_diagonal)
  expect_identical(status$max_off_diagonal, 20)
  expect_identical(status$tolerance, 0)
  status <- diagonal_status(diag(c(3, -1, 7)))
  expect_true(status$is_diagonal)
  expect_identical(status$max_off_diagonal, 0)
  expect_false(diagonal_status(matrix(c(1, NA, 0, 1), 2L))$is_diagonal)
})

test_that("the output estimate budgets the SD draws and coefficient factors of every structure", {

  old_limit <- getOption(.random_memory_option)
  on.exit(.random_memory_restore_option(old_limit), add = TRUE)
  .random_memory_restore_option(Inf)
  data <- data.frame(
    id = factor(c("a", "a")),
    f = factor(c("level1", "level2"), levels = paste0("level", 1:10))
  )
  compile <- function(structure){
    formula <- stats::as.formula(paste0("~ 1 + ", structure, "(0 + f | id)"))
    JAGS_formula(
      formula,
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = prior("normal", list(0, 1), list(0, Inf)),
          contrasts = c(f = "independent")
        )
      )
    )
  }
  n_draws <- 100L
  posterior_for <- function(result, cholesky){
    term <- result$formula_design$random_effects[[1L]]
    values <- stats::setNames(
      rep(1, length(unique(term$sd_parameter_names))),
      unique(term$sd_parameter_names)
    )
    if(cholesky){
      coordinates <- .bt_random_effect_cholesky_names(term, term$n_columns)
      values <- c(values, stats::setNames(
        as.vector(diag(term$n_columns)), as.vector(coordinates)
      ))
    }
    matrix(rep(values, each = n_draws), n_draws, dimnames = list(NULL, names(values)))
  }

  # 10 coefficient columns, 2 rows, 100 draws: the SD draw matrix of a DIAG
  # block and the Cholesky and correlation draws of a US block exceed what the
  # output alone admitted (9,600 bytes against a 10,000-byte ceiling)
  for(structure in c("diag", "us")){
    result <- compile(structure)
    posterior <- posterior_for(result, cholesky = identical(structure, "us"))
    for(diagonal_only in c(FALSE, TRUE)){
      .random_memory_restore_option(10000)
      expect_error(
        .random_memory_covariance(result, posterior, diagonal_only),
        "10 coefficient columns",
        fixed = TRUE,
        info = paste(structure, diagonal_only)
      )
      .random_memory_restore_option(1e6)
      measured <- .random_memory_measure(
        function() .random_memory_covariance(result, posterior, diagonal_only)
      )
      expect_gt(measured$admitted, 10000)
      if(measured$profile){
        expect_lte(measured$largest, measured$admitted)
      }
    }
  }

  # the components of the estimate
  estimate <- BayesTools:::.bt_random_effect_output_memory_estimate(
    "test", n_rows = 2L, n_draws = 100L, covariance = TRUE,
    n_columns = 10L, factor_structure = TRUE, known_group = TRUE
  )
  expect_identical(
    estimate$components,
    c(output_working = 3 * 8 * 100 * 2 * 2,
      sd_draws_working = 2 * 8 * 100 * 10,
      coefficient_factor_working = 2 * 8 * 100 * 10 * 10,
      coefficient_contraction_working = 4 * 8 * 2 * 10,
      known_group_working = 44 * 2 * 2)
  )
  diagonal <- BayesTools:::.bt_random_effect_output_memory_estimate(
    "test", n_rows = 2L, n_draws = 100L, covariance = TRUE, diagonal_only = TRUE,
    n_columns = 10L, known_group = TRUE
  )
  expect_identical(
    diagonal$components,
    c(output_working = 3 * 8 * 100 * 2,
      sd_draws_working = 2 * 8 * 100 * 10,
      known_group_working = 2 * 8 * 2 * 10)
  )
  expect_equal(estimate$peak_bytes, sum(estimate$components))
  expect_match(estimate$dimensions, "10 coefficient columns", fixed = TRUE)
})
