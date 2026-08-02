skip_if_not_test_profile("unit")

test_that("prior density ordinate validates input and has a stable schema", {

  normal_prior <- prior("normal", list(mean = 0, sd = 1))
  out <- prior_density_ordinate(normal_prior, 0)

  expect_s3_class(out, "prior_density_ordinate")
  expect_identical(
    names(unclass(out)),
    c(
      "schema_version", "value", "behavior", "log_density", "point_mass",
      "exact", "method", "reason", "provenance"
    )
  )
  expect_identical(out$schema_version, "1")
  expect_identical(out$value, 0)
  expect_identical(out$behavior, "regular")
  expect_equal(out$log_density, stats::dnorm(0, log = TRUE))
  expect_identical(out$point_mass, 0)
  expect_true(out$exact)
  expect_identical(out$method, "primitive")
  expect_null(out$reason)
  expect_true(out$method %in% c(
    "primitive", "point", "finite_mixture", "scalar_affine",
    "linear_normal", "named_transform", "unsupported_provenance"
  ))

  expect_error(
    prior_density_ordinate(normal_prior, NA_real_),
    "cannot contain NA/NaN",
    fixed = TRUE
  )
  expect_error(
    prior_density_ordinate(normal_prior, Inf),
    "must be finite",
    fixed = TRUE
  )
  expect_error(
    prior_density_ordinate(normal_prior, c(0, 1)),
    "must have length '1'",
    fixed = TRUE
  )
  expect_error(
    prior_density_ordinate(list(), 0),
    "must be a BayesTools prior or prior_linear_density object",
    fixed = TRUE
  )
})

test_that("primitive ordinates distinguish support and boundary limits", {

  truncated_normal <- prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
  expect_identical(
    prior_density_ordinate(truncated_normal, -1)$behavior,
    "zero"
  )
  expect_identical(
    prior_density_ordinate(truncated_normal, 0)$behavior,
    "regular"
  )

  finite_boundary <- prior_density_ordinate(
    prior("exp", list(rate = 2)),
    0
  )
  expect_identical(finite_boundary$behavior, "regular")
  expect_equal(finite_boundary$log_density, log(2))

  zero_boundary <- prior_density_ordinate(
    prior("gamma", list(shape = 2, rate = 1)),
    0
  )
  expect_identical(zero_boundary$behavior, "zero")
  expect_identical(zero_boundary$log_density, -Inf)

  infinite_boundary <- prior_density_ordinate(
    prior("gamma", list(shape = 0.5, rate = 1)),
    0
  )
  expect_identical(infinite_boundary$behavior, "infinite")
  expect_identical(infinite_boundary$log_density, Inf)

  expect_identical(
    prior_density_ordinate(prior("lognormal", list(0, 1)), 0)$behavior,
    "zero"
  )
  expect_identical(
    prior_density_ordinate(prior("beta", list(0.5, 2)), 0)$behavior,
    "infinite"
  )
  expect_identical(
    prior_density_ordinate(prior("beta", list(2, 1)), 1)$behavior,
    "regular"
  )

  uniform_prior <- prior("uniform", list(a = -2, b = 3))
  expect_identical(
    prior_density_ordinate(uniform_prior, -2)$behavior,
    "regular"
  )
  expect_identical(
    prior_density_ordinate(uniform_prior, 3)$behavior,
    "regular"
  )
})

test_that("nonlocal primitive zeros and exact discrete masses are structural", {

  expect_identical(
    prior_density_ordinate(
      prior("moment", list(mode = 0.5, location = 0.25)),
      0.25
    )$behavior,
    "zero"
  )
  expect_identical(
    prior_density_ordinate(
      prior("invmoment", list(mode = 0.5, df = 3, location = -0.25)),
      -0.25
    )$behavior,
    "zero"
  )

  point_at_zero <- prior_density_ordinate(prior("point", list(0)), 0)
  expect_identical(point_at_zero$behavior, "point_mass")
  expect_identical(point_at_zero$point_mass, 1)
  expect_identical(point_at_zero$log_density, -Inf)

  point_elsewhere <- prior_density_ordinate(prior("point", list(0)), 1)
  expect_identical(point_elsewhere$behavior, "zero")
  expect_identical(point_elsewhere$point_mass, 0)

  bernoulli <- prior_density_ordinate(
    prior("bernoulli", list(probability = 0.3)),
    1
  )
  expect_identical(bernoulli$behavior, "point_mass")
  expect_equal(bernoulli$point_mass, 0.3)

  none <- prior_density_ordinate(prior_none(), 0)
  expect_identical(none$behavior, "point_mass")
  expect_identical(none$point_mass, 1)
  expect_identical(none$method, "point")

  vector_prior <- prior("mnormal", list(mean = 0, sd = 1, K = 2))
  unsupported <- prior_density_ordinate(vector_prior, 0)
  expect_identical(unsupported$behavior, "unknown")
  expect_false(unsupported$exact)
})

test_that("finite mixtures combine continuous behavior and atoms exactly", {

  make_mixture <- function(priors){
    prior_mixture(
      priors,
      is_null = rep(FALSE, length(priors))
    )
  }

  regular <- make_mixture(list(
    prior("normal", list(0, 1), prior_weights = 1),
    prior("gamma", list(2, 1), prior_weights = 1)
  ))
  regular_out <- prior_density_ordinate(regular, 0)
  expect_identical(regular_out$behavior, "regular")
  expect_equal(regular_out$log_density, log(0.5) + stats::dnorm(0, log = TRUE))

  zero <- make_mixture(list(
    prior("gamma", list(2, 1), prior_weights = 1),
    prior("beta", list(2, 2), prior_weights = 1)
  ))
  expect_identical(prior_density_ordinate(zero, 0)$behavior, "zero")

  infinite <- make_mixture(list(
    prior("gamma", list(0.5, 1), prior_weights = 1),
    prior("normal", list(0, 1), prior_weights = 1)
  ))
  expect_identical(
    prior_density_ordinate(infinite, 0)$behavior,
    "infinite"
  )

  mixed_measure <- make_mixture(list(
    prior("point", list(0), prior_weights = 1),
    prior("normal", list(0, 1), prior_weights = 3)
  ))
  mixed_out <- prior_density_ordinate(mixed_measure, 0)
  expect_identical(mixed_out$behavior, "point_mass")
  expect_equal(mixed_out$point_mass, 0.25)
  expect_equal(
    mixed_out$log_density,
    log(0.75) + stats::dnorm(0, log = TRUE)
  )
  expect_identical(
    mixed_out$provenance$continuous_behavior,
    "regular"
  )

  point_and_infinite <- make_mixture(list(
    prior("point", list(0), prior_weights = 1),
    prior("gamma", list(0.5, 1), prior_weights = 1)
  ))
  point_and_infinite_out <- prior_density_ordinate(point_and_infinite, 0)
  expect_identical(point_and_infinite_out$behavior, "point_mass")
  expect_equal(point_and_infinite_out$point_mass, 0.5)
  expect_identical(
    point_and_infinite_out$provenance$continuous_behavior,
    "infinite"
  )

  zero_weight <- make_mixture(list(
    prior("normal", list(0, 1), prior_weights = 1),
    prior("gamma", list(0.5, 1), prior_weights = 1)
  ))
  attr(zero_weight, "prior_weights") <- c(1, 0)
  zero_weight_out <- prior_density_ordinate(zero_weight, 0)
  expect_identical(zero_weight_out$behavior, "regular")
  expect_equal(zero_weight_out$log_density, stats::dnorm(0, log = TRUE))

  spike_and_slab <- prior_spike_and_slab(
    prior("normal", list(0, 1)),
    prior_inclusion = prior("point", list(0.25))
  )
  spike_and_slab_out <- prior_density_ordinate(spike_and_slab, 0)
  expect_identical(spike_and_slab_out$behavior, "point_mass")
  expect_equal(spike_and_slab_out$point_mass, 0.75)
  expect_equal(
    spike_and_slab_out$log_density,
    log(0.25) + stats::dnorm(0, log = TRUE)
  )
})

test_that("finite-mixture structural precedence is deterministic", {

  result <- BayesTools:::.prior_density_ordinate_result
  undefined <- result(
    0, "undefined", NA_real_, method = "named_transform"
  )
  infinite <- result(
    0, "infinite", Inf, method = "primitive"
  )
  unknown <- result(
    0, "unknown", NA_real_, exact = FALSE,
    method = "unsupported_provenance"
  )

  expect_identical(
    BayesTools:::.prior_density_ordinate_combine(
      list(undefined, infinite),
      c(1, 1),
      0
    )$behavior,
    "undefined"
  )
  unknown_and_infinite <- BayesTools:::.prior_density_ordinate_combine(
    list(unknown, infinite),
    c(1, 1),
    0
  )
  expect_identical(unknown_and_infinite$behavior, "infinite")
  expect_identical(unknown_and_infinite$log_density, Inf)
  expect_true(unknown_and_infinite$exact)
})

test_that("nonzero scalar affine transformations retain behavior and Jacobian", {

  normal_prior <- prior("normal", list(mean = 1, sd = 2))
  positive <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = normal_prior),
    weights = c(beta = 2),
    n_grid = 128
  )
  negative <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = normal_prior),
    weights = c(beta = -2),
    n_grid = 128
  )

  positive_out <- prior_density_ordinate(positive, 2)
  negative_out <- prior_density_ordinate(negative, -2)
  expect_identical(positive_out$behavior, "regular")
  expect_identical(negative_out$behavior, "regular")
  expect_identical(positive_out$method, "scalar_affine")
  expect_identical(negative_out$method, "scalar_affine")
  expect_equal(
    positive_out$log_density,
    stats::dnorm(1, mean = 1, sd = 2, log = TRUE) - log(2)
  )
  expect_equal(negative_out$log_density, positive_out$log_density)

  point <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = prior("point", list(location = 0.1))),
    weights = c(beta = 3),
    n_grid = 128
  )
  point_out <- prior_density_ordinate(point, point$points$x)
  expect_identical(point_out$behavior, "point_mass")
  expect_identical(point_out$point_mass, 1)

  zero <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = normal_prior),
    weights = c(beta = 0),
    n_grid = 128
  )
  zero_out <- prior_density_ordinate(zero, 0)
  expect_identical(zero_out$behavior, "point_mass")
  expect_identical(zero_out$point_mass, 1)

  affine_mixture <- prior_mixture(
    list(
      prior("point", list(0), prior_weights = 1),
      prior("normal", list(0, 1), prior_weights = 1)
    ),
    is_null = c(TRUE, FALSE)
  )
  affine_mixture_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = affine_mixture),
    weights = c(beta = 2),
    n_grid = 128
  )
  affine_mixture_out <- prior_density_ordinate(affine_mixture_density, 0)
  expect_identical(affine_mixture_out$behavior, "point_mass")
  expect_equal(affine_mixture_out$point_mass, 0.5)
  expect_identical(
    affine_mixture_out$provenance$continuous_behavior,
    "regular"
  )
})

test_that("normal linear combinations are classified analytically", {

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(1, 2))
    ),
    weights = c(x = 2, y = -0.5),
    n_grid = 128
  )
  out <- prior_density_ordinate(density, 0)
  expect_identical(out$behavior, "regular")
  expect_identical(out$method, "linear_normal")
  expect_equal(
    out$log_density,
    stats::dnorm(0, mean = -0.5, sd = sqrt(5), log = TRUE)
  )
  expect_identical(out$provenance$weights, c(x = 2, y = -0.5))

  vector_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      beta = prior("mnormal", list(mean = 1, sd = 2, K = 2))
    ),
    weights = c("beta[1]" = 1, "beta[2]" = -0.5),
    n_grid = 128
  )
  vector_out <- prior_density_ordinate(vector_density, 0)
  expect_identical(vector_out$behavior, "regular")
  expect_identical(vector_out$method, "linear_normal")
  expect_equal(
    vector_out$log_density,
    stats::dnorm(0, mean = 0.5, sd = sqrt(5), log = TRUE)
  )

  shifted <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      shift = prior("point", list(1.5)),
      x = prior("normal", list(1, 2)),
      y = prior("normal", list(-1, 1))
    ),
    weights = c(shift = 2, x = -0.5, y = 2),
    n_grid = 128
  )
  shifted_out <- prior_density_ordinate(shifted, 0.5)
  expect_identical(shifted_out$behavior, "regular")
  expect_identical(shifted_out$method, "linear_normal")
  expect_equal(
    shifted_out$log_density,
    stats::dnorm(0.5, mean = 0.5, sd = sqrt(5), log = TRUE)
  )
})

test_that("density contexts retain structural ordinate provenance", {

  context <- BayesTools:::.prior_density_context(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(1, 2))
    ),
    column_names = c("x", "y"),
    n_grid = 128
  )
  density <- BayesTools:::.prior_density_from_context(
    context,
    weights = c(x = 2, y = -0.5)
  )
  out <- prior_density_ordinate(density, 0)
  expect_identical(out$behavior, "regular")
  expect_identical(out$method, "linear_normal")
  expect_identical(out$provenance$context$kind, "prior_density_context")
  expect_identical(
    out$provenance$context$standardized_weights,
    c(x = 2, y = -0.5)
  )

  row_density <- BayesTools:::.prior_density_from_context_rows(
    context,
    weights = matrix(
      c(1, 0, 2, 0),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("x", "y"))
    )
  )
  row_out <- prior_density_ordinate(row_density, 0)
  expect_identical(row_out$behavior, "regular")
  expect_identical(row_out$method, "finite_mixture")
  expect_identical(row_out$provenance$unique_rows, 2L)
  expect_identical(row_out$provenance$weight_dimensions, c(2L, 2L))
  expect_match(row_out$provenance$weights_hash, "^[0-9a-f]{8}$")
  expect_false("representative_weights" %in% names(row_out$provenance))

  model_context <- BayesTools:::.prior_density_model_mixture_context(
    prior_list = list(
      x = list(
        prior("point", list(0), prior_weights = 1),
        prior("normal", list(0, 1), prior_weights = 3)
      )
    ),
    column_names = "x",
    n_grid = 128
  )
  model_density <- BayesTools:::.prior_density_from_context(
    model_context,
    weights = c(x = 1)
  )
  model_out <- prior_density_ordinate(model_density, 0)
  expect_identical(model_out$behavior, "point_mass")
  expect_equal(model_out$point_mass, 0.25)
  expect_equal(
    model_out$log_density,
    log(0.75) + stats::dnorm(0, log = TRUE)
  )
  expect_identical(model_out$provenance$context, "model_mixture")
})

test_that("named monotone transformations classify interiors and boundaries", {

  normal_prior <- list(beta = prior("normal", list(0, 1)))
  exponential <- BayesTools:::.prior_linear_combination_density(
    prior_list = normal_prior,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "exp"
  )
  exp_interior <- prior_density_ordinate(exponential, 1)
  expect_identical(exp_interior$behavior, "regular")
  expect_identical(exp_interior$method, "named_transform")
  expect_equal(exp_interior$log_density, stats::dnorm(0, log = TRUE))
  expect_identical(
    prior_density_ordinate(exponential, 0)$behavior,
    "zero"
  )
  expect_identical(
    prior_density_ordinate(exponential, -1)$behavior,
    "zero"
  )

  hyperbolic <- BayesTools:::.prior_linear_combination_density(
    prior_list = normal_prior,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "tanh"
  )
  expect_equal(
    prior_density_ordinate(hyperbolic, 0)$log_density,
    stats::dnorm(0, log = TRUE)
  )
  expect_identical(
    prior_density_ordinate(hyperbolic, 1)$behavior,
    "zero"
  )
  expect_identical(
    prior_density_ordinate(hyperbolic, 2)$behavior,
    "zero"
  )

  linear <- BayesTools:::.prior_linear_combination_density(
    prior_list = normal_prior,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "lin",
    output_transformation_arguments = list(a = 1, b = -2)
  )
  expect_equal(
    prior_density_ordinate(linear, 1)$log_density,
    stats::dnorm(0, log = TRUE) - log(2)
  )

  positive_prior <- list(
    beta = prior(
      "lognormal",
      list(0, 1),
      truncation = list(lower = 0.01, upper = Inf)
    )
  )
  exponential_linear <- BayesTools:::.prior_linear_combination_density(
    prior_list = positive_prior,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "exp_lin",
    output_transformation_arguments = list(a = 0, b = 2)
  )
  expect_identical(
    prior_density_ordinate(exponential_linear, 1)$behavior,
    "regular"
  )
  expect_identical(
    prior_density_ordinate(exponential_linear, 0)$behavior,
    "zero"
  )

  constant_linear <- BayesTools:::.prior_linear_combination_density(
    prior_list = normal_prior,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "lin",
    output_transformation_arguments = list(a = 2, b = 0)
  )
  constant_linear_out <- prior_density_ordinate(constant_linear, 2)
  expect_identical(constant_linear_out$behavior, "point_mass")
  expect_identical(constant_linear_out$point_mass, 1)

  constant_exp_linear <- BayesTools:::.prior_linear_combination_density(
    prior_list = normal_prior,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "exp_lin",
    output_transformation_arguments = list(a = 2, b = 0)
  )
  constant_exp_linear_out <- prior_density_ordinate(
    constant_exp_linear,
    exp(2)
  )
  expect_identical(constant_exp_linear_out$behavior, "point_mass")
  expect_identical(constant_exp_linear_out$point_mass, 1)
})

test_that("composed named-transform boundary limits use source provenance", {

  lognormal_power <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("lognormal", list(0, 1))),
    weights = c(x = 2),
    source_transforms = c(x = "log"),
    n_grid = 128,
    output_transformation = "exp"
  )
  expect_identical(
    prior_density_ordinate(lognormal_power, 0)$behavior,
    "zero"
  )

  gamma_power <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("gamma", list(1, 1))),
    weights = c(x = 2),
    source_transforms = c(x = "log"),
    n_grid = 128,
    output_transformation = "exp"
  )
  expect_identical(
    prior_density_ordinate(gamma_power, 0)$behavior,
    "infinite"
  )

  lognormal_tanh <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("lognormal", list(0, 1))),
    weights = c(x = 1),
    source_transforms = c(x = "log"),
    n_grid = 128,
    output_transformation = "tanh"
  )
  expect_identical(
    prior_density_ordinate(lognormal_tanh, 1)$behavior,
    "zero"
  )

  student_tanh <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("t", list(0, 1, 3))),
    weights = c(x = 1),
    n_grid = 128
  )
  student_adaptive <- attr(student_tanh, "adaptive_evaluation", exact = TRUE)
  student_adaptive$arguments$output_transformation <- "tanh"
  attr(student_tanh, "adaptive_evaluation") <- student_adaptive
  expect_identical(
    prior_density_ordinate(student_tanh, 1)$behavior,
    "infinite"
  )

  positive_t <- prior(
    "t",
    list(location = 0, scale = 1, df = 3),
    truncation = list(lower = 0, upper = Inf)
  )
  dependent_boundary <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = positive_t),
    weights = c(x = 1),
    n_grid = 128
  )
  dependent_adaptive <- attr(
    dependent_boundary,
    "adaptive_evaluation",
    exact = TRUE
  )
  dependent_adaptive$arguments$output_transformation <- "exp_lin"
  dependent_adaptive$arguments$output_transformation_arguments <-
    list(a = 0, b = 1)
  attr(dependent_boundary, "adaptive_evaluation") <- dependent_adaptive
  dependent_out <- prior_density_ordinate(dependent_boundary, 0)
  expect_identical(dependent_out$behavior, "unknown")
  expect_false(dependent_out$exact)
})

test_that("unsupported transformations, convolutions, and products stay unknown", {

  custom_identity <- list(
    fun = function(x) x,
    inv = function(x) x,
    jac = function(x) rep(1, length(x))
  )
  custom <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = prior("normal", list(0, 1))),
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = custom_identity
  )
  custom_out <- prior_density_ordinate(custom, 0)
  expect_identical(custom_out$behavior, "unknown")
  expect_false(custom_out$exact)

  convolution <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("gamma", list(2, 1)),
      y = prior("gamma", list(2, 1))
    ),
    weights = c(x = 1, y = 1),
    n_grid = 128
  )
  convolution_out <- prior_density_ordinate(convolution, 1)
  expect_identical(convolution_out$behavior, "unknown")
  expect_false(convolution_out$exact)
  expect_true(is.finite(convolution_out$log_density))

  product_priors <- list(
    beta = prior("normal", list(0, 1)),
    sigma = prior("normal", list(0, 1))
  )
  attr(product_priors$beta, "multiply_by") <- "sigma"
  product <- BayesTools:::.prior_linear_combination_density(
    prior_list = product_priors,
    weights = c(beta = 1),
    n_grid = 128
  )
  expect_identical(prior_density_ordinate(product, 1)$behavior, "unknown")
  singular <- prior_density_ordinate(product, 0)
  expect_identical(singular$behavior, "infinite")
  expect_identical(singular$method, "unsupported_provenance")

  transformed_product <- BayesTools:::.prior_linear_combination_density(
    prior_list = product_priors,
    weights = c(beta = 1),
    n_grid = 128,
    output_transformation = "exp"
  )
  transformed_singular <- prior_density_ordinate(transformed_product, 1)
  expect_identical(transformed_singular$behavior, "infinite")
  expect_identical(transformed_singular$method, "named_transform")
})

test_that("invalid named transformation provenance is undefined", {

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(beta = prior("normal", list(0, 1))),
    weights = c(beta = 1),
    n_grid = 128
  )
  adaptive <- attr(density, "adaptive_evaluation", exact = TRUE)
  adaptive$arguments$output_transformation <- "exp_lin"
  adaptive$arguments$output_transformation_arguments <- list(a = 0, b = 1)
  attr(density, "adaptive_evaluation") <- adaptive

  out <- prior_density_ordinate(density, 1)
  expect_identical(out$behavior, "undefined")
  expect_true(out$exact)
  expect_identical(out$log_density, NA_real_)

  invalid_arguments <- adaptive
  invalid_arguments$arguments$output_transformation <- "lin"
  invalid_arguments$arguments$output_transformation_arguments <-
    list(a = Inf, b = 1)
  attr(density, "adaptive_evaluation") <- invalid_arguments
  invalid_out <- prior_density_ordinate(density, 1)
  expect_identical(invalid_out$behavior, "undefined")

  atom <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("point", list(0))),
    weights = c(x = 1),
    n_grid = 128
  )
  atom_adaptive <- attr(atom, "adaptive_evaluation", exact = TRUE)
  atom_adaptive$arguments$output_transformation <- "exp_lin"
  atom_adaptive$arguments$output_transformation_arguments <-
    list(a = 0, b = 1)
  attr(atom, "adaptive_evaluation") <- atom_adaptive
  atom_out <- prior_density_ordinate(atom, 1)
  expect_identical(atom_out$behavior, "undefined")
})

test_that("ordinary density underflow does not imply a structural zero", {

  normal_prior <- prior("normal", list(0, 1))
  expect_identical(pdf(normal_prior, 40), 0)

  out <- prior_density_ordinate(normal_prior, 40)
  expect_identical(out$behavior, "regular")
  expect_true(out$exact)
  expect_equal(out$log_density, stats::dnorm(40, log = TRUE))

  extreme <- prior_density_ordinate(normal_prior, 1e200)
  expect_identical(extreme$behavior, "regular")
  expect_identical(extreme$log_density, -Inf)
  expect_match(extreme$reason, "structurally regular")

  log_source <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("lognormal", list(0, 1))),
    weights = c(x = 1),
    source_transforms = c(x = "log"),
    n_grid = 128
  )
  log_source_extreme <- prior_density_ordinate(log_source, 1000)
  expect_identical(log_source_extreme$behavior, "regular")
  expect_identical(log_source_extreme$log_density, -Inf)
  expect_match(log_source_extreme$reason, "structurally regular")

  normal_sum <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(0, 1))
    ),
    weights = c(x = 1, y = 1),
    n_grid = 128
  )
  normal_sum_extreme <- prior_density_ordinate(normal_sum, 1e200)
  expect_identical(normal_sum_extreme$behavior, "regular")
  expect_identical(normal_sum_extreme$log_density, -Inf)
  expect_match(normal_sum_extreme$reason, "structurally regular")
})

test_that("legacy numerical metadata never establishes structural behavior", {

  stale <- structure(
    list(
      density = list(x = c(-1, 0, 1), y = c(1, 0, 1), mass = 1),
      points = data.frame(x = numeric(), p = numeric()),
      n_grid = 3L
    ),
    class = c("prior_linear_density", "prior_density")
  )
  attr(stale, "singular_density_points") <- 0
  attr(stale, "fft_clipping") <- list(clipped_value_count = 1L)
  attr(stale, "density_evaluator") <- function(value) rep(0, length(value))

  stale_zero <- prior_density_ordinate(stale, 0)
  expect_identical(stale_zero$behavior, "unknown")
  expect_false(stale_zero$exact)
  expect_identical(stale_zero$log_density, -Inf)

  stale_outside <- prior_density_ordinate(stale, 100)
  expect_identical(stale_outside$behavior, "unknown")

  structural <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(x = prior("normal", list(0, 1))),
    weights = c(x = 1),
    n_grid = 128
  )
  structural$density$y[] <- 0
  expect_identical(
    prior_density_ordinate(structural, 0)$behavior,
    "regular"
  )
  expect_identical(
    prior_density_ordinate(structural, 100)$behavior,
    "regular"
  )
})

test_that("ordinate provenance is compact and contains no closures", {

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(1, 2))
    ),
    weights = c(x = 2, y = -0.5),
    n_grid = 128
  )
  out <- prior_density_ordinate(density, 0)

  contains_function <- function(x){
    if(is.function(x)) return(TRUE)
    if(!is.list(x)) return(FALSE)
    any(vapply(x, contains_function, logical(1)))
  }

  expect_false(contains_function(unclass(out)))
  expect_lt(length(serialize(unclass(out), NULL)), 20000)
  expect_false("adaptive_evaluation" %in% names(out$provenance))
  expect_false("density_evaluator" %in% names(out$provenance))

  context <- BayesTools:::.prior_density_context(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(1, 2))
    ),
    column_names = c("x", "y"),
    n_grid = 128
  )
  row_density <- BayesTools:::.prior_density_from_context_rows(
    context,
    weights = matrix(
      c(1, 0, 2, 0, 3, 0),
      ncol = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("x", "y"))
    )
  )
  row_out <- prior_density_ordinate(row_density, 0)
  contains_matrix <- function(x){
    if(is.matrix(x)) return(TRUE)
    if(!is.list(x)) return(FALSE)
    any(vapply(x, contains_matrix, logical(1)))
  }
  forbidden_names <- function(x){
    if(!is.list(x)) return(character())
    c(
      intersect(
        names(x),
        c("grid", "density", "draws", "samples", "raw_weights",
          "representative_weights", "adaptive_evaluation",
          "density_evaluator")
      ),
      unlist(lapply(x, forbidden_names), use.names = FALSE)
    )
  }
  expect_false(contains_matrix(row_out$provenance))
  expect_length(forbidden_names(row_out$provenance), 0)
})
