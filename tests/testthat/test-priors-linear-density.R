skip_if_not_test_profile("unit")

test_that("linear density normalize warns only for large mass deviation", {

  quiet <- BayesTools:::.prior_linear_density_normalize(list(
    density = list(x = 0:1, y = c(1, 1), mass = 1 + 1e-8),
    points = NULL
  ), warn = TRUE)
  expect_equal(quiet$density$mass, 1, tolerance = 1e-12)

  # Intermediate partial masses stay quiet by default.
  expect_silent(
    BayesTools:::.prior_linear_density_normalize(list(
      density = list(x = 0:1, y = c(1, 1), mass = 0.5),
      points = NULL
    ))
  )

  expect_warning(
    noisy <- BayesTools:::.prior_linear_density_normalize(list(
      density = list(x = 0:1, y = c(1, 1), mass = 2.5),
      points = data.frame(x = 0.5, p = 0.5)
    ), warn = TRUE),
    "renormalizing to 1",
    fixed = TRUE
  )
  expect_equal(noisy$density$mass + sum(noisy$points$p), 1, tolerance = 1e-12)

  expect_error(
    BayesTools:::.prior_linear_density_normalize(list(
      density = list(x = 0:1, y = c(0, 0), mass = 0),
      points = NULL
    )),
    "zero total mass",
    fixed = TRUE
  )
})

test_that("linear prior density matches analytic normal sums", {

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(1, 2))
    ),
    weights = c(x = 2, y = -0.5),
    n_grid  = 1024
  )

  x <- density$density$x
  y <- density$density$y
  expected <- stats::dnorm(x, mean = -0.5, sd = sqrt(2^2 + 1^2))
  expected <- expected / (sum(expected) * (x[2] - x[1]))

  expect_equal(
    sum(abs(y - expected)) * (x[2] - x[1]),
    0,
    tolerance = 0.03
  )
  expect_equal(density$density$mass, 1, tolerance = 1e-8)
})

test_that("linear prior ordinates adapt across center and omitted tails", {

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = list(
      x = prior("normal", list(0, 1)),
      y = prior("normal", list(0, 1))
    ),
    weights   = c(x = 1, y = 1),
    n_grid    = 512,
    tail_prob = 1e-3
  )

  center <- BayesTools:::.prior_linear_density_height(density, 0)
  tail <- BayesTools:::.prior_linear_density_height(density, 8)
  expect_lt(
    abs(as.numeric(center) / stats::dnorm(0, sd = sqrt(2)) - 1),
    1e-4
  )
  expect_lt(
    abs(as.numeric(tail) / stats::dnorm(8, sd = sqrt(2)) - 1),
    1e-4
  )
  expect_true(isTRUE(attr(center, "adaptive_evaluation")$converged))
  expect_true(isTRUE(attr(tail, "adaptive_evaluation")$converged))
  expect_gt(attr(tail, "adaptive_evaluation")$refinements, 0)

  diagnostics <- attr(density, "numerical_diagnostics", exact = TRUE)
  expect_equal(diagnostics$tail_probability_per_source, 1e-3)
  expect_equal(
    diagnostics$intended_captured_probability_per_continuous_source,
    .998
  )
  expect_true(is.list(diagnostics$grid_normalization))
  expect_true(is.list(diagnostics$fft_clipping))
})

test_that("linear prior density handles multiply_by products and point mass", {

  priors <- list(
    beta = prior_mixture(
      list(
        prior("spike", list(0), prior_weights = 1),
        prior("normal", list(0, 1), prior_weights = 1)
      ),
      is_null = c(TRUE, FALSE)
    ),
    sigma = prior("cauchy", list(0, 1), list(0, 5))
  )
  attr(priors$beta, "multiply_by") <- "sigma"

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = priors,
    weights    = c(beta = 1),
    n_grid     = 512
  )

  expect_equal(BayesTools:::.prior_linear_density_point_mass(density, 0), 0.5, tolerance = 1e-8)
  expect_gt(BayesTools:::.prior_linear_density_height(density, 0), 0)
})

test_that("product densities count each mixed-measure component once", {

  make_dist <- function(inclusion){
    p <- prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior("spike", list(inclusion))
    )
    BayesTools:::.prior_linear_group_distribution(
      group = list(
        prior = p,
        weights = c(x = 1),
        indices = 1L
      ),
      dx = .02,
      tail_prob = 1e-4,
      source_transforms = c(x = NA_character_),
      n_grid = 256
    )
  }

  lhs <- make_dist(.5)
  rhs <- make_dist(.75)
  product <- BayesTools:::.prior_linear_density_product(
    lhs,
    rhs,
    n_grid = 256
  )

  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(product, 0),
    1 - .5 * .75,
    tolerance = 1e-12
  )
  expect_equal(product$density$mass, .5 * .75, tolerance = 1e-12)

  scaled <- BayesTools:::.prior_linear_density_product(
    lhs,
    BayesTools:::.prior_linear_density_point(.25),
    n_grid = 256
  )
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(scaled, 0),
    .5,
    tolerance = 1e-12
  )
  expect_equal(scaled$density$mass, .5, tolerance = 1e-12)
})

test_that("linear prior density treats zero-weight combinations as point priors", {

  priors <- list(beta = prior("normal", list(0, 1)))

  empty_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = priors,
    weights    = numeric(),
    n_grid     = 128
  )
  expect_equal(BayesTools:::.prior_linear_density_point_mass(empty_density, 0), 1)
  expect_null(empty_density$density)

  zero_density <- BayesTools:::.prior_linear_combination_density(
    prior_list = priors,
    weights    = c(beta = 0),
    n_grid     = 128
  )
  expect_equal(BayesTools:::.prior_linear_density_point_mass(zero_density, 0), 1)
  expect_null(zero_density$density)

  context <- BayesTools:::.prior_density_context(
    prior_list   = priors,
    column_names = "beta",
    n_grid       = 128
  )
  context_density <- BayesTools:::.prior_density_from_context(
    context,
    weights = c(beta = 0)
  )
  expect_equal(BayesTools:::.prior_linear_density_point_mass(context_density, 0), 1)
  expect_null(context_density$density)
})

test_that("linear prior density preserves every finite nonzero coefficient", {

  priors <- list(beta = prior("normal", list(0, 1)))
  tiny_weight <- .Machine$double.eps

  density <- BayesTools:::.prior_linear_combination_density(
    prior_list = priors,
    weights    = c(beta = tiny_weight),
    n_grid     = 256
  )

  expect_false(is.null(density$density))
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(density, 0),
    0
  )
  expect_equal(
    range(density$density$x) / tiny_weight,
    stats::qnorm(c(1e-4, 1 - 1e-4)),
    tolerance = 1e-8
  )

  support <- BayesTools:::.posterior_support_from_prior_list_weights(
    priors,
    c(beta = tiny_weight)
  )
  expect_equal(support$bounds, c(-Inf, Inf))
})


test_that("linear prior density rejects non-finite coefficient weights", {

  priors <- list(
    x = prior("normal", list(0, 1)),
    y = prior("normal", list(0, 1))
  )

  for(bad_weight in c(NA_real_, NaN)){
    expect_error(
      BayesTools:::.prior_linear_combination_density(
        prior_list = priors,
        weights    = c(x = 1, y = bad_weight),
        n_grid     = 128
      ),
      "The 'weights' argument cannot contain NA/NaN values.",
      fixed = TRUE
    )
  }

  for(bad_weight in c(Inf, -Inf)){
    expect_error(
      BayesTools:::.prior_linear_combination_density(
        prior_list = priors,
        weights    = c(x = 1, y = bad_weight),
        n_grid     = 128
      ),
      "The 'weights' argument must contain only finite values.",
      fixed = TRUE
    )
  }
})


test_that("linear prior density honors named scalar source transforms", {
  p <- prior("lognormal", list(0, 1))

  expect_equal(
    BayesTools:::.prior_linear_scalar_range(
      p,
      weight = 1,
      tail_prob = .025,
      source_transform = c(beta = "log")
    ),
    stats::qnorm(c(.025, .975)),
    tolerance = 1e-8
  )
})


test_that("row-wise prior densities mix row predictions, not averaged weights", {
  context <- BayesTools:::.prior_density_context(
    prior_list   = list(beta = prior("normal", list(0, 1))),
    column_names = "beta",
    n_grid       = 1024
  )

  row_density <- BayesTools:::.prior_density_from_context_rows(
    context,
    weights = matrix(c(1, 2), ncol = 1, dimnames = list(NULL, "beta"))
  )
  averaged_density <- BayesTools:::.prior_density_from_context(
    context,
    weights = c(beta = 1.5)
  )

  density_second_moment <- function(d){
    x <- d$density$x
    y <- d$density$y
    dx <- x[2] - x[1]
    sum(x^2 * y) * dx + sum(d$points$x^2 * d$points$p)
  }

  expect_equal(density_second_moment(row_density), mean(c(1^2, 2^2)), tolerance = .08)
  expect_equal(density_second_moment(averaged_density), 1.5^2, tolerance = .08)
})

test_that("row-wise prior densities preserve bitwise-distinct design rows", {

  context <- BayesTools:::.prior_density_context(
    prior_list = list(beta = prior("point", list(location = 1))),
    column_names = "beta",
    n_grid = 64
  )
  nearby <- 1 + .Machine$double.eps
  density <- BayesTools:::.prior_density_from_context_rows(
    context,
    weights = matrix(
      c(1, nearby),
      ncol = 1,
      dimnames = list(NULL, "beta")
    )
  )

  expect_identical(density$points$x, c(1, nearby))
  expect_equal(density$points$p, c(.5, .5))
})

test_that("plot_transformed_prior exposes transformed prior plotting as a public wrapper", {

  prior_list <- list(
    mu_intercept = prior("normal", list(0, 1)),
    mu_x = prior_mixture(
      list(
        prior("spike", list(0), prior_weights = 1),
        prior("normal", list(1, .25), prior_weights = 3)
      ),
      is_null = c(TRUE, FALSE)
    )
  )
  formula_scale <- list(mu = list(mu_x = list(mean = 5, sd = 2)))
  attr(formula_scale$mu, "log_intercept") <- FALSE

  plot <- plot_transformed_prior(
    prior_list    = prior_list,
    column_names  = c("mu_intercept", "mu_x"),
    formula_scale = formula_scale,
    parameter     = "mu_x",
    n_points      = 128,
    plot_type     = "ggplot",
    par_name      = "x"
  )

  expect_s3_class(plot, "ggplot")
  expect_true(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1))))
  expect_true(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomSegment"), logical(1))))

  device_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(device_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot_transformed_prior(
    prior_list    = prior_list,
    column_names  = c("mu_intercept", "mu_x"),
    formula_scale = formula_scale,
    parameter     = "mu_x",
    n_points      = 64,
    plot_type     = "base"
  ))
})

test_that("plot_transformed_prior returns NULL only for untransformed identity parameters", {

  prior_list <- list(
    mu_intercept = prior("normal", list(0, 1)),
    mu_x         = prior("normal", list(0, 1))
  )

  expect_null(plot_transformed_prior(
    prior_list   = prior_list,
    column_names = c("mu_intercept", "mu_x"),
    parameter    = "mu_x",
    plot_type    = "ggplot"
  ))

  transformed_plot <- plot_transformed_prior(
    prior_list     = prior_list,
    column_names   = c("mu_intercept", "mu_x"),
    parameter      = "mu_x",
    transformation = "exp",
    x_range        = c(-2, 2),
    n_points       = 64,
    plot_type      = "ggplot"
  )

  expect_s3_class(transformed_plot, "ggplot")
  expect_error(
    plot_transformed_prior(
      prior_list   = prior_list,
      column_names = c("mu_intercept", "mu_x"),
      parameter    = "mu_z"
    ),
    "Parameter 'mu_z' was not found in 'column_names'.",
    fixed = TRUE
  )
})
