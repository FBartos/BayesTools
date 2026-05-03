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
