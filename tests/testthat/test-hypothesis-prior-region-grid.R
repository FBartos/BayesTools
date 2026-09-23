skip_if_not_test_profile("unit")

test_that("direct prior regions split the continuous interpolant at the threshold", {

  density <- structure(list(
    density = list(x = c(-1, 0, 2), y = c(0, 2 / 3, 0), mass = .7),
    points = data.frame(x = .2, p = .3)
  ), class = c("prior_linear_density", "prior_density"))
  # The triangular density has probability 1/3 + (2/3 + .6) * .2/2 = .46
  # below .2. The atom contributes only when the comparison includes it.
  cases <- c(
    "theta < .2" = .7 * .46,
    "theta <= .2" = .7 * .46 + .3,
    "theta > .2" = .7 * .54,
    "theta >= .2" = .7 * .54 + .3,
    ".2 > theta" = .7 * .46,
    ".2 >= theta" = .7 * .46 + .3,
    ".2 < theta" = .7 * .54,
    ".2 <= theta" = .7 * .54 + .3,
    "theta < -2" = 0,
    "theta > -2" = 1,
    "theta < 3" = 1,
    "theta > 3" = 0,
    "theta < -1" = 0,
    "theta > 2" = 0,
    "theta <= 0" = .7 / 3,
    "theta >= 0" = .7 * 2 / 3 + .3
  )
  for(hypothesis in names(cases)){
    side <- hypothesis_parse(hypothesis)$statements[[1L]]$left
    expect_equal(
      .hypothesis_prior_density_prob(density, side, "theta"),
      unname(cases[[hypothesis]]), tolerance = 1e-14, info = hypothesis
    )
  }
})

test_that("normal and spike-and-slab region grids meet the unchanged refinement gate", {

  for(sd in c(1, .5)){
    scalar_prior <- prior("normal", list(mean = 0, sd = sd))
    density <- .prior_linear_combination_density(
      list(theta = scalar_prior), c(theta = 1), n_grid = 10000
    )
    side <- hypothesis_parse(paste0("theta < ", -.5 * sd))$statements[[1L]]$left
    expect_equal(
      .hypothesis_prior_density_prob(density, side, "theta"),
      stats::pnorm(-.5), tolerance = 2e-7
    )
  }

  mixture <- prior_mixture(list(
    prior("spike", list(location = 0)),
    prior("normal", list(mean = 0, sd = .5))
  ), is_null = c(TRUE, FALSE))
  density <- .prior_linear_combination_density(
    list(theta = mixture), c(theta = 1), n_grid = 10000
  )
  for(cut in c(0, .2)){
    for(op in c("<", "<=", ">", ">=")){
      side <- hypothesis_parse(paste("theta", op, cut))$statements[[1L]]$left
      lower <- op %in% c("<", "<=")
      continuous <- .5 * stats::pnorm(cut, sd = .5, lower.tail = lower)
      atom <- .5 * switch(op, "<" = 0 < cut, "<=" = 0 <= cut,
                          ">" = 0 > cut, ">=" = 0 >= cut)
      expect_equal(
        .hypothesis_prior_density_prob(density, side, "theta"),
        continuous + atom, tolerance = 2e-7, info = paste(op, cut)
      )
    }
  }
})

test_that("compound prior regions insert exact boundary knots on the interpolant", {

  # Triangular continuous part (trapezoid mass 1) plus an atom at .2. The
  # linear interpolant integrates exactly: P(-.5 < theta < 1) = .75,
  # P(theta > .2) = .54, P(theta^2 < .25) = .25 + .2916667 = .5416667.
  density <- structure(list(
    density = list(x = c(-1, 0, 2), y = c(0, 2 / 3, 0), mass = .7),
    points = data.frame(x = .2, p = .3)
  ), class = c("prior_linear_density", "prior_density"))
  cases <- c(
    "theta > -0.5 & theta < 1"  = .7 * .75 + .3,
    "!(theta <= 0.2)"           = .7 * .54,
    "!(theta < 0.2)"            = .7 * .54 + .3,
    "theta^2 < 0.25"            = .7 * (.25 + .5 * (2 / 3 + .5) / 2) + .3,
    "theta < -0.5 | theta > 1"  = .7 * .25,
    "abs(theta - 1) < 3"        = 1,
    "abs(theta - 5) < 1"        = 0
  )
  for(hypothesis in names(cases)){
    side <- hypothesis_parse(hypothesis)$statements[[1L]]$left
    expect_equal(
      .hypothesis_prior_density_prob(density, side, "theta"),
      unname(cases[[hypothesis]]), tolerance = 1e-10, info = hypothesis
    )
  }
})

test_that("region probabilities normalise the grid by its own trapezoid mass", {

  # Riemann normalisation (sum(y) * dx = 1) differs from the trapezoid mass
  # (2 / 3) of this grid; region masses are ratios of trapezoid integrals.
  density <- structure(list(
    density = list(x = c(0, 1, 2), y = c(1, 1, 1) / 3, mass = .6),
    points = data.frame(x = 1.5, p = .4)
  ), class = c("prior_linear_density", "prior_density"))
  cases <- c(
    "theta < 1"               = .6 * .5,
    "theta > 0.5"             = .6 * .75 + .4,
    "abs(theta - 1) < 0.5"    = .6 * .5,
    "theta < 0.5 | theta > 1" = .6 * .75 + .4
  )
  for(hypothesis in names(cases)){
    side <- hypothesis_parse(hypothesis)$statements[[1L]]$left
    expect_equal(
      .hypothesis_prior_density_prob(density, side, "theta"),
      unname(cases[[hypothesis]]), tolerance = 1e-10, info = hypothesis
    )
  }

  # Bounded supports with a positive endpoint ordinate (reference: analytic
  # distribution functions; tolerance 1e-4 is the refinement criterion).
  bounded <- list(
    list(prior = prior("normal", list(0, 1), list(0, Inf)),
         hypothesis = "theta > 0.5", exact = 2 * stats::pnorm(-.5)),
    list(prior = prior("normal", list(0, 1), list(0, Inf)),
         hypothesis = "theta > 0.5 & theta < 1",
         exact = 2 * (stats::pnorm(1) - stats::pnorm(.5))),
    list(prior = prior("uniform", list(0, 1)),
         hypothesis = "abs(theta - 0.5) < 0.25", exact = .5),
    list(prior = prior("exp", list(1)),
         hypothesis = "theta < 1", exact = stats::pexp(1)),
    list(prior = prior("exp", list(1)),
         hypothesis = "!(theta >= 1)", exact = stats::pexp(1))
  )
  for(case in bounded){
    grid <- .prior_linear_combination_density(
      list(theta = case$prior), c(theta = 1)
    )
    side <- hypothesis_parse(case$hypothesis)$statements[[1L]]$left
    expect_equal(
      .hypothesis_prior_density_prob(grid, side, "theta"),
      case$exact, tolerance = 1e-4, info = case$hypothesis
    )
  }
  grid <- .prior_linear_combination_density(
    list(theta = prior("normal", list(0, 1), list(0, Inf))), c(theta = 1)
  )
  side <- hypothesis_parse("theta > 0")$statements[[1L]]$left
  expect_equal(.hypothesis_prior_density_prob(grid, side, "theta"), 1,
               tolerance = 1e-12)
})

test_that("interval, union, and negated regions converge on deterministic prior grids", {

  # Reference: analytic normal probabilities; tolerance 1e-4 is the
  # documented grid-refinement criterion (all cases previously failed it).
  cases <- c(
    "theta > -0.1 & theta < 0.1" = stats::pnorm(.1) - stats::pnorm(-.1),
    "abs(theta) < 0.1"           = stats::pnorm(.1) - stats::pnorm(-.1),
    "theta > 0.13 & theta < 0.47" = stats::pnorm(.47) - stats::pnorm(.13),
    "theta > 1.5 & theta < 3"    = stats::pnorm(3) - stats::pnorm(1.5),
    "theta < -1 | theta > 1"     = 2 * stats::pnorm(-1),
    "theta^2 < 0.04"             = 2 * stats::pnorm(.2) - 1,
    "!(theta <= 0.123)"          = stats::pnorm(.123, lower.tail = FALSE),
    "exp(theta) > 1.1"           = stats::pnorm(log(1.1), lower.tail = FALSE)
  )
  for(n_grid in list(NULL, 128L)){
    grid <- do.call(.prior_linear_combination_density, c(
      list(prior_list = list(theta = prior("normal", list(0, 1))),
           weights = c(theta = 1)),
      if(!is.null(n_grid)) list(n_grid = n_grid)
    ))
    for(hypothesis in names(cases)){
      side <- hypothesis_parse(hypothesis)$statements[[1L]]$left
      expect_equal(
        .hypothesis_prior_density_prob(grid, side, "theta"),
        unname(cases[[hypothesis]]), tolerance = 1e-4,
        info = paste(hypothesis, if(is.null(n_grid)) "default" else n_grid)
      )
    }
  }

  context <- .prior_density_context(
    prior_list = list(a = prior("normal", list(0, 1)),
                      b = prior("normal", list(0, 1))),
    column_names = c("a", "b")
  )
  side <- hypothesis_parse("abs(theta) < 0.25")$statements[[1L]]$left
  expect_equal(
    .hypothesis_prior_density_prob(
      .prior_density_from_context(context, c(a = 1, b = 1)), side, "theta"
    ),
    stats::pnorm(.25, sd = sqrt(2)) - stats::pnorm(-.25, sd = sqrt(2)),
    tolerance = 1e-4
  )

  mixture <- prior_mixture(list(
    prior("spike", list(location = 0)),
    prior("normal", list(0, 1))
  ), is_null = c(TRUE, FALSE))
  grid <- .prior_linear_combination_density(
    list(theta = mixture), c(theta = 1)
  )
  interval <- hypothesis_parse("theta > -0.2 & theta < 0.2")$statements[[1L]]$left
  union <- hypothesis_parse("theta < -0.2 | theta > 0.2")$statements[[1L]]$left
  expect_equal(
    .hypothesis_prior_density_prob(grid, interval, "theta"),
    .5 + .5 * (stats::pnorm(.2) - stats::pnorm(-.2)), tolerance = 1e-4
  )
  expect_equal(
    .hypothesis_prior_density_prob(grid, union, "theta"),
    stats::pnorm(-.2), tolerance = 1e-4
  )

  set.seed(3)
  posterior <- structure(
    stats::rnorm(20000, .05, .1),
    class = c("marginal_posterior.simple", "marginal_posterior", "numeric"),
    prior_density = .prior_linear_combination_density(
      list(theta = prior("normal", list(0, 1))), c(theta = 1)
    ),
    posterior_atoms = posterior_atom_attribute()
  )
  out <- hypothesis_BF(posterior, hypothesis = "abs(theta) < 0.1",
                       parameter = "theta", columns = "all")
  prior_mass <- stats::pnorm(.1) - stats::pnorm(-.1)
  expect_equal(out[["prior"]], prior_mass / (1 - prior_mass), tolerance = 1e-4)
  expect_equal(out[["method"]], "prior-posterior odds")
})
