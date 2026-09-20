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
