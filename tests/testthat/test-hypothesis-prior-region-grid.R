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

# The grid evaluation alone, which region probabilities use when the prior
# density has no structural region probability or the region is not a union
# of intervals; densities with one are also checked through it.
grid_region_probability <- function(density, side, parameter = "theta"){
  .hypothesis_prior_density_grid_prob(
    density, .hypothesis_simple_parameter_comparison(side, parameter),
    .hypothesis_side_expression(side), parameter
  )
}

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
    expect_equal(
      grid_region_probability(density, side),
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
      expect_equal(
        grid_region_probability(density, side),
        continuous + atom, tolerance = 2e-7, info = paste("grid", op, cut)
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
    expect_equal(
      grid_region_probability(grid, side),
      case$exact, tolerance = 1e-4, info = paste("grid", case$hypothesis)
    )
  }
  grid <- .prior_linear_combination_density(
    list(theta = prior("normal", list(0, 1), list(0, Inf))), c(theta = 1)
  )
  side <- hypothesis_parse("theta > 0")$statements[[1L]]$left
  expect_equal(.hypothesis_prior_density_prob(grid, side, "theta"), 1,
               tolerance = 1e-12)
  expect_equal(grid_region_probability(grid, side), 1, tolerance = 1e-12)
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
      expect_equal(
        grid_region_probability(grid, side),
        unname(cases[[hypothesis]]), tolerance = 1e-4,
        info = paste("grid", hypothesis, if(is.null(n_grid)) "default" else n_grid)
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
  expect_equal(
    grid_region_probability(grid, interval),
    .5 + .5 * (stats::pnorm(.2) - stats::pnorm(-.2)), tolerance = 1e-4
  )
  expect_equal(
    grid_region_probability(grid, union),
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

# Independent references for the conditional-normal region probabilities
# (script .work/tmp/pr58-decisions/I/references.py of the BayesToolsVerse
# workspace): 30-digit mpmath. Gaussian convolutions X = G + w T integrate
# over the Gaussian variable g with the other term's exact CDF,
# P(l < X < u) = int phi(g) P(l - g < w T < u - g) dg, not over the package's
# integration variable; scale mixtures X = A + B S use mpmath tanh-sinh
# quadrature over s split at the conditional peaks. Probabilities are checked
# to 1e-8 absolute.
region_references <- list(
  normal_cauchy = c(
    "theta > 0" = .5,
    "theta > 0.5" = .34029874890903886948,
    "theta < -0.2" = .43429003358716276566,
    "theta > -0.1 & theta < 0.1" = .065979610218260501624,
    "theta > 3" = .032113597487249015036,
    "theta > 2 & theta < 4" = .04889811651364013153,
    "theta < -10" = .0080380756664539526223,
    "theta < -1 | theta > 1" = .41905203304630432791
  ),
  normal_t = c(
    "theta > 0" = .5,
    "theta > 0.5" = .31947220009838307904,
    "theta < -0.2" = .42553028174440189795,
    "theta > -0.1 & theta < 0.1" = .074802985908433690396,
    "theta > 3" = .0035329876235075740578,
    "theta > 2 & theta < 4" = .031341700342165305408,
    "theta < -10" = .00001830304247045283063,
    "theta < -1 | theta > 1" = .34899711873701399647
  )
)

region_probability <- function(density, hypothesis){
  side <- hypothesis_parse(hypothesis)$statements[[1L]]$left
  .hypothesis_prior_density_prob(density, side, "theta")
}

expect_region_probabilities <- function(density, references, tolerance = 1e-8){
  for(hypothesis in names(references)){
    probability <- region_probability(density, hypothesis)
    expect_lt(abs(probability - references[[hypothesis]]), tolerance,
              label = hypothesis)
  }
}

test_that("Gaussian-convolution region probabilities use the conditional-normal quadrature", {

  # Original-scale intercept b0 - (m / s) b1 of a scaled formula (m = 1,
  # s = 2) with a normal intercept and a Cauchy or t slope; these regions
  # previously stopped with grid non-convergence or were off by up to 7e-7.
  intercept <- prior("normal", list(0, 1))
  slopes <- list(normal_cauchy = prior("cauchy", list(0, .5)),
                 normal_t = prior("t", list(0, .5, 3)))
  for(case in names(slopes)){
    density <- .prior_linear_combination_density(
      list(b0 = intercept, b1 = slopes[[case]]), c(b0 = 1, b1 = -.5)
    )
    expect_region_probabilities(density, region_references[[case]])
  }

  # the same target through the formula-scale density context
  context <- .prior_density_context(
    prior_list = list(mu_intercept = intercept, mu_x = slopes$normal_cauchy),
    column_names = c("mu_intercept", "mu_x"),
    formula_scale = list(mu = list(mu_x = list(mean = 1, sd = 2)))
  )
  expect_region_probabilities(
    .prior_density_from_context(context, c(mu_intercept = 1)),
    region_references$normal_cauchy
  )

  # scale-disparate terms: regions beyond and around the narrow Gaussian peak
  disparate <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, .01)), b = prior("cauchy", list(0, .5))),
    c(a = 1, b = 1)
  )
  expect_region_probabilities(disparate, c(
    "theta > 5" = .031725642247176749714,
    "theta > 0.001" = .49936363541736656853,
    "theta < -3" = .052569014758854067472
  ))
  disparate <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, .01)), b = prior("cauchy", list(0, 1))),
    c(a = 1, b = 1)
  )
  expect_region_probabilities(disparate, c(
    "theta > -0.01 & theta < 0.01" = .0063653492281441294259
  ))

  # normal sums keep their closed form
  normal <- .prior_linear_combination_density(
    list(b0 = intercept, b1 = prior("normal", list(0, .5))), c(b0 = 1, b1 = -.5)
  )
  expect_equal(region_probability(normal, "theta > 0.5"),
               stats::pnorm(.5, sd = sqrt(1 + .25^2), lower.tail = FALSE),
               tolerance = 1e-14)
  expect_equal(region_probability(normal, "theta > 4 & theta < 5"),
               stats::pnorm(4, sd = sqrt(1 + .25^2), lower.tail = FALSE) -
                 stats::pnorm(5, sd = sqrt(1 + .25^2), lower.tail = FALSE),
               tolerance = 1e-12)

  # region odds end to end
  set.seed(1)
  posterior <- structure(
    stats::rnorm(20000, .3, .2),
    class = c("marginal_posterior.simple", "marginal_posterior", "numeric"),
    prior_density = .prior_linear_combination_density(
      list(b0 = intercept, b1 = slopes$normal_cauchy), c(b0 = 1, b1 = -.5)
    ),
    posterior_atoms = posterior_atom_attribute()
  )
  out <- hypothesis_BF(posterior, hypothesis = "theta > 0.5",
                       parameter = "theta", columns = "all")
  prior_mass <- region_references$normal_cauchy[["theta > 0.5"]]
  expect_lt(abs(out[["prior"]] - prior_mass / (1 - prior_mass)), 1e-7)
})

test_that("conditional-normal scale-mixture region probabilities match independent references", {

  # pure scale mixture N(0, 1) + N(0, 1) * gamma(3, 2): the refined grid
  # reported convergence at P(X > .5) = 0.381793, 2.7e-4 from the reference
  a <- prior("normal", list(0, 1))
  b <- prior("normal", list(0, 1))
  attr(b, "multiply_by") <- "s"
  density <- .prior_linear_combination_density(
    list(a = a, b = b, s = prior("gamma", list(3, 2))), c(a = 1, b = 1)
  )
  expect_region_probabilities(density, c(
    "theta > 0" = .5,
    "theta > 0.5" = .38151770904855646381,
    "theta < -0.2" = .45181835797143489515,
    "theta > -0.1 & theta < 0.1" = .048296755160832632211,
    "theta > 4" = .026126722864007111829
  ))

  # location peak of a narrow multiplied normal, N(0, .01) + N(3, .01) *
  # gamma(3, 2): regions around, straddling and beyond the peak s* = v / 3
  b <- prior("normal", list(3, .01))
  attr(b, "multiply_by") <- "s"
  density <- .prior_linear_combination_density(
    list(a = prior("normal", list(0, .01)), b = b, s = prior("gamma", list(3, 2))),
    c(a = 1, b = 1)
  )
  expect_region_probabilities(density, c(
    "theta > 1.5" = .91969144965196478572,
    "theta > 1.4 & theta < 1.6" = .024507286143457283953,
    "theta > 1.49 & theta < 1.51" = .0024525250666812009689,
    "theta < 0.3" = .0011518249956491194327,
    "theta > 10" = .038040843289549597714
  ))
})

test_that("mixture region probabilities sum their components and exact point masses", {

  intercept <- prior("normal", list(0, 1))
  slope <- prior("cauchy", list(0, .5))
  slope_ss <- prior_spike_and_slab(slope, prior_inclusion = prior("spike", list(.5)))
  intercept_ss <- prior_spike_and_slab(intercept, prior_inclusion = prior("spike", list(.3)))
  reference <- region_references$normal_cauchy
  normal_tail <- stats::pnorm(-.5)
  cauchy_tail <- stats::pcauchy(.5, scale = .25, lower.tail = FALSE)

  # spike-and-slab slope: the spike leaves the normal intercept
  density <- .prior_linear_combination_density(
    list(b0 = intercept, b1 = slope_ss), c(b0 = 1, b1 = -.5)
  )
  expect_region_probabilities(density, c(
    "theta > 0.5" = .5 * normal_tail + .5 * reference[["theta > 0.5"]],
    "theta > 2 & theta < 4" = .5 * (stats::pnorm(4) - stats::pnorm(2)) +
      .5 * reference[["theta > 2 & theta < 4"]]
  ))

  # spike-and-slab intercept: the spike leaves the scaled Cauchy slope
  density <- .prior_linear_combination_density(
    list(b0 = intercept_ss, b1 = slope), c(b0 = 1, b1 = -.5)
  )
  expect_region_probabilities(density, c(
    "theta > 0.5" = .7 * cauchy_tail + .3 * reference[["theta > 0.5"]]
  ))

  # both: the spike-spike component is a point mass .35 at zero, included
  # only by the inclusive relations
  density <- .prior_linear_combination_density(
    list(b0 = intercept_ss, b1 = slope_ss), c(b0 = 1, b1 = -.5)
  )
  expect_region_probabilities(density, c(
    "theta > 0" = .325,
    "theta >= 0" = .675,
    "!(theta < 0)" = .675,
    "theta > 0.5" = .35 * cauchy_tail + .15 * normal_tail +
      .15 * reference[["theta > 0.5"]]
  ))

  # model mixture: each model is a Gaussian convolution with its own
  # quadrature (the mixture grid cannot be refined within the grid limit)
  context <- .prior_density_build_context(
    list(a = list(prior("normal", list(0, .01), prior_weights = 1),
                  prior("t", list(0, 1, 3), prior_weights = 1)),
         b = list(prior("gamma", list(3, 2), prior_weights = 1),
                  prior("normal", list(0, 1), prior_weights = 1))),
    c("a", "b")
  )
  expect_region_probabilities(
    .prior_density_from_context(context, c(a = 1, b = 1)),
    c("theta < 0.3" = .2993480255493350221)
  )

  # distinct design rows are weighted by their counts
  context <- .prior_density_context(
    prior_list = list(b0 = intercept, b1 = slope), column_names = c("b0", "b1")
  )
  rows <- .prior_density_from_context_rows(
    context,
    rbind(c(b0 = 1, b1 = -.5), c(b0 = 1, b1 = 0), c(b0 = 1, b1 = -.5))
  )
  expect_region_probabilities(rows, c(
    "theta > 0.5" = (2 * reference[["theta > 0.5"]] + normal_tail) / 3
  ))
})

test_that("region probabilities outside a bounded structural support are exact", {

  # Two half-normal terms (a two-term convolution) and a gamma term multiplied
  # by a half-normal scale (a scale mixture) are supported on [0, Inf): a
  # region outside that support has probability 0 and a region containing it
  # probability 1, without a quadrature that could only integrate zeros.
  # References by hand-split integrals.
  half_normal <- prior("normal", list(0, 1), list(0, Inf))
  convolution <- .prior_linear_combination_density(
    list(a = half_normal, b = half_normal), c(a = 1, b = 1)
  )
  product_priors <- list(beta = prior("gamma", list(2, 2)), sigma = half_normal)
  attr(product_priors$beta, "multiply_by") <- "sigma"
  product <- .prior_linear_combination_density(product_priors, c(beta = 1))
  for(density in list(convolution, product)){
    expect_identical(as.numeric(region_probability(density, "theta < -0.5")), 0)
    expect_identical(as.numeric(region_probability(density, "theta > -0.1")), 1)
    expect_identical(as.numeric(region_probability(density, "theta < -1 | theta > -0.5")), 1)
  }
  expect_equal(
    as.numeric(region_probability(convolution, "theta < 1")),
    stats::integrate(function(a) 2 * stats::dnorm(a) * (2 * stats::pnorm(1 - a) - 1),
                     0, 1, rel.tol = 1e-12)$value,
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(region_probability(product, "theta < 1")),
    stats::integrate(function(s) 2 * stats::dnorm(s) * stats::pgamma(1 / s, 2, 2), 0, 1, rel.tol = 1e-12)$value +
      stats::integrate(function(s) 2 * stats::dnorm(s) * stats::pgamma(1 / s, 2, 2), 1, Inf, rel.tol = 1e-12)$value,
    tolerance = 1e-10
  )

  # in a mixture, components supported on [0, Inf) contribute exactly 0 to a
  # region below it: (N | T) + (0 | T) with T = N(0.5, 1)T(0, Inf)
  truncated <- prior("normal", list(.5, 1), list(0, Inf))
  mixture <- .prior_linear_combination_density(
    list(a = prior_mixture(list(prior("normal", list(0, 1)), truncated), is_null = c(FALSE, FALSE)),
         b = prior_mixture(list(prior("point", list(0)), truncated), is_null = c(TRUE, FALSE))),
    c(a = 1, b = 1)
  )
  normal_plus_truncated <- stats::integrate(
    function(t) stats::dnorm(t, .5) / stats::pnorm(.5) * stats::pnorm(-.5 - t), 0, Inf, rel.tol = 1e-12
  )$value
  expect_equal(
    as.numeric(region_probability(mixture, "theta < -0.5")),
    .25 * stats::pnorm(-.5) + .25 * normal_plus_truncated,
    tolerance = 1e-10
  )
})

test_that("log-identity and transformed region probabilities are exact", {

  # A log-intercept scale formula without scaled predictors reports
  # tau = exp(log(tau)): the region is the exact distribution function of
  # the declared invgamma(1, .15) prior, P(tau < x) = exp(-.15 / x).
  tau <- prior("invgamma", list(1, .15))
  density <- .prior_linear_combination_density(
    list(tau = tau), c(tau = 1),
    source_transforms = c(tau = "log"), output_transformation = "exp"
  )
  expect_region_probabilities(density, c(
    "theta > 0.1" = 1 - exp(-1.5),
    "theta < 0.5" = exp(-.3),
    "theta > 0.05 & theta < 0.3" = exp(-.5) - exp(-3)
  ), tolerance = 1e-14)

  # scaled log intercept with a lognormal prior: exp of a normal sum
  density <- .prior_linear_combination_density(
    list(tau = prior("lognormal", list(0, 1)), tau_x = prior("normal", list(0, .5))),
    c(tau = 1, tau_x = -.5),
    source_transforms = c(tau = "log", tau_x = NA), output_transformation = "exp"
  )
  expect_equal(region_probability(density, "theta > 0.5"),
               stats::pnorm(log(.5), sd = sqrt(1 + .25^2), lower.tail = FALSE),
               tolerance = 1e-14)

  # linear output transformations map the region (a negative slope swaps it)
  density <- .prior_linear_combination_density(
    list(b0 = prior("normal", list(0, 1)), b1 = prior("cauchy", list(0, .5))),
    c(b0 = 1, b1 = -.5),
    output_transformation = "lin", output_transformation_arguments = list(a = 1, b = -2)
  )
  expect_lt(abs(region_probability(density, "theta < 0") -
                  region_references$normal_cauchy[["theta > 0.5"]]), 1e-8)
})

test_that("region probabilities keep the grid only without a structural representation", {

  # a t30 plus gamma sum is a two-term convolution: its region probability is
  # the quadrature of the gamma density times the t30 distribution function
  region <- .hypothesis_prior_region(quote(theta < 1), "theta")
  convolution <- .prior_linear_combination_density(
    list(x = prior("t", list(0, 1, 30)), y = prior("gamma", list(3, 2))),
    c(x = 1, y = 1)
  )
  reference <- stats::integrate(function(y) stats::dgamma(y, 3, 2) * stats::pt(1 - y, 30),
                                0, 1, rel.tol = 1e-12)$value +
    stats::integrate(function(y) stats::dgamma(y, 3, 2) * stats::pt(1 - y, 30),
                     1, Inf, rel.tol = 1e-12)$value
  probability <- .prior_linear_density_region_probability(convolution, region)
  expect_identical(attr(probability, "numerical_diagnostics")$quadratures, 1L)
  expect_equal(as.numeric(probability), reference, tolerance = 1e-10)

  # adding a t5 term leaves no structural ordinate: the grid evaluation is
  # used unchanged
  skewed <- .prior_linear_combination_density(
    list(x = prior("t", list(0, 1, 30)), y = prior("gamma", list(3, 2)), z = prior("t", list(0, 1, 5))),
    c(x = 1, y = 1, z = 1)
  )
  expect_null(.prior_linear_density_region_probability(skewed, region))
  side <- hypothesis_parse("theta < 1")$statements[[1L]]$left
  expect_identical(
    .hypothesis_prior_density_prob(skewed, side, "theta"),
    .hypothesis_prior_density_grid_prob(
      skewed, .hypothesis_simple_parameter_comparison(side, "theta"),
      .hypothesis_side_expression(side), "theta"
    )
  )

  # the route follows the ordinate classification: structurally classified
  # combinations have an exact or quadrature region probability, and those
  # with an unknown ordinate keep the grid
  multiplied <- prior("normal", list(0, 1))
  attr(multiplied, "multiply_by") <- "s"
  densities <- list(
    .prior_linear_combination_density(
      list(a = prior("normal", list(0, 1)), b = prior("cauchy", list(0, .5))), c(a = 1, b = -.5)),
    .prior_linear_combination_density(
      list(a = prior("normal", list(0, 1)), b = prior("normal", list(1, 2))), c(a = 1, b = 2)),
    .prior_linear_combination_density(
      list(a = prior("normal", list(0, 1)), b = multiplied, s = prior("gamma", list(3, 2))),
      c(a = 1, b = 1)),
    skewed,
    convolution,
    .prior_linear_combination_density(
      list(a = prior("gamma", list(2, 1)), b = prior("gamma", list(3, 1))), c(a = 1, b = -1)),
    .prior_linear_combination_density(
      list(tau = prior("invgamma", list(1, .15)), tau_x = prior("normal", list(0, .5))),
      c(tau = 1, tau_x = -.5), source_transforms = c(tau = "log", tau_x = NA))
  )
  for(i in seq_along(densities)){
    ordinate <- prior_density_ordinate(densities[[i]], .3)
    probability <- .prior_linear_density_region_probability(densities[[i]], region)
    expect_identical(is.null(probability), identical(ordinate$behavior, "unknown"),
                     label = paste("density", i))
  }

  # regions that are not unions of intervals keep the grid
  expect_null(.hypothesis_prior_region(quote(abs(theta) < .1), "theta"))
  expect_null(.hypothesis_prior_region(quote(theta^2 < .25), "theta"))
  expect_null(.hypothesis_prior_region(quote(exp(theta) > 1.1), "theta"))
})

test_that("linear region relations become unions of intervals", {

  intervals <- function(condition){
    unname(.hypothesis_prior_region(str2lang(condition), "theta")$intervals)
  }
  expect_identical(intervals("theta > 0.5"), matrix(c(.5, Inf), 1L))
  expect_identical(intervals("0.5 > theta"), matrix(c(-Inf, .5), 1L))
  expect_identical(intervals("2 * theta < 1"), matrix(c(-Inf, .5), 1L))
  expect_identical(intervals("1 - theta <= 0"), matrix(c(1, Inf), 1L))
  expect_identical(intervals("!(theta <= 0.2)"), matrix(c(.2, Inf), 1L))
  expect_identical(intervals("theta > -0.1 & theta < 0.1"), matrix(c(-.1, .1), 1L))
  expect_identical(intervals("theta < -1 | theta > 1"), matrix(c(-Inf, 1, -1, Inf), 2L))
  expect_identical(intervals("theta < 0 | theta > 0"), matrix(c(-Inf, Inf), 1L))
  expect_identical(nrow(intervals("theta > 1 & theta < 0")), 0L)
  expect_identical(nrow(intervals("!(theta < 1 | theta > 0)")), 0L)
  expect_identical(intervals("(theta > -2 & theta < -1) | (theta > 3 & !(theta > 4))"),
                   matrix(c(-2, 3, -1, 4), 2L))
  expect_null(.hypothesis_prior_region(quote(theta > phi), "theta"))
})

test_that("region quadratures split at every region endpoint with the full budget per piece", {

  spec <- .prior_density_ordinate_gaussian_convolution_spec(
    list(b0 = prior("normal", list(0, 1)), b1 = prior("cauchy", list(0, .5))),
    c(b0 = 1, b1 = -.5), NULL
  )
  # each endpoint v has its peak at t* = -2 v, with local SD 2
  points <- .prior_conditional_normal_breakpoints(spec, c(-.1, .1))
  for(centre in c(.2, -.2)){
    for(peak_point in centre + c(0, -1, 1, -3, 3, -10, 10) * 2){
      expect_true(any(abs(points - peak_point) < 1e-12), label = format(peak_point))
    }
  }

  integral <- .prior_conditional_normal_region(spec, matrix(c(-.1, .1), 1L), 256L)
  integration <- integral$integration
  expect_identical(integration$budget, 256L)
  expect_true(integration$converged)
  expect_true(all(integration$piece_evaluations <= 256L))
  expect_identical(integration$evaluations, sum(integration$piece_evaluations))
  expect_gt(integration$evaluations, 256L)
  expect_lt(abs(integral$value -
                  region_references$normal_cauchy[["theta > -0.1 & theta < 0.1"]]), 1e-8)
})

# A conditional-normal scale mixture integrates every piece by quadrature
# (no piece is evaluated exactly), so mocked pieces determine the total.
scale_mixture_priors <- function(){
  multiplied <- prior("normal", list(0, 1))
  attr(multiplied, "multiply_by") <- "s"
  list(a = prior("normal", list(0, 1)), b = multiplied, s = prior("gamma", list(3, 2)))
}

test_that("Gaussian-convolution regions evaluate pieces outside the endpoint windows exactly", {

  # References: 40-digit mpmath, integrating over the Gaussian variable with
  # the other term's exact CDF (script
  # .work/tmp/pr58-decisions/I/references_r2.py of the workspace). The heavy
  # Cauchy tails (scale 1 and 5) previously stopped as "roundoff error" or
  # "probably divergent"; bounds at 1e4 missed up to 2e-6.
  normal <- prior("normal", list(0, 1))
  references <- list(
    "1" = c("theta > 0.5" = .364359843655328366816,
            "theta < -0.2" = .4444604638053888513949,
            "theta > 3" = .06148835450822882103696,
            "theta < -10 | theta > 10" = .03213115939557752564382),
    "5" = c("theta > 0.5" = .4440656746305845575946,
            "theta < -0.2" = .4774757698971268778466,
            "theta > 3" = .2315589935142253463258,
            "theta < -10 | theta > 10" = .1574045832897135562589)
  )
  for(scale in names(references)){
    density <- .prior_linear_combination_density(
      list(a = normal, b = prior("cauchy", list(0, as.numeric(scale)))), c(a = 1, b = -.5)
    )
    expect_region_probabilities(density, references[[scale]], tolerance = 1e-10)
  }

  t30 <- .prior_linear_combination_density(
    list(a = normal, b = prior("t", list(0, .5, 30))), c(a = 1, b = 1)
  )
  expect_region_probabilities(t30, c(
    "theta < 10000" = 1,
    "theta > -10000 & theta < 10000" = 1
  ), tolerance = 1e-10)
  expect_lt(abs(region_probability(t30, "theta > 10000") / 9.652759592340303478938e-109 - 1), 1e-6)
  t3 <- .prior_linear_combination_density(
    list(a = normal, b = prior("t", list(0, .5, 3))), c(a = 1, b = 1)
  )
  expect_region_probabilities(t3, c("theta < 10000" = .9999999999998621677691),
                              tolerance = 1e-10)
  expect_lt(abs(region_probability(t3, "theta > 10000") / 1.378322308848918731454e-13 - 1), 1e-8)
  gamma <- .prior_linear_combination_density(
    list(a = normal, b = prior("gamma", list(3, 2))), c(a = 1, b = 1)
  )
  expect_region_probabilities(gamma, c(
    "theta < 10000" = 1,
    "theta > -10000 & theta < 10000" = 1
  ), tolerance = 1e-10)
  expect_lt(abs(region_probability(gamma, "theta > 20") / 2.156584228136996354923e-14 - 1), 1e-8)
  # the probability beyond 1e4 underflows: an exactly zero total stops
  expect_error(
    region_probability(gamma, "theta > 10000"),
    "integration reported 'zero probability for a structurally positive region'",
    fixed = TRUE
  )

  # only the pieces inside the windows t* +- 10 of the bound .5 (t* = -1,
  # local SD 2) use quadrature; the others carry the 2 * Phi(-10) bound
  spec <- .prior_density_ordinate_gaussian_convolution_spec(
    list(a = normal, b = prior("cauchy", list(0, 5))), c(a = 1, b = -.5), NULL
  )
  integration <- .prior_conditional_normal_region(spec, matrix(c(.5, Inf), 1L), 4096L)$integration
  pieces <- cbind(utils::head(integration$breakpoints, -1L), integration$breakpoints[-1L])
  outside <- pieces[, 2L] <= -21 | pieces[, 1L] >= 19
  expect_identical(integration$exact_pieces, outside)
  expect_true(any(outside))
  expect_true(all(integration$piece_evaluations[outside] == 0L))
  masses <- stats::pcauchy(pieces[outside, 2L], scale = 5) - stats::pcauchy(pieces[outside, 1L], scale = 5)
  expect_equal(integration$piece_absolute_errors[outside], 2 * stats::pnorm(-10) * masses,
               tolerance = 1e-10)

  # scale mixtures integrate every piece
  expect_null(.prior_conditional_normal_region(
    .prior_conditional_normal_spec(
      scale_mixture_priors(),
      .prior_linear_split_multiply_groups(scale_mixture_priors(), c(a = 1, b = 1)), NULL
    ),
    matrix(c(.5, Inf), 1L), 4096L
  )$integration$exact_pieces)
})

test_that("rejected region quadratures stop instead of using the grid", {

  priors <- scale_mixture_priors()
  density <- .prior_linear_combination_density(priors, c(a = 1, b = 1))
  spec <- .prior_conditional_normal_spec(
    priors, .prior_linear_split_multiply_groups(priors, c(a = 1, b = 1)), NULL
  )
  n_pieces <- length(.prior_conditional_normal_breakpoints(spec, .5)) - 1L
  testthat::local_mocked_bindings(
    .prior_conditional_normal_piece = function(integrand, lower, upper, n_grid,
                                               relative, absolute){
      list(value = .1, abs.error = 1e-3, message = "roundoff error was detected",
           evaluations = 21L)
    },
    .package = "BayesTools"
  )
  side <- hypothesis_parse("theta > 0.5")$statements[[1L]]$left
  expect_error(
    .hypothesis_prior_density_prob(density, side, "theta"),
    paste0(
      "Conditional-normal prior probability was rejected by diagnostics: ",
      "integration reported 'roundoff error was detected' with absolute error ",
      format(n_pieces * 1e-3), ". Inspect the prior specification and the region bounds."
    ),
    fixed = TRUE
  )
})

test_that("quadrature region probabilities may exceed one only within their absolute error", {

  priors <- scale_mixture_priors()
  density <- .prior_linear_combination_density(priors, c(a = 1, b = 1))
  spec <- .prior_conditional_normal_spec(
    priors, .prior_linear_split_multiply_groups(priors, c(a = 1, b = 1)), NULL
  )
  n_pieces <- length(.prior_conditional_normal_breakpoints(spec, .5)) - 1L
  total <- 1 + 1e-12
  total_error <- 1e-11
  testthat::local_mocked_bindings(
    .prior_conditional_normal_piece = function(integrand, lower, upper, n_grid,
                                               relative, absolute){
      list(value = total / n_pieces, abs.error = total_error / n_pieces,
           message = "OK", evaluations = 21L)
    },
    .package = "BayesTools"
  )
  expect_identical(region_probability(density, "theta > 0.5"), 1)
  total <- 1 + 1e-6
  total_error <- 1e-9
  expect_error(
    region_probability(density, "theta > 0.5"),
    "Computed prior probability lies materially outside [0, 1].",
    fixed = TRUE
  )
})

test_that("point terms at negative locations do not warn in region probabilities", {

  # The structural route reads point-term offsets; a point term at a negative
  # location without a log source transformation must not evaluate its log
  # (which warned "NaNs produced" in every region test and ordinate).
  normal <- .prior_linear_combination_density(
    list(b = prior("normal", list(0, 1)), c = prior("point", list(-1))), c(b = 2, c = 1)
  )
  side <- hypothesis_parse("theta > 0")$statements[[1L]]$left
  expect_no_warning(probability <- .hypothesis_prior_density_prob(normal, side, "theta"))
  expect_equal(probability, stats::pnorm(.5, lower.tail = FALSE), tolerance = 1e-14)
  expect_no_warning(prior_density_ordinate(normal, .5))

  # Gaussian convolution shifted by the point term: P(b0 - .5 b1 - 1 > -.5)
  shifted <- .prior_linear_combination_density(
    list(b0 = prior("normal", list(0, 1)), b1 = prior("cauchy", list(0, .5)),
         c = prior("point", list(-1))),
    c(b0 = 1, b1 = -.5, c = 1)
  )
  side <- hypothesis_parse("theta > -0.5")$statements[[1L]]$left
  expect_no_warning(probability <- .hypothesis_prior_density_prob(shifted, side, "theta"))
  expect_lt(abs(probability - region_references$normal_cauchy[["theta > 0.5"]]), 1e-8)
})
