skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Harrell-Davis quantiles
# ============================================================================ #
#
# PURPOSE:
#   Tests for harrell_davis_quantile() in R/harrell-davis-quantile.R: exact
#   references, values pinned from Hmisc::hdquantile and mpmath, per-cell
#   accuracy of the weights, invariances, edge cases, shapes, and the
#   pre-registered accuracy study against an analytic band.
#
# DEPENDENCIES:
#   - No external packages required beyond testthat
#   - common-functions.R: hd_reference()
#
# SKIP CONDITIONS:
#   - The accuracy study with 10,000 draws is skipped on CRAN (skip_on_cran());
#     it runs whenever NOT_CRAN is true, which devtools::test() sets, so in
#     every profile of tools/test-profile.R (unit included). All other tests
#     are fast pure R tests.
#
# MODELS/FIXTURES:
#   - None required
#
# TAGS: @quantiles, @fast
# ============================================================================ #

test_that("harrell_davis_quantile matches exact rational references", {

  # m = 3, p = .5: Beta(2, 2), cell weights (7, 13, 7) / 27
  expect_equal(harrell_davis_quantile(c(0, 1, 10), .5), 83 / 27, tolerance = 1e-14)
  expect_equal(harrell_davis_quantile(c(1, 0, 0), .5), 7 / 27, tolerance = 1e-14)
  expect_equal(harrell_davis_quantile(c(0, 0, 1), .5), 7 / 27, tolerance = 1e-14)
  expect_equal(harrell_davis_quantile(c(1, 1, 0), .5), 20 / 27, tolerance = 1e-14)

  # m = 3, p = .25: Beta(1, 3), 1 - (1 - u)^3, cell weights (19, 7, 1) / 27;
  # p = .75: Beta(3, 1), u^3, cell weights (1, 7, 19) / 27 (a swap of the shape
  # parameters or a rank offset of the estimator changes these values)
  expect_equal(harrell_davis_quantile(c(0, 1, 10), .25), 17 / 27,  tolerance = 1e-14)
  expect_equal(harrell_davis_quantile(c(0, 1, 10), .75), 197 / 27, tolerance = 1e-14)
  expect_equal(
    harrell_davis_quantile(c(10, 0, 1), c(.75, .25, .5)),
    c(197, 17, 83) / 27,
    tolerance = 1e-14
  )

  # m = 2, p = .5: the mean of the two draws; p = .25 uses Beta(.75, 2.25)
  expect_equal(harrell_davis_quantile(c(3, 11), .5), 7, tolerance = 1e-14)
  w <- stats::pbeta(.5, .75, 2.25)
  expect_equal(harrell_davis_quantile(c(3, 11), .25), w * 3 + (1 - w) * 11, tolerance = 1e-14)
})

test_that("harrell_davis_quantile reproduces values pinned from Hmisc::hdquantile", {

  # Pinned from Hmisc 5.2.5 (Hmisc::hdquantile(x, p), one probability per call)
  # with the script .work/tmp/hd-bands/R1/make_hmisc_pins.R of the
  # BayesToolsVerse workspace; the draws are regenerated from the seeds:
  #   set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion"); x <- draw(m)
  # Hmisc is not a dependency of BayesTools.
  pins <- list(
    list(
      m = 1000, seed = 20261002L, draw = function(m) stats::rnorm(m),
      probs = c(.025, .5, .975, .140171, .30403, .67149, .605112, .800821),
      hd    = c(-1.953651007918050908, 0.039860609257444564, 1.944994697147360618,
                -1.044588490658620694, -0.491918146231295261, 0.479272378889210204,
                0.306852666356472126, 0.858442619671539120)
    ),
    list(
      m = 10000, seed = 20261003L, draw = function(m) stats::rexp(m),
      probs = c(.025, .5, .975, .389435, .978607, .926719, .889386, .242853),
      hd    = c(0.026413709306294737, 0.695523477091871678, 3.735276096088531617,
                0.488897614916946999, 3.905859207895304941, 2.551991889652609480,
                2.143688371667423098, 0.275931989174885484)
    )
  )

  for (pin in pins) {
    set.seed(pin$seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
    x <- pin$draw(pin$m)

    expect_equal(harrell_davis_quantile(x, pin$probs), pin$hd, tolerance = 1e-12)
    # one probability at a time and in a matrix with a second column
    expect_equal(
      vapply(pin$probs, function(p) harrell_davis_quantile(x, p), numeric(1)),
      pin$hd, tolerance = 1e-12
    )
    expect_equal(
      harrell_davis_quantile(cbind(x, -x), pin$probs)[, 1],
      pin$hd, tolerance = 1e-12
    )
  }
})

test_that("harrell_davis_quantile weights are accurate in every cell, also in the tails", {

  # Per-cell relative error of the weights against an independent computation on
  # the log scale (pbeta(log.p = TRUE), cells up to the mode by differences of
  # the lower tail and the cells beyond it by differences of the upper tail),
  # over the cells whose reference weight exceeds 1e-250. Measured maximum over
  # m = 100, 500, 2,000 and p = .001, .025, .5, .975, .999: 1.3e-13 (both
  # evaluate pbeta at the same double cell boundaries, so this is the rounding
  # of the two computations, not of the boundaries; against 80-digit mpmath
  # weights at the exact boundaries i / m both differ by up to 3e-13), so the
  # tolerance of 1e-11 keeps a margin of about 80. The naive
  # difference of distribution function values of the cells (as in
  # Hmisc::hdquantile) loses the small weights of the tail it computes by
  # differences of values near 1 and fails this bound (relative error 1 to 2 for
  # p = .025 and .5, 5e-8 for p = .975, at m = 500).
  m <- 500
  grid <- (0:m) / m
  floor_weight <- 1e-250
  tolerance <- 1e-11
  max_cell_error <- function(weights, reference) {
    keep <- is.finite(reference) & reference > floor_weight
    max(abs(weights[keep] - reference[keep]) / reference[keep])
  }
  for (p in c(.025, .5, .975)) {
    a <- p * (m + 1)
    b <- (1 - p) * (m + 1)
    weights <- BayesTools:::.harrell_davis_weights(p, m)

    log_lower <- stats::pbeta(grid, a, b, log.p = TRUE)
    log_upper <- stats::pbeta(grid, a, b, lower.tail = FALSE, log.p = TRUE)
    lower_ref <- exp(log_lower[-1]) * -expm1(log_lower[-(m + 1)] - log_lower[-1])
    upper_ref <- exp(log_upper[-(m + 1)]) * -expm1(log_upper[-1] - log_upper[-(m + 1)])
    reference <- ifelse(seq_len(m) <= which.max(weights), lower_ref, upper_ref)

    expect_lt(max_cell_error(weights, reference), tolerance)
    expect_equal(sum(weights), 1, tolerance = 1e-14)
    expect_true(all(weights >= 0))

    # the bound separates the tail-wise weights from the naive differences
    naive <- diff(stats::pbeta(grid, a, b))
    expect_gt(max_cell_error(naive, reference), tolerance)
  }
})

test_that("harrell_davis_quantile is accurate in the far tails of 100,000 draws", {

  # Binary draws y_(i) = 1(i > r) give HD(p) = sum_{i > r} w_i = 1 - I_{r/m}(a, b),
  # the Beta(a, b) mass above r / m (a = p (m + 1), b = (1 - p) (m + 1)); r draws
  # of one at the top, y_(i) = 1(i > m - r), at p' = 1 - p give the mass above
  # (m - r) / m of Beta(a', b'). Pins computed at 80 digits with mpmath (continued
  # fraction, cross-checked against mpmath.betainc and quadrature of the density)
  # by .work/tmp/hd-bands/R1b/mp_pins.py of the BayesToolsVerse workspace, for the
  # exact doubles p and 1 - p. Tolerances from the measured relative errors of the
  # double-precision results (.work/tmp/hd-bands/R1b/error_analysis.R): lower
  # tail at most 3.0e-15, tolerance 1e-13; upper tail at most 6.0e-13, which is
  # the rounding of the cell boundary (m - r) / m (up to 1.1e-16) times the Beta
  # density near it (about 6e3), tolerance 2e-11. Both keep a margin of about 30.
  m <- 100000
  pins <- list(
    list(r = 50,   p = .0005,      hd = 0.48120586272734142,   upper = FALSE),
    list(r = 1,    p = 1e-8,       hd = 2.1960687681330079e-4, upper = FALSE),
    list(r = 5000, p = .05,        hd = 0.49826346553611729,   upper = FALSE),
    list(r = 50,   p = 1 - .0005,  hd = 0.51879413727296992,   upper = TRUE),
    list(r = 1,    p = 1 - 1e-8,   hd = 0.99978039312208210,   upper = TRUE),
    list(r = 5000, p = 1 - .05,    hd = 0.50173653446385861,   upper = TRUE)
  )
  for (pin in pins) {
    draws <- if (pin$upper) as.numeric(seq_len(m) > m - pin$r) else as.numeric(seq_len(m) > pin$r)
    expect_equal(
      harrell_davis_quantile(draws, pin$p),
      pin$hd,
      tolerance = if (pin$upper) 2e-11 else 1e-13
    )
  }

  # shuffled draws in a matrix column give the same value
  set.seed(5)
  draws <- sample(as.numeric(seq_len(m) > 50))
  expect_equal(
    harrell_davis_quantile(unname(cbind(draws, 1 - draws)), .0005)[1, 1],
    0.48120586272734142,
    tolerance = 1e-13
  )
})

test_that("harrell_davis_quantile is equivariant and invariant as the estimator is", {

  set.seed(1)
  x     <- stats::rexp(300)
  probs <- c(.025, .1, .5, .9, .975)
  q     <- harrell_davis_quantile(x, probs)

  # affine (a > 0), reflection of the probabilities, and permutation
  expect_equal(harrell_davis_quantile(3 * x + 2, probs), 3 * q + 2, tolerance = 1e-12)
  expect_equal(harrell_davis_quantile(-3 * x + 2, 1 - probs), -3 * q + 2, tolerance = 1e-12)
  expect_equal(harrell_davis_quantile(-x, 1 - probs), -q, tolerance = 1e-12)
  expect_identical(harrell_davis_quantile(rev(x), probs), q)
  expect_identical(harrell_davis_quantile(sample(x), probs), q)

  # order and duplicates of the probabilities
  expect_identical(
    harrell_davis_quantile(x, c(.975, .5, .975, .025)),
    q[c(5, 3, 5, 1)]
  )

  # monotone in the probability
  grid <- seq(.001, .999, length.out = 300)
  expect_true(all(diff(harrell_davis_quantile(x, grid)) > 0))

  # constant draws are exact, also for scales that a rescaled form would round
  for (value in c(0, 1, 0.1, -3.3, 1e-300, 7e307, -1.5e200)) {
    expect_identical(
      harrell_davis_quantile(rep(value, 11), c(0, .025, .5, .975, 1)),
      rep(value, 5)
    )
  }
})

test_that("harrell_davis_quantile handles single draws, two draws, ties, and the end probabilities", {

  expect_identical(harrell_davis_quantile(5, c(0, .3, 1)), c(5, 5, 5))
  expect_identical(harrell_davis_quantile(matrix(c(5, 7), nrow = 1), c(0, .3, 1)), matrix(c(5, 5, 5, 7, 7, 7), nrow = 3))
  expect_identical(harrell_davis_quantile(Inf, .5), Inf)

  expect_identical(harrell_davis_quantile(c(4, 2), c(0, 1)), c(2, 4))
  expect_identical(harrell_davis_quantile(c(4, 2, 9, 7), c(0, 1)), c(2, 9))

  # ties: the estimator is a function of the sorted draws only
  x <- c(0, 0, 0, 1, 1, 1, 1, 5, 5, 9)
  expect_equal(
    harrell_davis_quantile(x, c(.1, .5, .9)),
    vapply(c(.1, .5, .9), function(p) hd_reference(x, p), numeric(1)),
    tolerance = 1e-14
  )
  # a single deviating draw in constant draws moves the estimate by its weight only
  y <- c(rep(2, 40), 3)
  expect_equal(harrell_davis_quantile(y, .5), 2 + BayesTools:::.harrell_davis_weights(.5, 41)[41], tolerance = 1e-14)
})

test_that("harrell_davis_quantile does not overflow for finite draws", {

  # the differences of the draws exceed the largest double
  expect_equal(harrell_davis_quantile(c(-1e308, 1e308), .5) / 1e308, 0, tolerance = 1e-14)
  w <- stats::pbeta(.5, .75, 2.25)
  expect_equal(
    harrell_davis_quantile(c(-1e308, 1e308), .25) / 1e308,
    -w + (1 - w),
    tolerance = 1e-14
  )
  big <- harrell_davis_quantile(c(-1.7e308, 0, 1.7e308), c(.025, .5, .975))
  expect_true(all(is.finite(big)))
  expect_equal(big[1], -big[3], tolerance = 1e-14)
  expect_lt(abs(big[2]) / 1.7e308, 1e-14)
  expect_true(big[1] > -1.7e308 && big[3] < 1.7e308)
})

test_that("harrell_davis_quantile falls back to the empirical quantile for infinite draws", {

  x <- c(1, 2, 3, 4, Inf)
  probs <- c(0, .025, .5, .975, 1)
  expect_identical(
    harrell_davis_quantile(x, probs),
    stats::quantile(x, probs, names = FALSE, type = 7)
  )
  x <- c(-Inf, -3, 0, 2, 5)
  expect_identical(
    harrell_davis_quantile(x, probs),
    stats::quantile(x, probs, names = FALSE, type = 7)
  )

  # only the column with the infinite draw falls back
  set.seed(1)
  finite <- stats::rnorm(50)
  x_matrix <- cbind(a = c(finite, Inf), b = c(finite, 0))
  q <- harrell_davis_quantile(x_matrix, c(.025, .5, .975))
  expect_identical(q[, "a"], stats::quantile(x_matrix[, "a"], c(.025, .5, .975), names = FALSE))
  expect_equal(q[, "b"], vapply(c(.025, .5, .975), function(p) hd_reference(x_matrix[, "b"], p), numeric(1)), tolerance = 1e-14)

  # infinite draws of both signs make the empirical interpolation NaN
  nan_matrix <- cbind(ok = c(1, 2), bad = c(-Inf, Inf))
  expect_error(
    harrell_davis_quantile(nan_matrix, .5),
    "column 2 of 'x'",
    class = "BayesTools_harrell_davis_undefined"
  )
  expect_error(harrell_davis_quantile(c(-Inf, Inf), .5), class = "BayesTools_harrell_davis")
  # ... unless the probabilities do not interpolate between them
  expect_identical(harrell_davis_quantile(c(-Inf, Inf), c(0, 1)), c(-Inf, Inf))

  # probabilities so close to 0 that the weights cannot be computed stop with the
  # classed error (the limit of the distribution function routine depends on the
  # number of draws and the R version: R 4.6 stops at 1e-320 with 1,000 draws but
  # not with 3 draws); a result is never NaN
  for (m in c(3, 1000)) {
    tiny <- tryCatch(
      harrell_davis_quantile(seq_len(m), 1e-320),
      BayesTools_harrell_davis_undefined = function(e) conditionMessage(e)
    )
    expect_true(
      (is.character(tiny) && grepl("Use values of 'probs' further from 0 and 1.", tiny, fixed = TRUE)) ||
        (is.numeric(tiny) && tiny >= 1 && tiny <= m)
    )
  }
  expect_equal(harrell_davis_quantile(c(1, 2, 3), 1e-12), 1, tolerance = 1e-9)
})

test_that("harrell_davis_quantile builds the weights only for columns with finite draws", {

  # a column of infinite draws is summarized by the empirical quantile whatever
  # the probability, also where the weights cannot be computed
  expect_identical(harrell_davis_quantile(rep(Inf, 1000), 1e-315), Inf)
  expect_identical(harrell_davis_quantile(rep(-Inf, 1000), 1 - 1e-12), -Inf)

  weights_built <- 0L
  testthat::local_mocked_bindings(
    .harrell_davis_weights = function(p, m) {
      weights_built <<- weights_built + 1L
      stop("the weights were built")
    },
    .package = "BayesTools"
  )

  infinite <- cbind(a = c(1, 2, Inf, 4), b = c(-Inf, 0, 1, 2))
  probs    <- c(.025, .5, .975)
  expect_identical(
    harrell_davis_quantile(infinite, probs),
    apply(infinite, 2, stats::quantile, probs = probs, names = FALSE, type = 7)
  )
  expect_identical(
    harrell_davis_quantile(infinite, c(0, 1)),
    matrix(c(1, Inf, -Inf, 2), nrow = 2, dimnames = list(NULL, c("a", "b")))
  )
  expect_identical(weights_built, 0L)

  # probabilities of 0 and 1 need no weights either
  expect_identical(harrell_davis_quantile(c(3, 1, 2), c(0, 1)), c(1, 3))
  expect_identical(weights_built, 0L)

  # a finite column builds them (once), after the infinite column was summarized
  expect_error(harrell_davis_quantile(cbind(infinite[, 1], c(1, 2, 3, 4)), probs), "the weights were built", fixed = TRUE)
  expect_identical(weights_built, 1L)
})

test_that("harrell_davis_quantile returns vectors and matrices with the documented shapes", {

  set.seed(2)
  x <- matrix(stats::rnorm(60), nrow = 20, dimnames = list(NULL, c("a", "b", "c")))
  q <- harrell_davis_quantile(x, c(.5, .025))
  expect_identical(dim(q), c(2L, 3L))
  expect_identical(colnames(q), c("a", "b", "c"))
  expect_null(rownames(q))
  expect_equal(q[, "b"], harrell_davis_quantile(x[, "b"], c(.5, .025)), tolerance = 1e-15)

  named <- harrell_davis_quantile(x, c(.5, .025, .975), names = TRUE)
  expect_identical(rownames(named), names(stats::quantile(1:2, c(.5, .025, .975))))
  expect_identical(
    names(harrell_davis_quantile(x[, 1], c(.5, .025), names = TRUE)),
    c("50%", "2.5%")
  )
  expect_null(names(harrell_davis_quantile(x[, 1], c(.5, .025))))

  # one column, no column, and no probabilities
  one <- harrell_davis_quantile(x[, 1, drop = FALSE], c(.5, .025))
  expect_identical(dim(one), c(2L, 1L))
  expect_identical(colnames(one), "a")
  expect_identical(dim(harrell_davis_quantile(matrix(numeric(), nrow = 5, ncol = 0), c(.5, .1))), c(2L, 0L))
  expect_identical(harrell_davis_quantile(x[, 1], numeric()), numeric())
  expect_identical(harrell_davis_quantile(x[, 1], NULL), numeric())
  expect_identical(dim(harrell_davis_quantile(x, numeric())), c(0L, 3L))

  # integer draws, classed draws, and names of the draws are not carried over
  expect_equal(harrell_davis_quantile(1:9, .5), 5, tolerance = 1e-14)
  expect_equal(harrell_davis_quantile(c(a = 1, b = 3, c = 2), .5), 2, tolerance = 1e-14)
  expect_null(names(harrell_davis_quantile(c(a = 1, b = 3, c = 2), .5)))
  expect_equal(harrell_davis_quantile(structure(c(3, 1, 2), class = "foo"), .5), 2, tolerance = 1e-14)
  expect_equal(harrell_davis_quantile(matrix(1:9, 3), .5), matrix(c(2, 5, 8), nrow = 1), tolerance = 1e-14)
})

test_that("harrell_davis_quantile validates its input", {

  expect_error(harrell_davis_quantile(numeric(), .5), "at least one draw", fixed = TRUE)
  expect_error(harrell_davis_quantile(matrix(numeric(), nrow = 0, ncol = 2), .5), "at least one draw", fixed = TRUE)
  expect_error(harrell_davis_quantile(c(1, NA, 3), .5), "The 'x' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(harrell_davis_quantile(c(1, NaN, 3), .5), "The 'x' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(harrell_davis_quantile(c("a", "b"), .5), "numeric vector or a numeric matrix", fixed = TRUE)
  expect_error(harrell_davis_quantile(factor(c(3, 1, 2)), .5), "numeric vector or a numeric matrix", fixed = TRUE)
  expect_error(harrell_davis_quantile(data.frame(a = 1:3), .5), "numeric vector or a numeric matrix", fixed = TRUE)
  expect_error(harrell_davis_quantile(array(1:8, c(2, 2, 2)), .5), "numeric vector or a numeric matrix", fixed = TRUE)

  expect_error(harrell_davis_quantile(1:3, -.1), "The 'probs' must be equal or higher than 0.", fixed = TRUE)
  expect_error(harrell_davis_quantile(1:3, 1.1), "The 'probs' must be equal or lower than 1.", fixed = TRUE)
  expect_error(harrell_davis_quantile(1:3, Inf), "The 'probs' must be equal or lower than 1.", fixed = TRUE)
  expect_error(harrell_davis_quantile(1:3, -Inf), "The 'probs' must be equal or higher than 0.", fixed = TRUE)
  expect_error(harrell_davis_quantile(1:3, NA_real_), "The 'probs' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(harrell_davis_quantile(1:3, c(.5, NaN)), "The 'probs' argument cannot contain NA/NaN values.", fixed = TRUE)
  expect_error(harrell_davis_quantile(1:3, "a"), "The 'probs' argument must be a numeric vector.", fixed = TRUE)
  expect_error(harrell_davis_quantile(1:3, .5, names = NA), "The 'names' argument cannot contain NA/NaN values.", fixed = TRUE)
})

# Pre-registered accuracy study (plan hd-bands): pointwise bands of lines
# y = b0 + b1 * x with b0 ~ N(0, 1) and b1 ~ N(0, .5^2), whose pointwise
# quantile is analytic, at p = .025 over a uniform grid of 41 points. The
# roughness of the error is the RMS of its second differences. Medians of the
# paired ratios over 200 seeded replicates: roughness of the empirical
# quantile over that of the HD estimate (at least 4 for 1,000 and 2 for
# 10,000 draws) and RMSE of the HD estimate over that of the empirical
# quantile (at most 1.02).
hd_accuracy_study <- function(m, minimum_roughness_ratio) {

  roughness <- function(error) sqrt(mean(diff(error, differences = 2)^2))
  rmse      <- function(error) sqrt(mean(error^2))
  grid      <- seq(-1, 1, length.out = 41)
  p         <- .025
  truth     <- stats::qnorm(p, 0, sqrt(1 + .25 * grid^2))

  set.seed(31, kind = "Mersenne-Twister", normal.kind = "Inversion")
  ratios <- vapply(seq_len(200), function(i) {
    b0  <- stats::rnorm(m)
    b1  <- stats::rnorm(m, 0, .5)
    y   <- outer(b0, rep(1, length(grid))) + outer(b1, grid)
    emp <- apply(y, 2, stats::quantile, probs = p, names = FALSE)
    hd  <- harrell_davis_quantile(y, p)[1, ]
    c(roughness = roughness(emp - truth) / roughness(hd - truth),
      rmse      = rmse(hd - truth) / rmse(emp - truth))
  }, numeric(2))

  expect_gte(stats::median(ratios["roughness", ]), minimum_roughness_ratio)
  expect_lte(stats::median(ratios["rmse", ]), 1.02)
}

test_that("harrell_davis_quantile attenuates the kinks of the empirical quantile without losing accuracy (1,000 draws)", {

  hd_accuracy_study(1000, 4)
})

test_that("harrell_davis_quantile attenuates the kinks of the empirical quantile without losing accuracy (10,000 draws)", {

  skip_on_cran()
  hd_accuracy_study(10000, 2)
})
