skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Prior densities of the levels of ordered terms in interactions
# ============================================================================ #
#
# A coefficient of an ordered term is theta_c = T_s * prod_k G_{k,s}[i_k(c)]:
# the total T_s of its slice s (the coordinates of the ordinary factors of the
# term, whose treatment or independent coding makes every coefficient touch one
# slice) times the cumulative allocation of every ordered factor k. The slices
# have independent totals and independent allocations unless an `id` shares an
# allocation. These tests check the selection of the slice of a level, the
# refusals of the combinations without a density route, and the independence
# of the terms of a sum.
#
# ============================================================================ #

N01 <- prior("normal", list(0, 1))
treatment <- prior_factor("normal", list(0, 1), contrast = "treatment")

.interaction_data <- function(){

  data.frame(
    f = ordered(rep(c("lo", "mid", "hi", "top"), 3), levels = c("lo", "mid", "hi", "top")),
    o = ordered(rep(c("p", "q", "r"), each = 4), levels = c("p", "q", "r")),
    g = factor(rep(c("A", "B", "C"), each = 4)),
    x = sin(seq_len(12))
  )
}

# The bound prior list of a formula.
.interaction_priors <- function(formula, priors){

  JAGS_formula(formula, "mu", .interaction_data(), priors)$prior_list
}

# The weights of a combination of the cells of a term: 'cells' are the names of
# the cells ("f=hi, g=B") and 'signs' their signs.
.cell_combination <- function(prior, cells, signs = rep(1, length(cells))){

  design <- BayesTools:::.factor_term_design_from_metadata(prior)
  rows <- design$design
  rownames(rows) <- design$cell_names
  weights <- colSums(rows[cells, , drop = FALSE] * signs)
  stats::setNames(weights, BayesTools:::.JAGS_prior_factor_names(
    BayesTools:::.prior_ordered_metadata(prior)$parameter_name, prior))
}

.interaction_share <- function(prior, cells, signs = rep(1, length(cells))){

  weights <- .cell_combination(prior, cells, signs)
  BayesTools:::.prior_ordered_linear_share(prior, weights, seq_along(weights))
}

# The structural route of a combination of the cells of a term, with extra
# terms of the combination when 'extra' (a prior list and its weights) is given.
.interaction_route <- function(prior, cells, signs = rep(1, length(cells)), extra = NULL){

  weights <- .cell_combination(prior, cells, signs)
  parameter <- BayesTools:::.prior_ordered_metadata(prior)$parameter_name
  prior_list <- c(stats::setNames(list(prior), parameter), extra$prior_list)
  BayesTools:::.prior_density_route_linear(prior_list, c(weights, extra$weights), NULL, 1024L)
}

.route_height <- function(route, value){

  ordinate <- BayesTools:::.prior_density_route_ordinate(route, value)
  expect_true(ordinate$exact)
  exp(ordinate$log_density)
}

.route_probability <- function(route, lower, upper = Inf){

  region <- list(intervals = BayesTools:::.prior_region_intervals(lower, upper),
                 indicator = function(x) x > lower & x < upper)
  result <- BayesTools:::.prior_density_route_region(route, region)
  expect_true(result$converged)
  result$probability
}

# Nested integrals over independent Beta variables, the references of the
# product of shares (stats::integrate(), independent of the tensor rules).
.nested_integral <- function(integrand, shapes, rel_tol = 1e-11){

  inner <- function(index, shares){
    if(index > length(shapes)) return(integrand(shares))
    stats::integrate(function(b) vapply(b, function(value){
      inner(index + 1L, c(shares, value))
    }, numeric(1)) * stats::dbeta(b, shapes[[index]][1L], shapes[[index]][2L]),
    0, 1, rel.tol = rel_tol, subdivisions = 500L)$value
  }
  inner(1L, numeric())
}

test_that("a level of an ordered interaction takes the allocation share of its slice", {

  allocation <- prior("dirichlet", list(c(2, 1, 1)))
  sliced <- .interaction_priors(~ f * g, list(intercept = N01, f = prior_ordered(N01),
    g = treatment, "f:g" = prior_ordered(N01, allocation = allocation)))[["mu_f__xXx__g"]]
  single <- JAGS_formula(~ f, "mu", .interaction_data(),
    list(intercept = N01, f = prior_ordered(N01, allocation = allocation)))$prior_list[["mu_f"]]

  for(slice in c("B", "C")) for(level in c("mid", "hi", "top")){
    share <- .interaction_share(sliced, paste0("f=", level, ", g=", slice))
    reference <- .interaction_share(single, level)
    expect_identical(share[c("type", "scale", "alpha")], reference[c("type", "scale", "alpha")],
                     info = paste(level, slice))
  }
  # the cumulative shares are Beta variables of the summed concentrations
  expect_identical(.interaction_share(sliced, "f=mid, g=B")$alpha, c(2, 2))
  expect_identical(.interaction_share(sliced, "f=hi, g=B")$alpha, c(3, 1))
  # the slices are different allocations, and the share records the slice used
  expect_identical(.interaction_share(sliced, "f=mid, g=B")$slice, 1L)
  expect_identical(.interaction_share(sliced, "f=mid, g=C")$slice, 2L)
  expect_false(identical(.interaction_share(sliced, "f=mid, g=B")$record$key,
                         .interaction_share(sliced, "f=mid, g=C")$record$key))
  # the top level is the total (a point share of 1) and a differences within a
  # slice is a subset of the increments
  expect_identical(.interaction_share(sliced, "f=top, g=C")[c("type", "scale")], list(type = "point", scale = 1))
  difference <- .interaction_share(sliced, c("f=top, g=B", "f=mid, g=B"), c(1, -1))
  expect_identical(difference[c("type", "scale", "alpha")], list(type = "beta", scale = 1, alpha = c(2, 2)))
  # a combination without any weight is the point 0
  none <- BayesTools:::.prior_ordered_linear_share(sliced, numeric(), integer())
  expect_identical(none, list(type = "point", scale = 0))
})

test_that("combinations without one allocation share are refused with the formula-measure condition", {

  refused <- function(share){
    condition <- tryCatch(share, error = function(e) e)
    expect_s3_class(condition, "BayesTools_formula_measure_unavailable")
    expect_identical(condition$reason, "unsupported_contribution_measure")
    condition
  }
  sliced <- .interaction_priors(~ f * g, list(intercept = N01, f = prior_ordered(N01),
    g = treatment, "f:g" = prior_ordered(N01)))[["mu_f__xXx__g"]]
  crossed <- .interaction_priors(~ f * o, list(intercept = N01, f = prior_ordered(N01),
    o = prior_ordered(N01), "f:o" = prior_ordered(N01)))[["mu_f__xXx__o"]]
  single <- .interaction_priors(~ f, list(intercept = N01, f = prior_ordered(N01)))[["mu_f"]]

  # coefficients of two slices: independent totals and allocations
  condition <- refused(.interaction_share(sliced, c("f=mid, g=B", "f=mid, g=C"), c(1, -1)))
  expect_match(conditionMessage(condition), "spans 2 slices", fixed = TRUE)
  # two ordered factors: the product of two shares
  condition <- refused(.interaction_share(crossed, "f=hi, o=q"))
  expect_match(conditionMessage(condition), "allocations of 2 ordered factors", fixed = TRUE)
  # a sum of levels that is not one allocation subset
  condition <- refused(.interaction_share(single, c("mid", "hi")))
  expect_match(conditionMessage(condition), "no beta reduction", fixed = TRUE)
  # one factor at its top level leaves one share
  expect_identical(.interaction_share(crossed, "f=top, o=q")[c("type", "alpha")],
                   list(type = "beta", alpha = c(1, 1)))
  expect_identical(.interaction_share(crossed, "f=hi, o=r")[c("type", "alpha")],
                   list(type = "beta", alpha = c(2, 1)))
})

test_that("fixed allocations never add increment weights across slices", {

  fixed <- c(.2, .3, .5)
  sliced <- .interaction_priors(~ f * g, list(intercept = N01, f = prior_ordered(N01),
    g = treatment, "f:g" = prior_ordered(N01, allocation = fixed)))[["mu_f__xXx__g"]]

  # within a slice the fixed weights contract: hi - mid = .5 - .2
  within <- .interaction_share(sliced, c("f=hi, g=B", "f=mid, g=B"), c(1, -1))
  expect_identical(within$type, "point")
  expect_equal(within$scale, .3, tolerance = 1e-14)
  expect_equal(.interaction_share(sliced, "f=mid, g=C")$scale, .2, tolerance = 1e-14)
  # across slices the same weights belong to different totals: +1 and -1 do not cancel
  condition <- tryCatch(.interaction_share(sliced, c("f=mid, g=B", "f=mid, g=C"), c(1, -1)),
    error = function(e) e)
  expect_s3_class(condition, "BayesTools_formula_measure_unavailable")
  # it is the difference of two independent fixed multiples of the totals,
  # N(0, .2 sqrt(2)), not the point 0 of cancelled weights
  difference <- .interaction_route(sliced, c("f=mid, g=B", "f=mid, g=C"), c(1, -1))
  expect_identical(difference$type, "normal")
  expect_equal(.route_height(difference, .1), stats::dnorm(.1, 0, .2 * sqrt(2)), tolerance = 1e-12)
  expect_equal(.route_probability(difference, .05), stats::pnorm(.05, 0, .2 * sqrt(2), lower.tail = FALSE),
               tolerance = 1e-12)
})

test_that("ordered terms that share a random allocation are not summed as independent terms", {

  priors <- function(id_x){
    .interaction_priors(~ f + x + f:x, list(intercept = N01, f = prior_ordered(N01, id = "shape"),
      x = N01, "f:x" = prior_ordered(N01, id = id_x)))
  }
  density_of <- function(prior_list, weights){
    BayesTools:::.prior_linear_combination_density(prior_list, weights)
  }
  # the middle level is the first increment, the top level all three increments
  middle <- function(term) stats::setNames(c(1, 0, 0), paste0(term, "[", 1:3, "]"))
  top <- function(term) stats::setNames(c(1, 1, 1), paste0(term, "[", 1:3, "]"))

  shared <- priors("shape")
  # the middle levels of two terms with one allocation: dependent
  both_middle <- c(middle("mu_f"), middle("mu_f__xXx__x"))
  expect_error(density_of(shared, both_middle), class = "BayesTools_formula_measure_unavailable")
  condition <- tryCatch(density_of(shared, both_middle), error = function(e) e)
  expect_match(conditionMessage(condition), "id:shape|factor:f", fixed = TRUE)
  # the route reports the same reason instead of an independent product
  route <- BayesTools:::.prior_density_route_linear(shared, both_middle, NULL)
  expect_identical(route$type, "unknown")
  expect_match(route$reason, "share the random allocation", fixed = TRUE)
  # the top levels are the totals (independent), whatever the allocation
  total <- density_of(shared, c(top("mu_f"), top("mu_f__xXx__x")))
  expect_s3_class(total, "prior_linear_density")
  expect_equal(prior_density_ordinate(total, 0)$log_density, stats::dnorm(0, 0, sqrt(2), log = TRUE),
               tolerance = 1e-12)
  # one level of one term and the total of the other are independent
  expect_s3_class(density_of(shared, c(middle("mu_f"), top("mu_f__xXx__x"))), "prior_linear_density")
  # different allocations are independent
  expect_s3_class(density_of(priors("other"), both_middle), "prior_linear_density")
})

test_that("the mixed-measure density of an interaction level has the atom of the spike and the exact curve", {

  allocation <- prior("dirichlet", list(c(2, 1, 1)))
  slab <- prior_spike_and_slab(N01, prior("point", list(.5)))
  priors <- .interaction_priors(~ f * g, list(intercept = N01, f = prior_ordered(slab),
    g = treatment, "f:g" = prior_ordered(slab, allocation = allocation)))
  displayed <- density(priors[["mu_f__xXx__g"]], x_range = c(-3, 3), n_points = 129L)
  expect_identical(attr(displayed, "method", exact = TRUE), "analytic_mixed_measure")

  # g = C, the middle level: 0.5 N(0, 1) x Beta(2, 2) beside the atom at 0 of mass 0.5
  level <- displayed[["mu_f__xXx__g[f=mid, g=C]"]]
  expect_identical(level$atoms, data.frame(location = 0, mass = .5))
  at <- level$continuous$x[c(20L, 45L, 80L, 110L)]
  reference <- vapply(at, function(x){
    .5 * stats::integrate(function(b) stats::dnorm(x / b) / b * stats::dbeta(b, 2, 2), 0, 1,
                          rel.tol = 1e-10)$value
  }, numeric(1))
  expect_equal(level$continuous$density[match(at, level$continuous$x)], reference, tolerance = 1e-6)
  # the top level of a slice is the spike-and-slab total
  expect_identical(displayed[["mu_f__xXx__g[f=top, g=B]"]]$atoms, data.frame(location = 0, mass = .5))
  # a reference cell of the treatment factor is the point 0
  expect_identical(displayed[["mu_f__xXx__g[f=top, g=A]"]]$atoms, data.frame(location = 0, mass = 1))
})

test_that("the difference of two slices is a product with independent shares and totals", {

  # T1 B1 - T2 B2, T ~ N(0, 1), B ~ Beta(2, 3): conditionally N(0, B1^2 + B2^2)
  allocation <- prior("dirichlet", list(c(2, 3)))
  sliced <- .interaction_priors(~ o * g, list(intercept = N01, o = prior_ordered(N01),
    g = treatment, "o:g" = prior_ordered(N01, allocation = allocation)))[["mu_o__xXx__g"]]
  route <- .interaction_route(sliced, c("o=q, g=B", "o=q, g=C"), c(1, -1))
  expect_identical(route$type, "ordered_product")
  expect_identical(route$spec$mode, "normal")
  expect_length(route$spec$structure$coordinates, 2L)

  # references: nested integrate() over the two shares
  density_at_zero <- .nested_integral(function(b) 1 / sqrt(2 * pi * (b[1L]^2 + b[2L]^2)),
                                      list(c(2, 3), c(2, 3)))
  expect_equal(density_at_zero, 0.7796905941, tolerance = 1e-9)
  expect_equal(.route_height(route, 0), density_at_zero, tolerance = 1e-6)
  tail_above <- .nested_integral(function(b) stats::pnorm(.1 / sqrt(b[1L]^2 + b[2L]^2), lower.tail = FALSE),
                                 list(c(2, 3), c(2, 3)))
  expect_equal(tail_above, 0.4231782899, tolerance = 1e-9)
  expect_equal(.route_probability(route, .1), tail_above, tolerance = 1e-6)
  expect_equal(.route_height(route, .7),
               .nested_integral(function(b) stats::dnorm(.7 / sqrt(b[1L]^2 + b[2L]^2)) / sqrt(b[1L]^2 + b[2L]^2),
                                list(c(2, 3), c(2, 3))), tolerance = 1e-6)
  # i.i.d. slices: the difference is symmetric
  expect_equal(.route_probability(route, 0), .5, tolerance = 1e-12)
  expect_equal(.route_probability(route, -Inf, -.1), tail_above, tolerance = 1e-6)
  # the diagnostics of the tensor rule
  integration <- BayesTools:::.prior_density_route_ordinate(route, .3)$provenance$integration
  expect_true(integration$converged)
  expect_gte(length(integration$rule_sizes), 2L)
})

test_that("a level of two ordered factors is the product of the total and two shares", {

  # T Bf Bo, T ~ N(0, 1), Bf ~ Beta(2, 2) (f mid), Bo ~ Beta(1, 2) (o q)
  crossed <- .interaction_priors(~ f * o, list(intercept = N01, f = prior_ordered(N01),
    o = prior_ordered(N01), "f:o" = prior_ordered(N01, allocation = list(
      f = prior("dirichlet", list(c(2, 1, 1))), o = prior("dirichlet", list(c(1, 2))))))
  )[["mu_f__xXx__o"]]
  route <- .interaction_route(crossed, "f=mid, o=q")
  expect_identical(route$type, "ordered_product")
  shapes <- list(c(2, 2), c(1, 2))
  for(value in c(-1.2, .05, .3, 1)){
    reference <- .nested_integral(function(b) stats::dnorm(value / prod(b)) / prod(b), shapes)
    expect_equal(.route_height(route, value), reference, tolerance = 1e-5, info = value)
  }
  expect_equal(.route_probability(route, .05),
               .nested_integral(function(b) stats::pnorm(.05 / prod(b), lower.tail = FALSE), shapes),
               tolerance = 1e-5)
  # the density at the offset: f_T(0) E[1 / (Bf Bo)] is infinite where a share has a = 1
  at_zero <- BayesTools:::.prior_density_route_ordinate(route, 0)
  expect_identical(at_zero$behavior, "infinite")
  # shares with a > 1 give the closed form f_T(0) E[1 / Bf] E[1 / Bo]
  smooth <- .interaction_priors(~ f * o, list(intercept = N01, f = prior_ordered(N01),
    o = prior_ordered(N01), "f:o" = prior_ordered(N01, allocation = list(
      f = prior("dirichlet", list(c(3, 1, 1))), o = prior("dirichlet", list(c(2.5, 2)))))))[["mu_f__xXx__o"]]
  expected <- stats::dnorm(0) * (3 + 2 - 1) / (3 - 1) * (2.5 + 2 - 1) / (2.5 - 1)
  expect_equal(.route_height(.interaction_route(smooth, "f=mid, o=q"), 0), expected, tolerance = 1e-12)
  # the top level of one factor leaves one share: the density of a single share product
  one <- .interaction_route(crossed, "f=top, o=q")
  expect_false(identical(one$type, "ordered_product"))
  expect_equal(.route_height(one, .3),
               .nested_integral(function(b) stats::dnorm(.3 / b[1L]) / b[1L], list(c(1, 2))), tolerance = 1e-6)
  # a leaf with a single coordinate agrees with the route of the share
  leaf <- BayesTools:::.prior_density_route_ordered_product(
    list(mu_f__xXx__o = crossed), .cell_combination(crossed, "f=top, o=q"), NULL, 1024L)
  expect_identical(leaf$type, "ordered_product")
  expect_equal(.route_height(leaf, .3), .route_height(one, .3), tolerance = 1e-6)
  expect_equal(.route_probability(leaf, .3), .route_probability(one, .3), tolerance = 1e-6)
})

test_that("combinations of increments that are not one cumulative share use the stick-breaking coordinates", {

  # f[mid] + f[hi] = T (2 G1 + G2), (G1, G2, G3) ~ Dirichlet(2, 1, 1.5)
  allocation <- prior("dirichlet", list(c(2, 1, 1.5)))
  dirichlet_density <- function(g1, g2){
    exp(lgamma(4.5) - lgamma(2) - lgamma(1) - lgamma(1.5)) * g1^(2 - 1) * g2^(1 - 1) * (1 - g1 - g2)^(1.5 - 1)
  }
  simplex_integral <- function(integrand){
    stats::integrate(function(g1) vapply(g1, function(a){
      stats::integrate(function(g2) integrand(a, g2) * dirichlet_density(a, g2), 0, 1 - a,
                       rel.tol = 1e-10, subdivisions = 400L)$value
    }, numeric(1)), 0, 1, rel.tol = 1e-10, subdivisions = 400L)$value
  }
  for(total in list(N01, prior("normal", list(0, 1), list(0, Inf)), prior("t", list(0, 1, 3)))){
    single <- .interaction_priors(~ f, list(intercept = N01,
      f = prior_ordered(total, allocation = allocation)))[["mu_f"]]
    route <- .interaction_route(single, c("mid", "hi"))
    expect_identical(route$type, "ordered_product")
    expect_length(route$spec$structure$coordinates, 2L)
    for(value in c(.2, .9, 1.7)){
      reference <- simplex_integral(function(g1, g2){
        scale <- 2 * g1 + g2
        BayesTools::mpdf(total, value / scale) / scale
      })
      expect_equal(.route_height(route, value), reference, tolerance = 1e-5,
                   info = paste(total$distribution, value))
    }
    reference <- simplex_integral(function(g1, g2) BayesTools::mccdf(total, .9 / (2 * g1 + g2)))
    expect_equal(.route_probability(route, .9), reference, tolerance = 1e-5, info = total$distribution)
  }
})

test_that("a shared allocation is one set of coordinates for all slices", {

  allocation <- prior("dirichlet", list(c(2, 1, 1.5)))
  shared <- .interaction_priors(~ f * g, list(intercept = N01, f = prior_ordered(N01), g = treatment,
    "f:g" = prior_ordered(N01, allocation = allocation, id = "a")))[["mu_f__xXx__g"]]

  # the same level in two slices: B (T1 - T2), B ~ Beta(2, 2.5), T1 - T2 ~ N(0, sqrt(2))
  route <- .interaction_route(shared, c("f=mid, g=B", "f=mid, g=C"), c(1, -1))
  expect_length(route$spec$structure$coordinates, 1L)
  for(value in c(.3, 1)){
    reference <- stats::integrate(function(b) stats::dnorm(value / (b * sqrt(2))) / (b * sqrt(2)) *
      stats::dbeta(b, 2, 2.5), 0, 1, rel.tol = 1e-12)$value
    expect_equal(.route_height(route, value), reference, tolerance = 1e-6)
  }
  # different levels of one allocation: nested subsets of the same Dirichlet
  # (G1 and G1 + G2), conditionally N(0, G1^2 + (G1 + G2)^2)
  nested <- .interaction_route(shared, c("f=mid, g=B", "f=hi, g=C"), c(1, -1))
  expect_length(nested$spec$structure$coordinates, 2L)
  dirichlet_density <- function(g1, g2){
    exp(lgamma(4.5) - lgamma(2) - lgamma(1) - lgamma(1.5)) * g1^(2 - 1) * (1 - g1 - g2)^(1.5 - 1)
  }
  for(value in c(.4, 1.1)){
    reference <- stats::integrate(function(g1) vapply(g1, function(a){
      stats::integrate(function(g2){
        scale <- sqrt(a^2 + (a + g2)^2)
        stats::dnorm(value / scale) / scale * dirichlet_density(a, g2)
      }, 0, 1 - a, rel.tol = 1e-10, subdivisions = 400L)$value
    }, numeric(1)), 0, 1, rel.tol = 1e-10, subdivisions = 400L)$value
    expect_equal(.route_height(nested, value), reference, tolerance = 1e-5)
  }
  # independent allocations are two sets of coordinates
  separate <- .interaction_priors(~ f * g, list(intercept = N01, f = prior_ordered(N01), g = treatment,
    "f:g" = prior_ordered(N01, allocation = allocation)))[["mu_f__xXx__g"]]
  expect_length(.interaction_route(separate, c("f=mid, g=B", "f=mid, g=C"), c(1, -1))$spec$structure$coordinates, 2L)
})

test_that("a spike-and-slab total has one indicator for all slices", {

  allocation <- prior("dirichlet", list(c(2, 3)))
  slab <- prior_spike_and_slab(N01, prior("point", list(.4)))
  sliced <- .interaction_priors(~ o * g, list(intercept = N01, o = prior_ordered(slab),
    g = treatment, "o:g" = prior_ordered(slab, allocation = allocation)))[["mu_o__xXx__g"]]
  route <- .interaction_route(sliced, c("o=q, g=B", "o=q, g=C"), c(1, -1))
  expect_identical(route$type, "mixture")
  # the atom at 0 has the mass of the indicator, 1 - pi, not (1 - pi)^2 of one indicator per slice
  atoms <- BayesTools:::.prior_ordered_route_points(route)
  expect_identical(atoms, data.frame(x = 0, p = .6))
  # the slab is the difference of two i.i.d. slabs
  slab_route <- .interaction_route(
    .interaction_priors(~ o * g, list(intercept = N01, o = prior_ordered(N01), g = treatment,
      "o:g" = prior_ordered(N01, allocation = allocation)))[["mu_o__xXx__g"]],
    c("o=q, g=B", "o=q, g=C"), c(1, -1))
  expect_equal(.route_probability(route, .1), .4 * .route_probability(slab_route, .1), tolerance = 1e-12)
  expect_equal(.route_probability(route, -Inf, Inf), 1)
  # the continuous part sums with the atom to one
  expect_equal(.route_probability(route, -Inf, 0) + .route_probability(route, 0, Inf), .4, tolerance = 1e-9)
  # the density of the continuous part
  expect_equal(.route_height(route, .3), .4 * .route_height(slab_route, .3), tolerance = 1e-12)
})

test_that("terms added to the ordered levels, and totals that are not normal", {

  allocation <- prior("dirichlet", list(c(2, 3)))
  sliced <- .interaction_priors(~ o * g, list(intercept = N01, o = prior_ordered(N01),
    g = treatment, "o:g" = prior_ordered(N01, allocation = allocation)))[["mu_o__xXx__g"]]
  intercept <- list(prior_list = list(mu_intercept = prior("normal", list(.5, 2))),
                    weights = c(mu_intercept = 1))
  # a normal term adds its variance to the conditional normal
  route <- .interaction_route(sliced, c("o=q, g=B", "o=q, g=C"), c(1, -1), extra = intercept)
  expect_identical(route$type, "ordered_product")
  reference <- .nested_integral(function(b) stats::dnorm(.8, .5, sqrt(b[1L]^2 + b[2L]^2 + 4)),
                                list(c(2, 3), c(2, 3)))
  expect_equal(.route_height(route, .8), reference, tolerance = 1e-6)
  # a point shifts it
  shifted <- .interaction_route(sliced, c("o=q, g=B", "o=q, g=C"), c(1, -1),
    extra = list(prior_list = list(mu_intercept = prior("point", list(1))), weights = c(mu_intercept = 1)))
  expect_equal(.route_height(shifted, 1.3), .route_height(.interaction_route(sliced,
    c("o=q, g=B", "o=q, g=C"), c(1, -1)), .3), tolerance = 1e-12)

  # several totals that are not normal, and a normal term added to a total that is not normal, have no route
  heavy <- .interaction_priors(~ o * g, list(intercept = N01, o = prior_ordered(N01), g = treatment,
    "o:g" = prior_ordered(prior("t", list(0, 1, 3)), allocation = allocation)))[["mu_o__xXx__g"]]
  refused <- .interaction_route(heavy, c("o=q, g=B", "o=q, g=C"), c(1, -1))
  expect_identical(refused$type, "unknown")
  expect_match(refused$reason, "closed-form law only when they are normal", fixed = TRUE)
  expect_identical(BayesTools:::.prior_density_route_ordinate(refused, .3)$behavior, "unknown")
  expect_error(BayesTools:::.prior_linear_combination_density(list(mu_o__xXx__g = heavy),
    .cell_combination(heavy, c("o=q, g=B", "o=q, g=C"), c(1, -1))),
    class = "BayesTools_formula_measure_unavailable")
  # one t total alone is the scale product
  single <- .interaction_route(heavy, "o=q, g=B")
  expect_false(identical(single$type, "unknown"))
  expect_equal(.route_height(single, .4),
               stats::integrate(function(b) BayesTools::mpdf(prior("t", list(0, 1, 3)), .4 / b) / b *
                                  stats::dbeta(b, 2, 3), 0, 1, rel.tol = 1e-11)$value, tolerance = 1e-6)
})

test_that("shares with a strong singularity, many coordinates, and point totals", {

  # a Beta exponent below 0.1 is a strong singularity: the ordinate is rejected, not guessed
  steep <- .interaction_priors(~ o * g, list(intercept = N01, o = prior_ordered(N01), g = treatment,
    "o:g" = prior_ordered(N01, allocation = prior("dirichlet", list(c(.05, 3))))))[["mu_o__xXx__g"]]
  route <- .interaction_route(steep, c("o=q, g=B", "o=q, g=C"), c(1, -1))
  ordinate <- BayesTools:::.prior_density_route_ordinate(route, .3)
  expect_false(ordinate$exact)
  expect_false(ordinate$provenance$integration$converged)
  expect_match(ordinate$provenance$integration$message, "strong singularity", fixed = TRUE)
  region <- BayesTools:::.prior_density_route_region(route, list(
    intervals = BayesTools:::.prior_region_intervals(.1, Inf), indicator = function(x) x > .1))
  expect_false(region$converged)
  expect_match(region$messages, "strong singularity", fixed = TRUE)
  # the prior density of the combination is unavailable (no grid is guessed)
  expect_error(BayesTools:::.prior_linear_combination_density(list(mu_o__xXx__g = steep),
    .cell_combination(steep, c("o=q, g=B", "o=q, g=C"), c(1, -1))),
    "strong singularity", class = "BayesTools_formula_measure_unavailable")

  # more allocation coordinates than the node cap admits are refused
  expect_identical(BayesTools:::.prior_ordered_product_sizes(5L), c(6L, 12L))
  expect_identical(BayesTools:::.prior_ordered_product_sizes(6L), 6L)
  data <- expand.grid(f = ordered(c("a", "b", "c", "d"), levels = c("a", "b", "c", "d")),
                      o = ordered(c("p", "q", "r", "s"), levels = c("p", "q", "r", "s")),
                      q = ordered(c("u", "v", "w", "x"), levels = c("u", "v", "w", "x")))
  three <- JAGS_formula(~ f * o * q, "mu", data, list(intercept = N01, f = prior_ordered(N01),
    o = prior_ordered(N01), q = prior_ordered(N01), "f:o" = prior_ordered(N01), "f:q" = prior_ordered(N01),
    "o:q" = prior_ordered(N01), "f:o:q" = prior_ordered(N01)))$prior_list[["mu_f__xXx__o__xXx__q"]]
  # three shares for the cell, six for a contrast of two cells
  cell <- .interaction_route(three, "f=c, o=r, q=w")
  expect_identical(cell$type, "ordered_product")
  expect_length(cell$spec$structure$coordinates, 3L)
  shapes <- list(c(2, 1), c(2, 1), c(2, 1))
  expect_equal(.route_height(cell, .2),
               .nested_integral(function(b) stats::dnorm(.2 / prod(b)) / prod(b), shapes, rel_tol = 1e-8),
               tolerance = 1e-5)
  contrast <- .interaction_route(three, c("f=b, o=q, q=v", "f=c, o=r, q=w"), c(1, -1))
  expect_identical(contrast$type, "unknown")
  expect_match(contrast$reason, "6 independent allocation coordinates", fixed = TRUE)

  # point totals: a function of the shares, affine in one share, nonlinear in several
  point <- .interaction_priors(~ f * o, list(intercept = N01, f = prior_ordered(N01), o = prior_ordered(N01),
    "f:o" = prior_ordered(prior("point", list(2)))))[["mu_f__xXx__o"]]
  affine <- .interaction_route(point, "f=top, o=q")
  expect_identical(affine$type, "scalar")
  expect_equal(.route_height(affine, 1), stats::dbeta(1 / 2, 1, 1) / 2, tolerance = 1e-12)
  nonlinear <- .interaction_route(point, "f=mid, o=q")
  expect_identical(nonlinear$type, "unknown")
  expect_match(nonlinear$reason, "nonlinear function of 2 allocation shares", fixed = TRUE)
})

test_that("a level of one slice is the total times a Beta share, as in a term without interactions", {

  # N(0, 1) x Beta(2, 3): the middle level of a three-level factor with Dirichlet(2, 3)
  allocation <- prior("dirichlet", list(c(2, 3)))
  sliced <- .interaction_priors(~ o * g, list(intercept = N01, o = prior_ordered(N01),
    g = treatment, "o:g" = prior_ordered(N01, allocation = allocation)))[["mu_o__xXx__g"]]
  single <- .interaction_priors(~ o, list(intercept = N01,
    o = prior_ordered(N01, allocation = allocation)))[["mu_o"]]
  cell <- .interaction_route(sliced, "o=q, g=B")
  level <- .interaction_route(single, "q")
  heights <- c(-1.2, -.3, .05, .4, 1)
  expected <- c(0.03717798, 0.60730055, 1.33854949, 0.44791123, 0.07066209)
  for(i in seq_along(heights)){
    reference <- stats::integrate(function(b) stats::dnorm(heights[[i]] / b) / b * stats::dbeta(b, 2, 3),
                                  0, 1, rel.tol = 1e-12)$value
    expect_equal(reference, expected[[i]], tolerance = 1e-7)
    expect_equal(.route_height(cell, heights[[i]]), reference, tolerance = 1e-6)
    # the slice changes nothing: the same ordinate as the term without the interaction
    expect_identical(.route_height(cell, heights[[i]]), .route_height(level, heights[[i]]))
  }
  expect_equal(.route_probability(cell, .4), 0.144692348, tolerance = 1e-7)
  expect_identical(.route_probability(cell, .4), .route_probability(level, .4))
  # the top level of a slice is its total
  top <- .interaction_route(sliced, "o=r, g=B")
  expect_equal(.route_height(top, .3), stats::dnorm(.3), tolerance = 1e-12)
  expect_equal(.route_probability(top, .3), stats::pnorm(.3, lower.tail = FALSE), tolerance = 1e-12)
})

# Totals with a bounded support: the conditional density f_T(x / c(B)) / |c(B)|
# jumps or has kinks in the shares where x / c(B) crosses a bound of T. The
# product T U V of a total with the cells of two uniform shares (f mid, o q)
# has closed-form references.
.bounded_cross_priors <- function(total){

  .interaction_priors(~ f * o, list(intercept = N01, f = prior_ordered(N01),
    o = prior_ordered(N01), "f:o" = prior_ordered(total, allocation = list(
      f = prior("dirichlet", list(c(1, .5, .5))), o = prior("dirichlet", list(c(1, 1)))))))
}

.bounded_cross_route <- function(total, cell = "f=mid, o=q"){

  .interaction_route(.bounded_cross_priors(total)[["mu_f__xXx__o"]], cell)
}

# The density of T U V at x < 1, T ~ U(1, 1 + width): (-log(x) log(1 + width) +
# log(1 + width)^2 / 2) / width.
.uniform_product_density <- function(x, width){

  (-log(x) * log1p(width) + log1p(width)^2 / 2) / width
}

test_that("a total with a bounded support is integrated with the breakpoints of its bounds", {

  # the case that rule doubling alone accepted at a relative error of 5e-3
  narrow <- .bounded_cross_route(prior("uniform", list(1, 1.01)))
  ordinate <- BayesTools:::.prior_density_route_ordinate(narrow, .72)
  expect_true(ordinate$exact)
  expect_equal(exp(ordinate$log_density), .uniform_product_density(.72, .01), tolerance = 1e-9)
  expect_equal(exp(ordinate$log_density), 0.331822869503, tolerance = 1e-9)
  # and through the public ordinate of the prior density of the level
  priors <- .bounded_cross_priors(prior("uniform", list(1, 1.01)))
  density <- BayesTools:::.prior_linear_combination_density(
    priors, .cell_combination(priors[["mu_f__xXx__o"]], "f=mid, o=q"), n_grid = 128L)
  public <- prior_density_ordinate(density, .72)
  expect_true(public$exact)
  expect_equal(exp(public$log_density), 0.331822869503, tolerance = 1e-9)

  # every value within the criterion of the closed form, for narrow to wide supports
  x <- seq(.01, .99, length.out = 99)
  for(width in c(1, .1, .01, .001)){
    route <- .bounded_cross_route(prior("uniform", list(1, 1 + width)))
    result <- BayesTools:::.prior_ordered_product_integrate_density(route$spec, x)
    expect_true(all(result$accepted), info = width)
    expect_lt(max(abs(result$value / .uniform_product_density(x, width) - 1)), 1e-7)
  }

  # T ~ U(1, 2) over both regimes of x (the density is log(2 / x)^2 / 2 above 1)
  wide <- .bounded_cross_route(prior("uniform", list(1, 2)))
  x <- c(seq(.05, .95, length.out = 11), seq(1.05, 1.95, length.out = 11))
  reference <- ifelse(x < 1, .uniform_product_density(x, 1), log(2 / x)^2 / 2)
  heights <- vapply(x, function(value) .route_height(wide, value), numeric(1))
  expect_lt(max(abs(heights / reference - 1)), 1e-7)

  # truncated normal totals with one or two finite bounds: a one-dimensional integral over T
  for(total in list(prior("normal", list(1, .5), list(0, 1.5)), prior("normal", list(0, 1), list(-Inf, 1)),
                    prior("normal", list(1, .3), list(.8, Inf)))){
    route <- .bounded_cross_route(total)
    x <- seq(.05, .95, length.out = 12)
    reference <- vapply(x, function(value){
      upper <- total$truncation$upper
      stats::integrate(function(t) BayesTools::mpdf(total, t) * (-log(value / t)) / t,
                       max(total$truncation$lower, value), if(is.finite(upper)) upper else 40,
                       rel.tol = 1e-12)$value
    }, numeric(1))
    heights <- vapply(x, function(value) .route_height(route, value), numeric(1))
    expect_lt(max(abs(heights / reference - 1)), 1e-7)
  }

  # regions: P(X > x) = E[S_T(x / W)] with W = U V
  survival <- function(y, width) pmin(1, pmax(0, (1 + width - y) / width))
  for(value in c(.2, .72, .95)){
    reference <- stats::integrate(function(w) -log(w) * survival(value / w, .01), value / 1.01, value,
                                  rel.tol = 1e-12)$value +
      stats::integrate(function(w) -log(w), value, 1, rel.tol = 1e-12)$value
    expect_equal(.route_probability(narrow, value), reference, tolerance = 1e-9, info = value)
    expect_equal(.route_probability(narrow, -Inf, value), 1 - reference, tolerance = 1e-9, info = value)
  }
})

test_that("a coefficient that changes sign and three coordinates are integrated with their breakpoints", {

  # T (G1 - G2), (G1, G2, G3) ~ Dirichlet(1, 1, 1), T ~ U(1, 2): G1 - G2 has the density 1 - |d|
  single <- .interaction_priors(~ f, list(intercept = N01, f = prior_ordered(prior("uniform", list(1, 2)),
    allocation = prior("dirichlet", list(c(1, 1, 1))))))[["mu_f"]]
  difference <- BayesTools:::.prior_density_route_linear(list(mu_f = single),
    c("mu_f[1]" = 1, "mu_f[2]" = -1), NULL, 1024L)
  expect_identical(difference$type, "ordered_product")
  x <- c(seq(-1.9, -.05, length.out = 10), seq(.05, 1.9, length.out = 10))
  reference <- vapply(x, function(value){
    stats::integrate(function(t) pmax(0, 1 - abs(value) / t) / t, max(1, abs(value)), 2, rel.tol = 1e-12)$value
  }, numeric(1))
  heights <- vapply(x, function(value) .route_height(difference, value), numeric(1))
  expect_lt(max(abs(heights / reference - 1)), 1e-7)
  tail <- stats::integrate(function(t) ifelse(.3 / t < 1, (1 - .3 / t)^2 / 2, 0), 1, 2, rel.tol = 1e-12)$value
  expect_equal(.route_probability(difference, .3), tail, tolerance = 1e-9)
  expect_equal(.route_probability(difference, -Inf, -.3), tail, tolerance = 1e-9)

  # the cell of three 3-level ordered factors is T times a product of three uniform shares (log W
  # is minus a Gamma(3) variable), T ~ U(1, 2)
  data <- expand.grid(f = ordered(c("a", "b", "c")), o = ordered(c("p", "q", "r")), q = ordered(c("u", "v", "w")))
  three <- JAGS_formula(~ f * o * q, "mu", data, list(intercept = N01, f = prior_ordered(N01),
    o = prior_ordered(N01), q = prior_ordered(N01), "f:o" = prior_ordered(N01), "f:q" = prior_ordered(N01),
    "o:q" = prior_ordered(N01), "f:o:q" = prior_ordered(prior("uniform", list(1, 2)))))$prior_list[[
      "mu_f__xXx__o__xXx__q"]]
  cell <- .interaction_route(three, "f=b, o=q, q=v")
  expect_length(cell$spec$structure$coordinates, 3L)
  x <- c(.02, .3, .9, 1.7)
  reference <- vapply(x, function(value){
    stats::integrate(function(t) log(t / value)^2 / 2 / t, max(1, value), 2, rel.tol = 1e-12)$value
  }, numeric(1))
  heights <- vapply(x, function(value) .route_height(cell, value), numeric(1))
  expect_lt(max(abs(heights / reference - 1)), 1e-7)
  expect_equal(.route_probability(cell, .3),
               stats::integrate(function(t) stats::pgamma(log(t / .3), 3), 1, 2, rel.tol = 1e-12)$value,
               tolerance = 1e-9)
})

test_that("a bounded total with more coordinates than the pieces admit is refused", {

  data <- expand.grid(f = ordered(c("a", "b", "c", "d")), o = ordered(c("p", "q", "r", "s")))
  crossed <- JAGS_formula(~ f * o, "mu", data, list(intercept = N01, f = prior_ordered(N01),
    o = prior_ordered(N01), "f:o" = prior_ordered(prior("uniform", list(1, 2)))))$prior_list[["mu_f__xXx__o"]]
  # four coordinates for the contrast of two cells: refused for a bounded total
  contrast <- .interaction_route(crossed, c("f=c, o=r", "f=b, o=q"), c(1, -1))
  expect_identical(contrast$type, "unknown")
  expect_match(contrast$reason, "no verified quadrature", fixed = TRUE)
  expect_identical(BayesTools:::.prior_density_route_ordinate(contrast, .3)$behavior, "unknown")
})

test_that("a combination with many allocation coordinates is refused before any grid is built", {

  # the contrast of the first and the last of 34 levels leaves 32 coordinates of one allocation
  levels <- sprintf("L%02d", 1:34)
  data <- data.frame(f = ordered(levels, levels = levels))
  priors <- JAGS_formula(~ f, "mu", data, list(intercept = N01, f = prior_ordered(N01)))$prior_list
  prior <- priors[["mu_f"]]
  weights <- stats::setNames(33:1, BayesTools:::.JAGS_prior_factor_names("mu_f", prior))
  route <- BayesTools:::.prior_density_route_linear(priors, weights, NULL, 64L)
  expect_identical(route$type, "unknown")
  expect_match(route$reason, "32 independent allocation coordinates", fixed = TRUE)
  # the range comes from the extrema of the reduced coefficients, never from a cube of 2^32 corners
  testthat::local_mocked_bindings(
    .prior_ordered_product_corners = function(structure) stop("a cube of the shares was enumerated"),
    .package = "BayesTools"
  )
  total_range <- BayesTools:::.prior_linear_group_range(
    list(prior = prior$total, weights = c(.ordered_total = 1), indices = 1L), 1e-4,
    c(.ordered_total = NA_character_))
  expect_equal(BayesTools:::.prior_ordered_product_range(prior, weights, seq_along(weights), 1e-4),
               33 * total_range, tolerance = 1e-12)
  elapsed <- system.time(condition <- tryCatch(
    BayesTools:::.prior_linear_combination_density(priors, weights, n_grid = 64L), error = identity
  ))[["elapsed"]]
  expect_s3_class(condition, "BayesTools_formula_measure_unavailable")
  expect_match(conditionMessage(condition), "32 independent allocation coordinates", fixed = TRUE)
  expect_lt(elapsed, 10)
})

test_that("a bounded total with integrable singularities at its bounds is rejected, never certified wrongly", {

  # Beta(.6, .6): the density is infinite at both bounds, so the pieces end in singularities
  # that Gauss-Legendre pieces resolve slowly
  total <- prior("beta", list(.6, .6))
  route <- .bounded_cross_route(total)
  x <- c(.05, .2, .5, .8, .95)
  reference <- vapply(x, function(value){
    stats::integrate(function(t) stats::dbeta(t, .6, .6) * (-log(value / t)) / t, value, 1,
                     rel.tol = 1e-9, subdivisions = 4000L, stop.on.error = FALSE)$value
  }, numeric(1))
  result <- BayesTools:::.prior_ordered_product_integrate_density(route$spec, x)
  # an accepted value is within the criterion of the reference; the others are rejected
  expect_true(all(abs(result$value / reference - 1)[result$accepted] < 1e-4))
  ordinate <- BayesTools:::.prior_density_route_ordinate(route, .5)
  expect_false(ordinate$exact)
  expect_false(ordinate$provenance$integration$converged)
  expect_match(ordinate$reason, "rejected by its diagnostics", fixed = TRUE)
  # the display grid of the combination is unavailable too: no height is read from a grid
  priors <- .interaction_priors(~ f * o, list(intercept = N01, f = prior_ordered(N01), o = prior_ordered(N01),
    "f:o" = prior_ordered(total, allocation = list(f = prior("dirichlet", list(c(1, .5, .5))),
                                                    o = prior("dirichlet", list(c(1, 1)))))))
  expect_error(BayesTools:::.prior_linear_combination_density(priors,
    .cell_combination(priors[["mu_f__xXx__o"]], "f=mid, o=q")),
    class = "BayesTools_formula_measure_unavailable")
})
