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
  # their design rows are refused when built as well, not deferred
  shared_context <- BayesTools:::.prior_density_build_context(shared, names(both_middle), n_grid = 256L)
  expect_error(BayesTools:::.prior_density_from_context_rows(shared_context, rbind(both_middle, both_middle)),
               "id:shape|factor:f", fixed = TRUE, class = "BayesTools_formula_measure_unavailable")
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
  # and so are design rows of that combination, not deferred
  steep_weights <- .cell_combination(steep, c("o=q, g=B", "o=q, g=C"), c(1, -1))
  steep_context <- BayesTools:::.prior_density_build_context(list(mu_o__xXx__g = steep), names(steep_weights),
                                                             n_grid = 256L)
  expect_error(BayesTools:::.prior_density_from_context_rows(steep_context, rbind(steep_weights, steep_weights)),
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

test_that("Beta totals with singular ends are integrated to the criterion", {

  # T ~ Beta(a, b) times U V: the density of X = T W, W = U V (density -log w), is
  # the integral over t of f_T(t) (-log(x / t)) / t, and the substitution z = (1 - t)^b
  # removes the singularity of the upper end for the independent reference
  density_reference <- function(x, a, b){
    stats::integrate(function(z){
      t <- 1 - z^(1 / b)
      t^(a - 2) * log(t / x) / (b * beta(a, b))
    }, 0, (1 - x)^b, rel.tol = 1e-11, abs.tol = 0, subdivisions = 4000L)$value
  }
  # P(X > x) = integral of f_T(t) P(W > x / t) with P(W > w) = 1 - w + w log w
  probability_reference <- function(x, a, b){
    stats::integrate(function(z){
      t <- 1 - z^(1 / b)
      w <- x / t
      t^(a - 1) / (b * beta(a, b)) * (1 - w + w * log(w))
    }, 0, (1 - x)^b, rel.tol = 1e-11, abs.tol = 0, subdivisions = 4000L)$value
  }

  # the case that rule doubling accepted at a relative error of 3.3e-4 (through the public ordinate)
  priors <- .bounded_cross_priors(prior("beta", list(1, .9)))
  density <- BayesTools:::.prior_linear_combination_density(
    priors, .cell_combination(priors[["mu_f__xXx__o"]], "f=mid, o=q"), n_grid = 128L)
  public <- prior_density_ordinate(density, 1e-4)
  expect_true(public$exact)
  expect_equal(exp(public$log_density), 39.5360676101493, tolerance = 1e-9)
  expect_equal(exp(public$log_density), density_reference(1e-4, 1, .9), tolerance = 1e-9)

  # singular lower ends (a < 1), upper ends (b < 1), both, and exponents between 1 and 2
  x <- c(.0001, .01, .1, .3, .72, .9, .99)
  for(shape in list(c(.6, .6), c(.9, 1), c(1, .6), c(.7, 1.3), c(1.5, .7), c(.5, .9), c(1.5, 1.5), c(2, 1.5))){
    route <- .bounded_cross_route(prior("beta", list(shape[[1L]], shape[[2L]])))
    expect_identical(route$type, "ordered_product")
    result <- BayesTools:::.prior_ordered_product_integrate_density(route$spec, x)
    expect_true(all(result$accepted), info = shape)
    reference <- vapply(x, density_reference, numeric(1), a = shape[[1L]], b = shape[[2L]])
    expect_lt(max(abs(result$value / reference - 1)), 1e-5)
    for(value in c(.01, .3, .95)){
      expect_equal(.route_probability(route, value), probability_reference(value, shape[[1L]], shape[[2L]]),
                   tolerance = 1e-8, info = c(shape, value))
    }
  }

  # a singular gamma total (shape < 1 at 0) and a truncated Beta whose bounds are regular
  gamma <- .bounded_cross_route(prior("gamma", list(.7, 1)))
  expect_lt(max(abs(vapply(c(.01, .5, 2, 5), function(value) .route_height(gamma, value), numeric(1)) /
    vapply(c(.01, .5, 2, 5), function(value){
      stats::integrate(function(t) stats::dgamma(t, .7, 1) * (-log(value / t)) / t, value, Inf,
                       rel.tol = 1e-11, abs.tol = 0)$value
    }, numeric(1)) - 1)), 1e-6)
})

test_that("totals with an unknown or very singular density at a bound are refused", {

  # exponent .3 < .5: the nodes next to the bound cannot be resolved to the criterion
  route <- .bounded_cross_route(prior("beta", list(.3, 1.5)))
  expect_identical(route$type, "unknown")
  expect_match(route$reason, "infinite with the exponent 0.3", fixed = TRUE)
  expect_error(BayesTools:::.prior_linear_combination_density(
    .bounded_cross_priors(prior("beta", list(.3, 1.5))),
    .cell_combination(.bounded_cross_priors(prior("beta", list(.3, 1.5)))[["mu_f__xXx__o"]], "f=mid, o=q")),
    class = "BayesTools_formula_measure_unavailable")
  # a family whose infinite density has no known exponent
  testthat::local_mocked_bindings(
    .prior_ordered_product_bound_exponent = function(prior, bound) NA_real_,
    .package = "BayesTools"
  )
  route <- .bounded_cross_route(prior("uniform", list(1, 2)))
  expect_identical(route$type, "unknown")
  expect_match(route$reason, "an unknown exponent", fixed = TRUE)
})

test_that("a non-monotone sequence of rules is not accepted", {

  route <- .bounded_cross_route(prior("uniform", list(1, 2)))
  sequences <- list(
    # differences 1.5e-3 then 9e-4: both small but the last is not at most half of the one before
    erratic = c(10, 10.3, 10.31, 10.3115, 10.3124, 10.3127, 10.3128),
    # differences 5e-3 then 5e-4: a convergent sequence
    convergent = c(10, 10.3, 10.31, 10.315, 10.3155, 10.31552, 10.315521)
  )
  rules_used <- function(values){
    testthat::local_mocked_bindings(
      .prior_ordered_product_nested = function(spec, values_, bounds, exponents, integrand, groups, n,
                                               nested, check = FALSE){
        index <<- index + 1L
        list(value = rep(values[[index]], length(groups)), bias = rep(0, length(groups)), evaluations = 1)
      },
      .package = "BayesTools"
    )
    index <- 0L
    result <- BayesTools:::.prior_ordered_product_piecewise(route$spec, matrix(.5, 1L, 1L),
      function(group, scale, shift) 0)
    list(accepted = result$accepted, rules = index)
  }
  erratic <- rules_used(sequences$erratic)
  convergent <- rules_used(sequences$convergent)
  # the convergent sequence is accepted at the fifth rule (difference 5e-4 after 5e-3)
  expect_true(convergent[["accepted"]])
  expect_identical(convergent[["rules"]], 5L)
  # the erratic one is not accepted at the fifth rule (9e-4 after 1.5e-3 is not halving)
  expect_gt(erratic[["rules"]], 5L)
})

test_that("values within the floating-point resolution of a singular bound are rejected, never certified wrongly", {

  # T ~ Beta(.5, .5) times U V near the bound x = 1 of the support: the density is ~ (1 - x)^1.5, the
  # quadrature nodes next to the breakpoint at distance 1 - x round onto it once 1 - x approaches the
  # floating-point resolution of the shares (the lost mass was 1.1e-3 of the density at 1 - 1e-10)
  reference <- function(x){
    delta <- 1 - x
    stats::integrate(function(s){
      2 * sqrt(delta) / pi * log1p(delta * (1 - s^2) / x) / (1 - delta * s^2)^1.5
    }, 0, 1, rel.tol = 1e-12, abs.tol = 0)$value
  }
  total <- prior("beta", list(.5, .5))
  priors <- .bounded_cross_priors(total)
  density <- BayesTools:::.prior_linear_combination_density(
    priors, .cell_combination(priors[["mu_f__xXx__o"]], "f=mid, o=q"), n_grid = 128L)

  for(delta in c(1e-2, 1e-4)){
    ordinate <- prior_density_ordinate(density, 1 - delta)
    expect_true(ordinate$exact, info = delta)
    expect_equal(exp(ordinate$log_density), reference(1 - delta), tolerance = 1e-5, info = delta)
  }
  for(delta in c(1e-6, 1e-8, 1e-10, 1e-12)){
    x <- 1 - delta
    ordinate <- prior_density_ordinate(density, x)
    # the density is bracketed analytically by 4 d^1.5 / (3 pi) and the same over x^2.5 (the
    # integral loses accuracy at 1 - 1e-12): it is either certified within the criterion or rejected
    bracket <- 4 * delta^1.5 / (3 * pi) * c(1, x^-2.5)
    if(ordinate$exact){
      value <- exp(ordinate$log_density)
      expect_true(value >= bracket[[1L]] * (1 - 1e-4) && value <= bracket[[2L]] * (1 + 1e-4), info = delta)
    }else{
      expect_false(ordinate$provenance$integration$converged)
      expect_match(ordinate$reason, "floating-point resolution", fixed = TRUE)
    }
  }
  # the unresolved mass at 1 - 1e-10 is not zeroed: it is rejected
  ordinate <- prior_density_ordinate(density, 1 - 1e-10)
  expect_false(ordinate$exact)
  expect_match(ordinate$reason, "not representable at this value", fixed = TRUE)
  # the height stops instead of returning the biased value
  expect_error(BayesTools:::.prior_linear_density_height(density, 1 - 1e-10), "rejected by diagnostics")
  # regions are not affected: the probabilities are continuous in the shares
  expect_equal(.route_probability(.bounded_cross_route(total), 1 - 1e-10),
               stats::integrate(function(t) stats::dbeta(t, .5, .5) * (1 - (1 - 1e-10) / t + (1 - 1e-10) / t *
                 log((1 - 1e-10) / t)), 1 - 1e-10, 1, rel.tol = 1e-12, subdivisions = 4000L,
                 stop.on.error = FALSE)$value, tolerance = 1e-4)
})

# The bound of the mass that the unresolved nodes of one inner integral leave out, with the standard
# arguments: the nodes 'pieces' with the density 'evaluated' of the scale density (its near-bound
# attributes included), the slopes of the shift and of the scale in the coordinate, and the cap on the
# inner integral and the numbers of unresolved outer nodes of the outer part.
.unresolved_bound <- function(spec, pieces, evaluated, scale, dscale, dshift, weight = 1, cap = Inf,
                              flagged = 0L, flagged_early = 0L){

  near <- attr(evaluated, "unresolved")
  near_bounds <- attr(evaluated, "near_bounds")
  evaluated[!is.finite(evaluated)] <- 0
  exponents <- BayesTools:::.prior_ordered_product_bound_exponents(spec$terms[[spec$random_terms[[1L]]]]$total)
  BayesTools:::.prior_ordered_product_unresolved_bias(
    spec, pieces, evaluated, pieces$unresolved | near, near_bounds, scale, rep(dscale, length(scale)),
    rep(dshift, length(scale)), as.numeric(names(exponents)), unname(exponents), weight,
    sum(pieces$w * evaluated), cap, flagged, flagged_early, 1L, 1L
  )
}

test_that("the mass left out by unresolved nodes is bounded by an envelope of the density, not by sampled values", {

  # the log envelope of the density of a total at a distance from its bound, against the density
  envelope <- function(total, bound, distance){
    exponent <- BayesTools:::.prior_ordered_product_bound_exponent(total, bound)
    exp(BayesTools:::.prior_ordered_product_density_envelope(total, rep(bound, length(distance)),
                                                              rep(exponent, length(distance)), distance))
  }
  distances <- 10^-(5:15)
  cases <- list(
    list(prior("beta", list(.5, .5)), 0, function(t) stats::dbeta(t, .5, .5)),
    list(prior("beta", list(.5, .5)), 1, function(t) stats::dbeta(t, .5, .5)),
    list(prior("beta", list(2.5, .7)), 0, function(t) stats::dbeta(t, 2.5, .7)),
    list(prior("beta", list(2.5, .7)), 1, function(t) stats::dbeta(t, .7, 2.5)),
    list(prior("beta", list(.6, 3.5)), 1, function(t) stats::dbeta(t, 3.5, .6)),
    list(prior("beta", list(.6, .8), truncation = list(lower = 0, upper = .6)), 0,
         function(t) stats::dbeta(t, .6, .8) / stats::pbeta(.6, .6, .8)),
    list(prior("gamma", list(.6, 2)), 0, function(t) stats::dgamma(t, .6, 2)),
    list(prior("gamma", list(.7, 1.5), truncation = list(lower = 0, upper = 5)), 0,
         function(t) stats::dgamma(t, .7, 1.5) / stats::pgamma(5, .7, 1.5))
  )
  for(case in cases){
    bound <- envelope(case[[1L]], case[[2L]], distances)
    # the density at the distance d from the bound (the Beta density at 1 - d is that of the
    # exchanged shapes at d)
    density <- case[[3L]](distances)
    expect_true(all(bound >= density), info = paste(case[[1L]]$distribution, case[[2L]]))
    # the envelope is the density up to the monotone factor over the distance: tight near the bound
    expect_lt(max((bound / density)[distances <= 1e-8]), 1 + 1e-6)
  }
  # no envelope without a known family, with another exponent, or at a distance that is not small
  expect_true(is.na(BayesTools:::.prior_ordered_product_density_envelope(
    prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf)), 0, 1, 1e-10)))
  expect_true(is.na(BayesTools:::.prior_ordered_product_density_envelope(
    prior("beta", list(.5, .5)), 0, .7, 1e-10)))
  expect_true(is.na(BayesTools:::.prior_ordered_product_density_envelope(
    prior("beta", list(.5, .5)), 0, .5, 1e-3)))
  expect_true(is.na(BayesTools:::.prior_ordered_product_density_envelope(
    prior("beta", list(.5, .5)), .3, .5, 1e-10)))

  # T ~ Beta(1e17, .5) times U V at x = .5: the whole density lies inside the unresolved interval,
  # where the sampled values of the resolved nodes are tiny (the estimate was 4e-10 of 0.693, the bound
  # 1e-23): an independent analytic bracket of f_X(.5) from E[log(T / x) / T; T > x]
  a <- 1e17
  x <- .5
  h <- 1e-5
  p_tail <- .5 / (a + .5) / h
  lower <- (1 - p_tail) * min(-log(x), log((1 - h) / x) / (1 - h))
  route <- .bounded_cross_route(prior("beta", list(a, .5)))
  exponents <- BayesTools:::.prior_ordered_product_bound_exponents(route$spec$terms[[1L]]$total)
  for(n in c(64L, 512L)){
    nested <- BayesTools:::.prior_ordered_product_nested(
      route$spec, matrix(x, 1L), as.numeric(names(exponents)), unname(exponents),
      function(group, scale, shift) BayesTools:::.prior_ordered_product_scale_density(route$spec, x, shift, scale),
      1L, n, TRUE, check = TRUE
    )
    expect_gte(nested$bias, lower - nested$value)
    expect_gt(nested$bias, 1e-4 * abs(nested$value))
  }
  ordinate <- BayesTools:::.prior_density_route_ordinate(route, x)
  expect_false(ordinate$exact)
  expect_match(ordinate$reason, "not representable at this value", fixed = TRUE)

  # T ~ Beta(1, .5) times U at x = .99: the contribution of the flagged nodes with the exact (logarithmic)
  # distance of each node from the breakpoint, which the sampled value of the resolved nodes underestimated
  e <- .5
  n <- 48L
  x <- .99
  spec <- .bounded_cross_route(prior("beta", list(1, e)))$spec
  pieces <- BayesTools:::.prior_ordered_product_split(matrix(c(0, x, 1), 1), matrix(c(NA, e, NA), 1),
                                                       matrix(c(NA, 1, NA), 1), 1, 1, n)
  evaluated <- BayesTools:::.prior_ordered_product_scale_density(spec, rep(x, length(pieces$v)), 0, pieces$v)
  unresolved <- pieces$unresolved | attr(evaluated, "unresolved")
  bound <- function(...) .unresolved_bound(spec, pieces, evaluated, pieces$v, 1, 0, ...)
  rule <- BayesTools:::.prior_ordered_beta_rule(n, 1, 1)
  left <- 12 * log(rule$u)
  right <- log1p(-rule$u)
  log_d <- pmax(left, right) + log1p(exp(-abs(left - right)))
  distance <- (1 - x) * exp(left - log_d)
  right_piece <- n + seq_len(n)
  flagged_mass <- sum((pieces$w[right_piece] * e * distance^(e - 1) * (x + distance)^-e)[unresolved[right_piece]])
  expect_gt(flagged_mass, 0)
  expect_gte(bound(), flagged_mass)
  # the nodes on the side of the breakpoint outside the support have the exact density 0
  expect_lt(bound(), 1.1 * flagged_mass)

  # the unresolved nodes of the last outer coordinate leave out the inner integral at their unrounded
  # value, bounded by the cap of the inner integral (or the computed integral plus the bound of its own
  # unresolved nodes, if larger); those of an earlier outer coordinate, and a missing cap, have no bound
  lost <- bound()
  inner <- sum(pieces$w * ifelse(is.finite(evaluated), evaluated, 0))
  expect_equal(bound(weight = .25, cap = 2, flagged = 1L), .25 * lost + .25 * 2 * (1 + 1e-6), tolerance = 1e-12)
  expect_equal(bound(weight = .25, cap = 2, flagged = 2L), .25 * lost + .25 * 2 * 2 * (1 + 1e-6), tolerance = 1e-12)
  expect_equal(bound(weight = .25, cap = 1e-12, flagged = 1L), .25 * lost + .25 * (inner + lost), tolerance = 1e-12)
  expect_identical(bound(weight = .25, cap = Inf, flagged = 1L), Inf)
  expect_identical(bound(weight = .25, cap = 2, flagged_early = 1L), Inf)
  # an unresolved node that has no singular edge to bound the mass at has no bound
  without_edge <- pieces
  without_edge$bound_left[] <- NA_real_
  without_edge$bound_right[] <- NA_real_
  expect_identical(.unresolved_bound(spec, without_edge, evaluated, pieces$v, 1, 0), Inf)
})

test_that("the bound counts the nodes on the side of a breakpoint inside the support, with the slope of the argument", {

  # T ~ Gamma(.6, 1) at its lower bound 0 with a shift that depends on the coordinate: the argument
  # (x - s(v)) / c(v) with s = .2 + ds v, c = 1 + .2 v is 0 at v* = .3 and positive on the side of v*
  # that makes the numerator positive; the exact density at a node is that of its exact distance
  # d = |ds| |v - v*| / c(v) from the bound (the distance in the coordinate comes from the logarithmic
  # form of the node, which the rounded node does not carry)
  shape <- .6
  spec <- .bounded_cross_route(prior("gamma", list(shape, 1)))$spec
  n <- 48L
  c_start <- 1
  dc <- .2
  for(case in list(list(ds = .5, x = .35, inside = "left"), list(ds = -.5, x = .05, inside = "right"))){
    pieces <- BayesTools:::.prior_ordered_product_split(matrix(c(0, .3, 1), 1), matrix(c(NA, shape, NA), 1),
                                                         matrix(c(NA, 0, NA), 1), 1, 1, n)
    scale <- c_start + dc * pieces$v
    shift <- .2 + case$ds * pieces$v
    evaluated <- BayesTools:::.prior_ordered_product_scale_density(spec, case$x, shift, scale)
    unresolved <- pieces$unresolved | attr(evaluated, "unresolved")
    bound <- .unresolved_bound(spec, pieces, evaluated, scale, dc, case$ds)
    # the exact flagged mass: the nodes of the piece inside the support (the left piece near its right
    # edge for ds > 0, the right piece near its left edge for ds < 0)
    piece <- if(case$inside == "left") seq_len(n) else n + seq_len(n)
    log_distance <- if(case$inside == "left") pieces$log_right[piece] else pieces$log_left[piece]
    distance <- abs(case$ds) * exp(log_distance) / scale[piece]
    exact <- stats::dgamma(distance, shape, 1) / scale[piece]
    flagged_mass <- sum((pieces$w[piece] * exact)[unresolved[piece]])
    expect_gt(flagged_mass, 0)
    expect_gte(bound, flagged_mass)
    # the computed value of a rounded node can exceed the exact density (its rounded distance is
    # smaller), and the larger of the two is bounded
    expect_lt(bound, 1.25 * flagged_mass)
  }
})

test_that("the supremum of the density near a regular bound is analytic, and unknown families have none", {

  # an independent grid maximum of the density in the interval of the radius around the bound
  radius <- 1e-3
  grid_maximum <- function(density, bound, lower, upper){
    points <- seq(max(bound - radius, lower), min(bound + radius, upper), length.out = 2001L)
    max(density(points[points > lower & points < upper]))
  }
  cases <- list(
    list(prior("uniform", list(1, 1 + 1e-9)), 1, 1, 1 + 1e-9, function(t) rep(1 / ((1 + 1e-9) - 1), length(t))),
    list(prior("beta", list(2.5, 1.5)), 1, 0, 1, function(t) stats::dbeta(t, 2.5, 1.5)),
    list(prior("beta", list(1, 1)), 0, 0, 1, function(t) stats::dbeta(t, 1, 1)),
    list(prior("beta", list(.6, .8), truncation = list(lower = .2, upper = .7)), .2, .2, .7,
         function(t) stats::dbeta(t, .6, .8) / (stats::pbeta(.7, .6, .8) - stats::pbeta(.2, .6, .8))),
    list(prior("beta", list(.6, .8), truncation = list(lower = .2, upper = .7)), .7, .2, .7,
         function(t) stats::dbeta(t, .6, .8) / (stats::pbeta(.7, .6, .8) - stats::pbeta(.2, .6, .8))),
    list(prior("gamma", list(2.5, 1.2), truncation = list(lower = .5, upper = 4)), .5, .5, 4,
         function(t) stats::dgamma(t, 2.5, 1.2) / (stats::pgamma(4, 2.5, 1.2) - stats::pgamma(.5, 2.5, 1.2))),
    list(prior("gamma", list(.7, 1.2), truncation = list(lower = .5, upper = 4)), 4, .5, 4,
         function(t) stats::dgamma(t, .7, 1.2) / (stats::pgamma(4, .7, 1.2) - stats::pgamma(.5, .7, 1.2))),
    list(prior("exp", list(2), truncation = list(lower = .1, upper = 3)), .1, .1, 3,
         function(t) stats::dexp(t, 2) / (stats::pexp(3, 2) - stats::pexp(.1, 2))),
    list(prior("normal", list(1, .5), truncation = list(lower = 0, upper = 3)), 3, 0, 3,
         function(t) stats::dnorm(t, 1, .5) / (stats::pnorm(3, 1, .5) - stats::pnorm(0, 1, .5))),
    list(prior("normal", list(1, .5), truncation = list(lower = 1 - 1e-4, upper = 3)), 1 - 1e-4, 1 - 1e-4, 3,
         function(t) stats::dnorm(t, 1, .5) / (stats::pnorm(3, 1, .5) - stats::pnorm(1 - 1e-4, 1, .5)))
  )
  for(case in cases){
    bound <- exp(BayesTools:::.prior_ordered_product_density_supremum(case[[1L]], case[[2L]], radius))
    reference <- grid_maximum(case[[5L]], case[[2L]], case[[3L]], case[[4L]])
    expect_gte(bound, reference * (1 - 1e-12))
    expect_lt(bound, 1.01 * reference)
  }
  # a family without an analytic supremum, and a density that is infinite at the bound, have none
  expect_true(is.na(BayesTools:::.prior_ordered_product_density_supremum(
    prior("lognormal", list(0, 1), truncation = list(lower = 0, upper = 2)), 2)))
  expect_true(is.na(BayesTools:::.prior_ordered_product_density_supremum(prior("beta", list(.5, .5)), 0)))
})

test_that("a bounded total of a narrow support is rejected where the rounding of the argument across its bounds matters", {

  # T ~ U(1, 1 + w) times U V: f_X(x) = (-log(x) log(1 + w) + log(1 + w)^2 / 2) / w for x <= 1, with the
  # width of the stored support. An argument that rounds across a bound changes the density from 1 / w to
  # 0, which is a relative 2.4e-4 of the density at w = 1e-12 (x = .7), 2.4e-3 at 1e-13 (x = .5) and 1.7e-2
  # at 1e-14 (x = .7): such values are rejected, the others certified within the criterion
  reference <- function(x, width){
    actual <- (1 + width) - 1
    (-log(x) * log1p(actual) + log1p(actual)^2 / 2) / actual
  }
  x <- c(.1, .3, .5, .7, .9, .99)
  for(width in c(1e-6, 1e-8, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14)){
    route <- .bounded_cross_route(prior("uniform", list(1, 1 + width)))
    result <- BayesTools:::.prior_ordered_product_integrate_density(route$spec, x)
    accepted <- result$accepted
    if(any(accepted)){
      expect_lt(max(abs(result$estimate / reference(x, width) - 1)[accepted]), 1e-4, label = width)
    }
    if(width >= 1e-9){
      expect_true(all(accepted), info = width)
    }
  }
  for(case in list(c(1e-12, .7), c(1e-13, .5), c(1e-14, .7))){
    route <- .bounded_cross_route(prior("uniform", list(1, 1 + case[[1L]])))
    result <- BayesTools:::.prior_ordered_product_integrate_density(route$spec, case[[2L]])
    expect_false(result$accepted, info = case)
    expect_true(result$unresolved, info = case)
    ordinate <- BayesTools:::.prior_density_route_ordinate(route, case[[2L]])
    expect_false(ordinate$exact, info = case)
    expect_match(ordinate$reason, "not representable at this value", fixed = TRUE)
  }
  # the same with a Beta total whose bound is regular
  route <- .bounded_cross_route(prior("beta", list(1, 1)))
  result <- BayesTools:::.prior_ordered_product_integrate_density(route$spec, c(.3, .7))
  expect_true(all(result$accepted))
  expect_equal(result$estimate, log(c(.3, .7))^2 / 2, tolerance = 1e-8)
})

test_that("the inner integral at the unresolved nodes of the last outer coordinate is capped, and not bounded otherwise", {

  # T ~ Beta(1, .5) times U and V: the inner integral over V at every U (also at the unrounded value of
  # an unresolved outer node) is below the cap sup(beta_p |c|) / |kappa| of the substitution to t
  route <- .bounded_cross_route(prior("beta", list(1, .5)))
  spec <- route$spec
  inner <- spec$structure$coordinates[[2L]]
  integral <- function(u, x){
    low <- BayesTools:::.prior_ordered_product_state(spec, cbind(u, 0))
    high <- BayesTools:::.prior_ordered_product_state(spec, cbind(u, 1))
    breakpoint <- (x - low$scale) / (high$scale - low$scale)
    pieces <- c(0, if(breakpoint > 0 && breakpoint < 1) breakpoint, 1)
    sum(vapply(seq_len(length(pieces) - 1L), function(i){
      stats::integrate(function(v){
        stats::dbeta(v, inner$a, inner$b) * BayesTools:::.prior_ordered_product_scale_density(
          spec, x, low$shift + (high$shift - low$shift) * v, low$scale + (high$scale - low$scale) * v)
      }, pieces[[i]], pieces[[i + 1L]], rel.tol = 1e-9, subdivisions = 1000L)$value
    }, numeric(1)))
  }
  for(x in c(.05, .3, .5, .9)){
    # (the adaptive reference integral resolves a breakpoint that is 1e-6 or more from the end)
    values <- unique(c(x + 10^-(2:6), x - 10^-(2:6), seq(.01, 1, length.out = 12)))
    values <- values[values > 0 & values <= 1]
    low <- BayesTools:::.prior_ordered_product_state(spec, cbind(values, 0))
    high <- BayesTools:::.prior_ordered_product_state(spec, cbind(values, 1))
    caps <- BayesTools:::.prior_ordered_product_inner_cap(inner, low, high, x)
    exact <- vapply(values, integral, numeric(1), x = x)
    expect_true(all(is.finite(caps)), info = x)
    expect_true(all(exact <= caps * (1 + 1e-6)), info = x)
  }
  # shares with a density that is infinite at an end, and a cancelling kappa, have no cap
  low <- list(scale = 1, shift = 0)
  high <- list(scale = 2, shift = 0)
  expect_identical(BayesTools:::.prior_ordered_product_inner_cap(list(a = 1, b = .1), low, high, .5), Inf)
  expect_identical(BayesTools:::.prior_ordered_product_inner_cap(list(a = 1, b = 1), low, high, 0), Inf)

  # T ~ Beta(1, .5), U uniform, V ~ Beta(1, .1): the inner integral grows as d^(e + b - 1) at the
  # breakpoint of the outer coordinate (the unrounded nodes carry at least 3.76e-10 of the density at
  # x = .5 with 64 nodes, which the rounded nodes did not), so the value has no bound and is rejected
  priors <- .interaction_priors(~ f * o, list(intercept = N01, f = prior_ordered(N01), o = prior_ordered(N01),
    "f:o" = prior_ordered(prior("beta", list(1, .5)), allocation = list(
      f = prior("dirichlet", list(c(1, .5, .5))), o = prior("dirichlet", list(c(1, .1)))))))
  singular <- .interaction_route(priors[["mu_f__xXx__o"]], "f=mid, o=q")
  expect_identical(singular$spec$structure$coordinates[[2L]]$b, .1)
  ordinate <- BayesTools:::.prior_density_route_ordinate(singular, .5)
  expect_false(ordinate$exact)
  expect_match(ordinate$reason, "not representable at this value", fixed = TRUE)
  # an unresolved node of a coordinate before the last one has no bound either
  expect_identical(
    BayesTools:::.prior_ordered_product_unresolved_bias(
      spec, list(row = integer(), w = numeric()), numeric(), logical(), matrix(FALSE, 0L, 2L), numeric(),
      numeric(), numeric(), c(0, 1), c(1, .5), 1, 1, 5, 0L, 1L, 1L, 1L),
    Inf)
})
