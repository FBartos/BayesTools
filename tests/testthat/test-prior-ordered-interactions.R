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
  expect_error(BayesTools:::.prior_linear_combination_density(
    list(mu_f__xXx__g = sliced),
    .cell_combination(sliced, c("f=mid, g=B", "f=mid, g=C"), c(1, -1))),
    class = "BayesTools_formula_measure_unavailable")
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
