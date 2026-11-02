skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Hypothesis tests on the totals and shares of ordered terms
# ============================================================================ #
#
# The totals (one per slice) and the shares of the random allocations of a
# prior_ordered() term are catalog quantities that hypothesis_BF() tests, by
# their labels (`f{total}`, `f[share: mid]`), their earlier labels
# (`f_ordered_total`, `f_ordered_allocation[mid]`), and their canonical names.
# The fits are synthetic, with posterior draws of known laws, so that each
# Bayes factor is compared with its analytic value.
#
# ============================================================================ #

N01 <- prior("normal", list(0, 1))

# A fit of draws from the priors of its ordered terms (totals, spike-and-slab
# indicators, gamma allocation nodes, and the increments they give) and of
# standard normal draws of the other coefficients. 'modify' receives the
# matrix of the ordered sources and the gamma nodes and returns it with the
# columns that the test fixes at known laws; the increments follow the
# modified columns.
.ordered_hypothesis_fit <- function(formula, priors, data, n = 20000L,
                                    modify = NULL, monitor_gamma = TRUE,
                                    formula_scale = NULL, seed = 1L){

  info  <- JAGS_formula(formula, "mu", data, priors, formula_scale = formula_scale)
  bound <- info$prior_list
  set.seed(seed)
  sources <- list()
  other   <- list()
  for(name in names(bound)){
    prior <- bound[[name]]
    if(!is.prior.ordered(prior)){
      columns <- if(is.prior.factor(prior)) .JAGS_prior_factor_names(name, prior) else name
      other[[name]] <- matrix(stats::rnorm(n * length(columns)), n, length(columns),
                              dimnames = list(NULL, columns))
      next
    }
    # a total with a parameter expression has no draws of its own: its column is
    # standard normal here, and the caller fixes it
    draws <- tryCatch(.prior_ordered_draws(prior, n), error = function(e) NULL)
    sources[[name]] <- if(is.null(draws)){
      total_names <- .prior_ordered_total_monitor_names(prior, name)
      matrix(stats::rnorm(n * length(total_names)), n, length(total_names),
             dimnames = list(NULL, total_names))
    }else{
      .prior_ordered_total_samples(draws, name)
    }
    for(record in .prior_ordered_dirichlet_records(prior)){
      eta_name <- .JAGS_prior_dirichlet_eta_name(record$node)
      if(!is.null(sources[[eta_name]])) next
      sources[[eta_name]] <- if(is.null(draws)){
        sapply(record$spec$alpha, function(alpha) stats::rgamma(n, alpha))
      }else{
        draws$allocation_samples[[record$key]] * stats::rexp(n)
      }
      colnames(sources[[eta_name]]) <- paste0(eta_name, "[", seq_len(record$dim), "]")
    }
  }
  sources <- do.call(cbind, unname(sources))
  if(!is.null(modify)) sources <- modify(sources)
  increments <- list()
  for(name in names(bound)[vapply(bound, is.prior.ordered, logical(1))]){
    increments[[name]] <- .bt_deterministic_node_evaluate(
      .bt_dnode_ordered_coefficients(.bt_ordered_spec(name, bound[[name]])),
      .bt_deterministic_lookup(sources)
    )
  }
  values <- do.call(cbind, unname(c(other, increments, list(sources))))
  if(!monitor_gamma) values <- values[, !startsWith(colnames(values), "prior_par_eta_"), drop = FALSE]
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(values)), sample = n,
         summary.pars = list(mutate = NULL), monitor = colnames(values)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list")     <- bound
  attr(fit, "formula_design") <- list(mu = info$formula_design)
  if(length(info$formula_scale) > 0L) attr(fit, "formula_scale") <- list(mu = info$formula_scale)

  attach_test_parameter_map(fit)
}

.ordered_hypothesis_data <- function(levels = c("lo", "mid", "hi")){

  data.frame(
    f = ordered(rep(levels, 3L), levels = levels),
    o = ordered(rep(c("p", "q", "r"), each = length(levels)), levels = c("p", "q", "r")),
    g = factor(rep(c("A", "B", "C"), length.out = 3L * length(levels))),
    x = sin(seq_len(3L * length(levels)))
  )
}

# The fit of ~ f with the total prior 'total' and a Dirichlet allocation.
.ordered_hypothesis_single <- function(total = N01, allocation = NULL, levels = c("lo", "mid", "hi"), ...){

  .ordered_hypothesis_fit(
    ~ f, list(intercept = N01, f = prior_ordered(total, allocation = allocation)),
    .ordered_hypothesis_data(levels), ...
  )
}

# The raw Bayes factor of a statement, with its warnings.
.ordered_hypothesis_bf <- function(fit, hypothesis, ...){

  warnings <- character()
  result <- withCallingHandlers(
    hypothesis_BF(fit, hypothesis, seed = 1, ...),
    warning = function(condition){
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )

  list(BF = attr(result, "raw_BF", exact = TRUE), result = result, warnings = warnings)
}

.ordered_gamma <- function(sources, ...){

  sources[, grep("^prior_par_eta_", colnames(sources))] <- cbind(...)

  sources
}

test_that("ordered totals and shares are catalog quantities with support and prior ordinates", {

  levels <- c("a", "b", "c", "d")
  alpha  <- c(2, 1, .5)
  fit <- .ordered_hypothesis_single(
    allocation = prior("dirichlet", list(alpha)), levels = levels
  )
  catalog <- parameter_catalog(fit)
  quantities <- catalog$quantities
  expect_identical(
    quantities$canonical_name[grepl("ordered", quantities$canonical_name)],
    c("mu_f_ordered_total", paste0("mu_f_ordered_allocation[", levels[-1L], "]"))
  )
  shares <- quantities[grepl("allocation", quantities$canonical_name), ]
  expect_identical(shares$status, rep("derived", 3L))
  expect_identical(shares$display_label, paste0("(mu) f[share: ", levels[-1L], "]"))
  expect_identical(shares$source_type, rep("composite", 3L))

  # the labels, the earlier labels, and the canonical names select one quantity
  for(level in levels[-1L]){
    canonical <- paste0("mu_f_ordered_allocation[", level, "]")
    for(alias in c(canonical, paste0("f[share: ", level, "]"), paste0("(mu) f[share: ", level, "]"),
                   paste0("f_ordered_allocation[", level, "]"), paste0("(mu) f_ordered_allocation[", level, "]"))){
      expect_identical(parameter_catalog_resolve(catalog, alias)$quantities$canonical_name, canonical,
                       info = alias)
    }
  }
  for(alias in c("mu_f_ordered_total", "f{total}", "(mu) f{total}", "f_ordered_total", "(mu) f_ordered_total")){
    expect_identical(parameter_catalog_resolve(catalog, alias)$quantities$canonical_name,
                     "mu_f_ordered_total", info = alias)
  }

  # exact supports: [0, 1] for a share, the support of the total prior
  expect_identical(unclass(shares$support[[1L]])$bounds, c(0, 1))
  expect_true(shares$support[[1L]]$exact)
  total <- quantities[quantities$canonical_name == "mu_f_ordered_total", ]
  expect_identical(unclass(total$support[[1L]])$bounds, c(-Inf, Inf))

  # the draws are the normalized gamma nodes
  eta <- as.matrix(fit$mcmc)[, paste0("prior_par_eta_mu_f_ordered_alloc_f_1[", 1:3, "]")]
  for(j in 1:3){
    selection <- parameter_catalog_resolve(catalog, shares$canonical_name[[j]])
    expect_equal(as.numeric(as.matrix(parameter_draws(fit, selection))), unname(eta[, j] / rowSums(eta)),
                 tolerance = 1e-14)
    # the prior ordinate of a share is the Beta marginal of the Dirichlet prior
    density <- parameter_prior_density(fit, selection)
    for(x in c(.2, .7)){
      expect_equal(exp(prior_density_ordinate(density, x)$log_density),
                   stats::dbeta(x, alpha[[j]], sum(alpha) - alpha[[j]]),
                   tolerance = 1e-4, info = paste(j, x))
    }
  }
  total_density <- parameter_prior_density(fit, parameter_catalog_resolve(catalog, "f{total}"))
  expect_equal(exp(prior_density_ordinate(total_density, .3)$log_density), stats::dnorm(.3), tolerance = 1e-6)

  # shares of the quantities of a fit whose gamma nodes are not monitored are
  # unavailable: only the total remains a quantity
  unmonitored <- .ordered_hypothesis_single(monitor_gamma = FALSE)
  expect_identical(
    parameter_catalog(unmonitored)$quantities$canonical_name[
      grepl("ordered", parameter_catalog(unmonitored)$quantities$canonical_name)],
    "mu_f_ordered_total"
  )
  expect_error(parameter_catalog_resolve(parameter_catalog(unmonitored), "f[share: mid]"),
               class = "BayesTools_parameter_not_found")
})

test_that("a share is tested against its known Beta posterior", {

  # the gamma draws (3, 1) normalize to a Beta(3, 1) share of the middle level
  # and a Beta(1, 3) share of the top level; the prior is Beta(1, 1)
  n <- 40000L
  fit <- .ordered_hypothesis_single(n = n, modify = function(sources){
    .ordered_gamma(sources, stats::rgamma(n, 3), stats::rgamma(n, 1))
  })
  # region: P(share > .5) is .875 a posteriori and .5 a priori, BF = 7. The
  # posterior odds of n draws have the relative error 1 / sqrt(n p (1 - p)),
  # 1.5% for n = 40000: the tolerance is four such errors.
  region <- .ordered_hypothesis_bf(fit, "f[share: mid] > 0.5")
  expect_equal(region$BF, (1 - stats::pbeta(.5, 3, 1)) / stats::pbeta(.5, 3, 1), tolerance = 0.06)
  # the labels, the earlier labels, and the canonical name test one quantity
  for(alias in c("(mu) f[share: mid] > 0.5", "f_ordered_allocation[mid] > 0.5", "mu_f_ordered_allocation[mid] > 0.5")){
    expect_equal(.ordered_hypothesis_bf(fit, alias)$BF, region$BF, tolerance = 1e-12, info = alias)
  }
  complement <- .ordered_hypothesis_bf(fit, "f[share: hi] < 0.5")
  expect_equal(complement$BF, region$BF, tolerance = 1e-12)
  # point: the prior ordinate over the posterior ordinate at .5
  point <- .ordered_hypothesis_bf(fit, "(mu) f[share: mid] = 0.5")
  expect_equal(point$BF, stats::dbeta(.5, 1, 1) / stats::dbeta(.5, 3, 1), tolerance = 0.05)
  # the share at its bounds follows the rules of the other bounded quantities
  expect_error(.ordered_hypothesis_bf(fit, "f[share: hi] > 0"),
               "Prior region mass for hypothesis '(mu) f[share: hi] > 0' is one", fixed = TRUE,
               class = "BayesTools_prior_region_mass_unavailable")
})

test_that("the total is tested against its known normal posterior", {

  # prior N(0, 1); posterior draws N(.8, .5)
  n <- 40000L
  fit <- .ordered_hypothesis_single(n = n, modify = function(sources){
    sources[, "mu_f_ordered_total"] <- stats::rnorm(n, .8, .5)
    sources
  })
  prior_density <- function(x) stats::dnorm(x, 0, 1)
  posterior_density <- function(x) stats::dnorm(x, .8, .5)
  point <- .ordered_hypothesis_bf(fit, "f{total} = 0")
  expect_equal(point$BF, prior_density(0) / posterior_density(0), tolerance = 0.06)
  expect_equal(.ordered_hypothesis_bf(fit, "(mu) f{total} = 0.8")$BF,
               prior_density(.8) / posterior_density(.8), tolerance = 0.02)
  # the posterior mass of the region is .945: four relative errors of its odds
  region <- .ordered_hypothesis_bf(fit, "f{total} > 0")
  expect_equal(region$BF, stats::pnorm(.8 / .5) / (1 - stats::pnorm(.8 / .5)), tolerance = 0.09)
  # earlier labels and the canonical name
  for(alias in c("f_ordered_total", "(mu) f_ordered_total", "mu_f_ordered_total")){
    expect_equal(.ordered_hypothesis_bf(fit, paste(alias, "= 0"))$BF, point$BF, tolerance = 1e-12, info = alias)
  }
  # the listed quantities
  quantities <- hypothesis_quantities(fit)
  expect_true(all(c("(mu) f{total}", "(mu) f[share: mid]", "(mu) f[share: hi]") %in% quantities$alias))
  expect_true(all(quantities$point_test[quantities$alias %in% c("(mu) f{total}", "(mu) f[share: mid]")]))
})

test_that("the total equals the top level of a cumulative fit", {

  n <- 40000L
  posterior <- function(sources){
    sources[, "mu_f_ordered_total"] <- stats::rnorm(n, .8, .5)
    sources
  }
  # a slab-only prior, unconditionally
  fit <- .ordered_hypothesis_single(n = n, modify = posterior)
  for(pair in list(c("f[hi] = 0", "f{total} = 0"), c("f[hi] < 0", "f{total} < 0"))){
    expect_equal(.ordered_hypothesis_bf(fit, pair[[1L]])$BF, .ordered_hypothesis_bf(fit, pair[[2L]])$BF,
                 tolerance = 1e-12, info = pair[[2L]])
  }
  # a spike-and-slab prior, within the inclusion event
  spike <- prior_spike_and_slab(N01, prior("point", list(.5)))
  fit <- .ordered_hypothesis_single(total = spike, n = n, modify = function(sources){
    sources[, "mu_f_ordered_total_indicator"] <- stats::rbinom(n, 1L, .5)
    sources[, "mu_f_ordered_total_variable"]  <- stats::rnorm(n, .8, .5)
    sources[, "mu_f_ordered_total"] <- sources[, "mu_f_ordered_total_indicator"] *
      sources[, "mu_f_ordered_total_variable"]
    sources
  })
  for(pair in list(c("f[hi] = 0", "f{total} = 0"), c("f[hi] < 0", "f{total} < 0"))){
    expect_equal(.ordered_hypothesis_bf(fit, pair[[1L]], conditional = TRUE)$BF,
                 .ordered_hypothesis_bf(fit, pair[[2L]], conditional = TRUE)$BF,
                 tolerance = 1e-12, info = pair[[2L]])
  }
})

test_that("spike-and-slab and mixture totals refuse point tests of the null without conditioning", {

  n <- 40000L
  included <- function(sources){
    sources[, "mu_f_ordered_total_indicator"] <- stats::rbinom(n, 1L, .5)
    sources[, "mu_f_ordered_total_variable"]  <- stats::rnorm(n, .8, .5)
    sources[, "mu_f_ordered_total"] <- sources[, "mu_f_ordered_total_indicator"] *
      sources[, "mu_f_ordered_total_variable"]
    sources
  }
  spike <- prior_spike_and_slab(N01, prior("point", list(.5)))
  fit <- .ordered_hypothesis_single(total = spike, n = n, modify = included)
  expect_error(.ordered_hypothesis_bf(fit, "f{total} = 0"), class = "BayesTools_point_mass_at_null")
  # the message names the quantity by its label
  condition <- tryCatch(.ordered_hypothesis_bf(fit, "f{total} = 0"), error = function(e) e)
  expect_match(conditionMessage(condition), "ordered quantity '(mu) f{total}'", fixed = TRUE)
  expect_false(grepl("mu_f_ordered_total", conditionMessage(condition), fixed = TRUE))
  # within the alternative the slab is tested against its own posterior N(.8, .5)
  conditional <- .ordered_hypothesis_bf(fit, "f{total} = 0.8", conditional = TRUE)
  expect_equal(conditional$BF, stats::dnorm(.8, 0, 1) / stats::dnorm(.8, .8, .5), tolerance = 0.06)
  # an unconditional region includes the models in which the total is zero: the
  # prior mass of 'total > 0' is .25 (inclusion .5 times the slab half), and the
  # posterior mass is the share of positive draws
  region <- .ordered_hypothesis_bf(fit, "f{total} > 0")
  expect_match(region$warnings, "full ensemble for (mu) f{total}", fixed = TRUE, all = FALSE)
  posterior_mass <- mean(as.matrix(fit$mcmc)[, "mu_f_ordered_total"] > 0)
  expect_equal(region$BF, (posterior_mass / (1 - posterior_mass)) / (.25 / .75), tolerance = 1e-6)
  conditional_region <- .ordered_hypothesis_bf(fit, "f{total} > 0", conditional = TRUE)
  expect_length(conditional_region$warnings, 0L)
  included_mass <- mean(as.matrix(fit$mcmc)[as.matrix(fit$mcmc)[, "mu_f_ordered_total_indicator"] == 1, "mu_f_ordered_total"] > 0)
  expect_equal(conditional_region$BF, (included_mass / (1 - included_mass)) / 1, tolerance = 1e-6)

  # a mixture of a point null and an alternative
  mixture <- prior_mixture(list(prior("point", list(0)), N01), is_null = c(TRUE, FALSE))
  fit <- .ordered_hypothesis_single(total = mixture, n = n, modify = function(sources){
    sources[, "mu_f_ordered_total_indicator"] <- stats::rbinom(n, 1L, .5) + 1L
    # the indicator selects the component: the point null (1) or the alternative (2)
    sources[, "mu_f_ordered_total"] <- ifelse(sources[, "mu_f_ordered_total_indicator"] == 1L, 0,
                                              stats::rnorm(n, .8, .5))
    sources
  })
  expect_error(.ordered_hypothesis_bf(fit, "f{total} = 0"), class = "BayesTools_point_mass_at_null")
  expect_equal(.ordered_hypothesis_bf(fit, "f{total} = 0.8", conditional = TRUE)$BF,
               stats::dnorm(.8, 0, 1) / stats::dnorm(.8, .8, .5), tolerance = 0.06)

  # a point total is fixed by the model; its shares are still tested
  fixed <- .ordered_hypothesis_single(total = prior("point", list(2.5)), n = 2000L)
  expect_error(.ordered_hypothesis_bf(fixed, "f{total} = 2.5"), "is fixed by the fitted model",
               class = "BayesTools_hypothesis_fixed")
  expect_error(.ordered_hypothesis_bf(fixed, "f{total} > 1"), class = "BayesTools_hypothesis_fixed")
  expect_true(is.finite(.ordered_hypothesis_bf(fixed, "f[share: mid] > 0.5")$BF))

  # a total with a parameter expression has no prior law: it is refused, its shares are not
  expression_total <- .ordered_hypothesis_single(
    total = prior("normal", list(0, expression(tau))), n = 2000L,
    modify = function(sources){
      sources[, "mu_f_ordered_total"] <- stats::rnorm(nrow(sources), .5, 1)
      sources
    })
  expect_error(.ordered_hypothesis_bf(expression_total, "f{total} > 0"), "parameter expression",
               class = "BayesTools_hypothesis_target")
  expect_error(.ordered_hypothesis_bf(expression_total, "f{total} = 0"), "parameter expression",
               class = "BayesTools_hypothesis_target")
  expect_true(is.finite(.ordered_hypothesis_bf(expression_total, "f[share: mid] > 0.5")$BF))
})

test_that("shares under a spike-and-slab total follow the prior outside the inclusion event", {

  # the shares of the draws without the total follow the prior (uniform), those
  # of the included draws a Beta(3, 1) for the middle level
  n <- 40000L
  spike <- prior_spike_and_slab(N01, prior("point", list(.5)))
  fit <- .ordered_hypothesis_single(total = spike, n = n, modify = function(sources){
    indicator <- stats::rbinom(n, 1L, .5)
    sources[, "mu_f_ordered_total_indicator"] <- indicator
    sources[, "mu_f_ordered_total_variable"]  <- stats::rnorm(n, .8, .5)
    sources[, "mu_f_ordered_total"] <- indicator * sources[, "mu_f_ordered_total_variable"]
    .ordered_gamma(sources, ifelse(indicator == 1L, stats::rgamma(n, 3), stats::rgamma(n, 1)),
                   stats::rgamma(n, 1))
  })
  selection <- parameter_catalog_resolve(parameter_catalog(fit), "f[share: mid]")
  indicator <- as.matrix(fit$mcmc)[, "mu_f_ordered_total_indicator"]
  unconditional <- parameter_mixed_posterior(fit, selection, transform_scaled = FALSE)
  conditional   <- parameter_mixed_posterior(fit, selection, conditional = TRUE, transform_scaled = FALSE)
  expect_length(unconditional, n)
  expect_identical(.bt_meta_get(conditional, "draw_index"), which(indicator == 1))
  expect_equal(mean(conditional), 3 / 4, tolerance = 0.01)
  expect_equal(mean(unconditional[indicator == 0]), 1 / 2, tolerance = 0.02)
  expect_length(posterior_metadata(unconditional, "atoms")$mass, 0L)
  # within the event: P(share > .5 | Beta(3, 1)) = .875 against the prior .5
  # (about 20000 included draws: four relative errors of the posterior odds)
  within <- .ordered_hypothesis_bf(fit, "f[share: mid] > 0.5", conditional = TRUE)
  expect_equal(within$BF, 7, tolerance = 0.09)
  expect_length(within$warnings, 0L)
  # on the full ensemble the posterior mass is the mixture of the prior and the Beta(3, 1)
  full <- .ordered_hypothesis_bf(fit, "f[share: mid] > 0.5")
  expect_match(full$warnings, "the total of its term is fixed by its null component", fixed = TRUE, all = FALSE)
  mass <- .5 * .5 + .5 * (1 - stats::pbeta(.5, 3, 1))
  expect_equal(full$BF, (mass / (1 - mass)) / 1, tolerance = 0.05)
})

test_that("sliced, crossed, and shared ordered terms have a quantity per total and share", {

  data <- .ordered_hypothesis_data()
  treatment <- prior_factor("normal", list(0, 1), contrast = "treatment")
  bf <- function(fit, hypothesis) .ordered_hypothesis_bf(fit, hypothesis)$BF

  # slices of an ordinary factor: a total and a set of shares per slice
  sliced <- .ordered_hypothesis_fit(~ f * g, list(intercept = N01, f = prior_ordered(N01), g = treatment,
                                                  "f:g" = prior_ordered(N01)), data)
  quantities <- parameter_catalog(sliced)$quantities
  expect_identical(
    quantities$canonical_name[grepl("ordered", quantities$canonical_name)],
    c("mu_f_ordered_total", "mu_f__xXx__g_ordered_total[1]", "mu_f__xXx__g_ordered_total[2]",
      "mu_f_ordered_allocation[mid]", "mu_f_ordered_allocation[hi]",
      "mu_f__xXx__g_ordered_allocation[1][mid]", "mu_f__xXx__g_ordered_allocation[1][hi]",
      "mu_f__xXx__g_ordered_allocation[2][mid]", "mu_f__xXx__g_ordered_allocation[2][hi]")
  )
  expect_identical(
    quantities$display_label[grepl("ordered", quantities$canonical_name)],
    c("(mu) f{total}", "(mu) f{total}:g{1}", "(mu) f{total}:g{2}",
      "(mu) f[share: mid]", "(mu) f[share: hi]",
      "(mu) f[share: mid]:g{1}", "(mu) f[share: hi]:g{1}",
      "(mu) f[share: mid]:g{2}", "(mu) f[share: hi]:g{2}")
  )
  # the statements by label, earlier label, and canonical name test the same quantity
  expect_equal(bf(sliced, "f{total}:g{2} > 0"), bf(sliced, "f:g_ordered_total[2] > 0"), tolerance = 1e-12)
  expect_equal(bf(sliced, "f{total}:g{2} > 0"), bf(sliced, "mu_f__xXx__g_ordered_total[2] > 0"), tolerance = 1e-12)
  expect_equal(bf(sliced, "(mu) f[share: hi]:g{1} < 0.5"), bf(sliced, "f:g_ordered_allocation[1][hi] < 0.5"),
               tolerance = 1e-12)
  expect_equal(bf(sliced, "f[share: hi]:g{1} < 0.5"), bf(sliced, "mu_f__xXx__g_ordered_allocation[1][hi] < 0.5"),
               tolerance = 1e-12)
  # a slice is a different quantity than the term without its slice
  expect_false(isTRUE(all.equal(bf(sliced, "f{total}:g{1} > 0"), bf(sliced, "f{total} > 0"))))

  # two ordered factors in one term: a share per factor
  crossed <- .ordered_hypothesis_fit(~ f * o, list(intercept = N01, f = prior_ordered(N01),
    o = prior_ordered(N01), "f:o" = prior_ordered(N01)), data)
  quantities <- parameter_catalog(crossed)$quantities
  expect_identical(
    quantities$display_label[grepl("ordered", quantities$canonical_name)],
    c("(mu) f{total}", "(mu) o{total}", "(mu) f:o{total}",
      "(mu) f[share: mid]", "(mu) f[share: hi]", "(mu) o[share: q]", "(mu) o[share: r]",
      "(mu) f[share: mid]:o", "(mu) f[share: hi]:o", "(mu) f:o[share: q]", "(mu) f:o[share: r]")
  )
  expect_equal(bf(crossed, "f[share: mid]:o > 0.5"), bf(crossed, "f:o_ordered_allocation_f[mid] > 0.5"),
               tolerance = 1e-12)
  expect_equal(bf(crossed, "f:o[share: q] > 0.5"), bf(crossed, "mu_f__xXx__o_ordered_allocation_o[q] > 0.5"),
               tolerance = 1e-12)
  expect_true(is.finite(bf(crossed, "f:o{total} > 0")))

  # a shared allocation has the shares of every term that uses it, over the same coordinates
  shared <- .ordered_hypothesis_fit(~ f + x + f:x, list(intercept = N01, f = prior_ordered(N01, id = "shape"),
    x = N01, "f:x" = prior_ordered(N01, id = "shape")), data)
  quantities <- parameter_catalog(shared)$quantities
  rows <- quantities[grepl("allocation", quantities$canonical_name), ]
  expect_identical(rows$display_label, c("(mu) f[share: mid]", "(mu) f[share: hi]",
                                         "(mu) f[share: mid]:x", "(mu) f[share: hi]:x"))
  expect_identical(rows$extraction_key[[1L]]$dependencies, rows$extraction_key[[3L]]$dependencies)
  expect_equal(bf(shared, "f[share: mid] > 0.5"), bf(shared, "f[share: mid]:x > 0.5"), tolerance = 1e-12)
  expect_true(is.finite(bf(shared, "f:x{total} > 0")))
})

test_that("statements on ordered quantities are rewritten without leaking their roots", {

  n <- 20000L
  fit <- .ordered_hypothesis_single(n = n, modify = function(sources){
    .ordered_gamma(sources, stats::rgamma(n, 3), stats::rgamma(n, 1))
  })
  roots <- c("mu_f_ordered_total{", "mu_f_ordered_allocation{", "ordered_total{", "ordered_allocation{")
  leaks <- function(text) any(vapply(roots, function(root) any(grepl(root, text, fixed = TRUE)), logical(1)))
  # the reported rows keep the statement as written
  result <- hypothesis_BF(fit, "f[share: mid] > 0.5", seed = 1)
  expect_identical(result$Alternative, "f[share: mid] > 0.5")
  expect_identical(rownames(result), "(mu) f[share: mid]")
  expect_false(leaks(c(result$Alternative, result$Null, rownames(result))))
  # a statement on two quantities names them by their labels
  condition <- tryCatch(hypothesis_BF(fit, "f[share: mid] > f[share: hi]"), error = function(e) e)
  expect_s3_class(condition, "BayesTools_hypothesis_ambiguous")
  expect_match(conditionMessage(condition), "(mu) f[share: mid], (mu) f[share: hi]", fixed = TRUE)
  expect_false(leaks(conditionMessage(condition)))
  condition <- tryCatch(hypothesis_BF(fit, "f{total} > 0 & f[share: hi] > 0.5"), error = function(e) e)
  expect_s3_class(condition, "BayesTools_hypothesis_ambiguous")
  expect_false(leaks(conditionMessage(condition)))
  # a refusal of the evaluation names the quantity by its label
  condition <- tryCatch(hypothesis_BF(fit, "f[share: hi] > 0"), error = function(e) e)
  expect_s3_class(condition, "BayesTools_prior_region_mass_unavailable")
  expect_false(leaks(conditionMessage(condition)))
  # conditioning needs an inclusion event, as for any other quantity
  expect_error(hypothesis_BF(fit, "f[share: mid] > 0.5", conditional = TRUE), "no inclusion gate")
  # the listed roots are rendered as statements
  quantities <- hypothesis_quantities(fit)
  expect_identical(quantities$parameter[quantities$alias == "(mu) f[share: mid]"], "mu_f_ordered_allocation{mid}")
  expect_identical(quantities$parameter[quantities$alias == "(mu) f{total}"], "mu_f_ordered_total")

  # levels with syntax-sensitive characters
  odd <- .ordered_hypothesis_single(levels = c("a", "b ]", "c, {d}"), n = n)
  odd_quantities <- parameter_catalog(odd)$quantities
  expect_identical(odd_quantities$display_label[grepl("allocation", odd_quantities$canonical_name)],
                   c("(mu) f[share: b ]]", "(mu) f[share: c, {d}]"))
  expect_true(is.finite(.ordered_hypothesis_bf(odd, "(mu) f[share: c, {d}] > 0.5")$BF))
})

test_that("the total is tested on the original scale only when its term is unscaled", {

  data <- .ordered_hypothesis_data()
  priors <- list(intercept = N01, f = prior_ordered(N01), x = N01, "f:x" = prior_ordered(N01))
  fit <- .ordered_hypothesis_fit(~ f + x + f:x, priors, data, formula_scale = list(x = TRUE))
  for(label in c("f{total}", "f:x{total}")){
    expect_true(is.finite(.ordered_hypothesis_bf(fit, paste(label, "> 0"), transform_scaled = FALSE)$BF))
    expect_error(.ordered_hypothesis_bf(fit, paste(label, "> 0"), transform_scaled = TRUE),
                 "defined on the fitted scale of the predictors", class = "BayesTools_hypothesis_target")
  }
  # the shares are unitless
  expect_equal(.ordered_hypothesis_bf(fit, "f[share: mid] > 0.5", transform_scaled = TRUE)$BF,
               .ordered_hypothesis_bf(fit, "f[share: mid] > 0.5", transform_scaled = FALSE)$BF,
               tolerance = 1e-12)
  # a term that the scaling leaves unchanged has the same total on both scales
  unscaled_priors <- list(intercept = N01, f = prior_ordered(N01), x = N01)
  unscaled <- .ordered_hypothesis_fit(~ f + x, unscaled_priors, data, formula_scale = list(x = TRUE))
  expect_equal(.ordered_hypothesis_bf(unscaled, "f{total} > 0", transform_scaled = TRUE)$BF,
               .ordered_hypothesis_bf(unscaled, "f{total} > 0", transform_scaled = FALSE)$BF,
               tolerance = 1e-12)
})
