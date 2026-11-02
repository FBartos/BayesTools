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

test_that("conditional ordered totals take the support of their restricted prior", {

  # spike at 0 and a Uniform(1, 2) slab: the draws of the inclusion event lie
  # in [1, 2], whose exact support the boundary correction of the posterior
  # ordinate needs (the unconditional support has the point 0 and is inexact)
  n <- 40000L
  spike <- prior_spike_and_slab(prior("uniform", list(1, 2)), prior("point", list(.5)))
  included <- function(slab){
    function(sources){
      sources[, "mu_f_ordered_total_indicator"] <- stats::rbinom(n, 1L, .5)
      sources[, "mu_f_ordered_total_variable"]  <- slab(n)
      sources[, "mu_f_ordered_total"] <- sources[, "mu_f_ordered_total_indicator"] *
        sources[, "mu_f_ordered_total_variable"]
      sources
    }
  }
  # the posterior is the prior: the Bayes factor is 1 at the boundaries and inside
  fit <- .ordered_hypothesis_single(total = spike, n = n, modify = included(function(n) stats::runif(n, 1, 2)))
  selection <- parameter_catalog_resolve(parameter_catalog(fit), "f{total}")
  expect_identical(unclass(selection$quantities$support[[1L]])$points, 0)
  expect_false(selection$quantities$support[[1L]]$exact)
  support <- posterior_metadata(
    parameter_mixed_posterior(fit, selection, conditional = TRUE, transform_scaled = FALSE), "support"
  )
  expect_identical(unclass(support)$bounds, c(1, 2))
  expect_length(unclass(support)$points, 0L)
  expect_true(support$exact)
  for(value in c(1, 1.5, 2)){
    expect_equal(.ordered_hypothesis_bf(fit, paste("f{total} =", value), conditional = TRUE)$BF, 1,
                 tolerance = 0.05, info = value)
  }
  # the posterior density 2 (2 - x) on [1, 2] against the prior density 1: the
  # Bayes factors are 1 / 2 at the bound 1 and 1 at 1.5 (the kernel error at
  # the reflected boundary is tolerated up to 8%)
  fit <- .ordered_hypothesis_single(total = spike, n = n,
                                    modify = included(function(n) 1 + stats::rbeta(n, 1, 2)))
  expect_equal(.ordered_hypothesis_bf(fit, "f{total} = 1", conditional = TRUE)$BF, 1 / 2, tolerance = 0.08)
  expect_equal(.ordered_hypothesis_bf(fit, "f{total} = 1.5", conditional = TRUE)$BF, 1, tolerance = 0.05)
  # a slab-only bounded prior (the unconditional support is exact)
  fit <- .ordered_hypothesis_single(total = prior("uniform", list(1, 2)), n = n)
  for(value in c(1, 1.5, 2)){
    expect_equal(.ordered_hypothesis_bf(fit, paste("f{total} =", value))$BF, 1, tolerance = 0.05, info = value)
  }
})

test_that("unconditional shares do not need the inclusion indicator of their total", {

  # a fit that did not monitor the indicator of the spike-and-slab total
  n <- 40000L
  spike <- prior_spike_and_slab(N01, prior("point", list(.5)))
  fit <- .ordered_hypothesis_single(total = spike, n = n, modify = function(sources){
    sources <- .ordered_gamma(sources, stats::rgamma(n, 3), stats::rgamma(n, 1))
    sources[, !grepl("_indicator$", colnames(sources)), drop = FALSE]
  })
  expect_false("mu_f_ordered_total_indicator" %in% colnames(as.matrix(fit$mcmc)))
  selection <- parameter_catalog_resolve(parameter_catalog(fit), "f[share: mid]")
  expect_identical(selection$quantities$status, "derived")
  expect_length(parameter_mixed_posterior(fit, selection, transform_scaled = FALSE), n)
  # the Beta(3, 1) posterior of the share against its Beta(1, 1) prior: BF 7
  unconditional <- .ordered_hypothesis_bf(fit, "f[share: mid] > 0.5")
  expect_equal(unconditional$BF, 7, tolerance = 0.06)
  # the inclusion event needs the indicator
  indicator_message <- paste0(
    "Ordered total indicator 'mu_f_ordered_total_indicator' is unavailable in 'draws'. ",
    "Include its monitored source coordinate."
  )
  expect_error(parameter_mixed_posterior(fit, selection, conditional = TRUE, transform_scaled = FALSE),
               indicator_message, fixed = TRUE)
  expect_error(.ordered_hypothesis_bf(fit, "f[share: mid] > 0.5", conditional = TRUE),
               indicator_message, fixed = TRUE)
})

test_that("planning errors of ordered quantities name their displayed labels", {

  # a normal total has no inclusion event to condition on
  fit <- .ordered_hypothesis_single(n = 1000L)
  event_message <- function(label){
    paste0("The inclusion event of '", label, "' is unavailable: the quantity has no ",
           "inclusion gate. Use 'conditional = FALSE'.")
  }
  for(statement in c("f[share: mid] > 0.5", "(mu) f[share: mid] > 0.5", "f_ordered_allocation[mid] > 0.5",
                     "mu_f_ordered_allocation[mid] > 0.5")){
    condition <- tryCatch(hypothesis_BF(fit, statement, conditional = TRUE), error = function(e) e)
    expect_identical(conditionMessage(condition), event_message("(mu) f[share: mid]"), info = statement)
  }
  for(statement in c("f{total} > 0", "f_ordered_total > 0", "mu_f_ordered_total > 0")){
    condition <- tryCatch(hypothesis_BF(fit, statement, conditional = TRUE), error = function(e) e)
    expect_identical(conditionMessage(condition), event_message("(mu) f{total}"), info = statement)
  }
  # sliced quantities
  data <- .ordered_hypothesis_data()
  sliced <- .ordered_hypothesis_fit(~ f * g, list(intercept = N01, f = prior_ordered(N01),
    g = prior_factor("normal", list(0, 1), contrast = "treatment"), "f:g" = prior_ordered(N01)), data, n = 1000L)
  for(statement in c("f{total}:g{2} > 0", "mu_f__xXx__g_ordered_total[2] > 0")){
    condition <- tryCatch(hypothesis_BF(sliced, statement, conditional = TRUE), error = function(e) e)
    expect_identical(conditionMessage(condition), event_message("(mu) f{total}:g{2}"), info = statement)
  }
  condition <- tryCatch(hypothesis_BF(sliced, "mu_f__xXx__g_ordered_allocation[2][hi] > 0.5", conditional = TRUE),
                        error = function(e) e)
  expect_identical(conditionMessage(condition), event_message("(mu) f[share: hi]:g{2}"))
  # levels with the characters of regular expressions and of their quoting (a
  # name is matched literally)
  for(level in c("a\\Eb", "a\\E(b", "a\\Qb", "a [b]", "a {b}", "a.b", "a|b", "a$1")){
    odd <- .ordered_hypothesis_single(levels = c("lo", level, "hi"), n = 1000L)
    for(statement in c(paste0("(mu) f[share: ", level, "] > 0.5"),
                       paste0("mu_f_ordered_allocation[", level, "] > 0.5"))){
      condition <- tryCatch(hypothesis_BF(odd, statement, conditional = TRUE), error = function(e) e)
      expect_identical(conditionMessage(condition), event_message(paste0("(mu) f[share: ", level, "]")),
                       info = paste(level, statement))
    }
  }
})

test_that("ordered totals with several bounded alternatives are tested within each component", {

  # a point null (component 1) and two bounded alternatives: the posterior
  # ordinate at a bound of a slab is that of the slab's own draws, reflected at
  # its own bounds, not the ordinate of the pooled alternatives
  n <- 80000L
  disjoint <- prior_mixture(list(prior("point", list(0), prior_weights = 5),
    prior("uniform", list(1, 2), prior_weights = 1), prior("uniform", list(3, 4), prior_weights = 3)),
    is_null = c(TRUE, FALSE, FALSE))
  overlapping <- prior_mixture(list(prior("point", list(0), prior_weights = 2),
    prior("uniform", list(1, 3), prior_weights = 1), prior("uniform", list(2, 4), prior_weights = 3)),
    is_null = c(TRUE, FALSE, FALSE))
  # posterior draws of the component probabilities 'probabilities'
  drawn <- function(probabilities, bounds){
    function(sources){
      component <- sample.int(3L, n, replace = TRUE, prob = probabilities)
      sources[, "mu_f_ordered_total_indicator"] <- component
      sources[, "mu_f_ordered_total"] <- 0
      for(i in 2:3){
        sources[component == i, "mu_f_ordered_total"] <- stats::runif(sum(component == i), bounds[[i - 1L]][[1L]], bounds[[i - 1L]][[2L]])
      }
      sources
    }
  }
  # Bayes factors of the null values: the prior density of a value inside the
  # slabs 'slabs' (weights, bounds) over the posterior density (the same
  # mixture with the component probabilities 'posterior'; both conditional on
  # the alternatives when 'conditional')
  density <- function(x, weights, bounds, conditional){
    weights <- weights[-1L] / (if(conditional) sum(weights[-1L]) else sum(weights))
    sum(vapply(seq_along(weights), function(i){
      weights[[i]] * stats::dunif(x, bounds[[i]][[1L]], bounds[[i]][[2L]])
    }, numeric(1)))
  }
  check <- function(total, bounds, prior_weights, posterior, values, conditional, label){
    fit <- .ordered_hypothesis_single(total = total, n = n, modify = drawn(posterior, bounds))
    for(value in values){
      expected <- density(value, prior_weights, bounds, conditional) / density(value, posterior, bounds, conditional)
      expect_equal(.ordered_hypothesis_bf(fit, paste("f{total} =", value), conditional = conditional)$BF,
                   expected, tolerance = 0.08, info = paste(label, conditional, value))
    }
  }
  disjoint_bounds <- list(c(1, 2), c(3, 4))
  overlapping_bounds <- list(c(1, 3), c(2, 4))
  for(conditional in c(TRUE, FALSE)){
    # the posterior is the prior: every Bayes factor is 1, at the bounds of the slabs too
    check(disjoint, disjoint_bounds, c(5, 1, 3), c(5, 1, 3), c(1, 1.5, 2, 3, 3.5, 4), conditional, "disjoint")
    check(overlapping, overlapping_bounds, c(2, 1, 3), c(2, 1, 3), c(1, 1.5, 2, 2.5, 3, 3.5, 4), conditional,
          "overlapping")
    # the posterior probabilities of the components differ from the prior ones
    check(disjoint, disjoint_bounds, c(5, 1, 3), c(2, 1, 1), c(1, 2, 3, 4), conditional, "disjoint reweighted")
    check(overlapping, overlapping_bounds, c(2, 1, 3), c(1, 2, 1), c(1, 2, 3, 4), conditional, "overlapping reweighted")
  }
})

# ---------------------------------------------------------------------------- #
# Savage-Dickey ordinates at the bounds of a single bounded alternative
# ---------------------------------------------------------------------------- #

# A fit of ~ x whose coefficient has the prior 'prior_x' and the draws 'values'
# (with the component indicator, spike-and-slab inclusion, and slab variable
# columns that the prior monitors, when given).
.ordered_hypothesis_coefficient_fit <- function(prior_x, values, indicator = NULL,
                                                inclusion = NULL, variable = NULL){

  info <- JAGS_formula(~ x, "mu", .ordered_hypothesis_data(), list(intercept = N01, x = prior_x))
  columns <- list(mu_intercept = stats::rnorm(length(values)), mu_x = values,
                  mu_x_indicator = indicator, mu_x_inclusion = inclusion, mu_x_variable = variable)
  columns <- do.call(cbind, columns[!vapply(columns, is.null, logical(1))])
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(columns)), sample = length(values),
         summary.pars = list(mutate = NULL), monitor = colnames(columns)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list")     <- info$prior_list
  attr(fit, "formula_design") <- list(mu = info$formula_design)

  attach_test_parameter_map(fit)
}

# Draws from a normal truncated to [lower, upper].
.ordered_hypothesis_truncated_normal <- function(n, mean, sd, lower, upper){

  stats::qnorm(stats::runif(n, stats::pnorm(lower, mean, sd), stats::pnorm(upper, mean, sd)), mean, sd)
}

.ordered_hypothesis_truncated_density <- function(x, mean, sd, lower, upper){

  stats::dnorm(x, mean, sd) / (stats::pnorm(upper, mean, sd) - stats::pnorm(lower, mean, sd))
}

test_that("a single bounded alternative is reflected at its bounds, also without conditioning", {

  # spike at 0 (inclusion probability 1/2) and a slab on [1, 2]: the posterior
  # ordinate at a bound of the slab is that of the slab's draws reflected at the
  # bound (the unconditional support of the spike and slab is not an interval)
  n <- 80000L
  set.seed(401)
  uniform_slab <- prior("uniform", list(1, 2))
  truncated_slab <- prior("normal", list(0, 1), list(1, 2))
  half <- prior("point", list(.5))
  indicator <- stats::rbinom(n, 1, .5)
  variable  <- stats::runif(n, 1, 2)
  # a posterior that differs from the prior: inclusion probability .7 and the
  # slab normal(0, 1.5) truncated to [1, 2]
  indicator_post <- stats::rbinom(n, 1, .7)
  variable_post  <- .ordered_hypothesis_truncated_normal(n, 0, 1.5, 1, 2)

  fits <- list(
    coefficient = list(
      same = .ordered_hypothesis_coefficient_fit(prior_spike_and_slab(uniform_slab, half),
        variable * indicator, indicator, rep(.5, n), variable),
      differs = .ordered_hypothesis_coefficient_fit(prior_spike_and_slab(truncated_slab, half),
        variable_post * indicator_post, indicator_post, rep(.5, n), variable_post)
    ),
    total = list(
      same = .ordered_hypothesis_single(total = prior_spike_and_slab(uniform_slab, half), n = n, seed = 719L),
      differs = .ordered_hypothesis_single(total = prior_spike_and_slab(truncated_slab, half), n = n, seed = 719L,
        modify = function(sources){
          sources[, "mu_f_ordered_total_indicator"] <- indicator_post
          sources[, "mu_f_ordered_total_variable"]  <- variable_post
          sources[, "mu_f_ordered_total"]           <- indicator_post * variable_post
          sources
        })
    )
  )
  names_of <- c(coefficient = "x", total = "f{total}")
  for(kind in names(fits)) for(conditional in c(TRUE, FALSE)) for(value in c(1, 1.5, 2)){
    statement <- paste(names_of[[kind]], "=", value)
    info <- paste(kind, conditional, value)
    # the posterior is the prior: the Bayes factor is 1 (it was 2 at the bounds)
    expect_equal(.ordered_hypothesis_bf(fits[[kind]]$same, statement, conditional = conditional)$BF, 1,
                 tolerance = 0.08, info = info)
    # the ratio of the prior and posterior densities of the slab, and of the
    # prior and posterior inclusion probabilities without conditioning
    expected <- .ordered_hypothesis_truncated_density(value, 0, 1, 1, 2) /
      .ordered_hypothesis_truncated_density(value, 0, 1.5, 1, 2) * if(conditional) 1 else .5 / .7
    expect_equal(.ordered_hypothesis_bf(fits[[kind]]$differs, statement, conditional = conditional)$BF,
                 expected, tolerance = 0.08, info = info)
  }

  # several alternatives with one support: the pooled draws are reflected at it
  set.seed(402)
  mixture <- prior_mixture(list(prior("point", list(0), prior_weights = 2),
    prior("uniform", list(1, 2), prior_weights = 1), prior("normal", list(0, 1), list(1, 2), prior_weights = 1)),
    is_null = c(TRUE, FALSE, FALSE))
  component <- sample.int(3L, n, replace = TRUE, prob = c(.5, .25, .25))
  values <- numeric(n)
  values[component == 2L] <- stats::runif(sum(component == 2L), 1, 2)
  values[component == 3L] <- .ordered_hypothesis_truncated_normal(sum(component == 3L), 0, 1, 1, 2)
  fit <- .ordered_hypothesis_coefficient_fit(mixture, values, indicator = component)
  for(conditional in c(TRUE, FALSE)) for(value in c(1, 1.5, 2)){
    expect_equal(.ordered_hypothesis_bf(fit, paste("x =", value), conditional = conditional)$BF, 1,
                 tolerance = 0.08, info = paste("two slabs", conditional, value))
  }

  # without an atom, or with an unbounded alternative, there is nothing to repair
  set.seed(403)
  slab_only <- .ordered_hypothesis_coefficient_fit(uniform_slab, stats::runif(n, 1, 2))
  indicator_u <- stats::rbinom(n, 1, .5)
  variable_u  <- stats::rnorm(n)
  unbounded <- .ordered_hypothesis_coefficient_fit(prior_spike_and_slab(N01, half),
    variable_u * indicator_u, indicator_u, rep(.5, n), variable_u)
  for(value in c(1, 1.5, 2)){
    expect_equal(.ordered_hypothesis_bf(slab_only, paste("x =", value))$BF, 1, tolerance = 0.08, info = value)
  }
  for(conditional in c(TRUE, FALSE)) for(value in c(-.5, .5)){
    expect_equal(.ordered_hypothesis_bf(unbounded, paste("x =", value), conditional = conditional)$BF, 1,
                 tolerance = 0.08, info = paste(conditional, value))
  }
})

test_that("the pooled ordinate is reflected at the common exact support of the continuous components", {

  point <- posterior_support_attribute(c(0, 0), points = 0, type = "points")
  unit <- posterior_support_attribute(c(1, 2))
  other <- posterior_support_attribute(c(1, 3))
  slab <- stats::runif(100, 1, 2)
  plan <- function(index, supports, samples = slab){
    BayesTools:::.Savage_Dickey_BF.component_plan(
      samples, index, BayesTools:::.posterior_components_new(index = index, supports = supports))
  }

  # one continuous component and several with one support are pooled at it
  expect_identical(plan(rep(2L, 100), list(point, unit))$pooled_bounds, c(1, 2))
  expect_identical(plan(rep(2:3, 50), list(point, unit, unit))$pooled_bounds, c(1, 2))
  # different supports are mixed per component, as before
  mixed <- plan(rep(2:3, 50), list(point, unit, other))
  expect_null(mixed$pooled_bounds)
  expect_identical(mixed$bounds, list(c(1, 2), c(1, 3)))
  # an inexact or missing component support leaves the pooled estimate as it was
  expect_null(plan(rep(2L, 100), list(point, posterior_support_attribute(c(1, 2), exact = FALSE))))
  expect_null(plan(rep(2L, 100), list(point, NULL)))
  expect_null(BayesTools:::.Savage_Dickey_BF.component_plan(slab, NULL, NULL))
  # a support that the draws leave is reported
  expect_match(plan(rep(2L, 100), list(point, posterior_support_attribute(c(1.5, 2))))$warning,
               "unusable for their posterior samples")
})

# ---------------------------------------------------------------------------- #
# Levels of ordered terms: the support of a level follows from its design
# ---------------------------------------------------------------------------- #

test_that("the support of an ordered level follows from the total and the allocation", {

  slab <- prior("uniform", list(1, 2))
  fit <- function(...){
    info <- JAGS_formula(~ f, "mu", .ordered_hypothesis_data(), list(intercept = N01, f = prior_ordered(...)))
    info$prior_list$mu_f
  }
  group <- function(prior, weights){
    columns <- paste0("mu_f[", seq_along(weights), "]")
    list(prior = prior, columns = columns, weights = weights[weights != 0], indices = which(weights != 0))
  }
  support <- function(prior, weights) BayesTools:::.posterior_support_group_ordered(group(prior, weights))

  # the top level of a cumulative contrast is the total: its support is exact
  dirichlet <- fit(slab)
  top <- support(dirichlet, c(1, 1))
  expect_identical(top$bounds, c(1, 2))
  expect_true(top$exact)
  # an inner level is the total times a share in (0, 1): a hull, not exact
  inner <- support(dirichlet, c(1, 0))
  expect_identical(inner$bounds, c(0, 2))
  expect_false(inner$exact)
  expect_identical(support(fit(prior("uniform", list(-1, 2))), c(1, 0))$bounds, c(-1, 2))
  expect_identical(support(fit(prior("normal", list(0, 1), list(0, Inf))), c(1, 0))$bounds, c(0, Inf))
  # a contrast of the levels is a hull of the weights
  contrast <- support(dirichlet, c(1, -1))
  expect_identical(contrast$bounds, c(-2, 2))
  expect_false(contrast$exact)
  # a fixed allocation fixes the share of every level
  fixed <- fit(slab, allocation = c(.25, .75))
  expect_identical(support(fixed, c(1, 0))$bounds, c(.25, .5))
  expect_true(support(fixed, c(1, 0))$exact)
  expect_identical(support(fixed, c(1, 1))$bounds, c(1, 2))
  # the support of a total that has none, or of a prior without bound metadata
  expect_null(BayesTools:::.posterior_support_group_ordered(
    list(prior = prior_ordered(slab), columns = c("a", "b"), weights = c(a = 1), indices = 1L)))

  # the totals of the slices of a crossed term are independent
  info <- JAGS_formula(~ f * g, "mu", .ordered_hypothesis_data(), list(intercept = N01, f = prior_ordered(slab),
    g = prior_factor("normal", list(0, 1), contrast = "treatment"), `f:g` = prior_ordered(slab)))
  crossed <- info$prior_list[["mu_f__xXx__g"]]
  crossed_support <- function(weights){
    BayesTools:::.posterior_support_group_ordered(list(prior = crossed, columns = paste0("c[", 1:4, "]"),
      weights = weights[weights != 0], indices = which(weights != 0)))
  }
  # the top level of one slice is that slice's total
  expect_identical(crossed_support(c(1, 1, 0, 0))$bounds, c(1, 2))
  expect_true(crossed_support(c(1, 1, 0, 0))$exact)
  # the top levels of both slices add
  expect_identical(crossed_support(c(1, 1, 1, 1))$bounds, c(2, 4))
  expect_true(crossed_support(c(1, 1, 1, 1))$exact)
  expect_false(crossed_support(c(1, 0, 1, 1))$exact)

  # the slices of a spike-and-slab total share its component: the sum of the
  # slices is exact for one active slice, and (conservatively) not for several
  spike_and_slab <- prior_spike_and_slab(prior("uniform", list(-1, 2)), prior("point", list(.5)))
  info <- JAGS_formula(~ f * g, "mu", .ordered_hypothesis_data(), list(intercept = N01, f = prior_ordered(spike_and_slab),
    g = prior_factor("normal", list(0, 1), contrast = "treatment"), `f:g` = prior_ordered(spike_and_slab)))
  shared <- info$prior_list[["mu_f__xXx__g"]]
  shared_support <- function(weights){
    BayesTools:::.posterior_support_group_ordered(list(prior = shared, columns = paste0("c[", 1:4, "]"),
      weights = weights[weights != 0], indices = which(weights != 0)))
  }
  expect_identical(shared_support(c(1, 1, 0, 0))$bounds, c(-1, 2))
  expect_true(shared_support(c(1, 1, 0, 0))$exact)
  expect_identical(shared_support(c(1, 1, 1, 1))$bounds, c(-2, 4))
  expect_false(shared_support(c(1, 1, 1, 1))$exact)
  # within a component of the total the slices are independent
  component <- shared
  component$total <- spike_and_slab[[which(attr(spike_and_slab, "components") == "alternative")]]
  expect_true(BayesTools:::.posterior_support_group_ordered(list(prior = component, columns = paste0("c[", 1:4, "]"),
    weights = c(1, 1, 1, 1), indices = 1:4))$exact)
})

test_that("the levels of an ordered term are reflected at the bounds of its total", {

  n <- 80000L
  slab <- prior("uniform", list(1, 2))

  # a slab-only total: the top level is the total, with its exact support
  fit <- .ordered_hypothesis_single(total = slab, n = n, seed = 719L)
  marginal <- marginal_posterior(as_mixed_posteriors(fit, "mu_f"), "mu_f", prior_samples = TRUE,
                                 use_formula = FALSE, n_samples = n)
  top <- BayesTools:::.posterior_support_get(marginal[["hi"]])
  expect_identical(top$bounds, c(1, 2))
  expect_true(top$exact)
  inner <- BayesTools:::.posterior_support_get(marginal[["mid"]])
  expect_identical(inner$bounds, c(0, 2))
  expect_false(inner$exact)
  for(value in c(1, 1.5, 2)){
    total <- .ordered_hypothesis_bf(fit, paste0("f{total} = ", value))$BF
    level <- .ordered_hypothesis_bf(fit, paste0("f[hi] = ", value))$BF
    # the posterior is the prior: the Bayes factor is 1 (it was 2 at the bounds)
    expect_equal(level, 1, tolerance = 0.08, info = value)
    expect_equal(level, total, tolerance = 1e-6, info = value)
  }
  # the inner levels have no exact bound and are tested as before
  for(value in c(.5, 1)){
    expect_equal(.ordered_hypothesis_bf(fit, paste0("f[mid] = ", value))$BF, 1, tolerance = 0.08, info = value)
  }

  # within the inclusion event of a spike-and-slab total, and without it: the
  # level of each component is reflected at the bounds of its total's component
  set.seed(404)
  fit <- .ordered_hypothesis_single(total = prior_spike_and_slab(slab, prior("point", list(.5))), n = n, seed = 719L)
  for(conditional in c(TRUE, FALSE)) for(value in c(1, 1.5, 2)){
    level <- .ordered_hypothesis_bf(fit, paste0("f[hi] = ", value), conditional = conditional)$BF
    expect_equal(level, 1, tolerance = 0.08, info = paste(conditional, value))
    expect_equal(level, .ordered_hypothesis_bf(fit, paste0("f{total} = ", value), conditional = conditional)$BF,
                 tolerance = 1e-6, info = paste(conditional, value))
  }
  # the inner levels have a hull support and are tested as before
  for(conditional in c(TRUE, FALSE)) for(value in c(.5, 1, 1.5)){
    expect_equal(.ordered_hypothesis_bf(fit, paste0("f[mid] = ", value), conditional = conditional)$BF, 1,
                 tolerance = 0.08, info = paste(conditional, value))
  }

  # several bounded alternatives: the level equals the total at the bounds of every slab
  disjoint <- prior_mixture(list(prior("point", list(0), prior_weights = 5),
    prior("uniform", list(1, 2), prior_weights = 1), prior("uniform", list(3, 4), prior_weights = 3)),
    is_null = c(TRUE, FALSE, FALSE))
  overlapping <- prior_mixture(list(prior("point", list(0), prior_weights = 2),
    prior("uniform", list(1, 3), prior_weights = 1), prior("uniform", list(2, 4), prior_weights = 3)),
    is_null = c(TRUE, FALSE, FALSE))
  for(total in list(disjoint = disjoint, overlapping = overlapping)){
    fit <- .ordered_hypothesis_single(total = total, n = n, seed = 719L)
    for(conditional in c(TRUE, FALSE)) for(value in 1:4){
      level <- .ordered_hypothesis_bf(fit, paste0("f[hi] = ", value), conditional = conditional)$BF
      expect_equal(level, 1, tolerance = 0.08, info = paste(conditional, value))
      expect_equal(level, .ordered_hypothesis_bf(fit, paste0("f{total} = ", value), conditional = conditional)$BF,
                   tolerance = 1e-6, info = paste(conditional, value))
    }
  }

  # unbounded totals are tested as before
  fit <- .ordered_hypothesis_single(total = N01, n = n, seed = 719L)
  for(value in c(-.5, .5)){
    level <- .ordered_hypothesis_bf(fit, paste0("f[hi] = ", value))$BF
    expect_equal(level, 1, tolerance = 0.08, info = value)
    expect_equal(level, .ordered_hypothesis_bf(fit, paste0("f{total} = ", value))$BF, tolerance = 1e-6, info = value)
  }

  # the levels of an ordinary factor with bounded coefficients are the
  # coefficients: exact bounds at both of them
  set.seed(405)
  info <- JAGS_formula(~ g, "mu", .ordered_hypothesis_data(), list(intercept = N01,
    g = prior_factor("uniform", list(1, 2), contrast = "treatment")))
  columns <- cbind(mu_intercept = stats::rnorm(n), `mu_g[1]` = stats::runif(n, 1, 2), `mu_g[2]` = stats::runif(n, 1, 2))
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(columns)), sample = n, summary.pars = list(mutate = NULL),
         monitor = colnames(columns)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list")     <- info$prior_list
  attr(fit, "formula_design") <- list(mu = info$formula_design)
  fit <- attach_test_parameter_map(fit)
  for(level in c("B", "C")) for(value in c(1, 2)){
    expect_equal(.ordered_hypothesis_bf(fit, paste0("g[", level, "] = ", value))$BF, 1, tolerance = 0.08,
                 info = paste(level, value))
  }
})

test_that("the levels of an ordered term with a mixture total carry the supports of its components", {

  n <- 4000L
  total <- prior_spike_and_slab(prior("uniform", list(1, 2)), prior("point", list(.5)))
  fit <- .ordered_hypothesis_single(total = total, n = n, seed = 719L)
  mixed <- as_mixed_posteriors(fit, "mu_f")
  marginal <- marginal_posterior(mixed, "mu_f", prior_samples = TRUE, use_formula = FALSE, n_samples = n)
  total_component <- BayesTools:::.bt_meta_get(mixed[["mu_f"]], "ordered_total_component")
  spike <- which(attr(total, "components") == "null")
  slab  <- which(attr(total, "components") == "alternative")

  for(level in c("mid", "hi")){
    components <- BayesTools:::.posterior_components_get(marginal[[level]])
    expect_s3_class(components, "BayesTools_posterior_components")
    # the component of a draw is that of the total
    expect_identical(as.integer(components$keys[components$index, "mu_f"]), as.integer(total_component))
    key <- function(component) match(component, components$keys[, "mu_f"])
    # the spike component is the point at zero, for every level
    expect_identical(components$supports[[key(spike)]]$bounds, c(0, 0))
    expect_true(components$supports[[key(spike)]]$exact)
    # the top level is the total of the slab (exact), an inner level a hull
    expect_identical(components$supports[[key(slab)]]$bounds, if(level == "hi") c(1, 2) else c(0, 2))
    expect_identical(components$supports[[key(slab)]]$exact, level == "hi")
  }

  # draws whose total component is unavailable stop instead of being pooled
  broken <- mixed
  broken[["mu_f"]] <- BayesTools:::.bt_meta_set(broken[["mu_f"]], "ordered_total_component", NULL)
  expect_error(marginal_posterior(broken, "mu_f", prior_samples = TRUE, use_formula = FALSE, n_samples = n),
               "Mixture component indices are unavailable")
})

test_that("the levels of a spike-and-slab factor are reflected at the bounds of the slab", {

  # the levels of an ordinary factor with a spike-and-slab prior are its
  # coefficients: the component supports give the reflection at both bounds
  n <- 80000L
  set.seed(406)
  prior_g <- prior_spike_and_slab(prior_factor("uniform", list(1, 2), contrast = "treatment"), prior("point", list(.5)))
  info <- JAGS_formula(~ g, "mu", .ordered_hypothesis_data(), list(intercept = N01, g = prior_g))
  indicator <- stats::rbinom(n, 1, .5)
  variable  <- cbind(stats::runif(n, 1, 2), stats::runif(n, 1, 2))
  columns <- cbind(mu_intercept = stats::rnorm(n), `mu_g[1]` = variable[, 1] * indicator,
                   `mu_g[2]` = variable[, 2] * indicator, mu_g_indicator = indicator, mu_g_inclusion = .5,
                   `mu_g_variable[1]` = variable[, 1], `mu_g_variable[2]` = variable[, 2])
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(columns)), sample = n, summary.pars = list(mutate = NULL),
         monitor = colnames(columns)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list")     <- info$prior_list
  attr(fit, "formula_design") <- list(mu = info$formula_design)
  fit <- attach_test_parameter_map(fit)
  for(conditional in c(TRUE, FALSE)) for(level in c("B", "C")) for(value in c(1, 1.5, 2)){
    expect_equal(.ordered_hypothesis_bf(fit, paste0("g[", level, "] = ", value), conditional = conditional)$BF, 1,
                 tolerance = 0.08, info = paste(conditional, level, value))
  }
})
