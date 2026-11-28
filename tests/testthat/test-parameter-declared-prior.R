# Direct random factor SDs, using the same compiler and fitted catalog as JAGS.
.declared_factor_sd_fit <- function(sd_prior, contrast){

  data <- data.frame(f = factor(rep(c("a", "b", "c"), 4L)),
                     id = factor(rep(1:4, each = 3L)))
  result <- JAGS_formula(
    ~ 1 + diag(1 + f | id), "mu", data, list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(
      sd = prior("gamma", list(2, 1)), terms = list(f = sd_prior),
      contrasts = c(f = contrast)
    ))
  )
  draws <- .generate_prior_sample_matrix(result$prior_list, 300, seed = 34)
  .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)),
    prior_list = result$prior_list, formula_design = list(mu = result$formula_design))
}

test_that("individual random factor SDs expose their scalar prior and variance law", {

  for(contrast in c("treatment", "independent")){
    sd_prior <- prior("uniform", list(1, 3))
    fit <- .declared_factor_sd_fit(sd_prior, contrast)
    catalog <- parameter_catalog(fit)
    sd <- parameter_catalog_resolve(catalog, "(mu) sd(f[b])")
    variance <- parameter_catalog_resolve(catalog, "(mu) var(f[b])")

    declared <- parameter_declared_prior(fit, sd)
    expect_false(is.prior.factor(declared))
    expect_equal(declared$parameters, sd_prior$parameters)
    expect_equal(lpdf(declared, 2), log(.5))
    expect_null(parameter_declared_prior(fit, variance))
    expect_identical(parameter_declared_prior(fit, variance, source = TRUE), declared)

    expect_equal(exp(prior_density_ordinate(parameter_prior_density(fit, sd), 2)$log_density), .5)
    expect_equal(exp(prior_density_ordinate(parameter_prior_density(fit, variance), 4)$log_density), .125)
    expect_equal(sd$quantities$support[[1L]]$bounds, c(1, 3))
    expect_equal(variance$quantities$support[[1L]]$bounds, c(1, 9))
  }
})

test_that("a factor SD keeps the declared spike-and-slab probability and conditional law", {

  inclusion <- prior("beta", list(1, 3))
  slab <- prior("uniform", list(1, 3))
  fit <- .declared_factor_sd_fit(prior_spike_and_slab(slab, inclusion), "independent")
  catalog <- parameter_catalog(fit)
  selection <- parameter_catalog_resolve(catalog, "(mu) sd(f[b])")
  scalar <- parameter_declared_prior(fit, selection)
  expect_true(is.prior.spike_and_slab(scalar))
  expect_false(.bt_prior_is_factor_family(scalar))
  expect_identical(.get_spike_and_slab_inclusion(scalar), inclusion)
  expect_identical(attr(scalar, "components"), c("alternative", "null"))
  law <- parameter_prior_density(fit, selection)
  expect_equal(law$points$x, 0)
  expect_equal(law$points$p, .75)
  expect_equal(exp(prior_density_ordinate(law, 2)$log_density), .125)

  conditional <- parameter_mixed_posterior(fit, selection, conditional = TRUE)
  conditional_law <- posterior_metadata(conditional, "prior_density")
  expect_equal(nrow(conditional_law$points), 0L)
  expect_equal(exp(prior_density_ordinate(conditional_law, 2)$log_density), .5)
  expect_equal(posterior_metadata(conditional, "support")$bounds, c(1, 3))
})

test_that("mixture members and labels remain attached to the selected SD marginal", {

  mixture <- prior_mixture(
    list(prior("point", list(.5), prior_weights = 2),
         prior("uniform", list(1, 3), prior_weights = 3)),
    components = c("fixed", "estimated")
  )
  fit <- .declared_factor_sd_fit(mixture, "treatment")
  selection <- parameter_catalog_resolve(parameter_catalog(fit), "(mu) sd(f[c])")
  scalar <- parameter_declared_prior(fit, selection)
  expect_true(is.prior.mixture(scalar))
  expect_false(.bt_prior_is_factor_family(scalar))
  expect_identical(attr(scalar, "components"), c("fixed", "estimated"))
  expect_equal(attr(scalar, "prior_weights"), c(2, 3))
  law <- parameter_prior_density(fit, selection)
  expect_equal(law$points$x, .5)
  expect_equal(law$points$p, .4)
  expect_equal(exp(prior_density_ordinate(law, 2)$log_density), .3)
})

test_that("one-column ordered and point-vector SDs keep their scalar laws", {

  data <- data.frame(f = ordered(rep(c("low", "high"), 4), levels = c("low", "high")),
                     id = factor(rep(1:4, each = 2L)))
  total <- prior("uniform", list(1, 3))
  for(ordered in c(TRUE, FALSE)){
    if(!ordered) data$f <- factor(data$f, ordered = FALSE)
    sd <- if(ordered) prior_ordered(total, allocation = 1) else prior("mpoint", list(.5, 1))
    result <- JAGS_formula(~ 1 + (0 + f || id), "mu", data,
      list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = sd)))
    draws <- .generate_prior_sample_matrix(result$prior_list, 12, seed = 4)
    fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), result$prior_list,
      formula_design = list(mu = result$formula_design))
    catalog <- parameter_catalog(fit)
    label <- if(ordered) "f" else "f[high]"
    sd <- parameter_catalog_resolve(catalog, paste0("(mu) sd(", label, ")"))
    variance <- parameter_catalog_resolve(catalog, paste0("(mu) var(", label, ")"))
    prior_sd <- parameter_declared_prior(fit, sd)
    law <- parameter_prior_density(fit, variance)
    if(ordered){
      expect_identical(prior_sd, total)
      expect_equal(exp(prior_density_ordinate(parameter_prior_density(fit, sd), 2)$log_density), .5)
      expect_equal(exp(prior_density_ordinate(law, 4)$log_density), .125)
    }else{
      expect_true(is.prior.point(prior_sd))
      expect_false(is.prior.vector(prior_sd))
      expect_equal(prior_sd$parameters$location, .5)
      expect_equal(law$points, data.frame(x = .25, p = 1))
    }
  }
})
