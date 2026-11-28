skip_if_not_test_profile("unit")

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

test_that("scalar SD marginals match the constructor and the fitted coordinate law", {

  simple <- prior("uniform", list(1, 3))
  spike <- prior_spike_and_slab(simple, prior("beta", list(1, 3)))
  mixture <- prior_mixture(list(prior("point", list(.5), prior_weights = 2), simple),
    components = c("fixed", "estimated"))
  for(contrast in c("treatment", "independent")){
    for(expected in list(simple, spike, mixture)){
      fit <- .declared_factor_sd_fit(expected, contrast)
      selection <- parameter_catalog_resolve(parameter_catalog(fit), "(mu) sd(f[b])")
      key <- selection$quantities$extraction_key[[1L]]
      scalar <- parameter_declared_prior(fit, selection)
      # The fitted prior additionally identifies its formula and random block.
      # The marginal must retain all distribution and mixture attributes.
      for(field in c("parameter", "random_sd", "random_factor", "random_grouping_factor",
                     "random_name", "random_structure")) attr(scalar, field) <- NULL
      expect_identical(scalar, expected)
      coordinate <- .bt_parameter_prior_density_coordinates(fit, selection$quantities,
        list(type = "coordinate", dependencies = key$source_parameter), n_grid = 512L,
        tail_prob = 1e-6)
      marginal <- parameter_prior_density(fit, selection, n_grid = 512L, tail_prob = 1e-6)
      expect_equal(marginal$points, coordinate$points)
      expect_equal(prior_density_ordinate(marginal, 2)$log_density,
        prior_density_ordinate(coordinate, 2)$log_density)
    }
  }
})

test_that("a two-column LKJ correlation exposes the prior of its fitted primitive", {

  data <- data.frame(x = c(1, 3, 4, 5, 1, 2), id = factor(rep(1:3, each = 2)))
  result <- JAGS_formula(~ 1 + x + us(1 + x | id), "mu", data,
    list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 1)),
      cor = prior_lkj(eta = 2))), formula_scale = TRUE)
  term <- result$formula_design$random_effects[[1L]]
  fit <- .prior_monitor_test_fit(result, c(JAGS_to_monitor(result$prior_list),
    term$correlation$primitive_names))
  fitted <- parameter_catalog_resolve(parameter_catalog(fit, transform_scaled = FALSE),
    "(mu) cor(intercept,x)")
  original <- parameter_catalog_resolve(parameter_catalog(fit), "(mu) cor(intercept,x)")
  expect_null(parameter_declared_prior(fit, fitted))
  expect_identical(parameter_declared_prior(fit, fitted, source = TRUE), prior("beta", list(2, 2)))
  expect_null(parameter_declared_prior(fit, original, source = TRUE))
})

test_that("a centered random slope leaves the original intercept SD's declared prior intact", {

  data <- data.frame(x = rep(c(-2, 2), 4), id = factor(rep(1:4, each = 2)))
  result <- JAGS_formula(~ 1 + x + us(1 + x | id), "mu", data,
    list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
    prior_random = prior_random(sd = prior("gamma", list(2, 3))), formula_scale = TRUE)
  term <- result$formula_design$random_effects[[1L]]
  fit <- .prior_monitor_test_fit(result, c(JAGS_to_monitor(result$prior_list),
    term$correlation$primitive_names))
  catalog <- parameter_catalog(fit)
  intercept <- parameter_catalog_resolve(catalog, "(mu) sd(intercept)")
  slope <- parameter_catalog_resolve(catalog, "(mu) sd(x)")
  expect_identical(intercept$quantities$extraction_key[[1]]$source_scale, 1)
  expect_identical(parameter_declared_prior(fit, intercept), result$prior_list[[term$sd_parameter_names[1]]])
  expect_null(parameter_declared_prior(fit, slope))
})

.declared_ordered_sd_fixture <- function(total, omit_indicator = FALSE, omit_total = FALSE){

  data <- data.frame(f = ordered(rep(c("low", "high"), 4), levels = c("low", "high")),
                     id = factor(rep(1:4, each = 2)))
  result <- JAGS_formula(~ 1 + (0 + f || id), "mu", data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(sd = prior_ordered(total, allocation = 1))))
  draws <- .generate_prior_sample_matrix(result$prior_list, 400, seed = 5)
  owner <- result$formula_design$random_effects[[1L]]$sd_parameter_names
  indicator <- paste0(prior_ordered_total_name(owner), "_indicator")
  if(omit_indicator) draws <- draws[, setdiff(colnames(draws), indicator), drop = FALSE]
  if(omit_total) draws <- draws[, setdiff(colnames(draws), prior_ordered_total_name(owner)), drop = FALSE]
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(draws)), result$prior_list,
    formula_design = list(mu = result$formula_design))
  list(fit = fit, draws = draws, owner = owner, indicator = indicator)
}

test_that("single-increment ordered SDs and variances use their total's point states", {

  totals <- list(
    spike = prior_spike_and_slab(prior("uniform", list(1, 3)), prior("beta", list(1, 3))),
    mixture = prior_mixture(list(prior("point", list(.5), prior_weights = 2),
      prior("uniform", list(1, 3)), prior("uniform", list(3, 5), prior_weights = 3)),
      is_null = c(TRUE, FALSE, FALSE)),
    fixed = prior("point", list(.5))
  )
  for(kind in names(totals)){
    fixture <- .declared_ordered_sd_fixture(totals[[kind]])
    fit <- fixture$fit
    draws <- fixture$draws
    owner <- fixture$owner
    indicator <- fixture$indicator
    active <- switch(kind, spike = draws[, indicator] == 1, mixture = draws[, indicator] != 1,
      fixed = rep(FALSE, nrow(draws)))
    atom <- if(kind == "spike") 0 else .5
    prior_atom_mass <- switch(kind, spike = .75, mixture = 2/6, fixed = 1)
    for(power in 1:2){
      label <- c("(mu) sd(f)", "(mu) var(f)")[power]
      selection <- parameter_catalog_resolve(parameter_catalog(fit), label)
      states <- parameter_gate_states(fit, selection)
      expect_true(states$known)
      expect_equal(states$atom, ifelse(active, NA_real_, atom^power))
      expect_equal(states$prior_atoms, data.frame(x = atom^power, mass = prior_atom_mass))
      law <- parameter_prior_density(fit, selection)
      expect_equal(law$points, data.frame(x = atom^power, p = prior_atom_mass))
      mixed <- parameter_mixed_posterior(fit, selection)
      atoms <- posterior_metadata(mixed, "atoms")
      expect_equal(as.numeric(atoms$locations), atom^power)
      expect_equal(atoms$mass, mean(!active))
      if(kind == "fixed"){
        expect_error(parameter_mixed_posterior(fit, selection, conditional = TRUE), "fixed")
        next
      }
      conditional <- parameter_mixed_posterior(fit, selection, conditional = TRUE)
      expect_equal(as.numeric(conditional), as.numeric(draws[active, owner])^power)
      expect_equal(posterior_metadata(conditional, "support")$bounds,
        c(1, if(kind == "mixture") 5 else 3)^power)
      conditional_law <- posterior_metadata(conditional, "prior_density")
      expect_equal(nrow(conditional_law$points), 0L)
      value <- 2^power
      height <- (if(kind == "mixture") .25 else 1) / 2
      jacobian <- power * 2^(power - 1L)
      expect_equal(exp(prior_density_ordinate(conditional_law, value)$log_density), height / jacobian)
      region <- hypothesis_BF(fit, paste(label, ">", value), seed = 1)
      prior_prob <- switch(kind, spike = .25/2, mixture = (3 + .5)/6)
      expect_equal(unname(attr(region, "raw_log_BF")),
        stats::qlogis(mean(draws[,owner] > 2)) - stats::qlogis(prior_prob), tolerance = 1e-8)
      point <- hypothesis_BF(fit, paste(label, "=", value), seed = 1)
      expect_true(is.finite(attr(point, "raw_log_BF")))
    }
  }
})

test_that("Bernoulli ordered SD totals declare their observed discrete branch and masses", {

  binary <- prior("bernoulli", list(.3))
  slab <- prior("uniform", list(1, 3))
  mixture <- prior_mixture(list(binary, slab), is_null = c(TRUE, FALSE))
  for(total in list(binary, mixture)){
    fixture <- .declared_ordered_sd_fixture(total)
    fit <- fixture$fit
    discrete <- if(is.prior.mixture(total)) fixture$draws[, fixture$indicator] == 1 else rep(TRUE, nrow(fixture$draws))
    weights <- if(is.prior.mixture(total)) c(.35, .15) else c(.7, .3)
    for(label in c("(mu) sd(f)", "(mu) var(f)")){
      selection <- parameter_catalog_resolve(parameter_catalog(fit), label)
      states <- parameter_gate_states(fit, selection)
      expect_true(states$known)
      expect_equal(states$atom, ifelse(discrete, fixture$draws[, fixture$owner], NA_real_))
      expect_equal(states$prior_atoms, data.frame(x = c(0, 1), mass = weights))
      atoms <- posterior_metadata(parameter_mixed_posterior(fit, selection), "atoms")
      expect_equal(atoms$mass, c(mean(states$atom %in% 0), mean(states$atom %in% 1)))
      if(is.prior.mixture(total)){
        conditional <- parameter_mixed_posterior(fit, selection, conditional = TRUE)
        expect_equal(nrow(posterior_metadata(conditional, "prior_density")$points), 0L)
        expect_equal(length(conditional), sum(!discrete))
      }else{
        expect_error(parameter_mixed_posterior(fit, selection, conditional = TRUE), "no inclusion gate")
      }
    }
  }
  # A binary alternative is atomic. The point-plan condition requires a
  # continuous alternative, as it does for a scalar point alternative.
  fit <- .declared_ordered_sd_fixture(prior_mixture(list(binary, slab), is_null = c(FALSE, TRUE)))$fit
  selection <- parameter_catalog_resolve(parameter_catalog(fit), "(mu) sd(f)")
  expect_error(parameter_mixed_posterior(fit, selection, conditional = TRUE), "not continuous")
  fit <- .declared_ordered_sd_fixture(mixture, omit_indicator = TRUE)$fit
  selection <- parameter_catalog_resolve(parameter_catalog(fit), "(mu) sd(f)")
  expect_false(parameter_gate_states(fit, selection)$known)
  expect_error(parameter_mixed_posterior(fit, selection, conditional = TRUE), "not monitored")
  for(total in list(binary, mixture)){
    fit <- .declared_ordered_sd_fixture(total, omit_total = TRUE)$fit
    selection <- parameter_catalog_resolve(parameter_catalog(fit), "(mu) sd(f)")
    expect_false(parameter_gate_states(fit, selection)$known)
    reason <- if(is.prior.mixture(total)) "ordered total coordinate 'mu__xREx__id_f_ordered_total' is not monitored" else "no inclusion gate"
    expect_error(parameter_mixed_posterior(fit, selection, conditional = TRUE), reason, fixed = TRUE)
  }
})

test_that("continuous ordered mixtures retain component metadata independently of gate states", {

  members <- list(prior("uniform", list(1, 2)), prior("uniform", list(3, 4)))
  for(omit in c(FALSE, TRUE)){
    fit <- .declared_ordered_sd_fixture(prior_mixture(members, components = c("a", "b")),
      omit_indicator = omit)$fit
    selected <- parameter_catalog_resolve(parameter_catalog(fit), "(mu) sd(f)")
    expect_null(parameter_gate_states(fit, selected))
    posterior <- parameter_mixed_posterior(fit, selected)
    expect_length(posterior_metadata(posterior, "atoms")$mass, 0)
    if(!omit) expect_length(.posterior_components_get(posterior)$supports, 2)

    fit <- .declared_ordered_sd_fixture(prior_mixture(members, is_null = c(TRUE, FALSE)),
      omit_indicator = omit)$fit
    selected <- parameter_catalog_resolve(parameter_catalog(fit), "(mu) sd(f)")
    states <- parameter_gate_states(fit, selected)
    expect_true(states$known)
    expect_true(all(states$continuous))
    if(omit){
      expect_null(states$event)
      expect_error(parameter_mixed_posterior(fit, selected, conditional = TRUE), "component indicator")
    }else{
      conditional <- parameter_mixed_posterior(fit, selected, conditional = TRUE)
      expect_equal(posterior_metadata(conditional, "support")$bounds, c(3, 4))
    }
  }
})
