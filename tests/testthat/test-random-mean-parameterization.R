skip_if_not_test_profile("unit")

.mean_parameterization_compile <- function(mode = "mean_centered",
                                           effect = prior("t", list(0, .7, 3)),
                                           sd = prior("cauchy", list(0, .4), list(0, Inf)),
                                           formula = ~ 1 + x + id(1 | study),
                                           marginalized = NULL, random = NULL){

  data <- data.frame(study = factor(rep(c("a", "b"), each = 4L)),
                     site = factor(rep(c("one", "two"), 4L)), x = rep(-1:2, 2L))
  if(is.null(random)){
    random <- prior_random(study = random_block(sd = sd, parameterization = mode,
      monitor = random_monitor(latent = TRUE, coefficients = TRUE)))
  }
  JAGS_formula(formula, "mu", data,
    prior_list = list(intercept = effect, x = prior("normal", list(0, .5))),
    prior_random = random, random_effects_compile = random_effects_compile(marginalized))
}

test_that("mean centering retains arbitrary effect priors and residual moderators", {

  for(case in c("student", "point", "log")){
    formula <- ~ 1 + x + id(1 | study)
    effect <- switch(case, student = prior("t", list(0, .7, 3)),
      point = prior("point", list(0)), log = prior("gamma", list(2, 1)))
    if(case == "log") attr(formula, "log(intercept)") <- TRUE
    translated <- .mean_parameterization_compile(effect = effect, formula = formula)
    centered <- .mean_parameterization_compile("centered", effect = effect, formula = formula)
    term <- translated$formula_design$random_effects[[1L]]
    plan <- term$mean_translation
    contribution <- if(case == "log") "log(mu_intercept)" else "mu_intercept"
    expect_identical(plan$fixed_intercept, "mu_intercept")
    expect_identical(plan$fixed_intercept_expression, contribution)
    expect_match(translated$formula_syntax,
      paste0(plan$location_name, "[g,1] ~ dnorm(", contribution, ","), fixed = TRUE)
    expect_match(translated$formula_syntax,
      paste0(plan$coefficient_name, "[g,1] <- ", plan$location_name, "[g,1] - ", contribution),
      fixed = TRUE)
    expect_match(translated$formula_syntax,
      paste0("mu[i] = mu_x * mu_data_x[i] + ", plan$location_name,
        "[", plan$group_map_name, "[i],1]"), fixed = TRUE)
    prior_fields <- function(priors) lapply(priors, function(prior){

      prior[c("distribution", "parameters", "truncation", "prior_weights")]
    })
    expect_equal(prior_fields(translated$prior_list), prior_fields(centered$prior_list))
    expect_identical(translated$add_parameters, centered$add_parameters)
    expect_identical(centered$formula_design$random_effects[[1L]]$parameterization_resolved,
      "centered")
  }
})

test_that("mean centering keeps semantic deviations and the canonical bridge measure", {

  translated <- .mean_parameterization_compile()
  centered <- .mean_parameterization_compile("centered")
  design <- translated$formula_design
  term <- design$random_effects[[1L]]
  z_names <- as.vector(.bt_random_effect_latent_names(term, term$n_groups, term$n_columns))
  coef_names <- as.vector(.bt_random_effect_coefficient_names(term, term$n_groups, term$n_columns))
  location_names <- paste0(term$mean_translation$location_name, "[", seq_len(term$n_groups), ",1]")
  map <- .bt_build_parameter_map(c(names(translated$prior_list), z_names, coef_names, location_names),
    prior_list = translated$prior_list, formula_design = list(mu = design))
  location <- match(location_names, map$coordinates$coordinate_name)
  expect_false(anyNA(location))
  expect_true(all(map$coordinates$internal[location]))
  expect_identical(map$coordinates$role[location], rep("random_mean_coordinate", 2L))
  expect_false(any(grepl("_xRE_MEANx", map$quantities$canonical_name, fixed = TRUE)))
  expect_false(any(grepl("_xRE_MEANx", map$aliases$alias, fixed = TRUE)))
  translated_bridge <- .bt_JAGS_formula_random_bridge_parameters(list(mu = design))
  expect_identical(translated_bridge,
    .bt_JAGS_formula_random_bridge_parameters(list(mu = centered$formula_design)))
  expect_identical(translated_bridge$parameters, z_names)

  # The actual compiler metadata, with deterministic parameter values rather
  # than a fabricated fit object, must reconstruct the original deviations.
  values <- matrix(c(.2, -.1, .6, -.3, .4), 1L,
    dimnames = list(NULL, c("mu_intercept", "mu_x", term$sd_parameter_names, z_names)))
  deviation <- .bt_try_random_effect_contribution_from_latent(term,
    term$model_matrix, term$group_map, values, translated$prior_list)
  expect_equal(as.numeric(deviation), .6 * c(-.3, .4)[term$group_map], tolerance = 1e-14)

  changed <- design
  changed$random_effects[[1L]]$mean_translation$fixed_intercept_expression <- "mu_x"
  condition <- tryCatch(.bt_JAGS_bridge_validate_formula_random_design("mu", design, changed),
    error = identity)
  expect_identical(conditionMessage(condition), paste0(
    "JAGS_bridgesampling() rebuilt formula random-effect design does not match the fitted design for parameter 'mu', ",
    "block 'study': random-effect parameterization differs. ",
    "Supply the same formula, data, scaling, and prior_random() metadata used to fit the model."))
  expect_null(conditionCall(condition))
})

test_that("mean centering requires one eligible explicitly selected block", {

  expect_error(prior_random(parameterization = "mean_centered"),
    "Mean-centered parameterization requires an explicitly named 'random_block()'.", fixed = TRUE)
  expect_error(.mean_parameterization_compile(sd = prior("point", list(0))),
    "Mean-centered parameterization is unavailable for random-effect block 'study': SD prior with an atom at zero.",
    fixed = TRUE)
  expect_error(.mean_parameterization_compile(marginalized = "study"),
    "Mean-centered parameterization is unavailable for marginalized random-effect blocks.", fixed = TRUE)
  expect_error(.mean_parameterization_compile(formula = ~ 1 + x + us(1 + x | study)),
    "Mean-centered parameterization is available only for a scalar random intercept.", fixed = TRUE)
  sd <- prior("normal", list(0, 1), list(0, Inf))
  expect_error(.mean_parameterization_compile(formula = ~ 1 + x + id(1 | study) + id(1 | site),
    random = prior_random(study = random_block(sd = sd, parameterization = "mean_centered"),
      site = random_block(sd = sd, parameterization = "mean_centered"))),
    "Mean-centered parameterization requires one explicitly selected random-effect block per formula.", fixed = TRUE)
})
