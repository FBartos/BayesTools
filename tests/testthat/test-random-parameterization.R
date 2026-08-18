skip_if_not_test_profile("unit")

.parameterization_sd_prior <- function(){

  prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
}

.parameterization_compile <- function(parameterization, data,
                                      formula = ~ 1 + random(
                                        1 | study,
                                        name = "study",
                                        covariance = "diag"
                                      ),
                                      sd = .parameterization_sd_prior()){

  JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      study = random_block(
        sd = sd,
        parameterization = parameterization
      )
    )
  )
}

test_that("random-effect parameterization validates public values", {

  expect_identical(prior_random()$parameterization, "noncentered")
  expect_identical(
    prior_random(parameterization = "centered")$parameterization,
    "centered"
  )
  expect_identical(
    prior_random(parameterization = "auto")$parameterization,
    "auto"
  )
  expect_null(random_block()$parameterization)
  expect_identical(
    random_block(parameterization = "centered")$parameterization,
    "centered"
  )

  expect_error(
    prior_random(parameterization = NULL),
    "parameterization"
  )
  expect_error(
    prior_random(parameterization = "adaptive"),
    "parameterization"
  )
  expect_error(
    prior_random(parameterization = c("centered", "auto")),
    "parameterization"
  )
  expect_error(
    random_block(parameterization = NA_character_),
    "parameterization"
  )
})

test_that("block parameterization inherits and overrides top-level policy", {

  specification <- prior_random(
    parameterization = "auto",
    study = random_block(parameterization = "centered"),
    site = random_block()
  )

  expect_identical(
    .bt_random_prior_for_block(specification, "study")$parameterization,
    "centered"
  )
  expect_identical(
    .bt_random_prior_for_block(specification, "site")$parameterization,
    "auto"
  )
  expect_identical(
    .bt_random_prior_for_block(specification, "unlisted")$parameterization,
    "auto"
  )
  expect_identical(specification$parameterization, "auto")
  expect_null(specification$blocks$site$parameterization)
})

test_that("mutated random parameterization metadata are rejected", {

  malformed <- prior_random(study = random_block())
  malformed$parameterization <- "adaptive"
  expect_error(.bt_check_prior_random(malformed), "parameterization")

  malformed <- prior_random(study = random_block())
  malformed$blocks$study$parameterization <- "adaptive"
  expect_error(.bt_check_prior_random(malformed), "parameterization")
})

test_that("random parameterization printing omits only the default", {

  expect_equal(
    utils::capture.output(print(prior_random())),
    "no random-effect priors specified"
  )
  expect_equal(
    utils::capture.output(print(prior_random(parameterization = "auto"))),
    c("settings", "  parameterization: auto")
  )
  expect_equal(
    utils::capture.output(print(random_block(parameterization = "centered"))),
    c("block", "  parameterization: centered")
  )

  specification <- prior_random(
    parameterization = "centered",
    study = random_block(parameterization = "noncentered")
  )
  expect_equal(utils::capture.output(print(specification)), c(
    "block: study",
    "  parameterization: noncentered",
    "settings",
    "  parameterization: centered"
  ))
})

test_that("centered and noncentered compilers preserve the same public scale", {

  data <- data.frame(
    study = factor(rep(c("a", "b"), each = 5L))
  )
  noncentered <- .parameterization_compile("noncentered", data)
  centered    <- .parameterization_compile("centered", data)
  noncentered_term <- noncentered$formula_design$random_effects[[1L]]
  centered_term    <- centered$formula_design$random_effects[[1L]]

  expect_identical(noncentered_term$parameterization_resolved, "noncentered")
  expect_identical(centered_term$parameterization_resolved, "centered")
  expect_equal(noncentered_term$sd_parameter_names,
               centered_term$sd_parameter_names)
  expect_equal(names(noncentered$prior_list), names(centered$prior_list))
  expect_match(
    noncentered$formula_syntax,
    "mu__xREx__study_xRE_Zx[i,j] ~ dnorm(0, 1)",
    fixed = TRUE
  )
  expect_match(
    centered$formula_syntax,
    "mu__xREx__study_xRE_COEFx[g,i] ~ dnorm",
    fixed = TRUE
  )
  expect_match(
    centered$formula_syntax,
    "mu__xREx__study_xRE_Zx[g,i] <- mu__xREx__study_xRE_COEFx[g,i] /",
    fixed = TRUE
  )
})

test_that("auto parameterization uses deterministic design diagnostics", {

  informative <- data.frame(
    study = factor(rep(c("a", "b"), each = 5L))
  )
  weak <- data.frame(
    study = factor(rep(c("a", "b"), each = 4L))
  )
  informative_result <- .parameterization_compile("auto", informative)
  weak_result        <- .parameterization_compile("auto", weak)
  informative_term <- informative_result$formula_design$random_effects[[1L]]
  weak_term        <- weak_result$formula_design$random_effects[[1L]]

  expect_identical(informative_term$parameterization_resolved, "centered")
  expect_identical(
    informative_term$parameterization_policy,
    BayesTools:::.bt_random_effect_auto_parameterization_policy()
  )
  expect_match(
    informative_term$parameterization_reason,
    "replication and conditioning"
  )
  expect_identical(weak_term$parameterization_resolved, "noncentered")
  expect_identical(
    weak_term$parameterization_policy,
    BayesTools:::.bt_random_effect_auto_parameterization_policy()
  )
  expect_identical(
    weak_term$parameterization_reason,
    "insufficient within-group information"
  )

  rank_deficient <- cbind(1, 1)
  expect_false(BayesTools:::.bt_random_effect_auto_centered_design(
    model_matrix = rank_deficient,
    group_map = rep(1L, nrow(rank_deficient))
  )$ok)

  unused_level <- data.frame(
    study = factor(rep(c("a", "b"), each = 5L), levels = c("a", "b", "c"))
  )
  unused_result <- .parameterization_compile("auto", unused_level)
  unused_term   <- unused_result$formula_design$random_effects[[1L]]
  expect_identical(unused_term$parameterization_resolved, "noncentered")
  expect_identical(
    unused_term$parameterization_reason,
    "one or more grouping levels are unobserved"
  )

  overflow_safe <- BayesTools:::.bt_random_effect_auto_centered_design(
    model_matrix = matrix(1e200, nrow = 5L, ncol = 1L),
    group_map = rep(1L, 5L)
  )
  expect_true(overflow_safe$ok)
})

test_that("centered parameterization rejects degenerate scale contracts", {

  data <- data.frame(study = factor(rep(c("a", "b"), each = 5L)))
  expect_error(
    .parameterization_compile(
      "centered",
      data,
      sd = prior("point", list(location = 0))
    ),
    "atom at zero",
    fixed = TRUE
  )

  block <- random_block(
    sd = .parameterization_sd_prior(),
    parameterization = "centered"
  )
  gated_binding <- BayesTools:::.bt_random_sd_binding(
    source = BayesTools:::.bt_prior_owned_sd_source("sd"),
    application = "block",
    factors = list(BayesTools:::.bt_random_variance_allocation_factor(
      weight_name = "allocation_weight",
      index = 1L,
      scale = "total_variance",
      n_targets = 2L,
      inclusion_name = "include_study"
    ))
  )
  eligibility <- BayesTools:::.bt_random_effect_centered_eligibility(
    block_prior = block,
    prior_list = list(sd = .parameterization_sd_prior()),
    sd_binding = gated_binding,
    row_indexed_external_sd = FALSE
  )
  expect_false(eligibility$ok)
  expect_identical(eligibility$reason, "variance-allocation inclusion gate")

  block$sd_source <- list(name = "external_sd")
  eligibility <- BayesTools:::.bt_random_effect_centered_eligibility(
    block_prior = block,
    prior_list = list(),
    sd_binding = NULL,
    row_indexed_external_sd = FALSE
  )
  expect_false(eligibility$ok)
  expect_match(eligibility$reason, "external SD source", fixed = TRUE)

  expect_true(BayesTools:::.bt_random_effect_prior_has_zero_atom(
    prior_spike_and_slab(.parameterization_sd_prior())
  ))
  expect_true(BayesTools:::.bt_random_effect_prior_has_zero_atom(
    prior_ordered(
      .parameterization_sd_prior(),
      allocation = c(1, 0)
    )
  ))

  external_block <- random_block(
    sd_source = random_sd_source("tau"),
    parameterization = "centered"
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | study),
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(study = external_block)
    ),
    "external SD source",
    fixed = TRUE
  )
})

test_that("centered scalar structures materialize covariance only internally", {

  data <- data.frame(
    study = factor(rep(c("a", "b"), each = 3L)),
    index = factor(rep(c("i1", "i2", "i3"), 2L)),
    time = rep(c(0, 0.5, 2), 2L)
  )
  formulas <- list(
    cs = ~ 1 + cs(index | study),
    hcs = ~ 1 + hcs(index | study),
    ar1 = ~ 1 + ar1(index | study),
    har = ~ 1 + har(index | study),
    car = ~ 1 + car(time | study)
  )

  for(structure in names(formulas)){
    result <- JAGS_formula(
      formula = formulas[[structure]],
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        study = random_block(
          sd = if(identical(structure, "car")){
            prior("point", list(location = 1))
          }else{
            .parameterization_sd_prior()
          },
          cor = prior("normal", list(0, 0.5)),
          monitor = random_monitor(correlation = FALSE),
          parameterization = "centered"
        )
      )
    )
    term <- result$formula_design$random_effects[[1L]]

    expect_identical(term$parameterization_resolved, "centered",
                     info = structure)
    if(identical(structure, "car")){
      expect_false(grepl("~ dmnorm.vcov", result$formula_syntax, fixed = TRUE),
                   info = structure)
      expect_false(grepl("_xRE_COVx", result$formula_syntax, fixed = TRUE),
                   info = structure)
      expect_false(grepl("_xRE_CORx_R", result$formula_syntax, fixed = TRUE),
                   info = structure)
      expect_match(
        result$formula_syntax,
        paste0(
          "mu__xREx__study_xRE_CAR_INNOV_VARx[2] <- ",
          "pexp(-2 * mu__xREx__study_xRE_CAR_LOG_PHIX[2], 1)"
        ),
        fixed = TRUE,
        info = structure
      )
      expect_match(
        result$formula_syntax,
        paste0(
          "mu__xREx__study_xRE_COEFx[g,i] ~ dnorm(",
          "mu__xREx__study_xRE_STDx[i] * ",
          "mu__xREx__study_xRE_CAR_PHIX[i]"
        ),
        fixed = TRUE,
        info = structure
      )
      expect_false(
        grepl(
          "pow(mu__xREx__study_rho",
          result$formula_syntax,
          fixed = TRUE
        ),
        info = structure
      )
    }else{
      expect_match(result$formula_syntax, "~ dmnorm.vcov", fixed = TRUE,
                   info = structure)
    }
    expect_false(any(grepl("_xRE_CORx_L", result$add_parameters, fixed = TRUE)),
                 info = structure)
    expect_false(any(grepl("_xRE_CORx_R", result$add_parameters, fixed = TRUE)),
                 info = structure)
    expect_true(term$correlation$rho_name %in% result$add_parameters,
                info = structure)
  }
})

test_that("marginalized random effects record backend-independent resolution", {

  block <- random_block(parameterization = "centered")
  resolved <- BayesTools:::.bt_random_effect_resolve_parameterization(
    block_prior = block,
    prior_list = list(),
    sd_binding = NULL,
    row_indexed_external_sd = FALSE,
    model_matrix = matrix(1, nrow = 2L),
    group_map = c(1L, 2L),
    compile_mode = "marginalized",
    block_name = "study"
  )

  expect_identical(resolved$requested, "centered")
  expect_identical(resolved$resolved, "marginalized")
})
