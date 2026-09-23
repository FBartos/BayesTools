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

test_that("block SD overrides replace SDs inherited through the other slot", {

  data <- data.frame(
    g = factor(rep(1:5, each = 4L)),
    h = factor(rep(1:4, 5L))
  )
  top_sd   <- .parameterization_sd_prior()
  block_sd <- prior("gamma", list(2, 2))
  sd_distributions <- function(specification){
    result <- JAGS_formula(
      formula = ~ 1 + (1 | g) + (1 | h),
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = specification
    )
    terms <- result$formula_design$random_effects
    stats::setNames(
      vapply(terms, function(term){
        result$prior_list[[term$sd_parameter_names[[1L]]]]$distribution
      }, character(1)),
      vapply(terms, `[[`, character(1), "block_name")
    )
  }

  covariance_override <- prior_random(
    sd = top_sd,
    g = random_block(covariance = random_covariance(sd = block_sd))
  )
  resolved <- .bt_random_prior_for_block(covariance_override, "g")
  expect_null(resolved$sd)
  expect_identical(resolved$covariance$sd, block_sd)
  expect_identical(
    sd_distributions(covariance_override),
    c(g = "gamma", h = "normal")
  )

  sd_override <- prior_random(
    covariance = random_covariance(sd = block_sd),
    g = random_block(sd = top_sd)
  )
  resolved <- .bt_random_prior_for_block(sd_override, "g")
  expect_identical(resolved$sd, top_sd)
  expect_null(resolved$covariance$sd)
  expect_identical(
    sd_distributions(sd_override),
    c(g = "normal", h = "gamma")
  )

  expect_error(
    sd_distributions(prior_random(
      sd = top_sd,
      g = random_block(sd = top_sd, covariance = random_covariance(sd = block_sd))
    )),
    "SD prior was supplied both",
    fixed = TRUE
  )
})

test_that("top-level correlation priors apply only to correlated blocks", {

  sd_prior <- .parameterization_sd_prior()
  lkj <- prior_lkj(eta = 2)
  rho <- prior("normal", list(0, 0.5))

  inherited <- .bt_random_prior_for_block(
    prior_random(sd = sd_prior, cor = lkj),
    "g1"
  )
  expect_true(.bt_random_block_has_inherited_correlation(inherited))
  for(case in list(list("us", 1L), list("diag", 3L), list("id", 2L))){
    resolved <- .bt_random_block_resolve_correlation(
      inherited, case[[1L]], case[[2L]]
    )
    expect_null(resolved$covariance$cor, info = case[[1L]])
    expect_false(.bt_random_block_has_inherited_correlation(resolved))
    expect_silent(.bt_validate_random_block_for_structure(
      resolved, case[[1L]], "g1"
    ))
  }
  resolved <- .bt_random_block_resolve_correlation(inherited, "us", 2L)
  expect_identical(resolved$covariance$cor, lkj)
  expect_false(.bt_random_block_has_inherited_correlation(resolved))
  expect_error(
    .bt_random_block_resolve_correlation(inherited, "cs", 3L, "g2"),
    "random-effect block 'g2' supplies an LKJ correlation prior, but structure 'cs' uses a scalar correlation prior",
    fixed = TRUE
  )

  scalar <- .bt_random_prior_for_block(
    prior_random(sd = sd_prior, cor = rho),
    "g1"
  )
  expect_silent(.bt_validate_random_block_for_structure(scalar, "us", "g1"))
  expect_null(.bt_random_block_resolve_correlation(scalar, "us", 1L)$covariance$cor)
  expect_null(.bt_random_block_resolve_correlation(scalar, "cs", 1L)$covariance$cor)
  expect_identical(
    .bt_random_block_resolve_correlation(scalar, "cs", 3L)$covariance$cor,
    rho
  )
  expect_error(
    .bt_random_block_resolve_correlation(scalar, "us", 2L, "g1"),
    "random-effect block 'g1' supplies a scalar correlation prior, but structure 'us' uses an LKJ correlation prior",
    fixed = TRUE
  )

  # A block-explicit correlation prior is not a default and is kept, so it is
  # still rejected where the block has no correlation parameter.
  explicit <- .bt_random_prior_for_block(
    prior_random(sd = sd_prior, g1 = random_block(cor = lkj)),
    "g1"
  )
  expect_false(.bt_random_block_has_inherited_correlation(explicit))
  expect_identical(
    .bt_random_block_resolve_correlation(explicit, "us", 1L)$covariance$cor,
    lkj
  )
  expect_error(
    .bt_validate_random_block_for_structure(explicit, "diag", "g1"),
    "structure 'diag' has no correlation parameter",
    fixed = TRUE
  )
})

test_that("block overrides can remove an inherited correlation prior", {

  set.seed(1)
  data <- data.frame(
    g1 = factor(rep(1:5, each = 6L)),
    g2 = factor(rep(1:6, 5L)),
    t  = factor(rep(1:3, 10L)),
    x  = stats::rnorm(30L)
  )
  sd_prior <- .parameterization_sd_prior()
  compile <- function(formula, specification){
    JAGS_formula(
      formula = formula,
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = specification
    )
  }
  cleared <- prior_random(
    sd = sd_prior,
    cor = prior_lkj(eta = 2),
    g1 = random_block(covariance = random_covariance(cor = NULL))
  )
  expect_null(.bt_random_prior_for_block(cleared, "g1")$covariance$cor)
  result <- compile(~ 1 + (1 | g1) + (1 + x | g2), cleared)
  terms <- result$formula_design$random_effects
  names(terms) <- vapply(terms, `[[`, character(1), "block_name")
  expect_null(terms$g1$correlation)
  expect_identical(terms$g2$correlation$eta, 2)

  # A top-level cor is not applied to a diag() block.
  result <- compile(
    ~ 1 + diag(1 | g1) + (1 + x | g2),
    prior_random(sd = sd_prior, cor = prior_lkj(eta = 2))
  )
  terms <- result$formula_design$random_effects
  names(terms) <- vapply(terms, `[[`, character(1), "block_name")
  expect_null(terms$g1$correlation)
  expect_identical(terms$g2$correlation$eta, 2)

  # Removing the inherited scalar prior also removes its inherited scale, so
  # the structure default (uniform on the raw correlation) applies.
  result <- compile(
    ~ 1 + cs(t | g2),
    prior_random(
      sd = sd_prior,
      covariance = random_covariance(
        cor = prior("normal", list(0, 0.5)),
        cor_scale = "logit"
      ),
      g2 = random_block(covariance = random_covariance(cor = NULL))
    )
  )
  rho_priors <- result$prior_list[grepl("_rho", names(result$prior_list))]
  expect_identical(names(rho_priors), "mu__xREx__g2_rho")
  expect_identical(rho_priors[[1L]]$distribution, "uniform")
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

test_that("centered eligibility includes structure-specific compiler contracts", {

  # Five groups with six replicates of four time indicators: the design
  # diagnostics alone favor centering.
  n_groups <- 5L
  car_matrix <- do.call(rbind, rep(list(diag(4L)[rep(1:4, times = 6L), ]), n_groups))
  car_groups <- rep(seq_len(n_groups), each = 24L)
  slope_matrix <- cbind(1, rep(c(-1, 1), length.out = 24L))
  slope_groups <- rep(1:3, each = 8L)
  resolve <- function(parameterization, prior_list, structure,
                      group_covariance = NULL, model_matrix = car_matrix,
                      group_map = car_groups){
    .bt_random_effect_resolve_parameterization(
      block_prior = random_block(parameterization = parameterization),
      prior_list = prior_list,
      sd_binding = NULL,
      row_indexed_external_sd = FALSE,
      model_matrix = model_matrix,
      group_map = group_map,
      compile_mode = "sampled",
      block_name = "id",
      structure = structure,
      group_covariance = group_covariance
    )
  }
  # The reasons mirror the centered CAR compiler checks of the SD support.
  zero_reason <- paste0(
    "centered CAR SD prior has an unrepresentable initial JAGS precision ",
    "at the lower centered SD support 0e+00; its support must be bounded ",
    "away from zero and infinity"
  )
  unbounded <- list(
    half_normal = list(sd = .parameterization_sd_prior()),
    gamma = list(sd = prior("gamma", list(2, 2))),
    tiny_point = list(sd = prior("point", list(location = 1e-200))),
    allocation = list()
  )
  car_reasons <- list(
    half_normal = zero_reason,
    gamma = zero_reason,
    tiny_point = paste0(
      "centered CAR SD prior has an unrepresentable initial JAGS precision ",
      "at the lower centered SD support ",
      format(1e-200, digits = 17, scientific = TRUE),
      "; its support must be bounded away from zero and infinity"
    ),
    allocation = paste0(
      "centered CAR requires one prior-owned block SD prior, which a ",
      "variance allocation does not provide"
    )
  )
  for(case in names(unbounded)){
    auto <- resolve("auto", unbounded[[case]], "car")
    expect_identical(auto$resolved, "noncentered", info = case)
    expect_identical(auto$reason, car_reasons[[case]], info = case)
    expect_error(
      resolve("centered", unbounded[[case]], "car"),
      paste0(
        "Centered parameterization is not available for random-effect ",
        "block 'id': ", car_reasons[[case]], "."
      ),
      fixed = TRUE
    )
  }
  bounded <- list(sd = prior("uniform", list(0.1, 2)))
  expect_identical(resolve("auto", bounded, "car")$resolved, "centered")
  expect_identical(resolve("centered", bounded, "car")$resolved, "centered")
  expect_identical(
    resolve("centered", list(sd = prior("point", list(location = 1))), "car")$resolved,
    "centered"
  )
  expect_identical(
    resolve("auto", unbounded$half_normal, "cs")$resolved,
    "centered"
  )

  kernel <- random_group_covariance(
    matrix(c(1, .2, .1, .2, 1, .15, .1, .15, 1), 3L, 3L,
           dimnames = rep(list(c("s1", "s2", "s3")), 2L))
  )
  covariance_reason <- paste0(
    "known group covariance with multiple random-effect columns requires ",
    "the exact noncentered parameterization"
  )
  auto <- resolve("auto", unbounded$half_normal, "us", kernel,
                  slope_matrix, slope_groups)
  expect_identical(auto$resolved, "noncentered")
  expect_identical(auto$reason, covariance_reason)
  expect_error(
    resolve("centered", unbounded$half_normal, "us", kernel,
            slope_matrix, slope_groups),
    covariance_reason,
    fixed = TRUE
  )
  expect_identical(
    resolve("auto", unbounded$half_normal, "us", NULL,
            slope_matrix, slope_groups)$resolved,
    "centered"
  )
  expect_identical(
    resolve("auto", unbounded$half_normal, "us", kernel,
            slope_matrix[, 1L, drop = FALSE], slope_groups)$resolved,
    "centered"
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

test_that("random-effect SD priors must have nonnegative support", {

  expect_error(
    prior_random(sd = prior("normal", list(0, 1))),
    "The 'sd' prior must have nonnegative support.",
    fixed = TRUE
  )
  expect_error(
    prior_random(sd = prior("normal", list(0, 1),
                            truncation = list(-1, Inf))),
    "The 'sd' prior must have nonnegative support.",
    fixed = TRUE
  )
  expect_error(
    random_block(sd = prior("cauchy", list(0, 1))),
    "The 'sd' prior must have nonnegative support.",
    fixed = TRUE
  )
  expect_error(
    random_covariance(sd = prior("normal", list(0, 1))),
    "The 'sd' prior must have nonnegative support.",
    fixed = TRUE
  )
  expect_error(
    prior_random(sd = prior("point", list(location = -1))),
    "The 'sd' prior must have nonnegative support.",
    fixed = TRUE
  )
  # A component of a mixture is checked like any other SD prior.
  expect_error(
    prior_random(sd = prior_spike_and_slab(
      prior_parameter = prior("normal", list(0, 1)),
      prior_inclusion = prior("beta", list(1, 1))
    )),
    "The 'sd' prior must have nonnegative support.",
    fixed = TRUE
  )

  # Families whose support is already nonnegative are accepted unchanged.
  expect_s3_class(
    prior_random(sd = prior("normal", list(0, 1),
                            truncation = list(0, Inf))),
    "prior_random"
  )
  expect_s3_class(prior_random(sd = prior("lognormal", list(0, 1))),
                  "prior_random")
  expect_s3_class(prior_random(sd = prior("invgamma", list(2, 1))),
                  "prior_random")
  expect_s3_class(prior_random(sd = prior("point", list(location = 0))),
                  "prior_random")
})

test_that("term-specific SD overrides are validated at construction", {

  # A bare prior and random_block(sd = ...) are equivalent term overrides;
  # both must be rejected before the backend could silently truncate them.
  expect_error(
    random_block(terms = list(intercept = prior("normal", list(0, 1)))),
    "The 'sd' prior must have nonnegative support.",
    fixed = TRUE
  )
  expect_error(
    random_block(terms = list(
      intercept = random_block(sd = prior("normal", list(0, 1)))
    )),
    "The 'sd' prior must have nonnegative support.",
    fixed = TRUE
  )
  expect_error(
    prior_random(study = random_block(terms = list(
      x = prior_mixture(list(prior("point", list(-1)), prior("gamma", list(2, 2))))
    ))),
    "The 'sd' prior must have nonnegative support.",
    fixed = TRUE
  )
  expect_s3_class(
    random_block(terms = list(
      intercept = .parameterization_sd_prior(),
      x = prior("gamma", list(2, 2))
    )),
    "random_block"
  )
})
