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

# four groups observed on a three-level factor twice, a row-level covariate,
# and a group-level covariate
.mean_design_data <- function(){

  data <- expand.grid(rep = 1:2, f = factor(c("a", "b", "c")), g = factor(1:4))
  data$x <- rep(c(-1, 0.5, 1.5), length.out = nrow(data)) + as.integer(data$g) / 10
  data$w <- c(-0.4, 0.1, 0.3, 0.6)[as.integer(data$g)]
  data
}

.mean_design_compile <- function(formula, random, data = .mean_design_data(), ...){

  unit <- prior("normal", list(0, 1))
  priors <- list(intercept = unit, x = unit, w = unit, v = unit,
                 f = prior_factor("normal", list(0, 1), contrast = "treatment"))
  terms <- c("intercept", attr(stats::terms(formula), "term.labels"))
  terms <- terms[!grepl("|", terms, fixed = TRUE)]
  JAGS_formula(formula, "mu", data, prior_list = priors[intersect(names(priors), terms)],
               prior_random = random, ...)
}

# The location maps reproduce the absorbed fixed columns on the fitted rows:
# Z A plus the stored residual = X_a for constant maps, the group values for
# group maps, and Z P = Z_p for a linked parent.
.expect_mean_translation_identity <- function(result){

  design <- result$formula_design
  X <- design$model_matrix
  blocks <- stats::setNames(design$random_effects,
                            vapply(design$random_effects, `[[`, character(1), "block_name"))
  for(term in blocks){
    spec <- term$mean_translation
    if(is.null(spec)) next
    Z <- term$model_matrix
    if(!is.null(spec$constant)){
      residual <- if(is.null(spec$constant$residual)) 0 else spec$constant$residual
      expect_lt(max(abs(Z %*% spec$constant$map + residual - X[, spec$constant$columns, drop = FALSE])), 1e-15)
      expect_identical(is.null(spec$residual_term), is.null(spec$constant$residual))
    }
    if(!is.null(spec$group)){
      expect_true(all(Z[, spec$group$column] == 1))
      expect_identical(unname(X[, spec$group$columns, drop = FALSE]),
                       spec$group$values[term$group_map, , drop = FALSE])
    }
    if(!is.null(spec$parent)){
      parent <- blocks[[spec$parent$block]]
      expect_identical(spec$parent$index[term$group_map], parent$group_map)
      expect_lt(max(abs(Z %*% spec$parent$map - parent$model_matrix)), 1e-12)
      expect_identical(parent$mean_translation$child, term$block_name)
    }
  }
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
    expect_identical(plan$absorbed_terms, "intercept")
    expect_identical(plan$constant$expressions, contribution)
    expect_match(translated$formula_syntax,
      paste0(plan$location_name, "[g,1] ~ dnorm(", contribution, ","), fixed = TRUE)
    expect_match(translated$formula_syntax,
      paste0(plan$coefficient_name, "[g,1] <- ", plan$location_name, "[g,1] - ", contribution),
      fixed = TRUE)
    expect_match(translated$formula_syntax,
      paste0("mu[i] = mu_x * mu_data_x[i] + ", plan$location_name,
        "[", plan$group_map_name, "[i],1]"), fixed = TRUE)
    # the scalar block reads no design columns
    expect_false(any(grepl("_xRE_DATAx", names(translated$data), fixed = TRUE)))
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
  changed$random_effects[[1L]]$mean_translation$constant$expressions <- "mu_x"
  condition <- tryCatch(.bt_JAGS_bridge_validate_formula_random_design("mu", design, changed),
    error = identity)
  expect_identical(conditionMessage(condition), paste0(
    "JAGS_bridgesampling() rebuilt formula random-effect design does not match the fitted design for parameter 'mu', ",
    "block 'study': random-effect parameterization differs. ",
    "Supply the same formula, data, scaling, and prior_random() metadata used to fit the model."))
  expect_null(conditionCall(condition))
})

test_that("mean centering absorbs the fixed terms that a multi-column block spans", {

  sd <- prior("normal", list(0, 1), list(0, Inf))
  block <- function(...) prior_random(g = random_block(sd = sd, parameterization = "mean_centered", ...))

  # indicator columns of f span the intercept and the treatment-coded f
  us <- .mean_design_compile(~ f + us(0 + f | g),
    block(cor = prior_lkj(eta = 1), contrasts = c(f = "independent")))
  spec <- us$formula_design$random_effects[[1L]]$mean_translation
  expect_identical(spec$absorbed_terms, c("intercept", "f"))
  expect_equal(spec$constant$map, rbind(c(1, 0, 0), c(1, 1, 0), c(1, 0, 1)))
  expect_identical(spec$constant$expressions, c("mu_intercept", "mu_f[1]", "mu_f[2]"))
  expect_match(us$formula_syntax, "mu__xREx__g_xRE_LOCx[2] <- mu_intercept + mu_f[1]", fixed = TRUE)
  # sequential conditionals through the LKJ factor
  expect_match(us$formula_syntax, paste0(
    "mu__xREx__g_xRE_MEANx[g,2] ~ dnorm(mu__xREx__g_xRE_LOCx[2] + mu__xREx__g_xRE_STDx[2] * ",
    "(mu__xREx__g_xRE_CORx_L[2,1] * mu__xREx__g_xRE_Zx[g,1]), ",
    "pow(mu__xREx__g_xRE_STDx[2] * mu__xREx__g_xRE_CORx_L[2,2], -2))"), fixed = TRUE)
  expect_match(us$formula_syntax, paste0(
    "mu[i] = inprod(mu__xREx__g_xRE_MEANx[mu__xREx__g_xRE_MAPx[i],1:3], ",
    "mu__xREx__g_xRE_DATAx[i,1:3])"), fixed = TRUE)
  expect_false(grepl("dmnorm", us$formula_syntax, fixed = TRUE))
  expect_false("mu_data_f" %in% names(us$data))
  expect_null(us$formula_design$jags_data_names[["f"]])
  .expect_mean_translation_identity(us)

  # a random intercept and slope span the intercept and the slope of x, also
  # when x is scaled (both parts read the same scaled column)
  for(scale in list(NULL, list(x = TRUE))){
    slope <- .mean_design_compile(~ x + (1 + x | g), block(cor = prior_lkj(eta = 1)), formula_scale = scale)
    spec <- slope$formula_design$random_effects[[1L]]$mean_translation
    expect_identical(spec$absorbed_terms, c("intercept", "x"))
    expect_equal(spec$constant$map, diag(2))
    .expect_mean_translation_identity(slope)
  }

  # independent and structured blocks: scalar conditionals, an indicator
  # layout reads the column of each row
  diag_block <- .mean_design_compile(~ f + diag(0 + f | g), block(contrasts = c(f = "independent")))
  expect_match(diag_block$formula_syntax,
    "mu__xREx__g_xRE_MEANx[g,3] ~ dnorm(mu__xREx__g_xRE_LOCx[3], pow(mu__xREx__g_xRE_STDx[3], -2))",
    fixed = TRUE)
  .expect_mean_translation_identity(diag_block)
  cs <- .mean_design_compile(~ f + cs(f | g), block())
  expect_match(cs$formula_syntax,
    "mu[i] = mu__xREx__g_xRE_MEANx[mu__xREx__g_xRE_MAPx[i],mu__xREx__g_xRE_COLx[i]]", fixed = TRUE)
  # the compound-symmetry prefix recurrence of the noncentered transform
  expect_match(cs$formula_syntax, paste0(
    "mu__xREx__g_xRE_MEANx[g,3] ~ dnorm(mu__xREx__g_xRE_LOCx[3] + mu__xREx__g_xRE_STDx[3] * ",
    "(mu__xREx__g_xRE_CS_PREFIXx[g,2]), pow(mu__xREx__g_xRE_STDx[3] * mu__xREx__g_xRE_CS_DIAGx[3], -2))"),
    fixed = TRUE)
  expect_false(grepl("_xRE_CORx_L", cs$formula_syntax, fixed = TRUE))
  ar1 <- .mean_design_compile(~ f + ar1(f | g), block())
  expect_match(ar1$formula_syntax, paste0(
    "mu__xREx__g_xRE_MEANx[g,2] ~ dnorm(mu__xREx__g_xRE_LOCx[2] + mu__xREx__g_xRE_STDx[2] * ",
    "(mu__xREx__g_rho * mu__xREx__g_xRE_UNIT_COEFx[g,1]), ",
    "pow(mu__xREx__g_xRE_STDx[2] * mu__xREx__g_xRE_AR_INNOVx, -2))"), fixed = TRUE)
  .expect_mean_translation_identity(ar1)
  .expect_mean_translation_identity(cs)

  # terms outside the span, and a term only partly in it (a random slope on
  # one indicator of f), stay in the predictor
  partial <- .mean_design_compile(~ x + f + (1 + x | g), block(cor = prior_lkj(eta = 1)))
  spec <- partial$formula_design$random_effects[[1L]]$mean_translation
  expect_identical(spec$absorbed_terms, c("intercept", "x"))
  expect_match(partial$formula_syntax, "mu[i] = inprod(mu_f, mu_data_f[i,]) + inprod(", fixed = TRUE)
  .expect_mean_translation_identity(partial)
  data <- .mean_design_data()
  data$b <- as.numeric(data$f == "b")
  indicator <- .mean_design_compile(~ f + (1 + b | g), block(cor = prior_lkj(eta = 1)), data = data)
  expect_identical(indicator$formula_design$random_effects[[1L]]$mean_translation$absorbed_terms, "intercept")

  # a map that reproduces its columns only up to rounding (mean-difference
  # random contrasts against treatment-coded fixed effects) keeps the residual
  # in the predictor
  meandif <- .mean_design_compile(~ f + (1 + f | g), block(cor = prior_lkj(eta = 1), contrasts = c(f = "meandif")))
  spec <- meandif$formula_design$random_effects[[1L]]$mean_translation
  expect_identical(spec$absorbed_terms, c("intercept", "f"))
  design <- meandif$formula_design
  residual <- design$model_matrix[, spec$constant$columns] - design$random_effects[[1L]]$model_matrix %*% spec$constant$map
  if(any(residual != 0)){
    expect_equal(meandif$data$mu__xREx__g_xRE_RESIDx, unname(residual))
    expect_match(meandif$formula_syntax,
      "+ inprod(mu__xREx__g_xRE_ACOEFx[1:3], mu__xREx__g_xRE_RESIDx[i,1:3])", fixed = TRUE)
  }else{
    expect_null(spec$constant$residual)
  }
  .expect_mean_translation_identity(meandif)
})

test_that("mean centering absorbs group-level terms and links nested blocks", {

  sd <- prior("normal", list(0, 1), list(0, Inf))
  centered <- function() random_block(sd = sd, parameterization = "mean_centered")

  # a group-level covariate enters the location of the random intercept
  moderator <- .mean_design_compile(~ x + w + (1 | g), prior_random(g = centered()))
  spec <- moderator$formula_design$random_effects[[1L]]$mean_translation
  expect_identical(spec$absorbed_terms, c("intercept", "w"))
  expect_identical(spec$group$expressions, "mu_w")
  expect_match(moderator$formula_syntax, paste0(
    "mu__xREx__g_xRE_LOCx[g,1] <- mu_intercept + inprod(mu__xREx__g_xRE_GLEVx[g,1:1], ",
    "mu__xREx__g_xRE_GCOEFx[1:1])"), fixed = TRUE)
  expect_equal(moderator$data$mu__xREx__g_xRE_GLEVx, matrix(c(-0.4, 0.1, 0.3, 0.6), ncol = 1L))
  expect_false("mu_data_w" %in% names(moderator$data))
  .expect_mean_translation_identity(moderator)

  # nested blocks: the coarser block takes the intercept, the finer one f
  # (constant within its groups) and the coarser block's locations
  nested <- .mean_design_compile(~ f + (1 | g/f),
    prior_random(sd = sd, g = centered(), f_g = centered()))
  blocks <- stats::setNames(nested$formula_design$random_effects,
    vapply(nested$formula_design$random_effects, `[[`, character(1), "block_name"))
  expect_identical(blocks$g$mean_translation$absorbed_terms, "intercept")
  expect_identical(blocks$g$mean_translation$child, "f_g")
  expect_identical(blocks$f_g$mean_translation$absorbed_terms, "f")
  expect_identical(blocks$f_g$mean_translation$parent$block, "g")
  expect_match(nested$formula_syntax, "mu[i] = mu__xREx__f_g_xRE_MEANx[mu__xREx__f_g_xRE_MAPx[i],1]", fixed = TRUE)
  expect_match(nested$formula_syntax,
    "mu__xREx__g_xRE_MEANx[mu__xREx__f_g_xRE_PARENTx[g],1]", fixed = TRUE)
  # the parent's rows enter only through its child's locations
  expect_false("mu__xREx__g_xRE_MAPx" %in% names(nested$data))
  .expect_mean_translation_identity(nested)

  # a parent whose map reproduces a nearly constant column only up to a
  # residual keeps that residual in the predictor when its child absorbs its
  # locations
  data <- .mean_design_data()
  data$v <- 1 + rep(c(-1e-10, 1e-10), length.out = nrow(data))
  residual <- .mean_design_compile(~ v + (1 | g/f),
    prior_random(sd = sd, g = centered(), f_g = centered()), data = data)
  blocks <- stats::setNames(residual$formula_design$random_effects,
    vapply(residual$formula_design$random_effects, `[[`, character(1), "block_name"))
  expect_identical(blocks$g$mean_translation$absorbed_terms, c("intercept", "v"))
  expect_identical(blocks$g$mean_translation$child, "f_g")
  expect_match(residual$formula_syntax, paste0(
    "mu[i] = mu__xREx__f_g_xRE_MEANx[mu__xREx__f_g_xRE_MAPx[i],1] + ",
    "inprod(mu__xREx__g_xRE_ACOEFx[1:2], mu__xREx__g_xRE_RESIDx[i,1:2])"), fixed = TRUE)
  .expect_mean_translation_identity(residual)

  # a declared but unobserved level has no group value and no parent: the
  # group-level term stays in the predictor and the blocks are not linked
  data <- .mean_design_data()
  data$g <- factor(data$g, levels = c(levels(data$g), "unused"))
  unused <- .mean_design_compile(~ x + w + (1 | g), prior_random(g = centered()), data = data)
  spec <- unused$formula_design$random_effects[[1L]]$mean_translation
  expect_identical(unused$formula_design$random_effects[[1L]]$n_groups, 5L)
  expect_identical(spec$absorbed_terms, "intercept")
  expect_null(spec$group)
  .expect_mean_translation_identity(unused)
})

test_that("mean centering requires an eligible block with something to translate", {

  expect_error(prior_random(parameterization = "mean_centered"),
    "Mean-centered parameterization requires an explicitly named 'random_block()'.", fixed = TRUE)
  expect_error(.mean_parameterization_compile(sd = prior("point", list(0))),
    "Mean-centered parameterization is unavailable for random-effect block 'study': SD prior with an atom at zero.",
    fixed = TRUE)
  expect_error(.mean_parameterization_compile(marginalized = "study"),
    "Mean-centered parameterization is unavailable for marginalized random-effect blocks.", fixed = TRUE)
  sd <- prior("normal", list(0, 1), list(0, Inf))
  # crossed blocks: the second has no fixed term left and no parent
  expect_error(.mean_parameterization_compile(formula = ~ 1 + x + id(1 | study) + id(1 | site),
    random = prior_random(study = random_block(sd = sd, parameterization = "mean_centered"),
      site = random_block(sd = sd, parameterization = "mean_centered"))),
    paste0("Mean-centered parameterization is unavailable for random-effect block 'site': its random-effect ",
      "design spans no fixed term that another mean-centered block does not already absorb, and no ",
      "mean-centered parent block."), fixed = TRUE)

  centered <- prior_random(g = random_block(sd = sd, parameterization = "mean_centered"))
  data <- .mean_design_data()
  # collinear random-effect columns
  data$x2 <- 2 * data$x
  expect_error(.mean_design_compile(~ x + (1 + x + x2 | g), centered, data = data),
    "Mean-centered parameterization is unavailable for random-effect block 'g': its random-effect design is rank deficient.",
    fixed = TRUE)
})
