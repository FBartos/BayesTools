.format_column       <- function(x, type, n_models){
  if(is.null(x) || length(x) == 0){
    return(x)
  }else{
    return(switch(
      type,
      "integer"         = round(x),
      "prior"           = .center_priors(x),
      "string_left"     = .string_left(x),
      "string"          = x,
      "hypothesis_label" = x,
      "estimate"        = format(round(x, digits = 3), nsmall = 3),
      "prior_prob"      = format(round(x, digits = 3), nsmall = 3),
      "post_prob"       = format(round(x, digits = 3), nsmall = 3),
      "probability"     = format(round(x, digits = 3), nsmall = 3),
      "marglik"         = format(round(x, digits = 2), nsmall = 2),
      "BF"              = .format_BF_column(x),
      "inclusion_BF"    = .format_BF_column(x),
      "BF_error"        = format(round(x, digits = 3), nsmall = 3),
      "n_models"        = paste0(round(x), "/", n_models),
      "ESS"             = round(x),
      "R_hat"           = format(round(x, digits = 3), nsmall = 3),
      "MCMC_error"      = format(round(x, digits = 5), nsmall = 5),
      "MCMC_SD_error"   = format(round(x, digits = 3), nsmall = 3),
      "min_ESS"             = round(x),
      "max_R_hat"           = format(round(x, digits = 3), nsmall = 3),
      "max_MCMC_error"      = format(round(x, digits = 5), nsmall = 5),
      "max_MCMC_SD_error"   = format(round(x, digits = 3), nsmall = 3)
    ))
  }
}
.format_column_names <- function(x, type, values){
  if(is.null(x) || length(x) == 0){
    return(x)
  }else{
    return(switch(
      type,
      "integer"         = x,
      "prior"           = .string_center(paste0("Prior ", x), values),
      "string_left"     = .string_left(x, values),
      "string"          = x,
      "hypothesis_label" = paste0(x, ":"),
      "estimate"        = x,
      "probability"     = x,
      "prior_prob"      = "Prior prob.",
      "post_prob"       = "Post. prob.",
      "marglik"         = "log(marglik)",
      "BF"              = if(is.null(attr(values, "name"))) "BF"           else attr(values, "name"),
      "inclusion_BF"    = if(is.null(attr(values, "name"))) "Inclusion BF" else attr(values, "name"),
      "BF_error"        = if(is.null(attr(values, "name"))) x              else attr(values, "name"),
      "n_models"        = "Models",
      "ESS"             = "ESS",
      "R_hat"           = "R-hat",
      "MCMC_error"      = if(is.null(attr(values, "name"))) "error(MCMC)" else attr(values, "name"),
      "MCMC_SD_error"   = "error(MCMC)/SD",
      "min_ESS"             = "min(ESS)",
      "max_R_hat"           = "max(R-hat)",
      "max_MCMC_error"      = "max[error(MCMC)]",
      "max_MCMC_SD_error"   = "max[error(MCMC)/SD]",
    ))
  }
}
.center_priors <- function(x){

  if(any(grepl("~", x) | grepl("=", x))){

    position_tilda  <- regexpr("~", x)
    position_equal  <- regexpr("=", x)
    from_right      <- sapply(seq_along(x), function(i){
      if(position_tilda[i] != -1){
        return(nchar(x[i]) - position_tilda[i])
      }else if(position_equal[i] != -1){
        return(nchar(x[i]) - position_equal[i])
      }else{
        return(0)
      }
    })
    add_to_right  <- ifelse(from_right == 0, 0, max(from_right) - from_right)
    x <- paste0(x, sapply(seq_along(x), function(i)paste0(rep(" ", add_to_right[i]), collapse = "")))

  }

  return(x)
}
.string_left   <- function(x, reference = x){

  if(length(x) > 0){

    add_to_right <- max(nchar(reference)) - nchar(x)
    add_to_right <- ifelse(add_to_right < 0, 0, add_to_right)
    x <- paste0(x, sapply(seq_along(x), function(i)paste0(rep(" ", add_to_right[i]), collapse = "")))

  }

  return(x)
}
.string_center <- function(x, reference = x){

  if(length(x) > 0){

    add_to_sides <- max(nchar(reference)) - nchar(x)
    add_to_sides <- ifelse(add_to_sides < 0, 0, add_to_sides)
    x <- paste0( sapply(seq_along(x), function(i)paste0(rep(" ", round(add_to_sides[i]/2)), collapse = "")), x, sapply(seq_along(x), function(i)paste0(rep(" ", round(add_to_sides[i]/2)), collapse = "")))

  }

  return(x)
}

.check_table_types      <- function(x, name, allow_NULL = FALSE){
  check_char(x, name, allow_values = c(
    "integer", "prior", "string_left", "string", "hypothesis_label", "estimate",
    "probability", "prior_prob", "post_prob",
    "marglik", "BF", "inclusion_BF", "BF_error", "n_models",
    "ESS", "R_hat", "MCMC_error", "MCMC_SD_error", "min_ESS", "max_R_hat", "max_MCMC_error", "max_MCMC_SD_error"),
    allow_NULL = allow_NULL)
}
.JAGS_estimates_diagnostic_columns <- function(){
  c("MCMC_error", "MCMC_SD_error", "ESS", "R_hat")
}
.JAGS_BF_diagnostic_columns <- function(){
  c("ESS", "MCMC_error", "BF_error_percent")
}
.JAGS_BF_diagnostic_column_types <- function(columns){
  unname(c(
    ESS              = "ESS",
    MCMC_error       = "MCMC_error",
    BF_error_percent = "BF_error"
  )[columns])
}
.normalize_diagnostic_columns <- function(columns, allowed_columns, name){

  if(is.null(columns) || length(columns) == 0){
    return(character())
  }

  if(is.logical(columns)){
    check_bool(columns, name, allow_NA = FALSE)
    if(columns){
      return(allowed_columns)
    }else{
      return(character())
    }
  }

  check_char(columns, name, check_length = 0, allow_NA = FALSE)

  shortcut_columns <- c("all", "none")
  if(any(columns %in% shortcut_columns)){
    if(length(columns) > 1){
      stop(paste0("The '", name, "' argument can use 'all' or 'none' only by itself."), call. = FALSE)
    }
    if(columns == "all"){
      return(allowed_columns)
    }else{
      return(character())
    }
  }

  check_char(columns, name, check_length = 0, allow_values = allowed_columns, allow_NA = FALSE)

  return(unique(columns))
}

# Inclusion rows summarise 0/1 inclusion-indicator draws. Their Mean is the
# posterior inclusion probability; the SD and quantiles of the draws are
# deterministic functions of that probability (not its uncertainty, which the
# MCMC error reports), so they are left blank. The rows are identified from
# metadata, never from their labels:
# - spike-and-slab and mixture indicators renamed by the estimates table
#   (`inclusion_columns`, recorded where the table creates them),
# - semantic random-effect inclusion quantities (summary priors with
#   `random_summary = "inclusion"`; random-SD spike-and-slab and mixture
#   components and variance-allocation gates),
# - raw variance-allocation gate indicators of the formula design, and
# - indicators of spike-and-slab totals of ordered-factor priors.
# `parameter_names` are the table's columns before formula renaming.
.runjags_summary_inclusion_rows <- function(parameter_names, prior_list,
                                            formula_design = NULL,
                                            inclusion_columns = character()){

  if(length(parameter_names) == 0L){
    return(logical())
  }

  random_inclusion <- vapply(
    prior_list[match(parameter_names, names(prior_list))],
    function(prior) !is.null(prior) &&
      identical(attr(prior, "random_summary", exact = TRUE), "inclusion"),
    logical(1)
  )

  ordered_total_indicators <- character()
  for(parameter in names(prior_list)){
    prior <- prior_list[[parameter]]
    if(is.prior.ordered(prior) && is.prior.spike_and_slab(prior$total)){
      ordered_total_indicators <- c(
        ordered_total_indicators,
        paste0(.prior_ordered_total_name(parameter), "_indicator")
      )
    }
  }

  indicator_columns <- c(
    inclusion_columns,
    .bt_random_variance_allocation_inclusion_indicator_names(formula_design),
    ordered_total_indicators
  )

  unname(random_inclusion | parameter_names %in% indicator_columns)
}

# 'structural' marks rows of structural constants (declared by their priors,
# such as the reference publication-weight bin): they have no MCMC
# diagnostics.
.runjags_summary_fast   <- function(model_samples, n_samples, n_chains, conditional, probs = c(0.025, 0.975), remove_diagnostics = FALSE,
                                    diagnostic_columns = .JAGS_estimates_diagnostic_columns(),
                                    inclusion = rep(FALSE, ncol(model_samples)),
                                    structural = rep(FALSE, ncol(model_samples))){

  diagnostic_columns <- .normalize_diagnostic_columns(diagnostic_columns, .JAGS_estimates_diagnostic_columns(), "diagnostic_columns")
  if(remove_diagnostics){
    diagnostic_columns <- character()
  }
  check_bool(inclusion, "inclusion", check_length = ncol(model_samples),
             allow_NULL = ncol(model_samples) == 0L)
  check_bool(structural, "structural", check_length = ncol(model_samples),
             allow_NULL = ncol(model_samples) == 0L)

  # compute quantiles dynamically
  quantile_cols <- lapply(probs, function(p) apply(model_samples, 2, stats::quantile, probs = p, na.rm = TRUE))
  names(quantile_cols) <- as.character(probs)

  # the chains needs to be kept merged for conditional summary (due to NAs in the chains)
  runjags_summary <- cbind.data.frame(
    "Mean"   = apply(model_samples, 2, mean,          na.rm = TRUE),
    "SD"     = apply(model_samples, 2, stats::sd,    na.rm = TRUE),
    as.data.frame(quantile_cols, check.names = FALSE)
  )

  # inclusion rows keep only the Mean (and the MCMC diagnostics)
  quantile_col_names <- as.character(probs)
  runjags_summary[inclusion, c("SD", quantile_col_names)] <- NA

  # don't produce fit diagnostics for conditional samples (different chain lengths etc...) or if remove_diagnostics is TRUE
  if(conditional || length(diagnostic_columns) == 0){
    return(runjags_summary)
  }

  # Split back the chains for diagnostics. Quantities that are undefined for
  # some product-space branches retain NA diagnostics while their estimates are
  # summarized from the branches on which they are defined.
  complete_columns <- colSums(is.na(model_samples)) == 0L
  runjags_diagnostics <- as.data.frame(matrix(
    NA_real_,
    nrow = ncol(model_samples),
    ncol = 4L,
    dimnames = list(
      colnames(model_samples),
      c("MCMC_error", "MCMC_SD_error", "ESS", "R_hat")
    )
  ))
  if(any(complete_columns)){
    diagnostic_samples <- model_samples[, complete_columns, drop = FALSE]
    model_samples_list <- split(
      as.data.frame(diagnostic_samples),
      rep(seq_len(n_chains), each = n_samples),
      drop = FALSE
    )
    model_samples_list <- coda::as.mcmc.list(lapply(
      model_samples_list,
      coda::as.mcmc
    ))
    mcmc_summary <- summary(model_samples_list, quantiles = NULL)$statistics

    # fix single parameter summaries
    if(is.null(dim(mcmc_summary))){
      mcmc_summary <- t(mcmc_summary)
    }

    complete_diagnostics <- cbind.data.frame(
      "MCMC_error"    = mcmc_summary[, "Time-series SE"],
      "MCMC_SD_error" = mcmc_summary[, "Time-series SE"] / mcmc_summary[, "SD"],
      "ESS"           = coda::effectiveSize(model_samples_list),
      "R_hat"         = if(n_chains > 1){
        coda::gelman.diag(
          model_samples_list,
          multivariate = FALSE,
          autoburnin = FALSE
        )$psrf[, 1]
      }else{
        NA
      }
    )
    runjags_diagnostics[complete_columns, ] <- complete_diagnostics
  }

  # remove incorrect NANs and NAs from the diagnostics
  runjags_diagnostics[is.nan(runjags_diagnostics[,"MCMC_SD_error"]),"MCMC_SD_error"] <- NA
  zero_ess <- !is.na(runjags_diagnostics[, "ESS"]) &
    runjags_diagnostics[, "ESS"] == 0
  runjags_diagnostics[zero_ess, "ESS"]                                               <- 0
  runjags_diagnostics[is.nan(runjags_diagnostics[,"R_hat"]),"R_hat"]                 <- NA

  # structural constants have no MCMC diagnostics
  runjags_diagnostics[structural, ] <- NA
  runjags_summary <- cbind.data.frame(runjags_summary, runjags_diagnostics[, diagnostic_columns, drop = FALSE])

  return(runjags_summary)
}
# Quantities whose draws can be undefined (NA) declare why through their
# 'undefined_draws' metadata (see parameter_draws() and the catalog's
# 'definedness'): "correlation" for an original-scale random-effect
# correlation (NA when an SD is zero), "allocation_active" for the variance
# share of a gated total-variance allocation (NA when no component is
# active), and "positive_definite" for raw Cholesky and LKJ coordinates (NA
# unless the correlation matrix is positive definite). Summaries of such quantities use
# the defined draws and report their share in a row footnote.
.bt_undefined_draws_reasons <- c(
  correlation       = "where the correlation is defined, i.e. both SDs are positive",
  allocation_active = "where the variance share is defined, i.e. an allocation component is active",
  ordered_parameterization = "where the model contains this ordered allocation parameterization",
  positive_definite = "where the correlation matrix is positive definite"
)

.bt_undefined_draws_footnote <- function(row, n_defined, n_draws, reason){

  if(!is.character(reason) || length(reason) != 1L || is.na(reason) ||
     !reason %in% names(.bt_undefined_draws_reasons)){
    stop(
      "Unknown 'undefined_draws' declaration for '", row, "'. Use one of ",
      paste0("\"", names(.bt_undefined_draws_reasons), "\"", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  stats::setNames(
    paste0(
      row, ": summarized over ", n_defined, " of ", n_draws, " draws ",
      .bt_undefined_draws_reasons[[reason]], "."
    ),
    row
  )
}

.runjags_conditional_warning <- function(parameters, n_samples, warning_limit = 500){
  if(n_samples == 0){
    return(sprintf("Conditional summary for %1$s parameter could not be computed due to no posterior samples.", paste0(parameters, collapse = ", ")))
  }else if(n_samples <= warning_limit){
    return(sprintf("Conditional summary for %1$s is based on %2$i samples.", paste0(parameters, collapse = ", "), n_samples))
  }else{
    return()
  }
}

# Transform mixed and marginal posterior samples (as used by the ensemble and
# marginal tables) from the standardized to the original predictor scale.
# Every column is mapped to the fitted coordinates it is a linear combination
# of through its 'quantities' draw metadata (a vector element without that
# metadata is the parameter it is named after); the coordinates are
# transformed through the fitted design of 'formula_scale', and each column is
# recombined from the transformed coordinates. The coordinates of combination
# columns (transformed factor levels) are recovered from the levels of the same
# element. Estimated marginal means of a formula are predictions at the stated
# predictor values; they do not depend on the predictor standardization and
# are returned unchanged.
.transform_scale_samples_list <- function(samples, formula_scale){

  if(is.null(formula_scale) || length(formula_scale) == 0){
    return(samples)
  }

  scaled_parameters <- names(formula_scale)
  elements <- list()
  for(name in names(samples)){
    x <- samples[[name]]
    if(inherits(x, "marginal_posterior.formula") ||
       !(is.numeric(x) || (is.list(x) && length(x) > 0L &&
                             all(vapply(x, is.numeric, logical(1)))))){
      next
    }
    element <- .transform_scale_element_columns(x, name, scaled_parameters)
    if(!is.null(element)){
      elements[[name]] <- element
    }
  }
  if(length(elements) == 0L){
    return(samples)
  }
  n_samples <- unique(vapply(elements, function(element){
    nrow(element$coordinates)
  }, integer(1)))
  if(length(n_samples) != 1L){
    stop("All posterior sample elements must contain the same number of samples.", call. = FALSE)
  }

  coordinate_values <- do.call(cbind, unname(lapply(elements, `[[`, "coordinates")))
  coordinate_values <- coordinate_values[
    ,
    !duplicated(colnames(coordinate_values)),
    drop = FALSE
  ]
  owners <- unlist(unname(lapply(elements, `[[`, "owners")))
  owners <- owners[!duplicated(names(owners))]
  .transform_scale_check_design_coordinates(
    owners        = owners,
    formula_scale = formula_scale
  )
  .transform_scale_check_random_coordinates(
    owners        = owners,
    formula_scale = formula_scale,
    declared      = unique(unlist(lapply(elements, `[[`, "declared"), use.names = FALSE))
  )
  transformed <- .apply_unscale_transform(coordinate_values, formula_scale)
  original_samples <- samples

  for(name in names(elements)){
    element <- elements[[name]]
    values <- transformed[, colnames(element$weights), drop = FALSE] %*%
      t(element$weights)
    x <- samples[[name]]
    if(is.list(x) && !is.numeric(x)){
      for(i in which(element$columns)){
        samples[[name]][[i]] <- .bt_draws_transform_values(x[[i]], function(level){
          values[, i]
        })
      }
    }else if(is.null(dim(x))){
      samples[[name]] <- .bt_draws_transform_values(x, function(level){
        values[, 1L]
      })
      samples[[name]] <- .bt_draws_set_original_scale_quantities(
        samples[[name]],
        .bt_draws_quantities(samples[[name]])$column
      )
    }else{
      samples[[name]][, element$columns] <- values[, element$columns, drop = FALSE]
      samples[[name]] <- .bt_meta_refresh(samples[[name]])
      # the transformed columns are their original-scale quantities
      samples[[name]] <- .bt_draws_set_original_scale_quantities(
        samples[[name]],
        colnames(x)[element$columns]
      )
    }
    parameter <- unique(owners[colnames(element$weights)])
    fixed_columns <- if(length(parameter)==1L) names(owners)[owners==parameter &
      !.formula_scale_matches_prefix(names(owners),parameter,"__xREx__") &
      !.formula_scale_matches_prefix(names(owners),parameter,"__xRE_ALLOCx") &
      !.formula_scale_matches_prefix(names(owners),parameter,"__xRE_SUMMARY__")]
    if(length(fixed_columns) && all(colnames(element$weights) %in% fixed_columns) &&
       !(is.list(samples[[name]]) && !is.numeric(samples[[name]]))){
      transform <- .bt_formula_coefficient_transform(fixed_columns,formula_scale[[parameter]],parameter)
      weights <- element$weights %*% transform$matrix[colnames(element$weights),,drop=FALSE]
      projections <- .bt_ordered_formula_projections(original_samples,weights,
        transform$source_transforms[transform$source_transforms!="identity"])
      samples[[name]] <- .bt_ordered_attach_linear_view(samples[[name]],projections,weights,original_samples)
    }
  }
  samples <- .bt_meta_set(samples,"transform_scaled",TRUE)
  samples <- .bt_meta_set(samples,"formula_scale",formula_scale)

  return(samples)
}

# The fitted fixed-effect coordinates transformed through the fitted design
# that 'formula_scale' carries must be coefficients of that design: a column of
# a term the design does not contain (e.g., a mixture of models whose formulas
# differ, transformed with the 'formula_scale' of a smaller model) has no
# original-scale value under it. 'owners' names the formula parameter of every
# coordinate. Random-effect and allocation coordinates have transforms of
# their own.
.transform_scale_check_design_coordinates <- function(owners, formula_scale){

  for(parameter in unique(owners)){
    spec <- attr(formula_scale[[parameter]], "unscale_design", exact = TRUE)
    if(is.null(spec)){
      next
    }
    coordinates <- names(owners)[owners == parameter]
    fixed <- !(
      .formula_scale_matches_prefix(coordinates, parameter, "__xREx__") |
        .formula_scale_matches_prefix(coordinates, parameter, "__xRE_ALLOCx") |
        .formula_scale_matches_prefix(coordinates, parameter, "__xRE_SUMMARY__")
    )
    outside <- setdiff(
      coordinates[fixed],
      .bt_formula_unscale_coefficient_names(spec, parameter)
    )
    if(length(outside) > 0L){
      .bt_formula_transform_stop(
        paste0(
          "Cannot transform ", paste0("'", outside, "'", collapse = ", "),
          " to the original predictor scale: the fitted design in ",
          "'formula_scale' of formula parameter '", parameter, "' does not ",
          "contain ", if(length(outside) > 1L) "these coefficients" else "this coefficient",
          ". Pass the 'formula_scale' of a fitted model whose formula contains ",
          "every term of the samples."
        ),
        parameter    = parameter,
        reason       = "coefficients_outside_design",
        coefficients = outside
      )
    }
  }

  invisible(TRUE)
}

# The random-effect SD coordinates of the samples are transformed with the
# random-effect structure that 'formula_scale' carries (its SD leaves), so they
# must be SDs of that structure, as fixed coordinates must be coefficients of
# its design: a coordinate the structure does not contain (a random slope the
# model of 'formula_scale' does not have) would stay on the fitted scale, and
# so would a coordinate of a block the structure leaves unchanged while the
# samples declare that its original-scale quantity is another one ('declared':
# the intercept SD of a block whose standardized slope the model of
# 'formula_scale' does not have). Both stop.
.transform_scale_check_random_coordinates <- function(owners, formula_scale,
                                                      declared = character()){

  for(parameter in unique(owners)){
    parameter_scale <- formula_scale[[parameter]]
    if(is.null(attr(parameter_scale, "unscale_design", exact = TRUE))){
      # standardization without the fitted design stops when transformed
      next
    }
    coordinates <- names(owners)[owners == parameter]
    random <- coordinates[
      .formula_scale_matches_prefix(coordinates, parameter, "__xREx__")
    ]
    if(length(random) == 0L){
      next
    }
    leaves <- attr(parameter_scale, "random_effect_sd_leaves", exact = TRUE)
    leaf_names <- unique(unlist(lapply(leaves, function(block){
      as.character(block$leaf_names)
    }), use.names = FALSE))
    outside <- setdiff(random, leaf_names)
    if(length(outside) > 0L){
      .bt_formula_transform_stop(
        paste0(
          "Cannot transform ", paste0("'", outside, "'", collapse = ", "),
          " to the original predictor scale: the random-effect structure in ",
          "'formula_scale' of formula parameter '", parameter, "' does not ",
          "contain ", if(length(outside) > 1L) "these random-effect SDs" else "this random-effect SD",
          ". Pass the 'formula_scale' of a fitted model whose random effects ",
          "contain every random-effect SD of the samples."
        ),
        parameter   = parameter,
        reason      = "random_effects_outside_structure",
        coordinates = outside
      )
    }
    groups <- .random_sd_column_unscale_groups(
      random_sd_cols = random,
      formula_scale  = parameter_scale,
      prefix         = parameter
    )
    transformed <- unlist(lapply(groups, `[[`, "leaf_names"), use.names = FALSE)
    unchanged <- setdiff(intersect(random, declared), transformed)
    if(length(unchanged) > 0L){
      .bt_formula_transform_stop(
        paste0(
          "Cannot transform ", paste0("'", unchanged, "'", collapse = ", "),
          " to the original predictor scale: the random-effect structure in ",
          "'formula_scale' of formula parameter '", parameter, "' leaves ",
          if(length(unchanged) > 1L) "these random-effect SDs" else "this random-effect SD",
          " unchanged, while the samples come from a model whose ",
          "original-scale SD differs from the fitted one. Pass the ",
          "'formula_scale' of the fitted model of the samples."
        ),
        parameter   = parameter,
        reason      = "random_effect_structure_differs",
        coordinates = unchanged
      )
    }
  }

  invisible(TRUE)
}

# The columns of one sample element as linear combinations of fitted
# coordinates: a list with the draws of those coordinates, the weight matrix
# (columns x coordinates), the columns that belong to a scaled formula, the
# formula parameter owning each coordinate ('owners', named by coordinate),
# and the coordinates of columns that declare a different original-scale
# quantity ('declared', from 'original_scale_quantities'). NULL when no column
# belongs to a scaled formula.
.transform_scale_element_columns <- function(x, name, scaled_parameters){

  levels <- is.list(x) && !is.numeric(x)
  values <- if(levels){
    do.call(cbind, lapply(x, as.numeric))
  }else if(is.null(dim(x))){
    matrix(as.numeric(x), ncol = 1L)
  }else{
    matrix(as.numeric(x), nrow = nrow(x))
  }
  quantities <- if(levels){
    level_quantities <- lapply(x, .bt_draws_quantities)
    if(any(vapply(level_quantities, is.null, logical(1)))){
      NULL
    }else{
      do.call(rbind, level_quantities)
    }
  }else{
    .bt_draws_quantities(x)
  }

  if(is.null(quantities)){
    formula_parameter <- .bt_meta_get(x, "formula_parameter")
    if(!any(formula_parameter %in% scaled_parameters)){
      return(NULL)
    }
    if(!levels && is.null(dim(x))){
      # a vector element without column metadata is its fitted coordinate
      colnames(values) <- name
      return(list(
        coordinates = values,
        weights     = matrix(1, nrow = 1L, ncol = 1L, dimnames = list(NULL, name)),
        columns     = TRUE,
        owners      = stats::setNames(formula_parameter[[1L]], name),
        declared    = character()
      ))
    }
    stop(
      "The posterior samples of '", name, "' do not identify their fitted ",
      "coordinates (draw metadata 'quantities'). Create them with ",
      "as_mixed_posteriors() or mix_posteriors().",
      call. = FALSE
    )
  }

  formula_parameters <- vapply(
    quantities$label_parts, `[[`, character(1), "formula_parameter"
  )
  scaled <- formula_parameters %in% scaled_parameters
  if(!any(scaled)){
    return(NULL)
  }
  # estimated marginal means are predictions, never coefficient combinations
  undefined <- scaled & vapply(quantities$label_parts, `[[`, logical(1), "marginal")
  if(any(undefined)){
    stop(
      "Cannot transform '",
      .bt_label(quantities$label_parts[undefined][[1L]], style = "table"),
      "' to the original predictor scale: it is not a combination of fitted ",
      "coefficients.",
      call. = FALSE
    )
  }
  coordinates <- unique(unlist(quantities$dependencies[scaled], use.names = FALSE))
  weights <- matrix(
    0,
    nrow = nrow(quantities),
    ncol = length(coordinates),
    dimnames = list(NULL, coordinates)
  )
  for(i in which(scaled)){
    weights[i, quantities$dependencies[[i]]] <- quantities$weights[[i]]
  }
  scaled_weights <- weights[scaled, , drop = FALSE]
  scaled_values <- values[, scaled, drop = FALSE]
  identity_columns <- nrow(scaled_weights) == ncol(scaled_weights) &&
    all(rowSums(scaled_weights != 0) == 1L) &&
    all(colSums(scaled_weights != 0) == 1L) &&
    all(scaled_weights[scaled_weights != 0] == 1)
  if(identity_columns){
    # every column is one fitted coordinate
    coordinate_values <- scaled_values[
      ,
      apply(scaled_weights != 0, 2L, which),
      drop = FALSE
    ]
  }else{
    # the levels of one factor term determine its coordinates when their
    # weights have full column rank
    decomposition <- qr(scaled_weights)
    if(decomposition$rank < ncol(scaled_weights)){
      stop(
        "Cannot transform '", name, "' to the original predictor scale: its ",
        "columns do not determine the fitted coefficients they combine.",
        call. = FALSE
      )
    }
    coordinate_values <- t(qr.coef(decomposition, t(scaled_values)))
  }
  colnames(coordinate_values) <- coordinates
  owners <- unlist(lapply(which(scaled), function(i){
    stats::setNames(
      rep(formula_parameters[[i]], length(quantities$dependencies[[i]])),
      quantities$dependencies[[i]]
    )
  }))

  original <- if(!levels) .bt_meta_get(x, "original_scale_quantities")
  declared <- if(is.null(original)){
    character()
  }else{
    unlist(
      quantities$dependencies[scaled & quantities$column %in% original$column],
      use.names = FALSE
    )
  }

  list(
    coordinates = coordinate_values,
    weights     = weights,
    columns     = scaled,
    owners      = owners[!duplicated(names(owners))],
    declared    = declared
  )
}
