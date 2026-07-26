#' Formula prediction with explicit random-effect target
#'
#' @description
#' Returns formula-scale predictions from fitted BayesTools formula metadata.
#' Unlike [JAGS_evaluate_formula()], this helper can return structured
#' prediction output for marginal random-effect targets, including marginal
#' random-effect covariance or simulated random-effect draws.
#'
#' @param fit model fitted with [JAGS_fit()] or posterior samples carrying
#' `formula_design` metadata.
#' @param parameter formula parameter name.
#' @param formula optional formula. If `NULL`, the fitted formula for
#' `parameter` is used.
#' @param data optional prediction data. If `NULL`, fitted source data are used.
#' @param fitted_rows optional integer vector mapping supplied prediction rows
#' to fitted observation indices. It is required for selected
#' posterior-indexed row sources. Reordering and duplicate indices are
#' supported; callback-computed row sources do not use it.
#' @param prior_list optional named prior list. If `NULL`, fitted priors are
#' used.
#' @param formula_target prediction target: `"conditional"`, `"fixed"`, or
#' `"marginal"`. The default includes fitted random-effect contributions.
#' @param blocks optional random-effect block names.
#' @param new_levels new-level policy for conditional or marginal random-effect
#' prediction. Use a `random_new_levels()` object or one of `"error"`, `"zero"`,
#' or `"sample"`.
#' @param marginal_method marginal representation used only with
#' `formula_target = "marginal"`: `"covariance"` returns fixed means plus
#' marginal random-effect covariance; `"sample"` returns materialized
#' random-effect draws. For known group covariance, `"sample"` jointly draws
#' the requested fitted levels from the corresponding covariance submatrix;
#' new levels remain unsupported.
#' @param seed optional random seed used when random-effect draws are simulated.
#' @param components whether to include component matrices where available.
#'
#' @return A list of class `BayesTools_formula_prediction` with fields `value`,
#' `mean`, `random`, `vcov`, `components`, and `metadata`.
#'
#' @seealso [JAGS_evaluate_formula()] [random_effects_marginal_vcov()]
#' @export
JAGS_predict_formula <- function(fit, parameter, formula = NULL, data = NULL,
                                 prior_list = NULL,
                                 formula_target = c("conditional", "fixed", "marginal"),
                                 blocks = NULL, new_levels = NULL,
                                 marginal_method = c("covariance", "sample"),
                                 seed = NULL, components = FALSE,
                                 fitted_rows = NULL){

  check_char(parameter, "parameter", allow_NA = FALSE)
  .bt_check_jags_node_name(parameter, "parameter")
  if(!is.null(fitted_rows) && is.null(data)){
    stop("'fitted_rows' can be supplied only with 'data'.", call. = FALSE)
  }
  if(!is.null(fitted_rows) && is.data.frame(data)){
    check_int(
      fitted_rows,
      "fitted_rows",
      lower = 1L,
      check_length = nrow(data),
      allow_NA = FALSE
    )
  }
  marginal_method_supplied <- !missing(marginal_method)
  formula_target <- match.arg(formula_target)
  marginal_method <- match.arg(marginal_method)
  if(!is.null(blocks)){
    check_char(blocks, "blocks", check_length = 0, allow_NA = FALSE)
    if(anyDuplicated(blocks)){
      stop("'blocks' must be unique.", call. = FALSE)
    }
  }
  if(!is.null(blocks) && identical(formula_target, "fixed")){
    stop(
      "'blocks' can be used only with formula_target = \"conditional\" or \"marginal\".",
      call. = FALSE
    )
  }
  if(!is.null(new_levels) && identical(formula_target, "fixed")){
    stop(
      "'new_levels' can be used only with formula_target = \"conditional\" or \"marginal\".",
      call. = FALSE
    )
  }
  if(isTRUE(marginal_method_supplied) &&
     !identical(formula_target, "marginal")){
    stop(
      "'marginal_method' can be used only with formula_target = \"marginal\".",
      call. = FALSE
    )
  }
  if(!is.null(new_levels)){
    new_levels <- .bt_random_new_levels_resolve(new_levels)
  }
  check_bool(components, "components", allow_NA = FALSE)
  if(!is.null(seed)){
    check_int(
      seed,
      "seed",
      lower = 0,
      upper = .Machine$integer.max,
      check_length = 1,
      allow_NA = FALSE
    )
  }

  fixed <- JAGS_evaluate_formula(
    fit = fit,
    formula = formula,
    parameter = parameter,
    data = data,
    prior_list = prior_list,
    formula_target = "fixed"
  )
  if(identical(formula_target, "fixed")){
    return(.bt_formula_prediction_object(
      value = fixed,
      mean = fixed,
      random = NULL,
      vcov = NULL,
      components = NULL,
      metadata = .bt_formula_prediction_metadata(
        parameter = parameter,
        formula_target = formula_target,
        marginal_method = NA_character_,
        blocks = blocks,
        new_levels = new_levels
      )
    ))
  }

  if(identical(formula_target, "conditional")){
    seed_state <- .bt_formula_prediction_seed(seed)
    on.exit(.bt_formula_prediction_restore_seed(seed_state), add = TRUE)
    value <- JAGS_evaluate_formula(
      fit = fit,
      formula = formula,
      parameter = parameter,
      data = data,
      prior_list = prior_list,
      formula_target = "conditional",
      blocks = blocks,
      new_levels = new_levels,
      fitted_rows = fitted_rows
    )
    random <- value - fixed
    return(.bt_formula_prediction_object(
      value = value,
      mean = fixed,
      random = random,
      vcov = NULL,
      components = if(isTRUE(components)) list(fixed = fixed, random = random) else NULL,
      metadata = .bt_formula_prediction_metadata(
        parameter = parameter,
        formula_target = formula_target,
        marginal_method = NA_character_,
        blocks = blocks,
        new_levels = new_levels
      )
    ))
  }

  posterior <- .bt_random_effect_marginal_covariance_posterior(
    fit = fit,
    posterior_samples = NULL
  )
  design <- .bt_random_effect_marginal_covariance_design(
    fit = fit,
    parameter = parameter
  )
  resolved_inputs <- .bt_JAGS_evaluate_formula_resolve_inputs(
    fit = fit,
    formula = formula,
    parameter = parameter,
    data = data,
    prior_list = prior_list,
    fitted_design = design
  )
  prior_list <- resolved_inputs$prior_list

  if(identical(marginal_method, "covariance")){
    vcov <- random_effects_marginal_vcov(
      fit = fit,
      parameter = design$parameter,
      data = data,
      posterior_samples = posterior,
      prior_list = prior_list,
      blocks = blocks,
      new_levels = new_levels,
      fitted_rows = fitted_rows
    )
    return(.bt_formula_prediction_object(
      value = fixed,
      mean = fixed,
      random = NULL,
      vcov = vcov,
      components = if(isTRUE(components)) list(fixed = fixed) else NULL,
      metadata = .bt_formula_prediction_metadata(
        parameter = parameter,
        formula_target = formula_target,
        marginal_method = marginal_method,
        blocks = blocks,
        new_levels = new_levels
      )
    ))
  }

  seed_state <- .bt_formula_prediction_seed(seed)
  on.exit(.bt_formula_prediction_restore_seed(seed_state), add = TRUE)
  random <- .bt_random_effects_marginal_sample(
    design = design,
    posterior = posterior,
    prior_list = prior_list,
    data = data,
    blocks = blocks,
    new_levels = new_levels,
    fitted_rows = fitted_rows
  )
  value <- fixed + random
  .bt_formula_prediction_object(
    value = value,
    mean = fixed,
    random = random,
    vcov = NULL,
    components = if(isTRUE(components)) list(fixed = fixed, random = random) else NULL,
    metadata = .bt_formula_prediction_metadata(
      parameter = parameter,
      formula_target = formula_target,
      marginal_method = marginal_method,
      blocks = blocks,
      new_levels = new_levels
    )
  )
}

.bt_formula_prediction_object <- function(value, mean, random, vcov,
                                          components, metadata){

  out <- list(
    value = value,
    mean = mean,
    random = random,
    vcov = vcov,
    components = components,
    metadata = metadata
  )
  class(out) <- c("BayesTools_formula_prediction", "list")

  out
}

.bt_formula_prediction_metadata <- function(parameter, formula_target,
                                            marginal_method, blocks,
                                            new_levels){

  list(
    parameter = parameter,
    formula_target = formula_target,
    marginal_method = marginal_method,
    blocks = blocks,
    new_levels = new_levels
  )
}

.bt_formula_prediction_seed <- function(seed){

  if(is.null(seed)){
    return(NULL)
  }
  state <- list(
    exists = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    value = if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }else{
      NULL
    }
  )
  set.seed(seed)
  state
}

.bt_formula_prediction_restore_seed <- function(state){

  if(is.null(state)){
    return(invisible(NULL))
  }
  if(isTRUE(state$exists)){
    assign(".Random.seed", state$value, envir = .GlobalEnv)
  }else if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
    rm(".Random.seed", envir = .GlobalEnv)
  }

  invisible(NULL)
}

.bt_random_effects_marginal_sample <- function(design, posterior, prior_list,
                                               data = NULL, blocks = NULL,
                                               new_levels = NULL,
                                               fitted_rows = NULL){

  selected <- .bt_random_effect_marginal_covariance_terms(
    design = design,
    blocks = blocks
  )
  output <- NULL
  for(random_term in selected$terms){
    block_new_levels <- .bt_random_effect_new_levels_policy(
      random_term = random_term,
      override = new_levels
    )
    block_data <- .bt_random_effect_marginal_covariance_block_data(
      design = design,
      random_term = random_term,
      data = data
    )
    new_row <- block_data$group_map > random_term$n_groups
    prediction_rows <- if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
      .bt_random_effect_prediction_fitted_rows(
        random_term = random_term,
        n_rows = nrow(block_data$model_matrix),
        data_supplied = block_data$data_supplied,
        fitted_rows = fitted_rows,
        new_row = new_row,
        context = "Marginal random-effect prediction"
      )
    }else{
      NULL
    }
    if(any(new_row) && !isTRUE(block_new_levels$allow)){
      new_groups <- block_data$group_levels[unique(block_data$group_map[new_row])]
      stop(
        "New random-effect level(s) for block '", random_term$block_name,
        "' require an explicit new-level policy: ",
        paste(new_groups, collapse = ", "),
        ". Use new_levels = \"zero\" or new_levels = \"sample\".",
        call. = FALSE
      )
    }
    rows <- seq_len(nrow(block_data$model_matrix))
    if(any(new_row) && identical(block_new_levels$method, "zero")){
      rows <- rows[!new_row]
    }
    contribution <- .bt_random_effect_group_contribution_sample(
      random_term = random_term,
      model_matrix = block_data$model_matrix,
      group_map = block_data$group_map,
      rows = rows,
      posterior = posterior,
      prior_list = prior_list,
      source_data = block_data$source_data,
      prediction_rows = prediction_rows
    )
    if(is.null(output)){
      output <- contribution
    }else{
      if(!identical(dim(output), dim(contribution))){
        stop(
          "Random-effect marginal sample block '", random_term$block_name,
          "' produced incompatible dimensions.",
          call. = FALSE
        )
      }
      output <- output + contribution
    }
  }

  output
}
