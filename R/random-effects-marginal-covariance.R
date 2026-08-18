#' Random-effect marginal variance-covariance samples
#'
#' @description
#' Constructs posterior-draw marginal covariance matrices implied by
#' BayesTools formula random effects. For each included random-effect block and
#' posterior draw, the helper builds the observation-level contribution
#' \eqn{Z_b G_b Z_b'} and sums contributions across blocks.
#'
#' @param fit a fitted object returned by [JAGS_fit()], an object carrying
#'   `formula_design` metadata, or a `BayesTools_formula_design` object.
#' @param parameter formula parameter name. Required when `fit` stores more than
#'   one formula design.
#' @param data optional data frame defining the observation rows. When `NULL`,
#'   the fitted formula design rows are used. New random-effect grouping levels
#'   require an explicit `new_levels` policy unless allowed by block metadata.
#' @param fitted_rows optional integer vector mapping supplied prediction rows
#'   to fitted observation indices. It is required whenever `data` is supplied
#'   and a selected block uses a posterior-indexed row source. Reordering and
#'   duplicate indices are supported. Callback-computed row sources do not use
#'   this mapping.
#' @param posterior_samples optional posterior sample matrix, data frame,
#'   `mcmc`, or `mcmc.list`. When `NULL`, samples are extracted from `fit`.
#' @param prior_list optional named prior list used to materialize fixed point
#'   hyperparameters. When `NULL`, priors are taken from `fit` or the formula
#'   design.
#' @param blocks optional character vector of random-effect block names to
#'   include. The default includes all formula random-effect blocks.
#' @param new_levels optional new-level policy. Use a `random_new_levels()`
#'   object or one of `"error"`, `"zero"`, or `"sample"`. `"zero"` assigns
#'   zero covariance contribution to unseen grouping levels. `"sample"` treats
#'   unseen grouping levels as independent draws from the block covariance `G`.
#' @param diagonal_only whether to return only observation-level marginal
#'   variances. This avoids allocating a `draw x row x row` array when callers
#'   need only `diag(Z G Z')`.
#' @param ... reserved for future extensions. Unused arguments are rejected.
#'
#' @return A list of class
#'   `BayesTools_random_effects_marginal_vcov` with fields:
#'   \describe{
#'     \item{samples}{A dense array with dimensions `draw x row x row`, or a
#'       `draw x row` matrix when `diagonal_only = TRUE`.}
#'     \item{metadata}{A list describing row ordering, included blocks,
#'       skipped blocks, structures, dimensions, and memory size. With
#'       `diagonal_only = TRUE`, it additionally identifies the `"diagonal"`
#'       variance representation, its sample dimensions, and the equivalent
#'       dense entry count. The default covariance metadata schema is
#'       unchanged.}
#'   }
#'
#' @seealso [JAGS_formula_design()] [JAGS_fit()]
#' @export
random_effects_marginal_vcov <- function(
    fit, parameter = NULL, data = NULL, posterior_samples = NULL,
    prior_list = NULL, blocks = NULL, new_levels = NULL,
    diagonal_only = FALSE, fitted_rows = NULL, ...){

  dots <- list(...)
  if(length(dots) > 0L){
    dot_names <- names(dots)
    if(is.null(dot_names)){
      dot_names <- rep("", length(dots))
    }
    unnamed <- !nzchar(dot_names)
    dot_names[unnamed] <- paste0("argument ", which(unnamed))
    stop(
      "Unused argument(s): ",
      paste(dot_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  check_char(parameter, "parameter", allow_NULL = TRUE, allow_NA = FALSE)
  check_char(blocks, "blocks", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_bool(diagonal_only, "diagonal_only", allow_NA = FALSE)
  if(!is.null(data) && !is.data.frame(data)){
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  if(!is.null(fitted_rows) && is.null(data)){
    stop("'fitted_rows' can be supplied only with 'data'.", call. = FALSE)
  }
  if(!is.null(fitted_rows)){
    check_int(
      fitted_rows,
      "fitted_rows",
      lower = 1L,
      check_length = nrow(data),
      allow_NA = FALSE
    )
  }
  if(!is.null(new_levels)){
    new_levels <- .bt_random_new_levels_resolve(new_levels)
  }

  design <- .bt_random_effect_marginal_covariance_design(
    fit = fit,
    parameter = parameter
  )
  posterior <- .bt_random_effect_marginal_covariance_posterior(
    fit = fit,
    posterior_samples = posterior_samples
  )
  prior_list <- .bt_random_effect_marginal_covariance_prior_list(
    prior_list = prior_list,
    fit = fit,
    design = design
  )

  .bt_random_effect_marginal_covariance_samples(
    design = design,
    posterior = posterior,
    prior_list = prior_list,
    data = data,
    blocks = blocks,
    new_levels = new_levels,
    fitted_rows = fitted_rows,
    diagonal_only = diagonal_only
  )
}


#' Random-effect marginal covariance factor states
#'
#' @description
#' Compiles the same structural covariance representation used by bridge
#' sampling and evaluates its draw-varying factor states for fitted rows. This
#' avoids constructing dense `draw x row x row` arrays when a downstream
#' likelihood can consume covariance factors directly.
#'
#' @param fit,parameter,posterior_samples,prior_list,blocks See
#'   [random_effects_marginal_vcov()].
#' @param row_blocks a list of integer row-index vectors that partitions the
#'   fitted rows without separating any structurally nonzero covariance
#'   contribution.
#' @param ... reserved for future extensions. Unused arguments are rejected.
#'
#' @return A list of class
#'   `BayesTools_random_effects_marginal_factor_states` with invariant
#'   `factor_plans`, one `factor_states` list per posterior draw, `row_blocks`,
#'   and structural `metadata`.
#'
#' @seealso [random_effects_marginal_vcov()] [JAGS_bridgesampling()]
#' @export
random_effects_marginal_factor_states <- function(
    fit, parameter = NULL, posterior_samples = NULL, prior_list = NULL,
    blocks = NULL, row_blocks, ...){

  dots <- list(...)
  if(length(dots) > 0L){
    dot_names <- names(dots)
    if(is.null(dot_names)){
      dot_names <- rep("", length(dots))
    }
    unnamed <- !nzchar(dot_names)
    dot_names[unnamed] <- paste0("argument ", which(unnamed))
    stop(
      "Unused argument(s): ",
      paste(dot_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  check_char(parameter, "parameter", allow_NULL = TRUE, allow_NA = FALSE)
  check_char(blocks, "blocks", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)

  design <- .bt_random_effect_marginal_covariance_design(
    fit = fit,
    parameter = parameter
  )
  posterior <- .bt_random_effect_marginal_covariance_posterior(
    fit = fit,
    posterior_samples = posterior_samples
  )
  prior_list <- .bt_random_effect_marginal_covariance_prior_list(
    prior_list = prior_list,
    fit = fit,
    design = design
  )
  selected <- .bt_random_effect_marginal_covariance_terms(
    design = design,
    blocks = blocks
  )
  random_effects <- selected$terms
  block_names <- vapply(random_effects, `[[`, character(1), "block_name")
  parameter_name <- design$parameter
  if(!is.character(parameter_name) || length(parameter_name) != 1L ||
     is.na(parameter_name) || !nzchar(parameter_name)){
    stop("Formula design parameter name is unavailable.", call. = FALSE)
  }

  row_blocks <- .bt_JAGS_bridge_marginal_random_row_blocks(
    row_blocks = row_blocks,
    n_rows = nrow(random_effects[[1L]]$model_matrix),
    parameter = parameter_name
  )
  if(is.null(row_blocks)){
    stop("'row_blocks' must be supplied.", call. = FALSE)
  }
  .bt_JAGS_bridge_validate_marginal_random_row_blocks(
    random_effects = random_effects,
    row_blocks = row_blocks,
    parameter = parameter_name
  )

  formula_design_list <- stats::setNames(list(design), parameter_name)
  formula_data_list <- stats::setNames(list(NULL), parameter_name)
  formula_prior_list <- stats::setNames(list(prior_list), parameter_name)
  marginal_random_spec <- stats::setNames(list(list(
    blocks = block_names,
    row_blocks = row_blocks,
    factor_state = TRUE
  )), parameter_name)
  evaluator <- .bt_JAGS_bridge_compile_marginal_random_evaluator(
    formula_design_list = formula_design_list,
    marginal_random_spec = marginal_random_spec,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = NULL
  )
  batch <- evaluator$factor_states(posterior)
  if(!is.null(batch)){
    batch <- batch[[parameter_name]]
    factor_plans  <- batch$factor_plans
    factor_states <- batch$factor_states
  }else{
    values <- lapply(seq_len(nrow(posterior)), function(row){
      samples <- posterior[row, ]
      formula_prior_parameters <- list()
      for(prior_name in names(prior_list)){
        prior <- prior_list[[prior_name]]
        if(is.prior.simplex(prior) &&
           identical(prior$distribution, "dirichlet")){
          K <- prior$parameters[["K"]]
          weights <- paste0(prior_name, "[", seq_len(K), "]")
          eta <- paste0(
            .JAGS_prior_dirichlet_eta_name(prior_name),
            "[", seq_len(K), "]"
          )
          if(all(weights %in% names(samples))){
            formula_prior_parameters[[prior_name]] <- unname(samples[weights])
          }else if(all(eta %in% names(samples))){
            eta_values <- unname(samples[eta])
            eta_sum    <- sum(eta_values)
            if(any(!is.finite(eta_values)) || any(eta_values < 0) ||
               !is.finite(eta_sum) || eta_sum <= 0){
              stop(
                "Dirichlet auxiliary coordinates for '", prior_name,
                "' must be finite, non-negative, and have a positive sum.",
                call. = FALSE
              )
            }
            formula_prior_parameters[[prior_name]] <-
              eta_values / eta_sum
          }
          if(all(weights %in% names(samples)) && all(eta %in% names(samples))){
            samples <- samples[setdiff(names(samples), eta)]
          }
          next
        }
        if(is.prior.point(prior)){
          formula_prior_parameters <- c(
            formula_prior_parameters,
            JAGS_marglik_parameters(samples, prior_list[prior_name])
          )
        }
      }
      evaluator$covariance(
        samples = samples,
        prior_parameters = list(),
        formula_prior_parameters = formula_prior_parameters,
        formula_parameters = list(),
        factor_covariance = FALSE,
        factor_state = TRUE
      )[[parameter_name]]
    })
    factor_plans  <- values[[1L]]$factor_plans
    factor_states <- lapply(values, `[[`, "factor_states")
  }
  structures <- stats::setNames(vapply(
    random_effects,
    .bt_random_effect_structure,
    character(1),
    context = "Random-effect marginal factor states"
  ), block_names)

  out <- list(
    factor_plans = factor_plans,
    factor_states = factor_states,
    row_blocks = row_blocks,
    metadata = list(
      parameter = parameter_name,
      n_draws = nrow(posterior),
      n_rows = nrow(random_effects[[1L]]$model_matrix),
      included_blocks = block_names,
      skipped_blocks = selected$skipped,
      structures = structures,
      representation = "factor_state"
    )
  )
  class(out) <- c(
    "BayesTools_random_effects_marginal_factor_states",
    "list"
  )
  out
}

#' Random-effect marginal variance factors
#'
#' @description
#' Extracts row-aligned design multipliers for one-column formula random-effect
#' blocks that are intended to be marginalized by a downstream likelihood. For
#' a block with SD parameter \eqn{\tau}, the marginal variance contribution is
#' \eqn{\tau^2} times `row_multiplier`. With known group covariance this factor
#' is built from the prepared group kernel, so the dense row-space covariance
#' factor is \eqn{Z K Z'}.
#'
#' @param formula_design a `BayesTools_formula_design` object returned by
#'   [JAGS_formula()] or [JAGS_formula_design()].
#' @param blocks optional character vector of random-effect block names. When
#'   `NULL`, only marginalized random-effect blocks are selected.
#' @param require_diagonal whether to require the row-space covariance factor to
#'   be diagonal. This protects downstream likelihoods that can only consume
#'   row-wise variance multipliers.
#' @param require_one_to_one whether each observation row must map to a distinct
#'   grouping level.
#'
#' @return A list of class
#'   `BayesTools_random_effects_marginal_variance_factors` with selected block
#'   metadata. Each block contains `row_multiplier`, `group_map`,
#'   `group_levels`, `model_matrix`, SD metadata, compile mode, and known group
#'   covariance metadata when present.
#'
#' @details Unsupported factor representations signal an error inheriting from
#'   `BayesTools_random_effects_marginal_variance_unavailable`. The condition
#'   carries stable `block_name` and `reason` fields so downstream packages can
#'   distinguish structural unavailability from invalid inputs. `reason` is one
#'   of `"row_indexed_sd"`, `"multiple_columns"`,
#'   `"non_diagonal_row_covariance"`, or `"repeated_groups"`.
#'
#' @seealso [random_effects_marginal_vcov()] [random_effects_compile()]
#' @export
random_effects_marginal_variance_factors <- function(
    formula_design,
    blocks = NULL,
    require_diagonal = TRUE,
    require_one_to_one = FALSE){

  if(!inherits(formula_design, "BayesTools_formula_design")){
    stop(
      "'formula_design' must be a BayesTools_formula_design object.",
      call. = FALSE
    )
  }
  check_char(blocks, "blocks", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_bool(require_diagonal, "require_diagonal", allow_NA = FALSE)
  check_bool(require_one_to_one, "require_one_to_one", allow_NA = FALSE)

  selected <- .bt_random_effect_marginal_variance_factor_terms(
    design = formula_design,
    blocks = blocks
  )
  random_terms <- selected$terms
  if(length(random_terms) == 0L){
    stop(
      if(is.null(blocks)){
        "No marginalized random-effect blocks were selected."
      }else{
        "No random-effect blocks were selected."
      },
      call. = FALSE
    )
  }

  block_metadata <- lapply(
    random_terms,
    .bt_random_effect_marginal_variance_factor_block,
    design = formula_design,
    require_diagonal = require_diagonal,
    require_one_to_one = require_one_to_one
  )
  names(block_metadata) <- vapply(random_terms, `[[`, character(1), "block_name")

  row_counts <- vapply(block_metadata, `[[`, integer(1), "n_rows")
  if(length(unique(row_counts)) != 1L){
    stop(
      "Selected random-effect blocks do not have a common row count.",
      call. = FALSE
    )
  }
  row_names <- lapply(block_metadata, `[[`, "row_names")
  common_row_names <- vapply(
    row_names,
    identical,
    logical(1),
    y = row_names[[1L]]
  )
  if(!all(common_row_names)){
    stop(
      "Selected random-effect blocks do not have common row names.",
      call. = FALSE
    )
  }

  out <- list(
    parameter = formula_design$parameter,
    n_rows = row_counts[[1L]],
    row_names = row_names[[1L]],
    included_blocks = names(block_metadata),
    skipped_blocks = selected$skipped,
    require_diagonal = require_diagonal,
    require_one_to_one = require_one_to_one,
    blocks = block_metadata
  )
  class(out) <- c(
    "BayesTools_random_effects_marginal_variance_factors",
    "list"
  )
  out
}

.bt_random_effect_marginal_variance_factor_terms <- function(design,
                                                             blocks = NULL){

  random_effects <- .bt_formula_design_random_effects(design)
  if(length(random_effects) == 0L){
    stop(
      "Formula design for parameter '", design$parameter,
      "' has no random-effect blocks.",
      call. = FALSE
    )
  }

  block_names <- vapply(random_effects, `[[`, character(1), "block_name")
  if(anyDuplicated(block_names)){
    stop("Formula random-effect block names must be unique.", call. = FALSE)
  }

  if(is.null(blocks)){
    modes <- .bt_random_effects_compile_modes_from_terms(random_effects)
    selected <- modes == "marginalized"
    skipped <- data.frame(
      block_name = block_names[!selected],
      reason = rep("not marginalized", sum(!selected)),
      stringsAsFactors = FALSE
    )
    return(list(terms = random_effects[selected], skipped = skipped))
  }
  if(anyDuplicated(blocks)){
    stop("'blocks' must be unique.", call. = FALSE)
  }

  unknown <- setdiff(blocks, block_names)
  if(length(unknown) > 0L){
    stop(
      "Unknown random-effect block(s): ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  selected <- random_effects[match(blocks, block_names)]
  skipped_names <- setdiff(block_names, blocks)
  skipped <- data.frame(
    block_name = skipped_names,
    reason = rep("not requested", length(skipped_names)),
    stringsAsFactors = FALSE
  )

  list(terms = selected, skipped = skipped)
}

.bt_random_effect_marginal_variance_factor_block <- function(random_term,
                                                             design,
                                                             require_diagonal,
                                                             require_one_to_one){

  block_data <- .bt_random_effect_marginal_covariance_block_data(
    design = design,
    random_term = random_term,
    data = NULL
  )

  if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    .bt_random_effect_marginal_variance_unavailable_stop(
      random_term = random_term,
      reason = "row_indexed_sd",
      "Random-effect marginal variance factors for block '",
      random_term$block_name,
      "' do not support row-indexed external SD sources."
    )
  }
  if(random_term$n_columns != 1L || ncol(block_data$model_matrix) != 1L){
    .bt_random_effect_marginal_variance_unavailable_stop(
      random_term = random_term,
      reason = "multiple_columns",
      "Random-effect marginal variance factors for block '",
      random_term$block_name,
      "' require one random-effect column."
    )
  }

  row_covariance <- .bt_random_effect_marginal_variance_base_covariance(
    random_term = random_term,
    model_matrix = block_data$model_matrix,
    group_map = block_data$group_map
  )
  dimnames(row_covariance) <- list(block_data$row_names, block_data$row_names)
  row_status <- .bt_random_effect_row_covariance_diagonal_status(row_covariance)
  one_to_one <- !anyDuplicated(block_data$group_map)

  if(isTRUE(require_diagonal) && !isTRUE(row_status$is_diagonal)){
    .bt_random_effect_marginal_variance_unavailable_stop(
      random_term = random_term,
      reason = "non_diagonal_row_covariance",
      "Random-effect marginal variance factors for block '",
      random_term$block_name,
      "' require diagonal row-space covariance, but the block has off-diagonal row covariance."
    )
  }
  if(isTRUE(require_one_to_one) && !isTRUE(one_to_one)){
    .bt_random_effect_marginal_variance_unavailable_stop(
      random_term = random_term,
      reason = "repeated_groups",
      "Random-effect marginal variance factors for block '",
      random_term$block_name,
      "' require a one-to-one row-to-group mapping, but grouping levels repeat."
    )
  }

  row_multiplier <- diag(row_covariance)
  names(row_multiplier) <- block_data$row_names
  column_names <- colnames(block_data$model_matrix)
  if(is.null(column_names)){
    column_names <- paste0("column", seq_len(ncol(block_data$model_matrix)))
  }

  list(
    block_name = random_term$block_name,
    grouping = random_term$group_label,
    structure = .bt_random_effect_structure(
      random_term,
      context = "Random-effect marginal variance factor metadata"
    ),
    compile_mode = .bt_random_effect_term_compile_mode(random_term),
    n_groups = length(block_data$group_levels),
    fitted_n_groups = random_term$n_groups,
    n_columns = random_term$n_columns,
    n_rows = nrow(block_data$model_matrix),
    row_names = block_data$row_names,
    row_order = seq_len(nrow(block_data$model_matrix)),
    group_levels = block_data$group_levels,
    group_map = block_data$group_map,
    column_names = column_names,
    model_matrix = block_data$model_matrix,
    sd_parameter_names = random_term$sd_parameter_names,
    sd_binding = random_term$sd_binding,
    homogeneous_sd = random_term$homogeneous_sd,
    group_covariance = .bt_random_effect_group_covariance_metadata(random_term),
    row_multiplier = row_multiplier,
    row_covariance_diagonal = row_status$is_diagonal,
    max_off_diagonal = row_status$max_off_diagonal,
    diagonal_tolerance = row_status$tolerance,
    one_to_one = one_to_one
  )
}

# Signal that a valid random-effect block cannot provide row variance factors.
.bt_random_effect_marginal_variance_unavailable_stop <- function(
    random_term, reason, ...){

  condition <- errorCondition(
    message = paste0(...),
    call = NULL,
    class = "BayesTools_random_effects_marginal_variance_unavailable",
    block_name = random_term$block_name,
    reason = reason
  )
  stop(condition)
}

.bt_random_effect_marginal_variance_base_covariance <- function(random_term,
                                                                model_matrix,
                                                                group_map){

  if(.bt_random_effect_has_known_group_covariance(random_term)){
    group_covariance <- .bt_random_effect_known_group_covariance(
      random_term,
      context = "Random-effect marginal variance factor"
    )
    if(any(group_map > length(group_covariance$levels))){
      stop(
        "Random-effect marginal variance factors for block '",
        random_term$block_name,
        "' cannot include new levels with known group covariance.",
        call. = FALSE
      )
    }
    return(
      group_covariance$kernel[group_map, group_map, drop = FALSE] *
        tcrossprod(model_matrix[, 1L])
    )
  }

  same_group <- outer(as.integer(group_map), as.integer(group_map), "==")
  storage.mode(same_group) <- "double"
  same_group * tcrossprod(model_matrix[, 1L])
}

.bt_random_effect_row_covariance_diagonal_status <- function(row_covariance){

  off_diagonal <- row_covariance
  diag(off_diagonal) <- 0
  max_off_diagonal <- max(abs(off_diagonal))

  list(
    is_diagonal = isTRUE(all(off_diagonal == 0)),
    max_off_diagonal = max_off_diagonal,
    tolerance = 0
  )
}

.bt_random_effect_marginal_covariance_samples <- function(design, posterior,
                                                          prior_list,
                                                          data = NULL,
                                                          blocks = NULL,
                                                          new_levels = NULL,
                                                          fitted_rows = NULL,
                                                          diagonal_only = FALSE,
                                                          parameters = NULL){

  selected <- .bt_random_effect_marginal_covariance_terms(
    design = design,
    blocks = blocks
  )
  random_terms <- selected$terms
  if(length(random_terms) == 0L){
    stop("No random-effect blocks were selected.", call. = FALSE)
  }

  samples <- NULL
  block_metadata <- vector("list", length(random_terms))
  for(block_i in seq_along(random_terms)){
    random_term <- random_terms[[block_i]]
    block_new_levels <- .bt_random_effect_new_levels_policy(
      random_term = random_term,
      override = new_levels
    )
    block_data <- .bt_random_effect_marginal_covariance_block_data(
      design = design,
      random_term = random_term,
      data = data
    )
    .bt_random_effect_check_memory(
      estimate = .bt_random_effect_output_memory_estimate(
        operation = "random-effect marginal covariance",
        n_rows = nrow(block_data$model_matrix),
        n_draws = nrow(posterior),
        covariance = TRUE,
        diagonal_only = diagonal_only
      ),
      block_name = random_term$block_name,
      alternative = if(isTRUE(diagonal_only)){
        paste0(
          "Reduce the number of rows, posterior draws, or selected blocks. ",
          "Raise the option (or set it to Inf) only after verifying the ",
          "operation's memory budget."
        )
      }else{
        paste0(
          "Use diagonal_only = TRUE when only marginal variances are needed, ",
          "or reduce the number of rows, posterior draws, or selected blocks. ",
          "Raise the option (or set it to Inf) only after verifying the ",
          "operation's memory budget."
        )
      }
    )
    prediction_rows <- if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
      .bt_random_effect_prediction_fitted_rows(
        random_term = random_term,
        n_rows = nrow(block_data$model_matrix),
        data_supplied = block_data$data_supplied,
        fitted_rows = fitted_rows,
        new_row = block_data$group_map > random_term$n_groups,
        context = "Random-effect marginal covariance"
      )
    }else{
      NULL
    }
    new_level_info <- .bt_random_effect_marginal_covariance_new_level_info(
      random_term = random_term,
      block_data = block_data,
      new_levels = block_new_levels
    )
    block_output <- .bt_random_effect_marginal_covariance_block_samples(
      random_term = random_term,
      model_matrix = block_data$model_matrix,
      group_map = block_data$group_map,
      source_data = block_data$source_data,
      prediction_rows = prediction_rows,
      posterior = posterior,
      prior_list = prior_list,
      parameters = parameters,
      diagonal_only = diagonal_only
    )
    if(any(new_level_info$row_mask) &&
       identical(block_new_levels$method, "zero")){
      if(isTRUE(diagonal_only)){
        block_output$samples[, new_level_info$row_mask] <- 0
      }else{
        block_output$samples[, new_level_info$row_mask, ] <- 0
        block_output$samples[, , new_level_info$row_mask] <- 0
      }
    }
    if(is.null(samples)){
      samples <- block_output$samples
    }else{
      .bt_random_effect_marginal_covariance_check_block_output(
        random_term = random_term,
        samples = samples,
        n_draws = nrow(posterior),
        row_names = block_data$row_names,
        diagonal_only = diagonal_only
      )
      samples <- samples + block_output$samples
    }
    block_metadata[[block_i]] <- .bt_random_effect_marginal_covariance_block_metadata(
      random_term = random_term,
      model_matrix = block_data$model_matrix,
      group_map = block_data$group_map,
      group_levels = block_data$group_levels,
      row_names = block_data$row_names,
      sample_dim = block_output$sample_dim,
      new_level_info = new_level_info,
      diagonal_only = diagonal_only
    )
  }
  names(block_metadata) <- vapply(random_terms, `[[`, character(1), "block_name")

  row_names <- if(isTRUE(diagonal_only)){
    colnames(samples)
  }else{
    dimnames(samples)[[2L]]
  }
  sample_entries <- prod(dim(samples))
  metadata <- list(
    parameter = design$parameter,
    n_draws = nrow(samples),
    n_rows = if(isTRUE(diagonal_only)) ncol(samples) else dim(samples)[2L],
    row_names = row_names,
    row_order = seq_along(row_names),
    data_source = if(is.null(data)) "fitted" else "data",
    dense = !isTRUE(diagonal_only),
    dense_entries = if(isTRUE(diagonal_only)) NA_real_ else sample_entries,
    estimated_size_bytes = as.numeric(utils::object.size(samples)),
    potentially_expensive = sample_entries > 1e7,
    included_blocks = names(block_metadata),
    skipped_blocks = selected$skipped,
    structures = vapply(block_metadata, `[[`, character(1), "structure"),
    blocks = block_metadata
  )
  if(isTRUE(diagonal_only)){
    metadata$representation <- "diagonal"
    metadata$quantity <- "variance"
    metadata$diagonal_only <- TRUE
    metadata$sample_dim <- dim(samples)
    metadata$sample_entries <- sample_entries
    metadata$equivalent_dense_entries <-
      as.numeric(nrow(samples)) * as.numeric(ncol(samples))^2
  }

  out <- list(samples = samples, metadata = metadata)
  class(out) <- c(
    "BayesTools_random_effects_marginal_vcov",
    "list"
  )
  out
}

.bt_random_effect_marginal_covariance_check_block_output <- function(
    random_term,
    samples,
    n_draws,
    row_names,
    diagonal_only = FALSE){

  expected_dim <- if(isTRUE(diagonal_only)){
    c(n_draws, length(row_names))
  }else{
    c(n_draws, length(row_names), length(row_names))
  }
  if(!identical(dim(samples), expected_dim)){
    stop(
      "Random-effect marginal covariance block '",
      random_term$block_name,
      "' produced dimensions ",
      paste(expected_dim, collapse = " x "),
      ", but previous block contributions have dimensions ",
      paste(dim(samples), collapse = " x "),
      ".",
      call. = FALSE
    )
  }

  actual_row_names <- if(isTRUE(diagonal_only)){
    colnames(samples)
  }else{
    dimnames(samples)[[2L]]
  }
  if(!identical(actual_row_names, row_names)){
    stop(
      "Random-effect marginal covariance block '",
      random_term$block_name,
      "' produced row names that do not match previous block contributions.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_random_effect_marginal_covariance_design <- function(fit,
                                                         parameter = NULL){

  if(inherits(fit, "BayesTools_formula_design")){
    if(is.null(parameter)){
      return(fit)
    }
    if(!identical(fit$parameter, parameter)){
      stop(
        "Formula design for parameter '", parameter, "' was not found.",
        call. = FALSE
      )
    }
    return(fit)
  }

  formula_design <- NULL
  if(is.list(fit) && length(fit) > 0L &&
     all(vapply(fit, inherits, logical(1), what = "BayesTools_formula_design"))){
    formula_design <- fit
  }else{
    formula_design <- try(JAGS_formula_design(fit), silent = TRUE)
    if(inherits(formula_design, "try-error")){
      formula_design <- NULL
    }
  }

  if(is.null(formula_design)){
    stop(
      "Random-effect marginal covariance construction needs fitted formula design metadata. ",
      "Use a fit produced by JAGS_fit() with formula_list, or pass a BayesTools_formula_design object.",
      call. = FALSE
    )
  }

  if(!is.null(parameter)){
    if(!parameter %in% names(formula_design)){
      matches <- vapply(
        formula_design,
        function(design) identical(design$parameter, parameter),
        logical(1)
      )
      if(!any(matches)){
        stop(
          "Formula design for parameter '", parameter, "' was not found.",
          call. = FALSE
        )
      }
      return(formula_design[[which(matches)[1L]]])
    }
    return(formula_design[[parameter]])
  }

  if(length(formula_design) != 1L){
    stop(
      "'parameter' is required when more than one formula design is stored.",
      call. = FALSE
    )
  }

  formula_design[[1L]]
}

.bt_random_effect_marginal_covariance_posterior <- function(
    fit,
    posterior_samples = NULL){

  if(is.null(posterior_samples)){
    if(inherits(fit, "BayesTools_formula_design") ||
       (is.list(fit) && length(fit) > 0L &&
        all(vapply(fit, inherits, logical(1), what = "BayesTools_formula_design")))){
      stop(
        "'posterior_samples' must be supplied when 'fit' is formula design metadata.",
        call. = FALSE
      )
    }
    posterior <- as.matrix(.fit_to_posterior(fit))
  }else if(is.matrix(posterior_samples) || is.data.frame(posterior_samples)){
    posterior <- as.matrix(posterior_samples)
  }else{
    posterior <- as.matrix(.fit_to_posterior(posterior_samples))
  }

  .bt_random_effect_marginal_covariance_validate_posterior(posterior)
}

.bt_random_effect_marginal_covariance_validate_posterior <- function(
    posterior, allow_zero_columns = FALSE){

  if(!is.matrix(posterior) || length(dim(posterior)) != 2L){
    stop("'posterior_samples' must be a two-dimensional sample matrix.", call. = FALSE)
  }
  if(!is.numeric(posterior)){
    stop("'posterior_samples' must be numeric.", call. = FALSE)
  }
  if(nrow(posterior) < 1L){
    stop("'posterior_samples' must contain at least one draw.", call. = FALSE)
  }
  zero_columns_allowed <- isTRUE(allow_zero_columns) && ncol(posterior) == 0L
  if(!zero_columns_allowed){
    posterior_names <- colnames(posterior)
    if(is.null(posterior_names) || length(posterior_names) != ncol(posterior) ||
       anyNA(posterior_names) || any(!nzchar(posterior_names))){
      stop("'posterior_samples' must have non-empty column names.", call. = FALSE)
    }
    if(anyDuplicated(posterior_names)){
      duplicated_names <- unique(posterior_names[duplicated(posterior_names)])
      stop(
        "'posterior_samples' column names must be unique. Duplicated column(s): ",
        paste0("'", duplicated_names[seq_len(min(4L, length(duplicated_names)))], "'",
               collapse = ", "),
        if(length(duplicated_names) > 4L) ", ..." else "",
        ".",
        call. = FALSE
      )
    }
  }

  if(is.null(attr(
    posterior,
    "BayesTools_random_effect_dirichlet_draw_cache",
    exact = TRUE
  ))){
    attr(
      posterior,
      "BayesTools_random_effect_dirichlet_draw_cache"
    ) <- new.env(parent = emptyenv())
  }

  posterior
}

.bt_random_effect_marginal_covariance_prior_list <- function(prior_list,
                                                             fit,
                                                             design){

  if(!is.null(prior_list)){
    check_list(prior_list, "prior_list", allow_NULL = FALSE)
    return(prior_list)
  }

  fit_prior_list <- attr(fit, "prior_list", exact = TRUE)
  if(is.list(fit_prior_list) && length(fit_prior_list) > 0L){
    return(fit_prior_list)
  }
  if(is.list(design$prior_list)){
    return(design$prior_list)
  }

  list()
}

.bt_random_effect_marginal_covariance_terms <- function(design,
                                                        blocks = NULL){

  random_effects <- .bt_formula_design_random_effects(design)
  if(length(random_effects) == 0L){
    stop(
      "Formula design for parameter '", design$parameter,
      "' has no random-effect blocks.",
      call. = FALSE
    )
  }

  block_names <- vapply(random_effects, `[[`, character(1), "block_name")
  if(anyDuplicated(block_names)){
    stop("Formula random-effect block names must be unique.", call. = FALSE)
  }

  if(is.null(blocks)){
    skipped <- data.frame(
      block_name = character(),
      reason = character(),
      stringsAsFactors = FALSE
    )
    return(list(terms = random_effects, skipped = skipped))
  }
  if(anyDuplicated(blocks)){
    stop("'blocks' must be unique.", call. = FALSE)
  }

  unknown <- setdiff(blocks, block_names)
  if(length(unknown) > 0L){
    stop(
      "Unknown random-effect block(s): ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  selected <- random_effects[match(blocks, block_names)]
  skipped_names <- setdiff(block_names, blocks)
  skipped <- data.frame(
    block_name = skipped_names,
    reason = rep("not requested", length(skipped_names)),
    stringsAsFactors = FALSE
  )

  list(terms = selected, skipped = skipped)
}

.bt_random_effect_marginal_covariance_block_data <- function(design,
                                                             random_term,
                                                             data = NULL){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect marginal covariance metadata"
  )

  if(is.null(data)){
    model_matrix <- random_term$model_matrix
    group_map <- random_term$group_map
    group_levels <- random_term$group_levels
    source_data <- design$source_data
  }else{
    random_data <- if(structure %in% c("cs", "hcs", "ar1", "car", "har")){
      data
    }else{
      .bt_random_effect_marginal_covariance_scaled_data(
        design = design,
        data = data
      )
    }
    prediction <- .bt_random_effect_prediction_data(
      random_term = random_term,
      data = random_data,
      group_data = data,
      allow_new_groups = .bt_random_effect_allows_new_levels(random_term),
      context = "random_effects_marginal_vcov()"
    )
    model_matrix <- prediction$model_matrix
    group_map <- prediction$group_map
    group_levels <- prediction$group_levels
    source_data <- data
  }

  .bt_random_effect_marginal_covariance_validate_block_data(
    random_term = random_term,
    model_matrix = model_matrix,
    group_map = group_map,
    group_levels = group_levels
  )

  row_names <- rownames(model_matrix)
  if(is.null(row_names)){
    row_names <- as.character(seq_len(nrow(model_matrix)))
  }

  list(
    model_matrix = model_matrix,
    group_map = group_map,
    group_levels = group_levels,
    data_supplied = !is.null(data),
    source_data = source_data,
    row_names = row_names
  )
}

.bt_random_effect_marginal_covariance_scaled_data <- function(design, data){

  if(is.null(design$formula_scale)){
    return(data)
  }

  fit <- structure(list(), class = "BayesTools_fit")
  attr(fit, "formula_scale") <- stats::setNames(
    list(design$formula_scale),
    design$parameter
  )

  .bt_apply_formula_scale_to_data(
    fit = fit,
    parameter = design$parameter,
    data = data,
    predictors_type = design$predictor_types
  )
}

.bt_random_effect_marginal_covariance_validate_block_data <- function(
    random_term,
    model_matrix,
    group_map,
    group_levels = random_term$group_levels){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect marginal covariance metadata"
  )
  if(!structure %in% c("id", "diag", "us", "cs", "hcs", "ar1", "car", "har")){
    stop(
      "Random-effect marginal covariance for block '",
      random_term$block_name,
      "' does not support covariance structure '",
      structure,
      "'.",
      call. = FALSE
    )
  }
  if(!is.matrix(model_matrix) || !is.numeric(model_matrix)){
    stop(
      "Random-effect marginal covariance metadata for block '",
      random_term$block_name,
      "' are missing numeric 'random_term$model_matrix'.",
      call. = FALSE
    )
  }
  n_columns <- .bt_random_effect_marginal_covariance_int_metadata(
    random_term$n_columns,
    field = "random_term$n_columns",
    random_term = random_term
  )
  fitted_n_groups <- .bt_random_effect_marginal_covariance_int_metadata(
    random_term$n_groups,
    field = "random_term$n_groups",
    random_term = random_term
  )
  if(ncol(model_matrix) != n_columns){
    stop(
      "Random-effect marginal covariance metadata for block '",
      random_term$block_name,
      "' have inconsistent model-matrix and column counts.",
      call. = FALSE
    )
  }
  if(nrow(model_matrix) < 1L){
    stop(
      "Random-effect marginal covariance metadata for block '",
      random_term$block_name,
      "' must contain at least one observation row.",
      call. = FALSE
    )
  }
  if(!is.numeric(group_map) || length(group_map) != nrow(model_matrix) ||
     anyNA(group_map) || any(!is.finite(group_map)) ||
     any(group_map < 1L) || any(group_map > .Machine$integer.max) ||
     any(group_map != floor(group_map)) ||
     any(group_map > length(group_levels))){
    stop(
      "Random-effect marginal covariance metadata for block '",
      random_term$block_name,
      "' are missing valid 'random_term$group_map'.",
      call. = FALSE
    )
  }
  if(is.null(random_term$group_levels) ||
     length(random_term$group_levels) != fitted_n_groups){
    stop(
      "Random-effect marginal covariance metadata for block '",
      random_term$block_name,
      "' are missing canonical 'random_term$group_levels'.",
      call. = FALSE
    )
  }
  if(is.null(group_levels) || length(group_levels) < fitted_n_groups){
    stop(
      "Random-effect marginal covariance metadata for block '",
      random_term$block_name,
      "' are missing row-level group levels.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_random_effect_marginal_covariance_new_level_info <- function(
    random_term,
    block_data,
    new_levels){

  new_row <- block_data$group_map > random_term$n_groups
  new_group_index <- if(length(block_data$group_levels) > random_term$n_groups){
    seq.int(random_term$n_groups + 1L, length(block_data$group_levels))
  }else{
    integer()
  }
  new_group_levels <- block_data$group_levels[new_group_index]
  if(length(new_group_levels) > 0L && !isTRUE(new_levels$allow)){
    stop(
      "New random-effect level(s) for block '", random_term$block_name,
      "' require an explicit new-level policy: ",
      paste(new_group_levels, collapse = ", "),
      ". Use new_levels = \"zero\" or new_levels = \"sample\".",
      call. = FALSE
    )
  }
  if(length(new_group_levels) > 0L &&
     .bt_random_effect_has_known_group_covariance(random_term)){
    stop(
      "New random-effect level(s) for block '",
      random_term$block_name,
      "' are not supported with known group covariance: ",
      paste(new_group_levels, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  list(
    policy = new_levels,
    group_levels = new_group_levels,
    group_index = new_group_index,
    row_mask = new_row,
    rows = which(new_row)
  )
}

.bt_random_effect_marginal_covariance_int_metadata <- function(x, field,
                                                               random_term){

  if(is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
     x == as.integer(x) && x >= 1L){
    return(as.integer(x))
  }

  stop(
    "Random-effect marginal covariance metadata",
    .bt_random_effect_metadata_block_detail(random_term),
    " are missing canonical '",
    field,
    "'.",
    call. = FALSE
  )
}

.bt_random_effect_marginal_covariance_block_samples <- function(
    random_term,
    model_matrix,
    group_map,
    source_data,
    prediction_rows,
    posterior,
    prior_list,
    parameters = NULL,
    diagonal_only = FALSE){

  if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    return(.bt_random_effect_marginal_covariance_row_indexed_block(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      source_data = source_data,
      prediction_rows = prediction_rows,
      posterior = posterior,
      prior_list = prior_list,
      parameters = parameters,
      diagonal_only = diagonal_only
    ))
  }

  n_columns <- ncol(model_matrix)
  sd_draws <- .bt_random_effect_sd_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior,
    prior_list = prior_list
  )
  if(is.null(sd_draws)){
    .bt_random_effect_marginal_covariance_missing_sd_stop(
      random_term = random_term,
      n_columns = n_columns
    )
  }
  .bt_random_effect_marginal_covariance_validate_draw_matrix(
    draws = sd_draws,
    n_draws = nrow(posterior),
    n_columns = n_columns,
    label = "SD",
    random_term = random_term,
    nonnegative = TRUE
  )

  if(.bt_random_effect_has_known_group_covariance(random_term)){
    return(.bt_random_effect_marginal_covariance_known_group_block(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      sd_draws = sd_draws,
      diagonal_only = diagonal_only
    ))
  }

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect marginal covariance metadata"
  )
  if(structure %in% c("diag", "id")){
    return(.bt_random_effect_marginal_covariance_independent(
      model_matrix = model_matrix,
      group_map = group_map,
      n_draws = nrow(posterior),
      column_scale_draws = sd_draws,
      diagonal_only = diagonal_only
    ))
  }
  if(identical(structure, "car")){
    return(.bt_random_effect_marginal_covariance_car(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      row_weights = NULL,
      column_weights = sd_draws,
      diagonal_only = diagonal_only
    ))
  }

  row_column <- .bt_random_effect_marginal_covariance_structured_columns(
    random_term = random_term,
    model_matrix = model_matrix
  )
  if(isTRUE(diagonal_only)){
    one_sparse <- .bt_random_effect_marginal_covariance_one_sparse_columns(
      model_matrix
    )
    if(!is.null(one_sparse)){
      row_sd <- sd_draws[, one_sparse$row_column, drop = FALSE] *
        matrix(
          one_sparse$row_value,
          nrow = nrow(posterior),
          ncol = nrow(model_matrix),
          byrow = TRUE
        )
      return(.bt_random_effect_marginal_variance_one_sparse(
        model_matrix = model_matrix,
        row_sd = row_sd
      ))
    }
  }
  if(!is.null(row_column)){
    return(.bt_random_effect_marginal_covariance_structured_one_hot(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      row_column = row_column,
      row_sd = sd_draws[, row_column, drop = FALSE],
      diagonal_only = diagonal_only
    ))
  }

  if(isTRUE(diagonal_only)){
    cholesky <- .bt_random_effect_cholesky_draws(
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior
    )
    if(is.null(cholesky)){
      .bt_random_effect_marginal_covariance_missing_correlation_stop(
        random_term = random_term,
        n_columns = n_columns,
        posterior = posterior
      )
    }
    return(.bt_random_effect_marginal_variance_cholesky(
      random_term = random_term,
      model_matrix = model_matrix,
      cholesky = cholesky,
      column_weights = sd_draws
    ))
  }

  correlation <- .bt_random_effect_marginal_covariance_correlation_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )

  list(
    samples = .bt_random_effect_marginal_covariance_expand(
      model_matrix = model_matrix,
      group_map = group_map,
      block_covariance = correlation,
      column_weights = sd_draws
    ),
    sample_dim = c(nrow(posterior), nrow(model_matrix), nrow(model_matrix))
  )
}

.bt_random_effect_marginal_covariance_known_group_block <- function(
    random_term,
    model_matrix,
    group_map,
    posterior,
    sd_draws,
    diagonal_only = FALSE){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  group_covariance <- .bt_random_effect_known_group_covariance(
    random_term,
    context = "Random-effect marginal covariance"
  )
  if(any(group_map > length(group_covariance$levels))){
    stop(
      "Random-effect marginal covariance for block '",
      random_term$block_name,
      "' cannot include new levels with known group covariance.",
      call. = FALSE
    )
  }

  row_names <- rownames(model_matrix)
  if(is.null(row_names)){
    row_names <- as.character(seq_len(n_rows))
  }
  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect marginal covariance metadata"
  )
  correlation <- if(structure %in% c("diag", "id") || n_columns == 1L){
    NULL
  }else{
    .bt_random_effect_marginal_covariance_correlation_draws(
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior
    )
  }
  group_kernel <- group_covariance$kernel[group_map, group_map, drop = FALSE]
  if(isTRUE(diagonal_only)){
    out <- matrix(
      NA_real_,
      nrow = n_draws,
      ncol = n_rows,
      dimnames = list(draw = NULL, row = row_names)
    )
  }else{
    out <- array(
      NA_real_,
      dim = c(n_draws, n_rows, n_rows),
      dimnames = list(draw = NULL, row = row_names, column = row_names)
    )
  }
  for(draw in seq_len(n_draws)){
    coefficient_covariance <- if(is.null(correlation)){
      diag(sd_draws[draw, ]^2, nrow = n_columns, ncol = n_columns)
    }else{
      correlation[draw, , ] * tcrossprod(sd_draws[draw, ])
    }
    design_covariance <- tcrossprod(
      model_matrix %*% coefficient_covariance,
      model_matrix
    )
    covariance <- group_kernel * design_covariance
    if(isTRUE(diagonal_only)){
      out[draw, ] <- diag(covariance)
    }else{
      out[draw, , ] <- covariance
    }
  }

  list(
    samples = out,
    sample_dim = if(isTRUE(diagonal_only)){
      c(n_draws, n_rows)
    }else{
      c(n_draws, n_rows, n_rows)
    }
  )
}

.bt_random_effect_marginal_covariance_row_indexed_block <- function(
    random_term,
    model_matrix,
    group_map,
    source_data,
    prediction_rows,
    posterior,
    prior_list,
    parameters = NULL,
    diagonal_only = FALSE){

  n_draws <- nrow(posterior)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)

  source_draws <- .bt_random_effect_row_indexed_source_draws(
    random_term = random_term,
    n_rows = n_rows,
    posterior = posterior,
    data = source_data,
    parameters = parameters,
    prediction_rows = prediction_rows,
    context = "Random-effect marginal covariance"
  )
  .bt_random_effect_marginal_covariance_validate_draw_matrix(
    draws = source_draws,
    n_draws = n_draws,
    n_columns = n_rows,
    label = "row-indexed SD source",
    random_term = random_term,
    nonnegative = TRUE
  )

  column_allocation <- .bt_random_effect_row_indexed_column_allocation_draws(
    random_term = random_term,
    posterior = posterior,
    prior_list = prior_list,
    n_columns = n_columns
  )
  if(is.null(column_allocation)){
    allocation <- .bt_random_effect_row_indexed_allocation_draws(
      random_term = random_term,
      posterior = posterior,
      prior_list = prior_list
    )
    .bt_random_effect_marginal_covariance_validate_draw_matrix(
      draws = matrix(allocation, ncol = 1L),
      n_draws = n_draws,
      n_columns = 1L,
      label = "row-indexed SD allocation",
      random_term = random_term,
      nonnegative = TRUE
    )
  }else{
    .bt_random_effect_marginal_covariance_validate_draw_matrix(
      draws = column_allocation,
      n_draws = n_draws,
      n_columns = n_columns,
      label = "row-indexed column SD allocation",
      random_term = random_term,
      nonnegative = TRUE
    )
  }

  if(is.null(column_allocation)){
    row_weights <- source_draws *
      matrix(allocation, nrow = n_draws, ncol = n_rows)
    column_weights <- NULL
  }else{
    row_weights <- source_draws
    column_weights <- column_allocation
  }

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect marginal covariance metadata"
  )
  if(structure %in% c("diag", "id")){
    return(.bt_random_effect_marginal_covariance_independent(
      model_matrix = model_matrix,
      group_map = group_map,
      n_draws = n_draws,
      column_scale_draws = column_weights,
      row_scale_draws = row_weights,
      diagonal_only = diagonal_only
    ))
  }
  if(identical(structure, "car")){
    return(.bt_random_effect_marginal_covariance_car(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      row_weights = row_weights,
      column_weights = column_weights,
      diagonal_only = diagonal_only
    ))
  }

  if(isTRUE(diagonal_only)){
    one_sparse <- .bt_random_effect_marginal_covariance_one_sparse_columns(
      model_matrix
    )
    if(!is.null(one_sparse)){
      row_sd <- if(is.null(column_allocation)){
        row_weights
      }else{
        source_draws *
          column_allocation[, one_sparse$row_column, drop = FALSE]
      }
      row_sd <- row_sd * matrix(
        one_sparse$row_value,
        nrow = n_draws,
        ncol = n_rows,
        byrow = TRUE
      )
      return(.bt_random_effect_marginal_variance_one_sparse(
        model_matrix = model_matrix,
        row_sd = row_sd
      ))
    }
  }

  row_column <- .bt_random_effect_marginal_covariance_structured_columns(
    random_term = random_term,
    model_matrix = model_matrix
  )
  if(!is.null(row_column)){
    row_sd <- if(is.null(column_allocation)){
      row_weights
    }else{
      source_draws * column_allocation[, row_column, drop = FALSE]
    }
    return(.bt_random_effect_marginal_covariance_structured_one_hot(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      row_column = row_column,
      row_sd = row_sd,
      diagonal_only = diagonal_only
    ))
  }

  if(isTRUE(diagonal_only)){
    cholesky <- .bt_random_effect_cholesky_draws(
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior
    )
    if(is.null(cholesky)){
      .bt_random_effect_marginal_covariance_missing_correlation_stop(
        random_term = random_term,
        n_columns = n_columns,
        posterior = posterior
      )
    }
    return(.bt_random_effect_marginal_variance_cholesky(
      random_term = random_term,
      model_matrix = model_matrix,
      cholesky = cholesky,
      row_weights = row_weights,
      column_weights = column_weights
    ))
  }

  correlation <- .bt_random_effect_marginal_covariance_correlation_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )

  list(
    samples = .bt_random_effect_marginal_covariance_expand(
      model_matrix = model_matrix,
      group_map = group_map,
      block_covariance = correlation,
      row_weights = row_weights,
      column_weights = column_weights
    ),
    sample_dim = c(n_draws, n_rows, n_rows)
  )
}

.bt_random_effect_marginal_covariance_one_sparse_columns <- function(
    model_matrix){

  if(any(!is.finite(model_matrix))){
    return(NULL)
  }
  nonzero <- model_matrix != 0
  if(any(rowSums(nonzero) != 1L)){
    return(NULL)
  }

  row_column <- max.col(nonzero, ties.method = "first")
  row_value <- model_matrix[cbind(seq_len(nrow(model_matrix)), row_column)]
  if(any(!is.finite(row_value))){
    return(NULL)
  }

  list(row_column = row_column, row_value = row_value)
}

.bt_random_effect_marginal_variance_one_sparse <- function(model_matrix,
                                                           row_sd){

  n_draws <- nrow(row_sd)
  n_rows  <- nrow(model_matrix)
  if(!is.matrix(row_sd) || !identical(dim(row_sd), c(n_draws, n_rows)) ||
     any(!is.finite(row_sd))){
    stop("One-sparse random-effect marginal SD draws are inconsistent.",
         call. = FALSE)
  }

  row_names <- rownames(model_matrix)
  if(is.null(row_names)){
    row_names <- as.character(seq_len(n_rows))
  }
  samples <- row_sd^2
  colnames(samples) <- row_names

  list(
    samples = samples,
    sample_dim = c(n_draws, n_rows)
  )
}

.bt_random_effect_marginal_variance_cholesky <- function(
    random_term,
    model_matrix,
    cholesky,
    row_weights = NULL,
    column_weights = NULL){

  .bt_random_effect_marginal_covariance_cholesky(
    random_term = random_term,
    model_matrix = model_matrix,
    group_map = NULL,
    cholesky = cholesky,
    row_weights = row_weights,
    column_weights = column_weights,
    diagonal_only = TRUE
  )
}

.bt_random_effect_marginal_covariance_car <- function(
    random_term,
    model_matrix,
    group_map,
    posterior,
    row_weights = NULL,
    column_weights = NULL,
    diagonal_only = FALSE){

  n_columns <- ncol(model_matrix)
  correlation <- .bt_random_effect_correlation_metadata(
    random_term = random_term,
    structure = "car",
    context = "Random-effect marginal covariance metadata"
  )
  .bt_random_effect_car_time_values(
    random_term = random_term,
    correlation = correlation,
    n_columns = n_columns,
    context = "Random-effect marginal covariance metadata"
  )
  cholesky <- .bt_random_effect_cholesky_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )
  if(is.null(cholesky)){
    .bt_random_effect_marginal_covariance_missing_correlation_stop(
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior
    )
  }

  .bt_random_effect_marginal_covariance_cholesky(
    random_term = random_term,
    model_matrix = model_matrix,
    group_map = group_map,
    cholesky = cholesky,
    row_weights = row_weights,
    column_weights = column_weights,
    diagonal_only = diagonal_only
  )
}

.bt_random_effect_marginal_covariance_cholesky <- function(
    random_term,
    model_matrix,
    group_map,
    cholesky,
    row_weights = NULL,
    column_weights = NULL,
    diagonal_only = FALSE){

  cholesky_dim <- dim(cholesky)
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  valid_cholesky <- is.array(cholesky) &&
    length(cholesky_dim) == 3L &&
    cholesky_dim[2L] == n_columns &&
    cholesky_dim[3L] == n_columns &&
    all(is.finite(cholesky))
  if(!isTRUE(valid_cholesky)){
    stop(
      "Random-effect marginal covariance Cholesky draws for block '",
      random_term$block_name,
      "' do not match the expected dimensions.",
      call. = FALSE
    )
  }
  n_draws <- cholesky_dim[1L]
  if(!is.null(row_weights) &&
     (!is.matrix(row_weights) || !is.numeric(row_weights) ||
      !identical(dim(row_weights), c(n_draws, n_rows)) ||
      any(!is.finite(row_weights)))){
    stop(
      "Random-effect marginal covariance row weights for block '",
      random_term$block_name,
      "' do not match the expected dimensions.",
      call. = FALSE
    )
  }
  if(!is.null(column_weights) &&
     (!is.matrix(column_weights) || !is.numeric(column_weights) ||
      !identical(dim(column_weights), c(n_draws, n_columns)) ||
      any(!is.finite(column_weights)))){
    stop(
      "Random-effect marginal covariance column weights for block '",
      random_term$block_name,
      "' do not match the expected dimensions.",
      call. = FALSE
    )
  }
  if(!isTRUE(diagonal_only) &&
     (!is.numeric(group_map) || length(group_map) != n_rows ||
      anyNA(group_map) || any(!is.finite(group_map)) ||
      any(group_map < 1L) || any(group_map > .Machine$integer.max) ||
      any(group_map != floor(group_map)))){
    stop(
      "Random-effect marginal covariance group map for block '",
      random_term$block_name,
      "' does not match the Cholesky expansion.",
      call. = FALSE
    )
  }

  row_names <- rownames(model_matrix)
  if(is.null(row_names)){
    row_names <- as.character(seq_len(n_rows))
  }
  if(isTRUE(diagonal_only)){
    samples <- matrix(
      0,
      nrow = n_draws,
      ncol = n_rows,
      dimnames = list(draw = NULL, row = row_names)
    )
    for(draw in seq_len(n_draws)){
      Z <- model_matrix
      if(!is.null(row_weights)){
        Z <- Z * row_weights[draw, ]
      }
      if(!is.null(column_weights)){
        Z <- sweep(
          Z,
          MARGIN = 2L,
          STATS = column_weights[draw, ],
          FUN = "*"
        )
      }
      samples[draw, ] <- rowSums((Z %*% cholesky[draw, , ])^2)
    }
    return(list(
      samples = samples,
      sample_dim = c(n_draws, n_rows)
    ))
  }

  samples <- array(
    0,
    dim = c(n_draws, n_rows, n_rows),
    dimnames = list(draw = NULL, row = row_names, column = row_names)
  )
  rows_by_group <- split(seq_len(n_rows), group_map)
  for(draw in seq_len(n_draws)){
    for(rows in rows_by_group){
      Z <- model_matrix[rows, , drop = FALSE]
      if(!is.null(row_weights)){
        Z <- Z * row_weights[draw, rows]
      }
      if(!is.null(column_weights)){
        Z <- sweep(
          Z,
          MARGIN = 2L,
          STATS = column_weights[draw, ],
          FUN = "*"
        )
      }
      ZL <- Z %*% cholesky[draw, , ]
      samples[draw, rows, rows] <- tcrossprod(ZL)
    }
  }
  list(
    samples = samples,
    sample_dim = c(n_draws, n_rows, n_rows)
  )
}

.bt_random_effect_marginal_covariance_structured_columns <- function(
    random_term, model_matrix){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect marginal covariance metadata"
  )
  if(!structure %in% c("cs", "hcs", "ar1", "car", "har")){
    return(NULL)
  }
  nonzero <- abs(model_matrix) > sqrt(.Machine$double.eps)
  if(any(rowSums(nonzero) != 1L)){
    return(NULL)
  }
  row_column <- max.col(nonzero, ties.method = "first")
  selected <- model_matrix[cbind(seq_len(nrow(model_matrix)), row_column)]
  if(any(abs(selected - 1) > sqrt(.Machine$double.eps))){
    return(NULL)
  }

  row_column
}

.bt_random_effect_marginal_covariance_structured_one_hot <- function(
    random_term, model_matrix, group_map, posterior, row_column, row_sd,
    diagonal_only = FALSE){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect marginal covariance metadata"
  )
  rho <- .bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = posterior,
    missing = "error",
    out_of_support = "error",
    context = "Random-effect marginal covariance metadata"
  )
  n_draws <- nrow(posterior)
  n_rows  <- nrow(model_matrix)
  if(!is.matrix(row_sd) || !identical(dim(row_sd), c(n_draws, n_rows)) ||
     any(!is.finite(row_sd)) || any(row_sd < 0)){
    stop(
      "Random-effect marginal covariance row SD draws for block '",
      random_term$block_name,
      "' do not match the structured one-hot design.",
      call. = FALSE
    )
  }

  row_names <- rownames(model_matrix)
  if(is.null(row_names)){
    row_names <- as.character(seq_len(n_rows))
  }
  if(isTRUE(diagonal_only)){
    out <- row_sd^2
    colnames(out) <- row_names
    return(list(
      samples = out,
      sample_dim = c(n_draws, n_rows)
    ))
  }

  exponent <- if(structure %in% c("cs", "hcs")){
    outer(row_column, row_column, "!=") * 1
  }else if(identical(structure, "car")){
    correlation <- .bt_random_effect_correlation_metadata(
      random_term = random_term,
      structure = structure,
      context = "Random-effect marginal covariance metadata"
    )
    time_values <- .bt_random_effect_car_time_values(
      random_term = random_term,
      correlation = correlation,
      n_columns = ncol(model_matrix),
      context = "Random-effect marginal covariance metadata"
    )
    abs(outer(time_values[row_column], time_values[row_column], "-"))
  }else{
    abs(outer(row_column, row_column, "-"))
  }
  same_group <- outer(group_map, group_map, "==")
  out <- array(
    0,
    dim = c(n_draws, n_rows, n_rows),
    dimnames = list(draw = NULL, row = row_names, column = row_names)
  )
  for(draw in seq_len(n_draws)){
    out[draw, , ] <- same_group * rho[draw]^exponent *
      tcrossprod(row_sd[draw, ])
  }

  list(
    samples = out,
    sample_dim = c(n_draws, n_rows, n_rows)
  )
}

.bt_random_effect_marginal_covariance_validate_draw_matrix <- function(
    draws,
    n_draws,
    n_columns,
    label,
    random_term,
    nonnegative = FALSE,
    context = "Random-effect marginal covariance"){

  if(!is.matrix(draws) || !is.numeric(draws) ||
     nrow(draws) != n_draws || ncol(draws) != n_columns){
    stop(
      context, " ",
      label,
      " draws for block '",
      random_term$block_name,
      "' do not match the expected dimensions.",
      call. = FALSE
    )
  }
  invalid <- !is.finite(draws)
  if(isTRUE(nonnegative)){
    invalid <- invalid | draws < 0
  }
  if(any(invalid)){
    stop(
      context, " ",
      label,
      " draws for block '",
      random_term$block_name,
      "' must be finite",
      if(isTRUE(nonnegative)) " and non-negative" else "",
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_random_effect_marginal_covariance_correlation_draws <- function(
    random_term,
    n_columns,
    posterior){

  cholesky <- .bt_random_effect_cholesky_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )
  if(is.null(cholesky)){
    .bt_random_effect_marginal_covariance_missing_correlation_stop(
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior
    )
  }

  .bt_random_effect_marginal_covariance_validate_correlation_cholesky(
    cholesky = cholesky,
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )

  out <- array(NA_real_, dim = dim(cholesky))
  for(draw in seq_len(nrow(posterior))){
    out[draw, , ] <- tcrossprod(cholesky[draw, , ])
  }
  if(any(!is.finite(out))){
    stop(
      "Random-effect marginal covariance correlation draws for block '",
      random_term$block_name,
      "' must define finite correlation matrices.",
      call. = FALSE
    )
  }

  out
}

.bt_random_effect_marginal_covariance_validate_correlation_cholesky <- function(
    cholesky,
    random_term,
    n_columns,
    posterior){

  if(!is.array(cholesky) || length(dim(cholesky)) != 3L ||
     !identical(dim(cholesky), c(nrow(posterior), n_columns, n_columns))){
    stop(
      "Random-effect marginal covariance correlation draws for block '",
      random_term$block_name,
      "' do not match the expected dimensions.",
      call. = FALSE
    )
  }
  if(any(!is.finite(cholesky))){
    stop(
      "Random-effect marginal covariance correlation draws for block '",
      random_term$block_name,
      "' must be finite.",
      call. = FALSE
    )
  }

  diagonal <- matrix(
    rowSums(matrix(
      cholesky^2,
      nrow = nrow(posterior) * n_columns,
      ncol = n_columns
    )),
    nrow = nrow(posterior),
    ncol = n_columns
  )
  if(any(!is.finite(diagonal))){
    stop(
      "Random-effect marginal covariance correlation draws for block '",
      random_term$block_name,
      "' must define finite correlation matrices.",
      call. = FALSE
    )
  }
  if(any(abs(diagonal - 1) > sqrt(.Machine$double.eps))){
    stop(
      "Random-effect marginal covariance correlation draws for block '",
      random_term$block_name,
      "' must define correlation matrices with unit diagonal.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# Expand independent random effects through their weighted design directly.
.bt_random_effect_marginal_covariance_independent <- function(
    model_matrix,
    group_map,
    n_draws,
    column_scale_draws = NULL,
    row_scale_draws = NULL,
    diagonal_only = FALSE){

  n_rows    <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  row_names <- rownames(model_matrix)
  if(is.null(row_names)){
    row_names <- as.character(seq_len(n_rows))
  }
  if(isTRUE(diagonal_only)){
    out <- matrix(
      0,
      nrow = n_draws,
      ncol = n_rows,
      dimnames = list(draw = NULL, row = row_names)
    )
    for(draw in seq_len(n_draws)){
      Z <- model_matrix
      if(!is.null(column_scale_draws)){
        Z <- sweep(
          Z,
          MARGIN = 2L,
          STATS = column_scale_draws[draw, ],
          FUN = "*"
        )
      }
      if(!is.null(row_scale_draws)){
        Z <- Z * row_scale_draws[draw, ]
      }
      out[draw, ] <- rowSums(Z^2)
    }
    return(list(
      samples = out,
      sample_dim = c(n_draws, n_rows)
    ))
  }
  out <- array(
    0,
    dim = c(n_draws, n_rows, n_rows),
    dimnames = list(draw = NULL, row = row_names, column = row_names)
  )
  rows_by_group <- split(seq_len(n_rows), group_map)
  for(draw in seq_len(n_draws)){
    for(rows in rows_by_group){
      Z <- model_matrix[rows, , drop = FALSE]
      if(!is.null(column_scale_draws)){
        Z <- sweep(
          Z,
          MARGIN = 2L,
          STATS = column_scale_draws[draw, ],
          FUN = "*"
        )
      }
      if(!is.null(row_scale_draws)){
        Z <- Z * row_scale_draws[draw, rows]
      }
      out[draw, rows, rows] <- tcrossprod(Z)
    }
  }

  list(
    samples = out,
    sample_dim = c(n_draws, n_rows, n_rows)
  )
}

.bt_random_effect_marginal_covariance_missing_sd_stop <- function(random_term,
                                                                  n_columns){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect marginal covariance metadata"
  )
  sd_names <- random_term$sd_parameter_names
  if(is.null(sd_names) || length(sd_names) != n_columns){
    stop(
      "Random-effect marginal covariance metadata for block '",
      random_term$block_name,
      "' with structure '",
      structure,
      "' are missing canonical 'random_term$sd_parameter_names'.",
      call. = FALSE
    )
  }

  expected <- .bt_random_effect_marginal_covariance_expected_sd_names(random_term)
  stop(
    "Random-effect marginal covariance for block '",
    random_term$block_name,
    "' with structure '",
    structure,
    "' cannot resolve SD draws. Expected posterior columns or point priors for: ",
    paste0("'", expected[seq_len(min(4L, length(expected)))], "'", collapse = ", "),
    if(length(expected) > 4L) ", ..." else "",
    ".",
    call. = FALSE
  )
}

.bt_random_effect_marginal_covariance_expected_sd_names <- function(random_term){

  expected <- character()
  sd_names <- random_term$sd_parameter_names
  if(is.character(sd_names)){
    expected <- c(expected, sd_names[!is.na(sd_names) & nzchar(sd_names)])
  }

  binding <- random_term$sd_binding
  if(!is.null(binding) && isTRUE(binding$true_allocation) &&
     length(binding$allocations) > 0L){
    allocation <- binding$allocations[[1L]]
    if(!is.null(allocation$source)){
      expected <- c(expected, .bt_random_sd_binding_source_name(allocation$source))
    }
    if(!is.null(allocation$weight_name)){
      expected <- c(
        expected,
        allocation$weight_name,
        paste0(.JAGS_prior_dirichlet_eta_name(allocation$weight_name), "[...]")
      )
    }
  }

  expected <- unique(expected[nzchar(expected)])
  if(length(expected) == 0L){
    expected <- "<unknown SD parameter>"
  }

  expected
}

.bt_random_effect_marginal_covariance_missing_correlation_stop <- function(
    random_term,
    n_columns,
    posterior){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect marginal covariance metadata"
  )

  if(structure %in% c("cs", "hcs", "ar1", "car", "har")){
    .bt_random_effect_rho_draws(
      random_term = random_term,
      posterior = posterior,
      missing = "error",
      out_of_support = "error",
      context = "Random-effect marginal covariance metadata"
    )
  }

  expected <- .bt_random_effect_cholesky_names(
    random_term = random_term,
    n_columns = n_columns
  )
  if(identical(structure, "us")){
    expected <- c(
      as.vector(expected),
      .bt_random_effect_lkj_primitive_names(
        random_term = random_term,
        n_columns = n_columns,
        context = "Random-effect marginal covariance metadata"
      )
    )
  }else{
    correlation <- random_term$correlation
    if(is.list(correlation)){
      expected <- c(expected, correlation$rho_name, correlation$sample_name)
    }
  }
  expected <- unique(expected[!is.na(expected) & nzchar(expected)])

  stop(
    "Random-effect marginal covariance for block '",
    random_term$block_name,
    "' with structure '",
    structure,
    "' cannot resolve correlation draws. Expected posterior column(s): ",
    paste0("'", expected[seq_len(min(4L, length(expected)))], "'", collapse = ", "),
    if(length(expected) > 4L) ", ..." else "",
    ".",
    call. = FALSE
  )
}

.bt_random_effect_marginal_covariance_expand <- function(
    model_matrix,
    group_map,
    block_covariance,
    row_weights = NULL,
    column_weights = NULL){

  n_draws <- dim(block_covariance)[1L]
  n_rows <- nrow(model_matrix)
  n_columns <- ncol(model_matrix)
  row_names <- rownames(model_matrix)
  if(is.null(row_names)){
    row_names <- as.character(seq_len(n_rows))
  }
  out <- array(0, dim = c(n_draws, n_rows, n_rows))
  dimnames(out) <- list(draw = NULL, row = row_names, column = row_names)

  rows_by_group <- split(seq_len(n_rows), group_map)
  n_draws_num <- as.numeric(n_draws)
  n_rows_num <- as.numeric(n_rows)
  group_info <- lapply(
    rows_by_group,
    function(rows){
      rows0 <- as.numeric(rows) - 1
      n_group_rows <- length(rows)
      list(
        rows = rows,
        model_matrix = model_matrix[rows, , drop = FALSE],
        array_offset = rep(rows0, times = n_group_rows) * n_draws_num +
          rep(rows0, each = n_group_rows) * n_draws_num * n_rows_num
      )
    }
  )
  intercept_only <- n_columns == 1L &&
    is.null(row_weights) &&
    is.null(column_weights) &&
    isTRUE(all(model_matrix[, 1L] == 1))

  for(draw in seq_len(n_draws)){
    G <- matrix(
      block_covariance[draw, , ],
      nrow = n_columns,
      ncol = n_columns
    )
    for(group in group_info){
      rows <- group$rows
      Z <- group$model_matrix
      if(!is.null(row_weights)){
        Z <- Z * row_weights[draw, rows]
      }
      if(!is.null(column_weights)){
        Z <- Z * matrix(
          column_weights[draw, ],
          nrow = length(rows),
          ncol = n_columns,
          byrow = TRUE
        )
      }
      index <- draw + group$array_offset
      if(intercept_only){
        out[index] <- out[index] + G[1L, 1L]
      }else if(n_columns == 1L){
        out[index] <- out[index] + G[1L, 1L] * as.vector(tcrossprod(Z[, 1L]))
      }else{
        out[index] <- out[index] + as.vector(Z %*% G %*% t(Z))
      }
    }
  }

  out
}

.bt_random_effect_marginal_covariance_block_metadata <- function(
    random_term,
    model_matrix,
    group_map,
    group_levels,
    row_names,
    sample_dim,
    new_level_info = NULL,
    diagonal_only = FALSE){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect marginal covariance metadata"
  )
  correlation <- .bt_random_effect_correlation_metadata(
    random_term = random_term,
    structure = structure,
    context = "Random-effect marginal covariance metadata"
  )

  metadata <- list(
    block_name = random_term$block_name,
    grouping = random_term$group_label,
    structure = structure,
    compile_mode = .bt_random_effect_term_compile_mode(random_term),
    n_groups = length(group_levels),
    fitted_n_groups = random_term$n_groups,
    n_columns = random_term$n_columns,
    n_rows = nrow(model_matrix),
    row_names = row_names,
    row_order = seq_len(nrow(model_matrix)),
    group_levels = group_levels,
    group_map = group_map,
    new_levels = if(is.list(new_level_info)) new_level_info$policy else NULL,
    new_group_levels = if(is.list(new_level_info)) new_level_info$group_levels else character(),
    new_level_rows = if(is.list(new_level_info)) new_level_info$rows else integer(),
    column_names = colnames(model_matrix),
    model_matrix = model_matrix,
    sd_parameter_names = random_term$sd_parameter_names,
    row_varying_sd = .bt_random_effect_has_row_indexed_external_sd(random_term),
    correlation_type = if(is.list(correlation)) correlation$type else NA_character_,
    group_covariance = .bt_random_effect_group_covariance_metadata(random_term),
    included = TRUE,
    skipped = FALSE,
    dense = !isTRUE(diagonal_only),
    dense_entries = if(isTRUE(diagonal_only)) NA_real_ else prod(sample_dim)
  )
  if(isTRUE(diagonal_only)){
    metadata$representation <- "diagonal"
    metadata$quantity <- "variance"
    metadata$diagonal_only <- TRUE
    metadata$sample_dim <- sample_dim
    metadata$sample_entries <- prod(sample_dim)
    metadata$equivalent_dense_entries <-
      as.numeric(sample_dim[[1L]]) * as.numeric(sample_dim[[2L]])^2
  }

  metadata
}
