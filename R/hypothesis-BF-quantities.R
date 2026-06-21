.as_hypothesis_quantities <- function(posterior, prior, parsed, parameter) {

  if(inherits(posterior, "marginal_inference")){
    return(.as_hypothesis_quantities_marginal_inference(
      posterior = posterior,
      parsed    = parsed,
      parameter = parameter
    ))
  }

  if(.hypothesis_inherits_marginal_posterior(posterior)){
    return(.as_hypothesis_quantities_marginal_posterior(
      posterior = posterior,
      parsed    = parsed,
      parameter = parameter
    ))
  }

  if(is.data.frame(posterior) || is.matrix(posterior)){
    if(is.null(prior)){
      stop("Prior draws are required for data-frame posterior inputs.",
           call. = FALSE)
    }
    return(list(.hypothesis_quantity_from_draws(
      posterior = as.data.frame(posterior, check.names = FALSE),
      prior     = as.data.frame(prior, check.names = FALSE),
      label     = "draws"
    )))
  }

  if(is.numeric(posterior)){
    if(is.null(prior)){
      stop("Prior draws are required for numeric posterior inputs.",
           call. = FALSE)
    }
    check_real(posterior, "posterior", check_length = 0, allow_NA = FALSE)
    if(is.null(parameter)){
      parameter <- .hypothesis_single_symbol(parsed)
    }
    if(is.null(parameter)){
      stop("The 'parameter' argument is required for numeric draws unless ",
           "the hypothesis contains exactly one named quantity.", call. = FALSE)
    }

    prior_info <- .hypothesis_prior_input_to_draws(
      prior     = prior,
      parameter = parameter,
      n         = max(length(posterior), 10000L)
    )
    posterior_df <- data.frame(posterior, check.names = FALSE)
    prior_df     <- data.frame(prior_info[["draws"]], check.names = FALSE)
    names(posterior_df) <- parameter
    names(prior_df)     <- parameter

    return(list(.hypothesis_quantity_from_draws(
      posterior     = posterior_df,
      prior         = prior_df,
      label         = parameter,
      parameter     = parameter,
      prior_density = prior_info[["density"]],
      prior_object  = prior_info[["prior_object"]]
    )))
  }

  stop("Unsupported posterior input for hypothesis_BF().", call. = FALSE)
}


.hypothesis_inherits_marginal_posterior <- function(x) {

  inherits(x, "marginal_posterior") ||
    any(grepl("^marginal_posterior\\.", class(x)))
}


.as_hypothesis_quantities_marginal_inference <- function(posterior, parsed,
                                                         parameter) {

  if(is.null(posterior[["conditional"]])){
    stop("'marginal_inference' object does not contain conditional draws.",
         call. = FALSE)
  }

  available <- names(posterior[["conditional"]])
  if(is.null(parameter)){
    level_refs <- .hypothesis_level_references(parsed)
    symbols <- .hypothesis_all_symbols(parsed)
    matches <- intersect(symbols, available)
    if(nrow(level_refs) > 0L &&
       length(unique(level_refs[["parameter"]])) == 1L &&
       unique(level_refs[["parameter"]]) %in% available){
      parameter <- unique(level_refs[["parameter"]])
    }else if(length(matches) == 1L){
      parameter <- matches
    }else if(length(matches) == 0L && length(available) == 1L){
      parameter <- available
    }else{
      stop("Specify 'parameter' for this marginal_inference object.",
           call. = FALSE)
    }
  }
  if(!parameter %in% available){
    stop("Parameter '", parameter, "' is not available in 'posterior'.",
         call. = FALSE)
  }

  return(.as_hypothesis_quantities_marginal_posterior(
    posterior = posterior[["conditional"]][[parameter]],
    parsed    = parsed,
    parameter = parameter
  ))
}


.as_hypothesis_quantities_marginal_posterior <- function(posterior, parsed,
                                                         parameter) {

  level_refs <- .hypothesis_level_references(parsed, parameter)
  symbol_parameter <- .hypothesis_single_symbol(parsed)
  if(is.null(parameter) && nrow(level_refs) > 0L &&
     length(unique(level_refs[["parameter"]])) == 1L){
    parameter <- unique(level_refs[["parameter"]])
  }
  if(is.null(parameter)){
    parameter <- symbol_parameter
  }
  if(is.null(parameter)){
    parameter <- attr(posterior, "parameter", exact = TRUE)
  }
  if(is.null(parameter)){
    stop("The 'parameter' argument is required for this marginal posterior.",
         call. = FALSE)
  }

  if(is.list(posterior)){
    level_refs <- .hypothesis_level_references(parsed, parameter)
    if(nrow(level_refs) > 0L){
      return(list(.hypothesis_quantity_from_marginal_posterior_levels(
        posterior  = posterior,
        parameter  = parameter,
        level_refs = level_refs,
        parsed     = parsed
      )))
    }

    out <- lapply(seq_along(posterior), function(i){
      level <- names(posterior)[i]
      if(is.null(level) || !nzchar(level)){
        level <- attr(posterior[[i]], "level", exact = TRUE)
      }
      if(is.null(level) || !nzchar(level)){
        level <- attr(posterior[[i]], "level_name", exact = TRUE)
      }
      label <- if(is.null(level) || !nzchar(level)){
        parameter
      }else{
        paste0(parameter, "[", level, "]")
      }
      .hypothesis_quantity_from_marginal_posterior(
        posterior = posterior[[i]],
        parameter = parameter,
        label     = label,
        parent    = posterior,
        index     = i
      )
    })
    return(out)
  }

  return(list(.hypothesis_quantity_from_marginal_posterior(
    posterior = posterior,
    parameter = parameter,
    label     = parameter
  )))
}


.hypothesis_level_references <- function(parsed, parameter = NULL) {

  symbols <- .hypothesis_all_symbols(parsed)
  refs <- regexec("^([^\\[]+)\\[([^\\]]+)\\]$", symbols, perl = TRUE)
  matches <- regmatches(symbols, refs)
  has_match <- vapply(matches, length, integer(1)) == 3L
  if(!any(has_match)){
    return(data.frame(
      symbol    = character(),
      parameter = character(),
      level     = character(),
      stringsAsFactors = FALSE
    ))
  }

  matches <- matches[has_match]
  out <- data.frame(
    symbol    = vapply(matches, `[[`, character(1), 1L),
    parameter = vapply(matches, `[[`, character(1), 2L),
    level     = trimws(vapply(matches, `[[`, character(1), 3L)),
    stringsAsFactors = FALSE
  )

  if(!is.null(parameter)){
    out <- out[out[["parameter"]] == parameter, , drop = FALSE]
  }

  return(unique(out))
}


.hypothesis_prior_input_to_draws <- function(prior, parameter, n) {

  prior_density <- NULL
  prior_object <- NULL
  if(inherits(prior, "prior")){
    prior_object <- prior
    prior_density <- tryCatch(
      .prior_linear_combination_density(
        prior_list = setNames(list(prior), parameter),
        weights    = setNames(1, parameter),
        n_grid     = .prior_linear_density_default_grid()
      ),
      error = function(e) NULL
    )
    prior <- rng(prior, n)
    if(is.matrix(prior) || is.data.frame(prior)){
      if(ncol(prior) != 1L){
        stop("Prior object must generate scalar draws for numeric posterior inputs.",
             call. = FALSE)
      }
      prior <- prior[, 1L]
    }
  }

  check_real(prior, "prior", check_length = 0, allow_NA = FALSE)

  return(list(
    draws        = prior,
    density      = prior_density,
    prior_object = prior_object
  ))
}


.hypothesis_quantity_from_draws <- function(posterior, prior, label,
                                            parameter = NULL,
                                            prior_density = NULL,
                                            prior_object = NULL) {

  out <- list(
    label              = label,
    parameter          = parameter,
    posterior_draws    = posterior,
    prior_draws        = prior,
    posterior_marginal = NULL,
    posterior_marginals = NULL,
    posterior_marginal_parent = NULL,
    posterior_marginal_index = NULL,
    posterior_marginal_indices = NULL,
    prior_densities    = NULL,
    prior_density      = prior_density,
    prior_object       = prior_object
  )
  class(out) <- "BayesTools_hypothesis_quantity"

  return(out)
}


.hypothesis_quantity_from_marginal_posterior_levels <- function(posterior,
                                                                parameter,
                                                                level_refs,
                                                                parsed) {

  levels <- unique(level_refs[["level"]])
  available <- names(posterior)
  missing <- setdiff(levels, available)
  if(length(missing) > 0L){
    stop("Hypothesis references unknown level '",
         paste(missing, collapse = "', '"), "' for parameter '", parameter, "'.",
         call. = FALSE)
  }

  posterior_draws <- lapply(levels, function(level) as.numeric(posterior[[level]]))
  n_draws <- vapply(posterior_draws, length, integer(1))
  if(length(unique(n_draws)) != 1L){
    stop("Level comparisons require equal-length posterior draws.",
         call. = FALSE)
  }

  .hypothesis_validate_level_conditionals(posterior, parameter, levels)

  posterior_df <- as.data.frame(posterior_draws, check.names = FALSE)
  names(posterior_df) <- paste0(parameter, "[", levels, "]")

  prior_df <- NULL
  if(.hypothesis_level_hypotheses_need_prior_draws(parsed, level_refs)){
    prior_df <- .hypothesis_prior_draws_from_marginal_levels(
      posterior = posterior,
      parameter = parameter,
      levels    = levels,
      n         = max(nrow(posterior_df), 10000L)
    )
  }

  posterior_marginals <- lapply(levels, function(level) {
    .hypothesis_marginal_child(posterior[[level]])
  })
  names(posterior_marginals) <- names(posterior_df)
  posterior_marginal_indices <- match(levels, names(posterior))
  names(posterior_marginal_indices) <- names(posterior_df)
  prior_densities <- lapply(levels, function(level) {
    attr(posterior[[level]], "prior_density", exact = TRUE)
  })
  names(prior_densities) <- names(posterior_df)

  out <- list(
    label              = parameter,
    parameter          = parameter,
    posterior_draws    = posterior_df,
    prior_draws        = prior_df,
    posterior_marginal = NULL,
    posterior_marginals = posterior_marginals,
    posterior_marginal_parent = posterior,
    posterior_marginal_index = NULL,
    posterior_marginal_indices = posterior_marginal_indices,
    prior_densities    = prior_densities,
    prior_density      = NULL,
    prior_object       = NULL
  )
  class(out) <- "BayesTools_hypothesis_quantity"

  return(out)
}


.hypothesis_marginal_child <- function(posterior) {

  if(.hypothesis_inherits_marginal_posterior(posterior) &&
     !inherits(posterior, "marginal_posterior")){
    class(posterior) <- unique(c(class(posterior), "marginal_posterior"))
  }

  return(posterior)
}


.hypothesis_validate_level_conditionals <- function(posterior, parameter,
                                                    levels) {

  keys <- vapply(levels, function(level) {
    .hypothesis_level_condition_key(posterior[[level]])
  }, character(1))
  if(all(keys == "<averaged>")){
    return(invisible(TRUE))
  }

  if(length(unique(keys)) > 1L){
    stop(
      "Level comparison for parameter '", parameter,
      "' uses different conditional posterior subsets. Use averaged marginals ",
      "or compare levels with identical conditionals.",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}


.hypothesis_level_condition_key <- function(level) {

  key <- attr(level, "condition_key", exact = TRUE)
  if(!is.null(key)){
    return(as.character(key))
  }

  conditional <- attr(level, "effective_conditional", exact = TRUE)
  conditional_rule <- attr(level, "effective_conditional_rule", exact = TRUE)
  if(is.null(conditional)){
    conditional <- attr(level, "conditional", exact = TRUE)
  }
  if(is.null(conditional_rule)){
    conditional_rule <- attr(level, "conditional_rule", exact = TRUE)
  }

  .hypothesis_conditional_key(conditional, conditional_rule)
}


.hypothesis_conditional_key <- function(conditional, conditional_rule = "AND") {

  if(is.null(conditional)){
    return("<averaged>")
  }
  if(is.null(conditional_rule)){
    conditional_rule <- "AND"
  }

  .condition_event_key(conditional, conditional_rule)
}


.hypothesis_context_condition_key <- function(context) {

  key <- context[["condition_key"]]
  if(!is.null(key)){
    return(as.character(key))
  }

  conditional <- context[["conditional"]]
  conditional_rule <- context[["conditional_rule"]]
  if(is.null(conditional)){
    return("<averaged>")
  }
  if(is.null(conditional_rule)){
    conditional_rule <- "AND"
  }

  .condition_event_key(conditional, conditional_rule)
}


.hypothesis_child_prior_context <- function(posterior, levels) {

  child_contexts <- lapply(levels, function(level) {
    attr(posterior[[level]], "prior_density_context", exact = TRUE)
  })
  has_context <- !vapply(child_contexts, is.null, logical(1))
  if(!any(has_context)){
    return(NULL)
  }

  if(!all(has_context)){
    stop(
      "Level comparisons require prior contexts for all conditional levels.",
      call. = FALSE
    )
  }
  valid_context <- vapply(child_contexts, .hypothesis_is_prior_density_context, logical(1))
  if(!all(valid_context)){
    stop("Invalid joint prior information for level-comparison hypotheses.",
         call. = FALSE)
  }

  keys <- vapply(child_contexts, .hypothesis_context_condition_key, character(1))
  if(length(unique(keys)) > 1L){
    stop(
      "Level comparison prior contexts use different conditional posterior subsets.",
      call. = FALSE
    )
  }

  child_contexts[[1]]
}


.hypothesis_level_hypotheses_need_prior_draws <- function(parsed, level_refs) {

  valid_symbols <- unique(level_refs[["symbol"]])
  for(hypothesis in parsed){
    if(!.hypothesis_sides_point_complement(hypothesis[["left"]],
                                           hypothesis[["right"]]) &&
       !.hypothesis_sides_point_complement(hypothesis[["right"]],
                                           hypothesis[["left"]])){
      return(TRUE)
    }

    point_side <- if(identical(hypothesis[["left"]][["type"]], "point")){
      hypothesis[["left"]]
    }else{
      hypothesis[["right"]]
    }
    symbol <- .hypothesis_direct_symbol(point_side[["expr"]])
    if(is.null(symbol) || !symbol %in% valid_symbols){
      return(TRUE)
    }
  }

  return(FALSE)
}


.hypothesis_prior_draws_from_marginal_levels <- function(posterior, parameter,
                                                         levels, n) {

  attr_prior_draws <- attr(posterior, "prior_draws", exact = TRUE)
  columns <- paste0(parameter, "[", levels, "]")
  if(!is.null(attr_prior_draws)){
    attr_prior_draws <- as.data.frame(attr_prior_draws, check.names = FALSE)
    missing <- setdiff(columns, names(attr_prior_draws))
    if(length(missing) == 0L){
      return(attr_prior_draws[, columns, drop = FALSE])
    }
  }

  child_prior_draws <- lapply(levels, function(level) {
    attr(posterior[[level]], "prior_draws", exact = TRUE)
  })
  if(all(!vapply(child_prior_draws, is.null, logical(1)))){
    n_draws <- vapply(child_prior_draws, length, integer(1))
    if(length(unique(n_draws)) != 1L){
      stop("Level comparisons require equal-length prior draws.",
           call. = FALSE)
    }
    out <- as.data.frame(child_prior_draws, check.names = FALSE)
    names(out) <- columns
    return(out)
  }

  context <- .hypothesis_child_prior_context(posterior, levels)
  if(is.null(context)){
    context <- attr(posterior, "prior_density_context", exact = TRUE)
  }
  if(is.null(context)){
    stop("Joint prior information is required for level-comparison hypotheses.",
         call. = FALSE)
  }
  if(!.hypothesis_is_prior_density_context(context)){
    stop("Invalid joint prior information for level-comparison hypotheses.",
         call. = FALSE)
  }

  level_weights <- lapply(levels, function(level) {
    weights <- attr(posterior[[level]], "linear_weights", exact = TRUE)
    if(is.null(weights)){
      stop("Linear prior weights are missing for level '", level, "'.",
           call. = FALSE)
    }
    .hypothesis_prepare_level_weights(weights)
  })
  names(level_weights) <- levels
  .hypothesis_validate_level_weights_context(level_weights, context)

  prior_matrix <- .hypothesis_prior_samples_from_context(context, n)
  row_i <- .hypothesis_level_weight_rows(level_weights, nrow(prior_matrix))

  prior_values <- lapply(levels, function(level){
    .hypothesis_apply_level_weights(prior_matrix, level_weights[[level]], row_i)
  })

  out <- as.data.frame(prior_values, check.names = FALSE)
  names(out) <- columns

  return(out)
}


.hypothesis_validate_level_weights_context <- function(level_weights, context) {

  for(level in names(level_weights)){
    weights <- .hypothesis_prepare_level_weights(level_weights[[level]])
    for(row_i in seq_len(nrow(weights))){
      .hypothesis_validate_context_weight_vector(context, weights[row_i, ])
    }
  }

  invisible(TRUE)
}


.hypothesis_validate_context_weight_vector <- function(context, weights) {

  nonzero_columns <- .hypothesis_nonzero_weight_columns(weights)
  if(length(nonzero_columns) == 0L){
    return(invisible(TRUE))
  }
  .hypothesis_check_weight_columns_available(nonzero_columns, context)

  if(inherits(context, "prior_density_context")){
    standardized_weights <- .prior_density_context_standardized_weights(
      context = context,
      weights  = weights
    )
    .prior_linear_weight_groups(context[["prior_list"]], standardized_weights)
    return(invisible(TRUE))
  }

  if(inherits(context, "prior_density_model_mixture_context")){
    for(model_i in seq_along(context[["model_weights"]])){
      model_prior_list <- .hypothesis_model_mixture_prior_list(context, model_i)
      .prior_linear_weight_groups(model_prior_list, weights)
    }
    return(invisible(TRUE))
  }

  if(inherits(context, "prior_density_conditional_context")){
    if(length(context[["prior_lists"]]) == 0L){
      stop("No prior models remain after applying the conditional event.",
           call. = FALSE)
    }
    for(prior_list in context[["prior_lists"]]){
      if(!is.null(context[["formula_scale"]]) &&
         length(context[["formula_scale"]]) > 0L){
        component_context <- .prior_density_context(
          prior_list    = prior_list,
          column_names  = context[["column_names"]],
          formula_scale = context[["formula_scale"]],
          n_grid        = context[["n_grid"]],
          tail_prob     = context[["tail_prob"]]
        )
        standardized_weights <- .prior_density_context_standardized_weights(
          context = component_context,
          weights  = weights
        )
        .prior_linear_weight_groups(prior_list, standardized_weights)
      }else{
        .prior_linear_weight_groups(prior_list, weights)
      }
    }
    return(invisible(TRUE))
  }

  stop("Unknown prior density context.", call. = FALSE)
}


.hypothesis_nonzero_weight_columns <- function(weights) {

  if(is.null(names(weights))){
    return(character())
  }

  names(weights)[is.finite(weights) &
                   abs(weights) > .prior_linear_density_zero_tol()]
}


.hypothesis_check_weight_columns_available <- function(nonzero_columns, context) {

  missing <- setdiff(nonzero_columns, context[["column_names"]])
  if(length(missing) > 0L){
    stop("Linear prior weights reference columns not available in the joint ",
         "prior context: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }

  invisible(TRUE)
}


.hypothesis_is_prior_density_context <- function(context) {

  inherits(context, "prior_density_context") ||
    inherits(context, "prior_density_model_mixture_context") ||
    inherits(context, "prior_density_conditional_context")
}


.hypothesis_prior_samples_from_context <- function(context, n) {

  if(inherits(context, "prior_density_context")){
    samples <- .generate_transformed_prior_samples(
      prior_list    = context[["prior_list"]],
      column_names  = context[["column_names"]],
      n_samples     = n,
      formula_scale = context[["formula_scale"]]
    )
    return(.hypothesis_complete_prior_sample_matrix(
      samples      = samples,
      column_names = context[["column_names"]],
      n            = n
    ))
  }

  if(inherits(context, "prior_density_model_mixture_context")){
    return(.hypothesis_prior_samples_from_model_mixture_context(context, n))
  }

  if(inherits(context, "prior_density_conditional_context")){
    return(.hypothesis_prior_samples_from_conditional_context(context, n))
  }

  stop("Unknown prior density context.", call. = FALSE)
}


.hypothesis_prior_samples_from_model_mixture_context <- function(context, n) {

  model_i <- sample(seq_along(context[["model_weights"]]), size = n,
                    replace = TRUE, prob = context[["model_weights"]])
  out <- .hypothesis_empty_prior_sample_matrix(context[["column_names"]], n)

  for(i in unique(model_i)){
    rows <- which(model_i == i)
    model_prior_list <- .hypothesis_model_mixture_prior_list(context, i)
    samples <- .generate_transformed_prior_samples(
      prior_list   = model_prior_list,
      column_names = context[["column_names"]],
      n_samples    = length(rows)
    )
    out[rows, ] <- .hypothesis_complete_prior_sample_matrix(
      samples      = samples,
      column_names = context[["column_names"]],
      n            = length(rows)
    )
  }

  return(out)
}


.hypothesis_model_mixture_prior_list <- function(context, model_i) {

  model_prior_list <- lapply(context[["prior_list"]], function(parameter_priors) {
    if(is.prior(parameter_priors)){
      return(parameter_priors)
    }
    parameter_priors[[model_i]]
  })
  names(model_prior_list) <- names(context[["prior_list"]])

  for(parameter in names(model_prior_list)){
    if(is.null(model_prior_list[[parameter]])){
      model_prior_list[[parameter]] <- prior("point", list(location = 0))
    }
  }

  return(model_prior_list)
}


.hypothesis_prior_samples_from_conditional_context <- function(context, n) {

  out <- .hypothesis_empty_prior_sample_matrix(context[["column_names"]], n)
  if(length(context[["prior_lists"]]) == 0L){
    stop("No prior models remain after applying the conditional event.", call. = FALSE)
  }

  model_i <- sample(seq_along(context[["model_weights"]]), size = n,
                    replace = TRUE, prob = context[["model_weights"]])
  for(i in unique(model_i)){
    rows <- which(model_i == i)
    samples <- .generate_transformed_prior_samples(
      prior_list    = context[["prior_lists"]][[i]],
      column_names  = context[["column_names"]],
      n_samples     = length(rows),
      formula_scale = context[["formula_scale"]]
    )
    out[rows, ] <- .hypothesis_complete_prior_sample_matrix(
      samples      = samples,
      column_names = context[["column_names"]],
      n            = length(rows)
    )
  }

  return(out)
}


.hypothesis_empty_prior_sample_matrix <- function(column_names, n) {

  matrix(
    0,
    nrow = n,
    ncol = length(column_names),
    dimnames = list(NULL, column_names)
  )
}


.hypothesis_complete_prior_sample_matrix <- function(samples, column_names, n) {

  samples <- as.matrix(samples)
  out <- .hypothesis_empty_prior_sample_matrix(column_names, n)
  columns <- intersect(column_names, colnames(samples))
  if(length(columns) > 0L){
    out[, columns] <- samples[, columns, drop = FALSE]
  }

  return(out)
}


.hypothesis_prepare_level_weights <- function(weights) {

  if(is.null(dim(weights))){
    weight_names <- names(weights)
    weights <- matrix(weights, nrow = 1L)
    colnames(weights) <- weight_names
  }else{
    weights <- as.matrix(weights)
  }
  if(is.null(colnames(weights))){
    stop("Linear prior weights must be named.", call. = FALSE)
  }

  weights
}


.hypothesis_level_weight_rows <- function(level_weights, n) {

  row_counts <- vapply(level_weights, nrow, integer(1))
  row_counts <- row_counts[row_counts > 1L]
  if(length(row_counts) == 0L){
    return(NULL)
  }
  if(length(unique(row_counts)) != 1L){
    stop("Level comparisons with row-varying prior weights require matching ",
         "weight rows across referenced levels.", call. = FALSE)
  }

  sample.int(row_counts[[1L]], size = n, replace = TRUE)
}


.hypothesis_apply_level_weights <- function(samples, weights, row_i = NULL) {

  samples <- as.matrix(samples)
  weights <- .hypothesis_prepare_level_weights(weights)
  nonzero_columns <- colnames(weights)[colSums(abs(weights), na.rm = TRUE) >
                                          .prior_linear_density_zero_tol()]
  missing <- setdiff(nonzero_columns, colnames(samples))
  if(length(missing) > 0L){
    stop("Linear prior weights reference columns not available in the joint ",
         "prior context: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }

  columns <- intersect(colnames(weights), colnames(samples))
  if(length(columns) == 0L){
    return(rep(0, nrow(samples)))
  }
  if(nrow(weights) == 1L){
    return(as.numeric(samples[, columns, drop = FALSE] %*%
                        as.numeric(weights[1L, columns])))
  }

  if(is.null(row_i)){
    row_i <- sample.int(nrow(weights), size = nrow(samples), replace = TRUE)
  }
  if(any(row_i > nrow(weights))){
    stop("Shared level-weight rows exceed available weight rows.",
         call. = FALSE)
  }

  rowSums(samples[, columns, drop = FALSE] *
            weights[row_i, columns, drop = FALSE])
}


.hypothesis_quantity_from_marginal_posterior <- function(posterior, parameter,
                                                         label,
                                                         parent = NULL,
                                                         index = NULL) {

  if(.hypothesis_inherits_marginal_posterior(posterior) &&
     !inherits(posterior, "marginal_posterior")){
    class(posterior) <- unique(c(class(posterior), "marginal_posterior"))
  }

  posterior_df <- data.frame(as.numeric(posterior), check.names = FALSE)
  names(posterior_df) <- parameter

  out <- list(
    label              = label,
    parameter          = parameter,
    posterior_draws    = posterior_df,
    prior_draws        = NULL,
    posterior_marginal = posterior,
    posterior_marginals = NULL,
    posterior_marginal_parent = parent,
    posterior_marginal_index = index,
    posterior_marginal_indices = NULL,
    prior_densities    = NULL,
    prior_density      = attr(posterior, "prior_density", exact = TRUE),
    prior_object       = NULL
  )
  class(out) <- "BayesTools_hypothesis_quantity"

  return(out)
}
