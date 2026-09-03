#' Random-effect marginal covariance update plan
#'
#' @description
#' Resolves a selected public random-effect quantity to the exact scalar
#' dependence of its marginal covariance contribution. The result is compiled
#' exclusively from the fitted parameter map and formula random-effect
#' metadata; posterior draws are neither inspected nor used to infer an update
#' form.
#'
#' @param fit a fitted object carrying a BayesTools parameter map and formula
#'   design metadata.
#' @param selection a one-quantity selection returned by
#'   [parameter_catalog_resolve()].
#'
#' @return A list of class
#'   `BayesTools_random_effects_marginal_update_plan`. Field `family` is one of
#'   `"affine"`, `"factor"`, `"markov"`, or `"unsupported"`. An affine plan
#'   additionally records the scalar covariance coefficient transform and the
#'   affected random-effect blocks.
#'
#' @details
#' `family = "affine"` means that, conditional on every non-selected posterior
#' coordinate, the row covariance has the exact form `A + h(value) * B`, where
#' `value` is the coordinate declared by `coefficient_input` (`"source"` or
#' `"quantity"`) and `h` is recorded in `coefficient_transform`.
#' `family = "factor"` identifies an exact candidate-dependent covariance
#' factor, while `family = "markov"` identifies a compiled structured Markov
#' state. Unsupported plans carry a structural reason. The accessor never
#' estimates affineness from evaluated covariance matrices.
#'
#' @seealso [parameter_catalog()] [random_effects_marginal_factor_states()]
#' @export
random_effects_marginal_update_plan <- function(fit, selection){

  catalog <- parameter_catalog(fit)
  .bt_validate_parameter_selection(selection, catalog = catalog)
  if(nrow(selection$quantities) != 1L){
    stop("'selection' must contain exactly one parameter quantity.",
         call. = FALSE)
  }

  quantity <- selection$quantities[1L, , drop = FALSE]
  key <- quantity$extraction_key[[1L]]
  if(!identical(quantity$provider, "BayesTools") ||
     !identical(key$type, "random_summary")){
    return(.bt_random_effect_marginal_update_unavailable(
      quantity = quantity,
      reason = "not_random_effect",
      message = "The selected quantity is not a BayesTools random-effect quantity."
    ))
  }

  parameter <- quantity$formula_parameter
  design <- .bt_random_effect_marginal_covariance_design(
    fit = fit,
    parameter = parameter
  )
  plan <- .bt_random_effect_marginal_update_compile(
    fit = fit,
    selection = selection,
    quantity = quantity,
    key = key
  )
  plan$invariant_covariance <-
    .bt_random_effect_marginal_update_invariant_covariance(
      design = design,
      plan = plan,
      key = key
    )
  plan$quantity_id       <- quantity$quantity_id
  plan$canonical_name    <- quantity$canonical_name
  plan$formula_parameter <- parameter
  plan$source_parameter  <- key$source_parameter
  plan$source_transform  <- key$source_transform
  plan$dependencies      <- key$dependencies
  class(plan) <- c(
    "BayesTools_random_effects_marginal_update_plan",
    "list"
  )
  plan
}


#' Random-effect marginal covariance update grid
#'
#' @description
#' Compiles the draw- and candidate-varying coefficient state for an exact
#' non-affine random-effect covariance update. The output remains compact in
#' the number of posterior draws and grid values; it does not materialize one
#' random-effect state for every draw-by-grid combination.
#'
#' @param fit a fitted object carrying formula random-effect metadata.
#' @param update a plan returned by
#'   [random_effects_marginal_update_plan()].
#' @param values finite candidate values on the source-coordinate scale
#'   declared by `update`.
#' @param posterior_samples optional posterior sample object. By default the
#'   posterior is obtained from `fit`.
#' @param prior_list optional formula prior list used to reconstruct random SD
#'   bindings. By default it is obtained from `fit`.
#'
#' @return A list of class
#'   `BayesTools_random_effects_marginal_update_grid`. Factor updates contain
#'   posterior coefficient scales, correlation Cholesky factors, and the
#'   selected candidate scales. Markov updates contain posterior coefficient
#'   scales plus candidate correlation Cholesky factors, transitions, and
#'   innovation variances.
#'
#' @details
#' This accessor supports plans whose `family` is `"factor"` or `"markov"`.
#' Every state is evaluated by the same compiled formula machinery used by
#' [random_effects_marginal_factor_states()]. No covariance form is inferred
#' from posterior samples or candidate covariance matrices.
#'
#' @seealso [random_effects_marginal_update_plan()]
#'   [random_effects_marginal_factor_states()]
#' @export
random_effects_marginal_update_grid <- function(
    fit, update, values, posterior_samples = NULL, prior_list = NULL){

  if(!inherits(
    update,
    "BayesTools_random_effects_marginal_update_plan"
  ) || !update$family %in% c("factor", "markov")){
    stop(
      "'update' must be a factor or Markov random-effect marginal update plan.",
      call. = FALSE
    )
  }
  if(!is.numeric(values) || length(values) < 1L ||
     anyNA(values) || any(!is.finite(values))){
    stop("'values' must contain finite numeric source-coordinate values.",
         call. = FALSE)
  }
  if(!is.character(update$blocks) || length(update$blocks) != 1L ||
     is.na(update$blocks) || !nzchar(update$blocks)){
    stop("Non-affine random-effect updates must identify one fitted block.",
         call. = FALSE)
  }
  source <- update$source_parameter
  if(!is.character(source) || length(source) != 1L ||
     is.na(source) || !nzchar(source)){
    stop("Non-affine random-effect updates must identify one source coordinate.",
         call. = FALSE)
  }

  compiled <- .bt_random_effect_marginal_update_grid_evaluator(
    fit = fit,
    update = update,
    posterior_samples = posterior_samples,
    prior_list = prior_list
  )
  posterior <- compiled$posterior
  if(!source %in% colnames(posterior)){
    stop(
      "Random-effect update source coordinate '", source,
      "' is unavailable in 'posterior_samples'.",
      call. = FALSE
    )
  }
  evaluator <- compiled$evaluator
  parameter <- update$formula_parameter
  block     <- update$blocks[[1L]]
  scale <- evaluator$coefficient_scales(
    posterior = posterior,
    parameter = parameter,
    block = block
  )
  if(is.null(scale)){
    stop("Random-effect coefficient scales are unavailable for this update.",
         call. = FALSE)
  }

  if(identical(update$family, "factor")){
    component <- update$component_index
    if(!is.numeric(component) || length(component) != 1L ||
       is.na(component) || component < 1L ||
       component > ncol(scale) || component != floor(component)){
      stop("Factor updates must identify one valid coefficient component.",
           call. = FALSE)
    }
    component <- as.integer(component)
    candidate_scale <- matrix(
      NA_real_,
      nrow = length(values),
      ncol = nrow(posterior)
    )
    for(value_i in seq_along(values)){
      candidate <- posterior
      candidate[, source] <- values[[value_i]]
      current <- evaluator$coefficient_scales(
        posterior = candidate,
        parameter = parameter,
        block = block
      )
      if(is.null(current) || !identical(dim(current), dim(scale))){
        stop("Candidate random-effect coefficient scales are unavailable.",
             call. = FALSE)
      }
      candidate_scale[value_i, ] <- current[, component]
    }
    cholesky <- evaluator$coefficient_cholesky(
      posterior = posterior,
      parameter = parameter,
      block = block
    )
    if(is.null(cholesky)){
      stop("Factor updates require coefficient-correlation Cholesky states.",
           call. = FALSE)
    }
    out <- list(
      family = "factor",
      block = block,
      component_index = component,
      coefficient_scale = unname(scale),
      coefficient_cholesky = unname(cholesky),
      candidate_scale = unname(candidate_scale)
    )
  }else{
    candidate <- posterior[rep(1L, length(values)), , drop = FALSE]
    candidate[, source] <- values
    cholesky <- evaluator$coefficient_cholesky(
      posterior = candidate,
      parameter = parameter,
      block = block
    )
    if(is.null(cholesky)){
      stop("Markov updates require candidate correlation Cholesky states.",
           call. = FALSE)
    }
    n_columns <- ncol(scale)
    if(n_columns < 2L || !identical(
      dim(cholesky),
      c(length(values), n_columns, n_columns)
    )){
      stop("Candidate Markov coefficient states have inconsistent dimensions.",
           call. = FALSE)
    }
    transition <- matrix(NA_real_, nrow = length(values),
                         ncol = n_columns - 1L)
    innovation <- matrix(NA_real_, nrow = length(values),
                         ncol = n_columns - 1L)
    for(value_i in seq_along(values)){
      current <- matrix(cholesky[value_i, , ], n_columns, n_columns)
      transition[value_i, ] <-
        current[cbind(2:n_columns, seq_len(n_columns - 1L))] /
        diag(current)[seq_len(n_columns - 1L)]
      innovation[value_i, ] <- diag(current)[2:n_columns]^2
    }
    out <- list(
      family = "markov",
      block = block,
      coefficient_scale = unname(scale),
      candidate_cholesky = unname(cholesky),
      candidate_transition = unname(transition),
      candidate_innovation_variance = unname(innovation)
    )
  }
  class(out) <- c(
    "BayesTools_random_effects_marginal_update_grid",
    "list"
  )
  out
}


.bt_random_effect_marginal_update_grid_evaluator <- function(
    fit, update, posterior_samples, prior_list){

  parameter <- update$formula_parameter
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
    blocks = update$blocks
  )
  random_effects <- selected$terms
  if(length(random_effects) != 1L){
    stop("Non-affine random-effect updates must resolve to one fitted block.",
         call. = FALSE)
  }
  n_rows <- nrow(random_effects[[1L]]$model_matrix)
  row_blocks <- list(seq_len(n_rows))
  .bt_JAGS_bridge_validate_marginal_random_row_blocks(
    random_effects = random_effects,
    row_blocks = row_blocks,
    parameter = parameter
  )
  evaluator <- .bt_JAGS_bridge_compile_marginal_random_evaluator(
    formula_design_list = stats::setNames(list(design), parameter),
    marginal_random_spec = stats::setNames(list(list(
      blocks = update$blocks,
      row_blocks = row_blocks,
      factor_state = TRUE
    )), parameter),
    formula_data_list = stats::setNames(list(NULL), parameter),
    formula_prior_list = stats::setNames(list(prior_list), parameter),
    model_data = NULL
  )

  list(evaluator = evaluator, posterior = posterior)
}


.bt_random_effect_marginal_update_invariant_covariance <- function(
    design, plan, key){

  if(!identical(plan$family, "affine") ||
     !identical(plan$update, "scale") ||
     !identical(key$source_type, "identity") ||
     !is.character(key$source_parameter) ||
     length(key$source_parameter) != 1L ||
     is.na(key$source_parameter) || !nzchar(key$source_parameter)){
    return(NULL)
  }

  random_terms <- .bt_formula_design_random_effects(design)
  block_names <- vapply(random_terms, `[[`, character(1), "block_name")
  if(!setequal(plan$blocks, block_names) ||
     any(vapply(random_terms, function(random_term){
       random_term$n_columns != 1L &&
         !.bt_random_effect_structure(
           random_term,
           context = "Random-effect marginal covariance update plan"
         ) %in% c("id", "diag")
     }, logical(1))) ||
     any(vapply(random_terms, function(random_term){
       !identical(
         unique(random_term$sd_parameter_names),
         key$source_parameter
       )
     }, logical(1)))){
    return(NULL)
  }

  block_metadata <- lapply(random_terms, function(random_term){
    block_data <- .bt_random_effect_marginal_covariance_block_data(
      design = design,
      random_term = random_term,
      data = NULL
    )
    if(random_term$n_columns == 1L){
      covariance <- .bt_random_effect_marginal_variance_base_covariance(
        random_term = random_term,
        model_matrix = block_data$model_matrix,
        group_map = block_data$group_map
      )
    }else{
      group_kernel <- if(.bt_random_effect_has_known_group_covariance(
        random_term
      )){
        known <- .bt_random_effect_known_group_covariance(
          random_term,
          context = "Random-effect marginal covariance update plan"
        )
        known$kernel[block_data$group_map, block_data$group_map, drop = FALSE]
      }else{
        outer(block_data$group_map, block_data$group_map, "==") * 1
      }
      covariance <- group_kernel * tcrossprod(block_data$model_matrix)
    }
    list(
      n_rows = nrow(block_data$model_matrix),
      row_covariance = covariance
    )
  })
  row_counts <- vapply(block_metadata, `[[`, integer(1), "n_rows")
  if(length(unique(row_counts)) != 1L){
    stop(
      "Random-effect blocks do not have a common fitted row count.",
      call. = FALSE
    )
  }
  basis <- Reduce(
    `+`,
    lapply(block_metadata, `[[`, "row_covariance")
  )

  list(
    reference_coefficient = 0,
    base_covariance = matrix(0, nrow = row_counts[[1L]],
                             ncol = row_counts[[1L]]),
    update_covariance = unname(basis)
  )
}


.bt_random_effect_marginal_update_compile <- function(fit, selection, quantity,
                                                       key){

  evaluator <- key$evaluator
  if(evaluator %in% c("allocation_sd", "allocation_var")){
    random_term <- if(nzchar(key$random_block)){
      .bt_parameter_catalog_find_random_term(fit, key)
    }else{
      NULL
    }
    allocation <- .bt_parameter_catalog_find_allocation(
      fit = fit,
      key = key,
      random_term = random_term
    )
    if(identical(
      .bt_random_effect_allocation_scale_metadata(
        allocation,
        context = "Random-effect marginal covariance update plan"
      ),
      "total_variance"
    ) && length(.bt_random_effect_summary_allocation_gate_names(
      allocation
    )) > 0L){
      return(.bt_random_effect_marginal_update_unavailable(
        quantity = quantity,
        reason = "gated_realized_allocation",
        message = paste0(
          "The selected realized allocation aggregate has no single ",
          "unconditional scalar covariance update across inclusion states."
        )
      ))
    }
    direct_source <- key$source_type %in% c(
      "identity", "one_to_one_transform"
    ) && is.character(key$source_parameter) &&
      length(key$source_parameter) == 1L &&
      !is.na(key$source_parameter) && nzchar(key$source_parameter)
    coefficient_input <- if(direct_source) "source" else "quantity"
    coefficient_transform <- if(direct_source ||
                                    identical(evaluator, "allocation_sd")){
      list(type = "square")
    }else{
      list(type = "identity")
    }
    return(.bt_random_effect_marginal_update_affine(
      update = "scale",
      blocks = .bt_random_effect_marginal_update_allocation_blocks(
        fit,
        key,
        allocation
      ),
      coefficient_transform = coefficient_transform,
      coefficient_input = coefficient_input
    ))
  }

  if(identical(evaluator, "allocation") &&
     quantity$quantity %in% c("var_prop", "var_mult", "sd_mult") &&
     key$source_type %in% c("identity", "one_to_one_transform")){
    random_term <- if(nzchar(key$random_block)){
      .bt_parameter_catalog_find_random_term(fit, key)
    }else{
      NULL
    }
    allocation <- .bt_parameter_catalog_find_allocation(
      fit = fit,
      key = key,
      random_term = random_term
    )
    return(.bt_random_effect_marginal_update_affine(
      update = "allocation",
      blocks = .bt_random_effect_marginal_update_allocation_blocks(
        fit,
        key,
        allocation
      ),
      coefficient_transform = list(type = "identity"),
      coefficient_input = "source",
      allocation = list(
        weight_name = allocation$weight_name,
        index = key$index,
        n_targets = allocation$n_targets
      )
    ))
  }

  if(evaluator %in% c("sd", "sd_variance")){
    random_term <- .bt_parameter_catalog_find_random_term(fit, key)
    return(.bt_random_effect_marginal_update_sd(
      fit = fit,
      quantity = quantity,
      key = key,
      random_term = random_term
    ))
  }

  if(identical(evaluator, "rho")){
    random_term <- .bt_parameter_catalog_find_random_term(fit, key)
    structure <- .bt_random_effect_structure(
      random_term,
      context = "Random-effect marginal covariance update plan"
    )
    if(structure %in% c("cs", "hcs") ||
       (structure %in% c("ar1", "har") &&
        identical(random_term$n_columns, 2L))){
      return(.bt_random_effect_marginal_update_affine(
        update = "correlation",
        blocks = random_term$block_name,
        coefficient_transform = .bt_random_effect_marginal_update_transform(
          fit,
          selection
        ),
        coefficient_input = "source",
        component_index = NA_integer_,
        structure = structure
      ))
    }
    if(structure %in% c("ar1", "car", "har")){
      return(list(
        family = "markov",
        update = "correlation",
        blocks = random_term$block_name,
        structure = structure
      ))
    }
  }

  if(identical(evaluator, "correlation") &&
     identical(key$source_type, "one_to_one_transform")){
    random_term <- .bt_parameter_catalog_find_random_term(fit, key)
    if(identical(random_term$n_columns, 2L)){
      return(.bt_random_effect_marginal_update_affine(
        update = "correlation",
        blocks = random_term$block_name,
        coefficient_transform = .bt_random_effect_marginal_update_transform(
          fit,
          selection
        ),
        coefficient_input = "source",
        component_index = key$index,
        structure = "us"
      ))
    }
  }

  .bt_random_effect_marginal_update_unavailable(
    quantity = quantity,
    reason = "non_scalar_covariance_path",
    message = paste0(
      "The selected random-effect quantity has no declared exact scalar ",
      "marginal covariance update path."
    )
  )
}


.bt_random_effect_marginal_update_sd <- function(fit, quantity, key,
                                                  random_term){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Random-effect marginal covariance update plan"
  )
  if(isTRUE(key$allocation_derived)){
    allocations <- random_term$sd_binding$allocations
    allocations <- Filter(function(allocation){
      identical(
        .bt_random_effect_allocation_target_metadata(allocation),
        "sd_component"
      )
    }, allocations)
    if(length(allocations) != 1L){
      return(.bt_random_effect_marginal_update_unavailable(
        quantity = quantity,
        reason = "ambiguous_allocation_scale",
        message = paste0(
          "The selected allocation-derived SD does not identify one exact ",
          "variance allocation."
        )
      ))
    }
    return(.bt_random_effect_marginal_update_affine(
      update = "scale",
      blocks = .bt_random_effect_marginal_update_allocation_blocks(
        fit,
        key,
        allocations[[1L]]
      ),
      coefficient_transform = if(identical(key$evaluator, "sd_variance")){
        list(type = "identity")
      }else{
        list(type = "square")
      },
      coefficient_input = "quantity",
      component_index = key$index,
      structure = structure
    ))
  }

  n_columns <- random_term$n_columns
  homogeneous <- length(unique(random_term$sd_parameter_names)) == 1L
  affine <- n_columns == 1L || structure %in% c("id", "diag") || homogeneous
  if(isTRUE(affine)){
    return(.bt_random_effect_marginal_update_affine(
      update = if(n_columns == 1L || homogeneous) "scale" else "column_scale",
      blocks = random_term$block_name,
      coefficient_transform = list(type = "square"),
      coefficient_input = "source",
      component_index = key$index,
      structure = structure
    ))
  }

  if(structure %in% c("us", "hcs", "har") ||
     structure %in% c("cs", "ar1", "car")){
    return(list(
      family = "factor",
      update = "column_scale",
      blocks = random_term$block_name,
      component_index = key$index,
      structure = structure
    ))
  }

  .bt_random_effect_marginal_update_unavailable(
    quantity = quantity,
    reason = "unsupported_sd_structure",
    message = paste0(
      "The selected random-effect SD has no declared marginal covariance ",
      "update path."
    )
  )
}


.bt_random_effect_marginal_update_affine <- function(
    update, blocks, coefficient_transform, coefficient_input,
    component_index = NA_integer_, structure = "", allocation = NULL){

  list(
    family = "affine",
    update = update,
    blocks = unname(blocks),
    coefficient_transform = coefficient_transform,
    coefficient_input = coefficient_input,
    component_index = component_index,
    structure = structure,
    allocation = allocation
  )
}


.bt_random_effect_marginal_update_allocation_blocks <- function(
    fit, key, allocation){

  design <- .bt_random_effect_marginal_covariance_design(
    fit = fit,
    parameter = key$formula_parameter
  )
  random_terms <- .bt_formula_design_random_effects(design)
  block_names <- vapply(random_terms, `[[`, character(1), "block_name")
  allocations <- design$random_allocations
  if(is.null(allocations)){
    allocations <- list()
  }
  expand <- function(current, trail = character()){
    label <- current$label
    if(!is.character(label) || length(label) != 1L || is.na(label) ||
       !nzchar(label) || label %in% trail){
      stop(
        "Variance-allocation metadata contain an invalid dependency graph.",
        call. = FALSE
      )
    }
    terms <- unname(current$terms)
    if(!is.character(terms) || length(terms) == 0L || anyNA(terms) ||
       any(!nzchar(terms))){
      stop(
        "Variance-allocation metadata do not identify valid components.",
        call. = FALSE
      )
    }
    if(all(terms %in% block_names)){
      return(terms)
    }
    components <- current$component_labels
    if(!is.character(components) || length(components) != length(terms) ||
       anyNA(components) || any(!nzchar(components))){
      components <- names(current$terms)
    }
    if(!is.character(components) || length(components) != length(terms) ||
       anyNA(components) || any(!nzchar(components))){
      stop(
        "Variance-allocation metadata do not identify valid components.",
        call. = FALSE
      )
    }
    out <- character()
    for(i in seq_along(terms)){
      children <- Filter(function(candidate){
        parent <- candidate$parent
        is.list(parent) && identical(parent$allocation, label) &&
          identical(parent$component, components[[i]])
      }, allocations)
      if(length(children) == 1L){
        out <- c(out, expand(children[[1L]], c(trail, label)))
      }else if(length(children) == 0L && terms[[i]] %in% block_names){
        out <- c(out, terms[[i]])
      }else{
        stop(
          "Variance-allocation metadata do not resolve to fitted random-effect blocks.",
          call. = FALSE
        )
      }
    }
    out
  }

  blocks <- expand(allocation)
  if(length(blocks) == 0L || anyDuplicated(blocks) ||
     !all(blocks %in% block_names)){
    stop(
      "Variance-allocation metadata do not identify unique random-effect blocks.",
      call. = FALSE
    )
  }
  blocks
}


.bt_random_effect_marginal_update_transform <- function(fit, selection){

  parameter_transform(fit, selection)
}


.bt_random_effect_marginal_update_unavailable <- function(quantity, reason,
                                                           message){

  out <- list(
    family = "unsupported",
    update = "",
    blocks = character(),
    reason = reason,
    message = message,
    quantity_id = quantity$quantity_id,
    canonical_name = quantity$canonical_name,
    formula_parameter = quantity$formula_parameter
  )
  class(out) <- c(
    "BayesTools_random_effects_marginal_update_plan",
    "list"
  )
  out
}
