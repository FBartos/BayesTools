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
      blocks = .bt_random_effect_marginal_update_allocation_blocks(allocation),
      coefficient_transform = coefficient_transform,
      coefficient_input = coefficient_input
    ))
  }

  if(identical(evaluator, "allocation") &&
     quantity$quantity %in% c("var_prop", "var_ratio", "sd_ratio") &&
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
      blocks = .bt_random_effect_marginal_update_allocation_blocks(allocation),
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


.bt_random_effect_marginal_update_allocation_blocks <- function(allocation){

  blocks <- unname(allocation$terms)
  if(!is.character(blocks) || length(blocks) == 0L || anyNA(blocks) ||
     any(!nzchar(blocks)) || anyDuplicated(blocks)){
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
