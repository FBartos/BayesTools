# Internal random-effect summary helpers.

.bt_random_effect_summary_samples <- function(model_samples, prior_list,
                                              formula_design = NULL,
                                              coordinates = NULL,
                                              mode = c("standard", "full", "raw", "none"),
                                              formula_scale = NULL){

  mode <- match.arg(mode)
  if(is.null(coordinates)){
    coordinates <- .bt_build_parameter_coordinates(
      columns = colnames(model_samples),
      prior_list = prior_list,
      formula_design = formula_design,
      formula_scale = formula_scale
    )
  }
  .bt_validate_parameter_coordinates(coordinates)
  if(is.null(formula_design) || length(formula_design) == 0L || identical(mode, "raw")){
    return(list(model_samples = model_samples, prior_list = prior_list))
  }

  random_design <- .bt_random_effect_summary_designs(formula_design)
  if(length(random_design) == 0L){
    return(list(model_samples = model_samples, prior_list = prior_list))
  }

  model_samples <- as.matrix(model_samples)

  if(!identical(mode, "none")){
    derived <- .bt_random_effect_summary_derived_samples(
      model_samples = model_samples,
      prior_list = prior_list,
      random_design = random_design,
      mode = mode,
      formula_scale = formula_scale
    )
    if(ncol(derived$model_samples) > 0L){
      model_samples <- cbind(model_samples, derived$model_samples)
      prior_list <- c(prior_list, derived$prior_list)
    }
  }

  .bt_random_effect_summary_remove_raw(
    model_samples = model_samples,
    prior_list = prior_list,
    coordinates = coordinates
  )
}

.bt_random_effect_summary_designs <- function(formula_design){

  if(inherits(formula_design, "BayesTools_formula_design")){
    formula_design <- list(formula_design)
  }
  if(!is.list(formula_design)){
    return(list())
  }

  formula_design[vapply(formula_design, .bt_formula_design_has_any_random_effects, logical(1))]
}

.bt_random_effect_summary_derived_samples <- function(model_samples, prior_list,
                                                      random_design, mode,
                                                      formula_scale = NULL){

  columns <- list()
  summary_priors <- list()
  used_names <- character()

  add_summary <- function(name, values, parameter, type, label,
                          component_label = NULL,
                          block = NULL, grouping = NULL,
                          structure = NULL, effect_label = NULL,
                          allocation = NULL,
                          allocation_metadata = NULL,
                          allocation_index = NULL,
                          component = NULL){
    values <- .bt_random_effect_summary_validate_values(
      values = values,
      name = name,
      label = label,
      n_draws = nrow(model_samples)
    )
    name <- .bt_random_effect_summary_unique_name(name, used_names)
    used_names <<- c(used_names, name)
    columns[[name]] <<- values
    summary_priors[[name]] <<- .bt_random_effect_summary_prior(
      parameter = parameter,
      type = type,
      label = label,
      component_label = component_label,
      block = block,
      grouping = grouping,
      structure = structure,
      effect_label = effect_label,
      allocation = allocation,
      allocation_metadata = allocation_metadata,
      allocation_index = allocation_index,
      component = component
    )
    invisible(NULL)
  }

  seen_allocations <- character()
  add_allocation_summary <- function(allocation, parameter, block = NULL,
                                     grouping = NULL, random_term = NULL){
    if(is.null(allocation) || allocation$weight_name %in% seen_allocations){
      return(invisible(NULL))
    }
    scale_values <- .bt_random_effect_summary_allocation_scale_samples(
      allocation = allocation,
      model_samples = model_samples,
      prior_list = prior_list
    )
    allocation_owner <- .bt_random_effect_allocation_public_name(allocation)
    if(!is.null(scale_values)){
      sd_quantity <- .bt_random_effect_allocation_sd_quantity(allocation)
      add_summary(
        name = .bt_random_effect_summary_name(
          parameter = parameter,
          type = sd_quantity,
          parts = allocation$label
        ),
        values = scale_values,
        parameter = parameter,
        type = sd_quantity,
        label = .bt_random_effect_semantic_name(
          parameter = "",
          owner = allocation_owner,
          quantity = sd_quantity,
          formula_prefix = FALSE
        ),
        allocation = allocation$label,
        allocation_metadata = allocation
      )
      var_quantity <- .bt_random_effect_allocation_var_quantity(allocation)
      add_summary(
        name = .bt_random_effect_summary_name(
          parameter = parameter,
          type = var_quantity,
          parts = allocation$label
        ),
        values = scale_values^2,
        parameter = parameter,
        type = var_quantity,
        label = .bt_random_effect_semantic_name(
          parameter = "",
          owner = allocation_owner,
          quantity = var_quantity,
          formula_prefix = FALSE
        ),
        allocation = allocation$label,
        allocation_metadata = allocation
      )
    }
    allocation_summary <- .bt_random_effect_summary_allocation_samples(
      allocation = allocation,
      random_term = random_term,
      model_samples = model_samples,
      prior_list = prior_list,
      include_multipliers = mode %in% c("standard", "full")
    )
    allocation_target <- .bt_random_effect_summary_allocation_target(allocation)
    for(i in seq_along(allocation_summary$names)){
      add_summary(
        name = allocation_summary$names[i],
        values = allocation_summary$values[, i],
        parameter = parameter,
        type = allocation_summary$types[i],
        label = allocation_summary$labels[i],
        component_label = NULL,
        block = if(identical(allocation_target, "sd_component")) block else NULL,
        grouping = if(identical(allocation_target, "sd_component")) grouping else NULL,
        structure = if(identical(allocation_target, "sd_component") && !is.null(random_term)) {
          .bt_random_effect_summary_term_structure(random_term)
        }else{
          NULL
        },
        effect_label = if(identical(allocation_target, "sd_component") && !is.null(random_term)) .bt_random_effect_public_name(random_term) else NULL,
        allocation = allocation$label,
        allocation_metadata = allocation,
        allocation_index = allocation_summary$indices[i],
        component = allocation_summary$components[i]
      )
    }
    inclusion_summary <- .bt_random_effect_summary_allocation_inclusion_samples(
      allocation = allocation,
      model_samples = model_samples
    )
    for(i in seq_along(inclusion_summary$names)){
      add_summary(
        name = inclusion_summary$names[i],
        values = inclusion_summary$values[, i],
        parameter = parameter,
        type = inclusion_summary$types[i],
        label = inclusion_summary$labels[i],
        component_label = inclusion_summary$labels[i],
        allocation = allocation$label,
        allocation_metadata = allocation,
        allocation_index = inclusion_summary$indices[i],
        component = inclusion_summary$components[i]
      )
    }
    seen_allocations <<- c(seen_allocations, allocation$weight_name)
    invisible(NULL)
  }

  for(design in random_design){
    parameter <- design$parameter
    for(random_term in design$random_effects){
      display_structure <- .bt_random_effect_summary_term_structure(random_term)
      sd_summary <- .bt_random_effect_summary_sd_samples(
        random_term = random_term,
        model_samples = model_samples,
        prior_list = prior_list,
        parameter = parameter,
        formula_scale = formula_scale
      )
      correlation_model_samples <- .bt_random_effect_summary_complete_scaled_samples(
        random_term = random_term,
        model_samples = model_samples,
        prior_list = prior_list,
        parameter = parameter,
        formula_scale = formula_scale
      )
      for(i in seq_along(sd_summary$names)){
        add_summary(
          name = .bt_random_effect_summary_name(
            parameter = parameter,
            type = "sd",
            parts = c(random_term$block_name, sd_summary$components[i])
          ),
          values = sd_summary$values[, i],
          parameter = parameter,
          type = "sd",
          label = .bt_random_effect_sd_summary_label(
            component = sd_summary$components[i],
            random_term = random_term
          ),
          component_label = .bt_random_effect_sd_component_summary_label(
            component = sd_summary$components[i]
          ),
          block = random_term$block_name,
          grouping = random_term$group_label,
          structure = display_structure,
          effect_label = .bt_random_effect_public_name(random_term),
          component = sd_summary$components[i]
        )
      }

      inclusion_summary <- .bt_random_effect_summary_inclusion_samples(
        random_term = random_term,
        model_samples = model_samples,
        prior_list = prior_list,
        parameter = parameter
      )
      for(i in seq_along(inclusion_summary$names)){
        add_summary(
          name = inclusion_summary$names[i],
          values = inclusion_summary$values[, i],
          parameter = parameter,
          type = inclusion_summary$types[i],
          label = inclusion_summary$labels[i],
          component_label = inclusion_summary$component_labels[i],
          block = random_term$block_name,
          grouping = random_term$group_label,
          structure = display_structure,
          effect_label = .bt_random_effect_public_name(random_term),
          component = inclusion_summary$components[i]
        )
      }

      rho <- .bt_random_effect_summary_rho_samples(random_term, model_samples)
      if(!is.null(rho)){
        add_summary(
          name = .bt_random_effect_summary_name(
            parameter = parameter,
            type = "cor",
            parts = random_term$block_name
          ),
          values = rho,
          parameter = parameter,
          type = "cor",
          label = .bt_random_effect_semantic_name(
            parameter = "",
            owner = .bt_random_effect_public_name(random_term),
            quantity = "cor",
            formula_prefix = FALSE
          ),
          component_label = "cor",
          block = random_term$block_name,
          grouping = random_term$group_label,
          structure = display_structure,
          effect_label = .bt_random_effect_public_name(random_term)
        )
      }

      correlation_summary <- .bt_random_effect_summary_correlation_samples(
        random_term = random_term,
        model_samples = correlation_model_samples
      )
      for(i in seq_along(correlation_summary$labels)){
        add_summary(
          name = .bt_random_effect_summary_name(
            parameter = parameter,
            type = "cor",
            parts = c(random_term$block_name, correlation_summary$parts[[i]])
          ),
          values = correlation_summary$values[, i],
          parameter = parameter,
          type = "cor",
          label = .bt_random_effect_semantic_name(
            parameter = "",
            owner = .bt_random_effect_public_name(random_term),
            quantity = "cor",
            arguments = correlation_summary$parts[[i]],
            formula_prefix = FALSE
          ),
          component_label = paste0("cor(", correlation_summary$labels[i], ")"),
          block = random_term$block_name,
          grouping = random_term$group_label,
          structure = display_structure,
          effect_label = .bt_random_effect_public_name(random_term),
          component = correlation_summary$labels[i]
        )
      }

      if(!is.null(random_term$sd_binding) &&
         length(random_term$sd_binding$allocations) > 0L){
        for(allocation in random_term$sd_binding$allocations){
          add_allocation_summary(
            allocation = allocation,
            parameter = parameter,
            block = random_term$block_name,
            grouping = random_term$group_label,
            random_term = random_term
          )
        }
      }
    }
    for(allocation in design$random_allocations){
      add_allocation_summary(
        allocation = allocation,
        parameter = parameter
      )
    }
  }

  if(length(columns) == 0L){
    summary_matrix <- matrix(nrow = nrow(model_samples), ncol = 0L)
  }else{
    summary_matrix <- do.call(cbind, columns)
    summary_order <- .bt_random_effect_summary_column_order(
      summary_priors = summary_priors,
      random_design = random_design
    )
    summary_matrix <- summary_matrix[, summary_order, drop = FALSE]
    summary_priors <- summary_priors[summary_order]
  }
  list(model_samples = summary_matrix, prior_list = summary_priors)
}

.bt_random_effect_summary_column_order <- function(summary_priors,
                                                   random_design){

  parameter_names <- names(summary_priors)
  if(length(parameter_names) == 0L){
    return(parameter_names)
  }

  summary_type <- .bt_random_effect_summary_prior_attributes(
    summary_priors,
    "random_summary"
  )
  summary_parameter <- .bt_random_effect_summary_prior_attributes(
    summary_priors,
    "parameter"
  )
  summary_allocation <- .bt_random_effect_summary_prior_attributes(
    summary_priors,
    "random_allocation"
  )
  is_allocation_summary <- summary_type %in% c(
    "var_prop",
    "var_ratio",
    "sd_ratio"
  ) | (summary_type == "inclusion" & nzchar(summary_allocation))
  used <- stats::setNames(rep(FALSE, length(parameter_names)), parameter_names)
  ordered <- character()

  add_matches <- function(matches){
    matches[is.na(matches)] <- FALSE
    selected <- parameter_names[matches & !used]
    if(length(selected) > 0L){
      ordered <<- c(ordered, selected)
      used[selected] <<- TRUE
    }
    invisible(NULL)
  }

  for(design in random_design){
    parameter <- design$parameter
    parameter_match <- summary_parameter == parameter

    for(allocation in design$random_allocations){
      add_matches(
        parameter_match &
          summary_type %in% c("sd_total", "var_total", "sd_common", "var_common") &
          summary_allocation == allocation$label
      )
    }
    add_matches(parameter_match & summary_type %in%
      c("sd_total", "var_total", "sd_common", "var_common"))

    for(random_term in design$random_effects){
      add_matches(
        parameter_match &
          !is_allocation_summary &
          !summary_type %in% c("sd_total", "var_total", "sd_common", "var_common") &
          .bt_random_effect_summary_column_matches_term(
            summary_priors,
            random_term
          )
      )
    }

    add_matches(
      parameter_match &
        !is_allocation_summary &
        !summary_type %in% c("sd_total", "var_total", "sd_common", "var_common")
    )
    for(allocation in design$random_allocations){
      add_matches(
        parameter_match &
          is_allocation_summary &
          summary_allocation == allocation$label
      )
    }
    add_matches(parameter_match & is_allocation_summary)
  }

  add_matches(rep(TRUE, length(parameter_names)))
  ordered
}

.bt_random_effect_summary_prior_attributes <- function(summary_priors,
                                                       attribute){

  vapply(summary_priors, function(prior){
    value <- attr(prior, attribute, exact = TRUE)
    if(is.null(value) || length(value) == 0L || is.na(value[1L])){
      return("")
    }
    as.character(value[1L])
  }, character(1))
}

.bt_random_effect_summary_column_matches_term <- function(summary_priors,
                                                          random_term){

  term_names <- unique(c(
    random_term$block_name,
    .bt_random_effect_public_name(random_term),
    random_term$group_label,
    .bt_random_effect_summary_group_label(random_term)
  ))
  term_names <- term_names[!is.na(term_names) & nzchar(term_names)]

  if(length(term_names) == 0L){
    return(rep(FALSE, length(summary_priors)))
  }

  random_factor <- .bt_random_effect_summary_prior_attributes(
    summary_priors,
    "random_factor"
  )
  random_name <- .bt_random_effect_summary_prior_attributes(
    summary_priors,
    "random_name"
  )
  random_grouping <- .bt_random_effect_summary_prior_attributes(
    summary_priors,
    "random_grouping_factor"
  )

  random_factor %in% term_names |
    random_name %in% term_names |
    random_grouping %in% term_names
}
