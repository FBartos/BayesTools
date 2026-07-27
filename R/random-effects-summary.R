# Internal random-effect summary helpers.

.bt_random_effect_summary_samples <- function(model_samples, prior_list,
                                              formula_design = NULL,
                                              parameter_registry = NULL,
                                              mode = c("standard", "full", "raw", "none"),
                                              formula_scale = NULL){

  mode <- match.arg(mode)
  if(is.null(parameter_registry)){
    parameter_registry <- .bt_build_parameter_registry(
      columns = colnames(model_samples),
      prior_list = prior_list,
      formula_design = formula_design,
      formula_scale = formula_scale
    )
  }
  .bt_validate_parameter_registry(parameter_registry)
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
    parameter_registry = parameter_registry
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

  for(prior_name in names(prior_list)){
    prior <- prior_list[[prior_name]]
    if(isTRUE(attr(prior, "random_sd_total"))){
      allocation <- attr(prior, "random_allocation")
      values <- .bt_random_effect_parameter_draws(prior_name, model_samples, prior_list)
      if(is.null(values)){
        .bt_random_effect_summary_missing_sd_total_stop(allocation)
      }
      add_summary(
        name = .bt_random_effect_summary_name(
          parameter = attr(prior, "parameter"),
          type = "sd_total",
          parts = allocation
        ),
        values = values,
        parameter = attr(prior, "parameter"),
        type = "sd_total",
        label = paste0("sd_total(", allocation, ")"),
        allocation = allocation
      )
    }
  }

  seen_allocations <- character()
  add_allocation_summary <- function(allocation, parameter, block = NULL,
                                     grouping = NULL, random_term = NULL){
    if(is.null(allocation) || allocation$weight_name %in% seen_allocations){
      return(invisible(NULL))
    }
    allocation_summary <- .bt_random_effect_summary_allocation_samples(
      allocation = allocation,
      random_term = random_term,
      model_samples = model_samples,
      prior_list = prior_list,
      include_multipliers = identical(mode, "full")
    )
    allocation_target <- .bt_random_effect_summary_allocation_target(allocation)
    for(i in seq_along(allocation_summary$names)){
      add_summary(
        name = allocation_summary$names[i],
        values = allocation_summary$values[, i],
        parameter = parameter,
        type = allocation_summary$types[i],
        label = allocation_summary$labels[i],
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
      display_group <- .bt_random_effect_summary_group_label(random_term)
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
            group = display_group,
            random_term = random_term
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
            type = "rho",
            parts = random_term$block_name
          ),
          values = rho,
          parameter = parameter,
          type = "rho",
          label = paste0("rho(", display_group, ")"),
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
          label = paste0("cor(", correlation_summary$labels[i], " | ", display_group, ")"),
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
    "var_frac",
    "var_ratio",
    "sd_multiplier"
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
          summary_type == "sd_total" &
          summary_allocation == allocation$label
      )
    }
    add_matches(parameter_match & summary_type == "sd_total")

    for(random_term in design$random_effects){
      add_matches(
        parameter_match &
          !is_allocation_summary &
          summary_type != "sd_total" &
          .bt_random_effect_summary_column_matches_term(
            summary_priors,
            random_term
          )
      )
    }

    add_matches(
      parameter_match &
        !is_allocation_summary &
        summary_type != "sd_total"
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
