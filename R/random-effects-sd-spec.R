# Internal random-effect SD specification helpers.

.bt_random_effect_term_indexes <- function(model_matrix, has_intercept){

  terms_indexes <- attr(model_matrix, "assign")
  if(isTRUE(has_intercept)){
    terms_indexes <- terms_indexes + 1
    terms_indexes[1] <- 0
  }

  terms_indexes
}

.bt_random_effect_sd_spec <- function(parameter, parameter_suffix,
                                      prior_list, random_term,
                                      grouping_factor,
                                      allocation_info, model_matrix,
                                      model_terms, model_terms_type,
                                      predictors_type, data,
                                      random_structure, has_intercept,
                                      homogeneous_sd){

  if(!is.null(allocation_info)){
    return(.bt_random_effect_allocated_sd_spec(
      parameter = parameter,
      parameter_suffix = parameter_suffix,
      allocation_info = allocation_info,
      model_matrix = model_matrix,
      model_terms = model_terms,
      model_terms_type = model_terms_type,
      predictors_type = predictors_type,
      data = data,
      random_structure = random_structure,
      has_intercept = has_intercept,
      homogeneous_sd = homogeneous_sd,
      structured_index = random_term$structured_index
    ))
  }

  n_columns <- ncol(model_matrix)
  terms_indexes <- .bt_random_effect_term_indexes(model_matrix, has_intercept)
  sd_parameter_names <- rep(NA_character_, n_columns)
  syntax <- character()
  new_prior_list <- list()
  random_scale_terms <- character()

  if(isTRUE(homogeneous_sd)){
    prior_name <- paste0(parameter_suffix, "_sd")
    .bt_random_effect_check_scalar_sd_prior(prior_list[["sd"]], "sd")
    new_prior_list[[prior_name]] <- .bt_random_effect_set_prior_metadata(
      prior_list[["sd"]],
      block = random_term$block_name,
      grouping = grouping_factor,
      type = "sd",
      structure = random_structure,
      name = .bt_random_effect_public_name(random_term)
    )
    random_scale_terms[[prior_name]] <- "sd"
    sd_parameter_names <- rep(paste0(parameter, "_sd"), n_columns)
    syntax <- .bt_random_effect_sd_assignment_syntax(parameter, sd_parameter_names)
  }else{
    if(isTRUE(has_intercept)){
      prior_name <- paste0(parameter_suffix, "_intercept")
      .bt_random_effect_check_scalar_sd_prior(prior_list[["intercept"]], "intercept")
      new_prior_list[[prior_name]] <- .bt_random_effect_set_prior_metadata(
        prior_list[["intercept"]],
        block = random_term$block_name,
        grouping = grouping_factor,
        type = "sd",
        structure = random_structure,
        name = .bt_random_effect_public_name(random_term)
      )
      random_scale_terms[[prior_name]] <- "intercept"
      sd_parameter_names[1L] <- paste0(parameter, "_intercept")
      syntax <- c(
        syntax,
        .bt_random_effect_sd_assignment_syntax(
          parameter,
          sd_parameter_names[1L],
          column_index = 1L
        )
      )
    }

    for(term_index in unique(terms_indexes[terms_indexes > 0])){
      term_spec <- .bt_random_effect_sd_term_spec(
        term_index = term_index,
        parameter = parameter,
        parameter_suffix = parameter_suffix,
        prior_list = prior_list,
        random_term = random_term,
        grouping_factor = grouping_factor,
        model_terms = model_terms,
        model_terms_type = model_terms_type,
        terms_indexes = terms_indexes,
        predictors_type = predictors_type,
        data = data,
        random_structure = random_structure
      )
      sd_parameter_names[term_spec$columns] <- term_spec$sd_parameter_names
      syntax <- c(syntax, term_spec$syntax)
      new_prior_list[[term_spec$prior_name]] <- term_spec$prior
      random_scale_terms[[term_spec$prior_name]] <- term_spec$scale_term
    }
  }

  sd_leaves <- .bt_random_effect_resolved_sd_leaves(
    parameter = parameter,
    model_matrix = model_matrix,
    model_terms = model_terms,
    model_terms_type = model_terms_type,
    predictors_type = predictors_type,
    data = data,
    random_structure = random_structure,
    has_intercept = has_intercept,
    homogeneous_sd = homogeneous_sd,
    structured_index = random_term$structured_index
  )
  if(!identical(unname(sd_leaves$leaf_names_by_column), unname(sd_parameter_names))){
    stop("Random-effect SD leaf metadata are inconsistent with generated SD parameters.", call. = FALSE)
  }

  list(
    terms_indexes = terms_indexes,
    sd_parameter_names = sd_parameter_names,
    sd_leaves = sd_leaves,
    random_scale_terms = random_scale_terms,
    add_parameters = character(),
    syntax = syntax,
    prior_list = new_prior_list,
    allocation_info = NULL
  )
}

.bt_random_effect_sd_term_spec <- function(term_index, parameter,
                                           parameter_suffix, prior_list,
                                           random_term, grouping_factor,
                                           model_terms, model_terms_type,
                                           terms_indexes, predictors_type,
                                           data, random_structure){

  model_term <- model_terms[term_index]
  columns <- which(terms_indexes == term_index)
  this_prior <- prior_list[[model_term]]

  attr(this_prior, "interaction") <- grepl("__xXx__", model_term)
  if(.is_prior_interaction(this_prior)){
    attr(this_prior, "interaction_terms") <- strsplit(model_term, "__xXx__")[[1]]
  }

  if(!is.null(attr(this_prior, "multiply_by"))){
    stop("'multiply_by' attribute is inadmissible for random effects", call. = FALSE)
  }

  if(identical(model_terms_type[[term_index]], "continuous")){
    .bt_random_effect_check_scalar_sd_prior(this_prior, model_term)
    sd_parameter_names <- rep(paste0(parameter, "_", model_term), length(columns))
    syntax <- .bt_random_effect_sd_assignment_syntax(
      parameter,
      sd_parameter_names,
      column_index = columns
    )
  }else if(identical(model_terms_type[[term_index]], "factor")){
    factor_spec <- .bt_random_effect_factor_sd_term_spec(
      prior = this_prior,
      parameter = parameter,
      model_term = model_term,
      columns = columns,
      terms_indexes = terms_indexes,
      term_index = term_index,
      predictors_type = predictors_type,
      data = data,
      random_structure = random_structure
    )
    this_prior <- factor_spec$prior
    sd_parameter_names <- factor_spec$sd_parameter_names
    syntax <- .bt_random_effect_sd_assignment_syntax(
      parameter,
      sd_parameter_names,
      column_index = columns
    )
  }else{
    stop("Unrecognized model term.", call. = FALSE)
  }

  this_prior <- .bt_random_effect_set_prior_metadata(
    this_prior,
    block = random_term$block_name,
    grouping = grouping_factor,
    type = "sd",
    structure = random_structure,
    name = .bt_random_effect_public_name(random_term)
  )
  this_prior <- .bt_random_effect_forward_component_metadata(this_prior)

  prior_name <- paste0(parameter_suffix, "_", model_term)
  list(
    columns = columns,
    sd_parameter_names = sd_parameter_names,
    syntax = syntax,
    prior_name = prior_name,
    prior = this_prior,
    scale_term = model_term
  )
}

.bt_random_effect_check_scalar_sd_prior <- function(prior, term){

  if(.bt_random_prior_has_non_scalar_family(prior)){
    stop(
      "Random-effect SD prior for '", term,
      "' must be an ordinary scalar prior.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_random_effect_check_structured_sd_priors <- function(prior_list,
                                                        model_terms,
                                                        random_structure,
                                                        homogeneous_sd){

  if(!random_structure %in% c("cs", "hcs", "ar1", "car", "har")){
    return(invisible(TRUE))
  }

  terms <- if(isTRUE(homogeneous_sd)) "sd" else model_terms
  for(term in terms){
    .bt_random_effect_check_scalar_sd_prior(prior_list[[term]], term)
  }

  invisible(TRUE)
}

.bt_random_effect_factor_sd_term_spec <- function(prior, parameter,
                                                  model_term, columns,
                                                  terms_indexes,
                                                  term_index,
                                                  predictors_type, data,
                                                  random_structure){

  contrast <- .bt_random_effect_factor_term_contrast(model_term, predictors_type, data)
  if(random_structure %in% c("cs", "hcs", "ar1", "car", "har")){
    .bt_random_effect_check_scalar_sd_prior(prior, model_term)
    prior_type <- "prior.independent"
    attr(prior, "levels") <- length(columns)
    attr(prior, "level_names") <- .bt_random_effect_factor_prior_level_names(
      model_term = model_term,
      predictors_type = predictors_type,
      data = data
    )
    sd_parameter_names <- paste0(parameter, "_", model_term, "[", seq_along(columns), "]")
    prior <- .bt_random_effect_set_factor_prior_class(prior, prior_type)

  }else if(contrast %in% c("contr.treatment", "contr.independent")){
    if(identical(contrast, "contr.treatment")){
      prior_type <- "prior.treatment"
      attr(prior, "levels") <- length(columns) + 1L
    }else{
      prior_type <- "prior.independent"
      attr(prior, "levels") <- length(columns)
    }

    attr(prior, "level_names") <- .bt_random_effect_factor_prior_level_names(
      model_term = model_term,
      predictors_type = predictors_type,
      data = data
    )

    sd_parameter_names <- paste0(parameter, "_", model_term, "[", seq_along(columns), "]")
    prior <- .bt_random_effect_set_factor_prior_class(prior, prior_type)

  }else if(contrast %in% c("contr.orthonormal", "contr.meandif")){
    sd_parameter_names <- rep(paste0(parameter, "_", model_term), length(columns))
  }else{
    stop("Unsupported factor contrasts for the random effects.", call. = FALSE)
  }

  list(prior = prior, sd_parameter_names = sd_parameter_names)
}

.bt_random_effect_set_factor_prior_class <- function(prior, prior_type){

  if(is.prior.simple(prior)){
    class(prior) <- c(class(prior), "prior.factor", prior_type)
  }else if(is.prior.spike_and_slab(prior) || is.prior.mixture(prior)){
    for(p in seq_along(prior)){
      class(prior[[p]]) <- c(class(prior[[p]]), "prior.factor", prior_type)
    }
    if(is.prior.spike_and_slab(prior)){
      class(prior) <- c(
        class(prior)[!class(prior) %in% c("prior.simple_spike_and_slab")],
        "prior.factor_spike_and_slab",
        prior_type
      )
    }else{
      class(prior) <- c(
        class(prior)[!class(prior) %in% c("prior.simple_mixture")],
        "prior.factor_mixture",
        prior_type
      )
    }
  }

  prior
}

.bt_random_effect_forward_component_metadata <- function(prior){

  if(is.prior.spike_and_slab(prior) || is.prior.mixture(prior)){
    for(p in seq_along(prior)){
      attr(prior[[p]], "levels") <- attr(prior, "levels")
      attr(prior[[p]], "level_names") <- attr(prior, "level_names")
      attr(prior[[p]], "interaction") <- attr(prior, "interaction")
      attr(prior[[p]], "interaction_terms") <- attr(prior, "interaction_terms")
    }
  }

  prior
}

.bt_random_effect_allocated_sd_spec <- function(parameter, parameter_suffix,
                                                allocation_info, model_matrix,
                                                model_terms, model_terms_type,
                                                predictors_type, data,
                                                random_structure,
                                                has_intercept,
                                                homogeneous_sd,
                                                structured_index = NULL){

  terms_indexes <- .bt_random_effect_term_indexes(model_matrix, has_intercept)
  if(identical(allocation_info$components, "sd")){
    return(.bt_random_effect_allocated_sd_leaves_spec(
      parameter = parameter,
      parameter_suffix = parameter_suffix,
      allocation_info = allocation_info,
      model_matrix = model_matrix,
      model_terms = model_terms,
      model_terms_type = model_terms_type,
      predictors_type = predictors_type,
      data = data,
      random_structure = random_structure,
      has_intercept = has_intercept,
      homogeneous_sd = homogeneous_sd,
      structured_index = structured_index
    ))
  }

  if(!(isTRUE(homogeneous_sd) || ncol(model_matrix) == 1L)){
    stop(
      "Block-level variance allocation supports random-effect blocks with one SD component. ",
      "Use random intercepts, single-column random slopes, homogeneous SD structures, or components = \"sd\" for heterogeneous blocks.",
      call. = FALSE
    )
  }

  sd_leaves <- .bt_random_effect_resolved_sd_leaves(
    parameter = parameter,
    model_matrix = model_matrix,
    model_terms = model_terms,
    model_terms_type = model_terms_type,
    predictors_type = predictors_type,
    data = data,
    random_structure = random_structure,
    has_intercept = has_intercept,
    homogeneous_sd = homogeneous_sd,
    structured_index = structured_index
  )

  if(isTRUE(homogeneous_sd)){
    sd_term <- "sd"
  }else if(isTRUE(has_intercept)){
    sd_term <- "intercept"
  }else{
    allocated_term_indexes <- unique(terms_indexes[terms_indexes > 0])
    if(length(allocated_term_indexes) != 1L){
      stop(
        "Variance allocation for non-intercept random-effect blocks requires exactly one random-effect term.",
        call. = FALSE
      )
    }
    sd_term <- model_terms[allocated_term_indexes]
  }

  sd_parameter <- if(ncol(model_matrix) == 1L){
    sd_leaves$leaf_names_by_column[[1L]]
  }else{
    paste0(parameter, "_", sd_term)
  }
  sd_suffix <- paste0(parameter_suffix, "_", sd_term)
  n_columns <- ncol(model_matrix)
  sd_parameter_names <- rep(sd_parameter, n_columns)

  if(!identical(unname(sd_leaves$leaf_names_by_column), unname(sd_parameter_names))){
    stop("Random-effect SD leaf metadata are inconsistent with generated SD parameters.", call. = FALSE)
  }

  syntax <- paste0(sd_parameter, " = ", allocation_info$source_expression)
  syntax <- c(
    syntax,
    .bt_random_effect_sd_assignment_syntax(parameter, sd_parameter_names)
  )

  random_scale_terms <- stats::setNames(sd_term, sd_suffix)

  list(
    terms_indexes = terms_indexes,
    sd_parameter_names = sd_parameter_names,
    sd_leaves = sd_leaves,
    random_scale_terms = random_scale_terms,
    add_parameters = sd_parameter,
    syntax = syntax,
    prior_list = list(),
    allocation_info = allocation_info
  )
}

.bt_random_effect_resolved_sd_leaves <- function(parameter, model_matrix,
                                                 model_terms,
                                                 model_terms_type,
                                                 predictors_type, data,
                                                 random_structure,
                                                 has_intercept,
                                                 homogeneous_sd,
                                                 structured_index = NULL){

  n_columns <- ncol(model_matrix)
  terms_indexes <- .bt_random_effect_term_indexes(model_matrix, has_intercept)
  leaf_terms <- rep(NA_character_, n_columns)
  leaf_names <- rep(NA_character_, n_columns)

  if(isTRUE(homogeneous_sd)){
    leaf_terms[] <- "sd"
    leaf_names[] <- paste0(parameter, "_sd")
  }else{
    if(isTRUE(has_intercept)){
      leaf_terms[1L] <- "intercept"
      leaf_names[1L] <- paste0(parameter, "_intercept")
    }
    for(i in unique(terms_indexes[terms_indexes > 0])){
      columns <- which(terms_indexes == i)
      model_term <- model_terms[i]
      if(identical(model_terms_type[[i]], "continuous")){
        leaf_terms[columns] <- model_term
        leaf_names[columns] <- paste0(parameter, "_", model_term)
      }else if(identical(model_terms_type[[i]], "factor")){
        contrast <- .bt_random_effect_factor_term_contrast(model_term, predictors_type, data)
        if(random_structure %in% c("cs", "hcs", "ar1", "car", "har") ||
           contrast %in% c("contr.treatment", "contr.independent")){
          leaf_term_labels <- .bt_random_effect_factor_sd_leaf_terms(
            model_term = model_term,
            columns = columns,
            predictors_type = predictors_type,
            data = data,
            contrast = contrast,
            random_structure = random_structure
          )
          for(j in seq_along(columns)){
            leaf_terms[columns[j]] <- leaf_term_labels[j]
            leaf_names[columns[j]] <- paste0(parameter, "_", model_term, "[", j, "]")
          }
        }else if(contrast %in% c("contr.orthonormal", "contr.meandif")){
          leaf_terms[columns] <- model_term
          leaf_names[columns] <- paste0(parameter, "_", model_term)
        }else{
          stop("Unsupported factor contrasts for the random effects.", call. = FALSE)
        }
      }else{
        stop("Unrecognized model term.", call. = FALSE)
      }
    }
  }

  if(any(is.na(leaf_names))){
    stop("Random-effect SD leaf metadata could not be resolved.", call. = FALSE)
  }
  unique_leaf_names <- unique(leaf_names)
  unique_leaf_terms <- vapply(unique_leaf_names, function(leaf_name){
    leaf_terms[match(leaf_name, leaf_names)]
  }, character(1))
  leaf_index_by_column <- match(leaf_names, unique_leaf_names)

  out <- list(
    terms_indexes = terms_indexes,
    leaf_names_by_column = leaf_names,
    leaf_terms_by_column = leaf_terms,
    leaf_names = unique_leaf_names,
    leaf_terms = unique_leaf_terms,
    leaf_index_by_column = leaf_index_by_column,
    column_names = colnames(model_matrix),
    column_index = seq_len(n_columns),
    n_leaves = length(unique_leaf_names)
  )
  class(out) <- c("BayesTools_random_effect_sd_leaves", "list")

  out
}

.bt_random_effect_factor_prior_level_names <- function(model_term,
                                                       predictors_type,
                                                       data){

  is_interaction <- grepl("__xXx__", model_term, fixed = TRUE)
  components <- .bt_random_effect_term_components(model_term)
  factor_components <- components[
    components %in% names(predictors_type) &
      predictors_type[components] == "factor"
  ]

  if(length(factor_components) == 0L){
    return(NULL)
  }

  if(is_interaction){
    level_names <- lapply(factor_components, function(factor_component){
      levels(data[, factor_component])
    })
    names(level_names) <- factor_components
    return(level_names)
  }

  levels(data[, factor_components[1L]])
}

.bt_random_effect_factor_sd_leaf_terms <- function(model_term, columns,
                                                   predictors_type, data,
                                                   contrast,
                                                   random_structure){

  n_parameters <- length(columns)
  level_names <- .bt_random_effect_factor_prior_level_names(
    model_term = model_term,
    predictors_type = predictors_type,
    data = data
  )

  if(is.null(level_names)){
    return(paste0(model_term, "[", seq_len(n_parameters), "]"))
  }

  if(is.list(level_names)){
    candidate_level_names <- list(level_names)
    treatment_level_names <- lapply(level_names, function(level_name){
      level_name[-1L]
    })
    if(all(lengths(treatment_level_names) > 0L)){
      candidate_level_names <- c(candidate_level_names, list(treatment_level_names))
    }
    if(identical(contrast, "contr.treatment") &&
       !random_structure %in% c("cs", "hcs", "ar1", "car", "har")){
      candidate_level_names <- rev(candidate_level_names)
    }

    for(candidate in candidate_level_names){
      if(prod(lengths(candidate)) != n_parameters){
        next
      }
      return(.format_factor_level_parameter_names(
        parameter = model_term,
        level_names = candidate,
        n_parameters = n_parameters
      ))
    }

    return(paste0(model_term, "[", seq_len(n_parameters), "]"))
  }

  level_names <- as.character(level_names)
  if(length(level_names) == n_parameters){
    parameter_levels <- level_names
  }else if(length(level_names) - 1L == n_parameters){
    parameter_levels <- level_names[-1L]
  }else{
    parameter_levels <- seq_len(n_parameters)
  }

  paste0(model_term, "[", parameter_levels, "]")
}

.bt_random_effect_sd_assignment_syntax <- function(parameter, sd_parameter_names,
                                                  column_index = seq_along(sd_parameter_names)){

  paste0(
    parameter,
    "_xRE_STDx[",
    column_index,
    "] = ",
    sd_parameter_names
  )
}

.bt_random_effect_allocated_sd_leaves_spec <- function(parameter,
                                                       parameter_suffix,
                                                       allocation_info,
                                                       model_matrix,
                                                       model_terms,
                                                       model_terms_type,
                                                       predictors_type,
                                                       data,
                                                       random_structure,
                                                       has_intercept,
                                                       homogeneous_sd,
                                                       structured_index = NULL){

  leaves <- .bt_random_effect_resolved_sd_leaves(
    parameter = parameter,
    model_matrix = model_matrix,
    model_terms = model_terms,
    model_terms_type = model_terms_type,
    predictors_type = predictors_type,
    data = data,
    random_structure = random_structure,
    has_intercept = has_intercept,
    homogeneous_sd = homogeneous_sd,
    structured_index = structured_index
  )
  K <- length(leaves$leaf_names)
  if(K < 2L){
    stop(
      "SD-component variance allocation requires at least two resolved SD components.",
      call. = FALSE
    )
  }

  allocation_prior <- .bt_random_variance_allocation_prior(allocation_info, K)
  allocation_prior <- .bt_random_effect_set_allocation_metadata(
    allocation_prior,
    allocation = allocation_info$label,
    terms = allocation_info$terms
  )

  syntax <- character()
  for(leaf_i in seq_len(K)){
    expression <- .bt_random_variance_allocation_expression(
      source_name = allocation_info$source_node,
      weight_name = allocation_info$weight_name,
      index = leaf_i,
      scale = allocation_info$scale,
      n_targets = K
    )
    syntax <- c(syntax, paste0(leaves$leaf_names[leaf_i], " = ", expression))
  }
  syntax <- c(
    syntax,
    .bt_random_effect_sd_assignment_syntax(parameter, leaves$leaf_names_by_column)
  )

  factor <- .bt_random_variance_allocation_factor(
    weight_name = allocation_info$weight_name,
    index = NA_integer_,
    scale = allocation_info$scale,
    n_targets = K
  )
  allocation_info$n_targets <- K
  allocation_info$leaf_names <- leaves$leaf_names
  allocation_info$leaf_terms <- leaves$leaf_terms
  allocation_info$leaf_index_by_column <- leaves$leaf_index_by_column
  allocation_info$factor_template <- factor
  allocation_info$factors <- allocation_info$parent_factors

  random_scale_terms <- stats::setNames(
    leaves$leaf_terms,
    sub(paste0("^", parameter), parameter_suffix, leaves$leaf_names)
  )
  prior_list <- list(allocation_prior)
  names(prior_list) <- allocation_info$weight_suffix

  list(
    terms_indexes = leaves$terms_indexes,
    sd_parameter_names = leaves$leaf_names_by_column,
    sd_leaves = leaves,
    random_scale_terms = random_scale_terms,
    add_parameters = leaves$leaf_names,
    syntax = syntax,
    prior_list = prior_list,
    allocation_info = allocation_info
  )
}
