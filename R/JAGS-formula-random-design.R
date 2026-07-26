.bt_random_effects_interface <- function(random_effects, prior_random = NULL){

  if(length(random_effects) == 0L){
    return("none")
  }
  if(!is.null(prior_random)){
    return("prior_random")
  }

  stop("Formula random effects require 'prior_random'.", call. = FALSE)
}

.bt_random_effect_dense_complexity <- function(structure, n_groups,
                                               n_columns, n_rows,
                                               monitor_policy,
                                               centered = FALSE){

  is_unstructured <- identical(structure, "us") && n_columns > 1L
  is_structured   <- structure %in% c("cs", "hcs", "ar1", "car", "har") &&
    n_columns > 1L

  cholesky_products <- if(is_structured && isTRUE(centered)){
    n_columns * (n_columns - 1) * (n_columns + 1) / 6
  }else if(is_unstructured){
    n_columns^2
  }else{
    0
  }
  transform_products <- if(is_unstructured){
    n_groups * n_columns^2 + n_rows * n_columns
  }else if(is_structured && isTRUE(centered)){
    n_groups * n_columns^2 + n_rows
  }else if(is_structured){
    n_groups * n_columns + n_rows
  }else{
    n_groups * n_columns + n_rows * n_columns
  }
  monitored_values <- 0
  if(isTRUE(monitor_policy$latent)){
    monitored_values <- monitored_values + n_groups * n_columns
  }
  if(isTRUE(monitor_policy$coefficients)){
    monitored_values <- monitored_values + n_groups * n_columns
  }
  if(is_unstructured){
    monitored_values <- monitored_values + n_columns^2
    if(isTRUE(monitor_policy$correlation)){
      monitored_values <- monitored_values + n_columns^2
    }
    monitored_values <- monitored_values + n_columns * (n_columns - 1) / 2
  }else if(is_structured && isTRUE(centered)){
    monitored_values <- monitored_values + 2 * n_columns^2
  }else if(is_structured){
    monitored_values <- monitored_values + 1L
  }

  c(
    cholesky_products = cholesky_products,
    transform_products = transform_products,
    monitored_values = monitored_values
  )
}

# Estimate emitted syntax and monitored nodes for a group-local block.
.bt_random_effect_group_local_complexity <- function(
    layout, n_rows, monitor_policy, row_indexed_external_sd = FALSE,
    column_allocation = FALSE){

  if(!inherits(layout, "BayesTools_random_effect_structured_local_layout")){
    stop("'layout' must be a structured local layout.", call. = FALSE)
  }
  check_int(n_rows, "n_rows", lower = 1, allow_NA = FALSE)
  check_bool(
    row_indexed_external_sd,
    "row_indexed_external_sd",
    allow_NA = FALSE
  )
  check_bool(column_allocation, "column_allocation", allow_NA = FALSE)

  n_local           <- layout$n_local
  correlation_nodes <- if(layout$global_n_columns > 1L) 1L else 0L
  prefix_nodes      <- if(layout$structure %in% c("cs", "hcs")){
    sum(pmax(lengths(layout$group_columns) - 1L, 0L))
  }else{
    0L
  }
  syntax_nodes     <-
    n_local * (2L + !isTRUE(row_indexed_external_sd)) +
    prefix_nodes + n_rows +
    if(isTRUE(column_allocation)) layout$global_n_columns else 0L
  monitored_values <- correlation_nodes
  if(isTRUE(monitor_policy$latent)){
    monitored_values <- monitored_values + n_local
  }
  if(isTRUE(monitor_policy$coefficients)){
    monitored_values <- monitored_values + n_local
  }
  if(isTRUE(column_allocation)){
    monitored_values <- monitored_values + layout$global_n_columns
  }

  c(
    syntax_nodes = syntax_nodes,
    monitored_values = monitored_values
  )
}

# Stop before a group-local block exceeds configurable compiler limits.
.bt_random_effect_check_group_local_complexity <- function(
    random_term, layout, n_rows, monitor_policy,
    row_indexed_external_sd = FALSE, column_allocation = FALSE){

  multiplier <- getOption(
    "BayesTools.random_effects_complexity_multiplier",
    1
  )
  if(!is.numeric(multiplier) || length(multiplier) != 1L || is.na(multiplier) ||
     multiplier <= 0){
    stop(
      "Option 'BayesTools.random_effects_complexity_multiplier' must be a positive numeric scalar.",
      call. = FALSE
    )
  }
  if(is.infinite(multiplier)){
    return(invisible(NULL))
  }

  estimates <- .bt_random_effect_group_local_complexity(
    layout = layout,
    n_rows = n_rows,
    monitor_policy = monitor_policy,
    row_indexed_external_sd = row_indexed_external_sd,
    column_allocation = column_allocation
  )
  limits   <- multiplier * c(
    syntax_nodes = 2e5,
    monitored_values = 2e4
  )
  exceeded <- estimates > limits
  if(!any(exceeded)){
    return(invisible(estimates))
  }

  details <- paste0(
    names(estimates)[exceeded], "=",
    format(estimates[exceeded], scientific = FALSE, trim = TRUE),
    " (limit ",
    format(limits[exceeded], scientific = FALSE, trim = TRUE),
    ")"
  )
  stop(
    "Random-effect block '", random_term$block_name,
    "' (", toupper(layout$structure), "; ", layout$n_local,
    " active group-index cells) exceeds the current group-local JAGS compiler limits: ",
    paste(details, collapse = ", "), ". This representation emits syntax ",
    "and monitored nodes for each active cell and can exhaust memory before ",
    "sampling. Reduce the structured index dimension or active cells, or ",
    "disable latent monitoring when coefficient reconstruction and bridge ",
    "sampling are not needed.",
    call. = FALSE
  )
}

.bt_random_effect_check_dense_complexity <- function(random_term, structure,
                                                      n_groups, n_columns,
                                                      n_rows, monitor_policy,
                                                      centered = FALSE){

  status <- .bt_random_effect_dense_complexity_status(
    structure = structure,
    n_groups = n_groups,
    n_columns = n_columns,
    n_rows = n_rows,
    monitor_policy = monitor_policy,
    centered = centered
  )
  if(all(is.infinite(status$limits))){
    return(invisible(NULL))
  }
  if(!any(status$exceeded)){
    return(invisible(status$estimates))
  }

  estimates <- status$estimates
  limits    <- status$limits
  exceeded  <- status$exceeded
  details <- paste0(
    names(estimates)[exceeded], "=",
    format(estimates[exceeded], scientific = FALSE, trim = TRUE),
    " (limit ",
    format(limits[exceeded], scientific = FALSE, trim = TRUE),
    ")"
  )
  stop(
    "Random-effect block '", random_term$block_name,
    "' (", toupper(structure), "; ", n_groups, " groups x ",
    n_columns, " columns) exceeds the current dense JAGS compiler limits: ",
    paste(details, collapse = ", "), ". This representation can exhaust ",
    "memory before sampling. Reduce the structured index dimension or request ",
    "parameterization = \"noncentered\" (or \"auto\") so an eligible block ",
    "can use the group-local structured compiler.",
    call. = FALSE
  )
}

.bt_random_effect_dense_complexity_status <- function(structure, n_groups,
                                                       n_columns, n_rows,
                                                       monitor_policy,
                                                       centered = FALSE){

  multiplier <- getOption(
    "BayesTools.random_effects_complexity_multiplier",
    1
  )
  if(!is.numeric(multiplier) || length(multiplier) != 1L || is.na(multiplier) ||
     multiplier <= 0){
    stop(
      "Option 'BayesTools.random_effects_complexity_multiplier' must be a positive numeric scalar.",
      call. = FALSE
    )
  }
  if(is.infinite(multiplier)){
    return(list(
      estimates = .bt_random_effect_dense_complexity(
        structure = structure,
        n_groups = n_groups,
        n_columns = n_columns,
        n_rows = n_rows,
        monitor_policy = monitor_policy,
        centered = centered
      ),
      limits = stats::setNames(rep(Inf, 3L), c(
        "cholesky_products", "transform_products", "monitored_values"
      )),
      exceeded = stats::setNames(rep(FALSE, 3L), c(
        "cholesky_products", "transform_products", "monitored_values"
      ))
    ))
  }

  estimates <- .bt_random_effect_dense_complexity(
    structure = structure,
    n_groups = n_groups,
    n_columns = n_columns,
    n_rows = n_rows,
    monitor_policy = monitor_policy,
    centered = centered
  )
  limits <- multiplier * c(
    cholesky_products = 2e5,
    transform_products = 5e6,
    monitored_values = 2e4
  )
  exceeded <- estimates > limits

  list(estimates = estimates, limits = limits, exceeded = exceeded)
}

.JAGS_random_effect_formula <- function(formula, parameter, data,
                                        prior_random = NULL,
                                        sd_binding_context = NULL,
                                        group_data = data,
                                        compile_mode = c("sampled", "marginalized")){

  if(is.null(prior_random)){
    stop("Formula random effects require 'prior_random'.", call. = FALSE)
  }
  compile_mode <- match.arg(compile_mode)
  sampled_random_effect <- identical(compile_mode, "sampled")

  random_term <- .bt_as_random_effect_term(formula)
  .bt_validate_random_effect_term_supported(random_term)

  # extract the grouping factor information
  grouping_factor <- random_term$group_label
  grouping_metadata <- .bt_random_group_metadata(random_term, group_data)
  grouping_factor_levels <- grouping_metadata$levels
  grouping_mapping <- grouping_metadata$map

  formula <- random_term$term_formula
  random_structure <- .bt_random_term_structure(random_term, prior_random)
  structured_formula <- .bt_random_effect_normalize_structured_formula(
    formula = formula,
    data = data,
    structure = random_structure
  )
  formula <- structured_formula$formula
  data <- structured_formula$data
  random_term$term_formula <- formula
  random_term$structured_index <- structured_formula$index

  # obtain predictors characteristics factors (copy from formula)
  formula_terms    <- stats::terms(formula)
  has_intercept    <- attr(formula_terms, "intercept") == 1
  predictors       <- as.character(attr(formula_terms, "variables"))[-1]
  if(any(!predictors %in% colnames(data)))
    stop(paste0("The ", paste0("'", predictors[!predictors %in% colnames(data)], "'", collapse = ", ")," predictor variable is missing in the data set."))
  predictors_type  <- sapply(predictors, function(predictor){
    if(is.factor(data[[predictor]]) | is.character(data[[predictor]])){
      return("factor")
    }else{
      return("continuous")
    }
  })
  if(length(predictors) == 0L){
    predictors_type <- stats::setNames(character(), character())
  }
  model_terms      <- c(if(has_intercept) "intercept", attr(formula_terms, "term.labels"))
  model_terms_type <- sapply(model_terms, function(model_term){
    model_term <- strsplit(model_term, ":")[[1]]
    if(length(model_term) == 1 && model_term == "intercept"){
      return("continuous")
    }else if(any(predictors_type[model_term] == "factor")){
      return("factor")
    }else{
      return("continuous")
    }
  })

  homogeneous_sd <- .bt_random_effect_homogeneous_sd(random_term, random_structure)
  sd_binding <- .bt_random_sd_binding_for_block(
    sd_binding_context,
    random_term$block_name
  )

  block_prior <- .bt_random_prior_for_block(prior_random, random_term$block_name)
  if(random_structure %in% c("cs", "hcs", "ar1", "car", "har") &&
     length(block_prior$contrasts) > 0L){
    stop(
      "Random-effect block '", random_term$block_name,
      "' uses structure '", random_structure,
      "', whose level basis is defined by the covariance structure; ",
      "random_block(contrasts = ...) is not supported for this block.",
      call. = FALSE
    )
  }
  if(is.null(sd_binding) && !is.null(block_prior$sd_source)){
    sd_binding <- .bt_random_sd_binding(
      source = block_prior$sd_source,
      application = "block",
      factors = list(),
      true_allocation = FALSE,
      allocations = list()
    )
  }
  bound_sd <- !is.null(sd_binding)
  .bt_validate_random_block_for_structure(
    block_prior,
    structure = random_structure,
    block_name = random_term$block_name
  )
  row_indexed_external_sd <- .bt_random_sd_binding_has_row_external_source(sd_binding)
  if(bound_sd){
    if(!is.null(block_prior$terms)){
      stop(
        "Random-effect block '", random_term$block_name,
        "' cannot use term-specific SD overrides while it is controlled by an SD binding.",
        call. = FALSE
      )
    }
    prior_list <- list()
    original_prior_names <- character()
  }else{
    prior_list <- .bt_random_prior_terms_to_prior_list(block_prior, model_terms, homogeneous_sd)
    original_prior_names <- names(prior_list)
  }
  monitor_policy <- block_prior$monitor
  if(isTRUE(sampled_random_effect) &&
     isTRUE(row_indexed_external_sd) && isTRUE(monitor_policy$coefficients)){
    stop(
      "Group-level coefficient monitoring is not supported for random-effect block '",
      random_term$block_name,
      "' with a row-indexed external SD source.",
      call. = FALSE
    )
  }
  if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
    monitor_policy$latent <- TRUE
  }

  if(!bound_sd){
    if(homogeneous_sd){
      check_list(prior_list, "prior_list", check_names = "sd", allow_other = TRUE, all_objects = TRUE)
    }else{
      check_list(prior_list, "prior_list", check_names = model_terms, allow_other = TRUE, all_objects = TRUE)
    }
    .bt_random_effect_check_structured_sd_priors(
      prior_list = prior_list,
      model_terms = model_terms,
      random_structure = random_structure,
      homogeneous_sd = homogeneous_sd
    )
    prior_list <- .bt_random_effect_force_nonnegative_priors(prior_list)
  }
  data <- .bt_random_effect_apply_factor_prior_contrasts(
    data = data,
    predictors_type = predictors_type,
    model_terms = model_terms,
    model_terms_type = model_terms_type,
    prior_list = prior_list,
    contrast_overrides = block_prior$contrasts
  )

  # get the design matrix. For no-intercept random formulas, add an intercept
  # while constructing the matrix and drop it afterwards. This preserves
  # BayesTools factor-prior contrasts instead of forcing raw level indicators.
  random_design <- .bt_random_effect_design_matrix(
    formula,
    data,
    preserve_no_intercept_contrasts = !random_structure %in% c("cs", "hcs", "ar1", "car", "har"),
    structure = random_structure,
    block_name = random_term$block_name
  )
  model_frame <- random_design$model_frame
  model_matrix <- random_design$model_matrix
  car_metadata <- random_design$car
  raw_column_names <- colnames(model_matrix)
  random_factor_predictors <- names(predictors_type)[predictors_type == "factor"]
  random_xlevels <- lapply(random_factor_predictors, function(predictor){
    if(predictor %in% names(model_frame) && is.factor(model_frame[[predictor]])){
      levels(model_frame[[predictor]])
    }else{
      NULL
    }
  })
  names(random_xlevels) <- random_factor_predictors
  random_xlevels <- random_xlevels[!vapply(random_xlevels, is.null, logical(1))]
  structure_owned_basis <- random_structure %in% c(
    "cs", "hcs", "ar1", "car", "har"
  )
  random_contrast_matrices <- if(isTRUE(structure_owned_basis)){
    lapply(random_xlevels, function(level_names){
      out <- diag(length(level_names))
      dimnames(out) <- list(level_names, level_names)
      out
    })
  }else{
    .bt_concrete_factor_contrasts(
      model_frame,
      names(random_xlevels),
      context = paste0(
        "Random-effect design for block '",
        random_term$block_name,
        "'"
      )
    )
  }

  # check whether intercept is unique parameter
  if(sum(grepl("intercept", names(prior_list))) > 1)
    stop("only the intercept parameter can contain 'intercept' in its name.")
  # check whether any reserved term is in usage
  .bt_validate_random_effect_reserved_name(
    names(prior_list),
    context = "naming variables or prior distributions"
  )

  # replace interaction signs (due to JAGS incompatibility)
  colnames(model_matrix)  <- gsub(":", "__xXx__", colnames(model_matrix))
  column_names            <- colnames(model_matrix)
  names(prior_list)       <- gsub(":", "__xXx__", names(prior_list))
  names(model_terms_type) <- gsub(":", "__xXx__", names(model_terms_type))
  model_terms             <- gsub(":", "__xXx__", model_terms)

  # prepare syntax & data based on the formula
  parameter_suffix <- paste0("_xREx__", random_term$block_name) # priors should not be named with parameter name (done on exit from formula)
  parameter        <- paste0(parameter, "_", parameter_suffix) # variables should be named already here
  random_syntax    <- NULL
  JAGS_data        <- list()
  new_prior_list   <- list()
  random_scale_terms <- character()
  add_parameters <- character()
  jags_modules <- character()
  required_packages <- character()
  correlation_metadata <- NULL

  ### in essence, the following prepares constructors that:
  # 1) samples standardized random effects xRE_Zx[ids, predictors] from a multivariate normal distribution
  # 2) create a vector of by-parameter standard deviation of the random effects xRE_STDx[predictors]
  # 3) multiplies the standardized random effects by parameter-specific standard deviations to create xRE_COEFx[ids, predictors] matrix
  # 4) computes the per observation formula output based on indexing the by-id COEF and selecting the observation variables
  # 5) appends the per-observation output to the higher order formula (done in the formula call itself)

  n_id  <- length(grouping_factor_levels)
  n_par <- ncol(model_matrix)
  if(n_par < 1L){
    stop("Random-effect term '", random_term$block_name, "' does not generate any design columns.", call. = FALSE)
  }

  # Resolve the canonical SD leaves and bindings before choosing a sampled
  # parameterization. Centered eligibility depends on their zero/external
  # support, not only on the user-facing block prior.
  sd_spec <- .bt_random_effect_sd_spec(
    parameter = parameter,
    parameter_suffix = parameter_suffix,
    prior_list = prior_list,
    random_term = random_term,
    grouping_factor = grouping_factor,
    sd_binding = sd_binding,
    model_matrix = model_matrix,
    model_terms = model_terms,
    model_terms_type = model_terms_type,
    predictors_type = predictors_type,
    data = data,
    random_structure = random_structure,
    has_intercept = has_intercept,
    homogeneous_sd = homogeneous_sd
  )
  terms_indexes <- sd_spec$terms_indexes
  sd_parameter_names <- sd_spec$sd_parameter_names
  sd_leaves <- sd_spec$sd_leaves
  random_scale_terms <- c(random_scale_terms, sd_spec$random_scale_terms)
  add_parameters <- c(add_parameters, sd_spec$add_parameters)
  random_syntax <- c(random_syntax, sd_spec$syntax)
  new_prior_list <- c(new_prior_list, sd_spec$prior_list)
  sd_binding <- sd_spec$sd_binding
  row_indexed_external_sd <- .bt_random_sd_binding_has_row_external_source(sd_binding)

  parameterization <- .bt_random_effect_resolve_parameterization(
    block_prior = block_prior,
    prior_list = sd_spec$prior_list,
    sd_binding = sd_binding,
    row_indexed_external_sd = row_indexed_external_sd,
    model_matrix = model_matrix,
    group_map = grouping_mapping,
    n_groups = n_id,
    compile_mode = compile_mode,
    block_name = random_term$block_name
  )
  dense_status <- .bt_random_effect_dense_complexity_status(
    structure = random_structure,
    n_groups = n_id,
    n_columns = n_par,
    n_rows = nrow(model_matrix),
    monitor_policy = monitor_policy,
    centered = identical(parameterization$resolved, "centered")
  )
  structured_layout <- NULL
  if(isTRUE(sampled_random_effect) &&
     random_structure %in% c("cs", "hcs", "ar1", "car", "har")){
    structured_layout <- .bt_random_effect_structured_local_layout(
      model_matrix = model_matrix,
      group_map = grouping_mapping,
      structure = random_structure,
      parameter_stem = parameter,
      n_groups = n_id,
      exact_indicator = isTRUE(random_design$exact_indicator),
      column_coordinates = if(identical(random_structure, "car")){
        car_metadata$time_values
      }else{
        seq_len(n_par)
      }
    )
  }
  dense_latent_cells <- n_id * n_par
  sparse_advantage <- !is.null(structured_layout) &&
    dense_latent_cells - structured_layout$n_local >= 100L &&
    4 * structured_layout$n_local <= dense_latent_cells
  group_local <- !is.null(structured_layout) &&
    isTRUE(structured_layout$all_groups_observed) &&
    identical(parameterization$resolved, "noncentered") &&
    !isTRUE(monitor_policy$coefficients) &&
    (any(dense_status$exceeded) || isTRUE(sparse_advantage)) &&
    structured_layout$n_local < dense_latent_cells
  latent_layout <- if(isTRUE(group_local)) structured_layout else NULL
  if(isTRUE(sampled_random_effect) && isTRUE(group_local)){
    .bt_random_effect_check_group_local_complexity(
      random_term = random_term,
      layout = latent_layout,
      n_rows = nrow(model_matrix),
      monitor_policy = monitor_policy,
      row_indexed_external_sd = row_indexed_external_sd,
      column_allocation = isTRUE(row_indexed_external_sd) &&
        identical(sd_binding$application, "column")
    )
  }else if(isTRUE(sampled_random_effect)){
    .bt_random_effect_check_dense_complexity(
      random_term = random_term,
      structure = random_structure,
      n_groups = n_id,
      n_columns = n_par,
      n_rows = nrow(model_matrix),
      monitor_policy = monitor_policy,
      centered = identical(parameterization$resolved, "centered")
    )
  }
  group_covariance <- .bt_random_effect_prepare_known_group_covariance(
    random_term = random_term,
    group_levels = grouping_factor_levels,
    n_columns = n_par,
    model_matrix = model_matrix,
    random_structure = random_structure,
    compile_mode = compile_mode,
    row_indexed_external_sd = row_indexed_external_sd
  )
  # step 1:
  if(isTRUE(sampled_random_effect) &&
     identical(parameterization$resolved, "noncentered") &&
     !isTRUE(group_local) && !is.null(group_covariance)){
    group_mean_name <- paste0(parameter, "_xRE_GROUP_MUx")
    group_precision_name <- paste0(parameter, "_xRE_GROUP_PRECx")
    group_latent_name <- paste0(parameter, "_xRE_GROUP_Zx")
    random_syntax <- c(random_syntax, paste0(
      " ", group_latent_name, "[1:", n_id, "] ~ dmnorm(",
      group_mean_name, "[1:", n_id, "], ",
      group_precision_name, "[1:", n_id, ",1:", n_id, "])\n",
      " for(i in 1:", n_id, "){\n",
      "   ", paste0(parameter, "_xRE_Zx"), "[i,1] = ", group_latent_name, "[i]\n",
      " }\n"
    ))
    JAGS_data[[group_mean_name]] <- rep(0, n_id)
    JAGS_data[[group_precision_name]] <- group_covariance$precision
  }else if(isTRUE(sampled_random_effect) &&
           identical(parameterization$resolved, "noncentered") &&
           !isTRUE(group_local)){
    random_syntax <- c(random_syntax, paste0(
      " for(i in 1:",n_id,"){\n",
      "   for(j in 1:", n_par, "){\n",
      "     ", paste0(parameter, "_xRE_Zx"), "[i,j] ~ dnorm(0, 1)\n",
      "   }\n",
      " }\n"
    ))
  }

  # step 2 was resolved before step 1 so parameterization sees canonical scales.
  if(!is.null(group_covariance) && isTRUE(row_indexed_external_sd)){
    stop(
      "Known group covariance for random-effect block '",
      random_term$block_name,
      "' does not support row-indexed external SD sources.",
      call. = FALSE
    )
  }

  # step 3
  if(random_structure == "us" && n_par == 1L){
    block_prior <- .bt_random_prior_for_block(prior_random, random_term$block_name)
    if(!is.null(block_prior$covariance$cor)){
      stop(
        "Single-column random-effect structure 'us' has no correlation parameter; remove the 'cor' prior.",
        call. = FALSE
      )
    }
  }
  if(random_structure %in% c("diag", "id") ||
     (random_structure == "us" && n_par == 1L)){
    if(isTRUE(sampled_random_effect) &&
       identical(parameterization$resolved, "centered") &&
       !is.null(group_covariance)){
      group_precision_name <- paste0(parameter, "_xRE_GROUP_PRECx")
      centered_precision_name <- paste0(parameter, "_xRE_GROUP_CENTER_PRECx")
      JAGS_data[[group_precision_name]] <- group_covariance$precision
      random_syntax <- c(random_syntax, paste0(
        " for(i in 1:", n_id, "){\n",
        "   for(j in 1:", n_id, "){\n",
        "     ", centered_precision_name, "[i,j] <- ",
        group_precision_name, "[i,j] / pow(", parameter,
        "_xRE_STDx[1], 2)\n",
        "   }\n",
        " }\n",
        " ", parameter, "_xRE_COEFx[1:", n_id,
        ",1] ~ dmnorm(rep(0, ", n_id, "), ", centered_precision_name,
        "[1:", n_id, ",1:", n_id, "])\n",
        " for(i in 1:", n_id, "){\n",
        "   ", parameter, "_xRE_Zx[i,1] <- ", parameter,
        "_xRE_COEFx[i,1] / ", parameter, "_xRE_STDx[1]\n",
        " }\n"
      ))
    }else if(isTRUE(sampled_random_effect) &&
             identical(parameterization$resolved, "centered")){
      random_syntax <- c(random_syntax, .bt_JAGS_centered_independent_random(
        parameter = parameter,
        K = n_par,
        n_groups = n_id,
        sd_name = paste0(parameter, "_xRE_STDx")
      ))
    }else if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
      random_syntax <- c(random_syntax, paste0(
        " for(i in 1:",n_par,"){\n",
        "   ",paste0(parameter, "_xRE_UNIT_COEFx"),"[1:",n_id,",i] = ",paste0(parameter, "_xRE_Zx"),"[1:",n_id,",i]\n",
        " }\n"
      ))
    }else if(isTRUE(sampled_random_effect)){
      random_syntax <- c(random_syntax, paste0(
        " for(i in 1:",n_par,"){\n",
        "   ",paste0(parameter, "_xRE_COEFx"),"[1:",n_id,",i] = ",paste0(parameter, "_xRE_Zx"),"[1:",n_id,",i] * ",paste0(parameter, "_xRE_STDx"),"[i]\n",
        " }\n"
      ))
    }
  }else if(random_structure == "us"){
    block_prior <- .bt_random_prior_for_block(prior_random, random_term$block_name)
    if(n_par > 1L && is.null(block_prior$covariance$cor)){
      stop(
        "Random-effect block '", random_term$block_name,
        "' with structure 'us' requires an LKJ correlation prior. ",
        "Supply 'cor = prior_lkj(eta = ...)'.",
        call. = FALSE
      )
    }
    lkj_prior <- .bt_random_block_lkj_prior(block_prior)
    lkj_module <- JAGS_lkj_corr_cholesky(
      name = paste0(parameter, "_xRE_CORx"),
      K = n_par,
      eta = lkj_prior$eta,
      include_correlation = isTRUE(monitor_policy$correlation) && isTRUE(lkj_prior$include_correlation),
      include_primitives = isTRUE(monitor_policy$lkj_primitives) || isTRUE(lkj_prior$include_primitives)
    )
    random_syntax <- c(random_syntax, lkj_module$syntax)
    add_parameters <- c(add_parameters, lkj_module$monitor, lkj_module$primitive_names)
    jags_modules <- c(jags_modules, lkj_module$jags_module)
    required_packages <- c(required_packages, lkj_module$required_packages)
    correlation_metadata <- list(
      type = "lkj",
      eta = lkj_prior$eta,
      primitive_names = lkj_module$primitive_names,
      primitive_bounds = lkj_module$primitive_bounds,
      cholesky_name = lkj_module$cholesky_name,
      correlation_name = lkj_module$correlation_name
    )
    if(isTRUE(sampled_random_effect) &&
       identical(parameterization$resolved, "centered")){
      correlation_expression <- function(row, column){
        terms <- paste0(
          lkj_module$cholesky_name, "[", row, ",", seq_len(n_par), "] * ",
          lkj_module$cholesky_name, "[", column, ",", seq_len(n_par), "]"
        )
        paste(terms, collapse = " + ")
      }
      random_syntax <- c(random_syntax, .bt_JAGS_centered_correlated_random(
        parameter = parameter,
        K = n_par,
        n_groups = n_id,
        sd_name = paste0(parameter, "_xRE_STDx"),
        correlation_expression = correlation_expression,
        cholesky_name = lkj_module$cholesky_name
      ))
    }else if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
      random_syntax <- c(random_syntax, paste0(
        " for(g in 1:",n_id,"){\n",
        "   for(i in 1:",n_par,"){\n",
        "     ",paste0(parameter, "_xRE_UNIT_COEFx"),"[g,i] = inprod(", lkj_module$cholesky_name, "[i,1:", n_par, "], ", paste0(parameter, "_xRE_Zx"), "[g,1:", n_par, "])\n",
        "   }\n",
        " }\n"
      ))
    }else if(isTRUE(sampled_random_effect)){
      random_syntax <- c(random_syntax, paste0(
        " for(g in 1:",n_id,"){\n",
        "   for(i in 1:",n_par,"){\n",
        "     ",paste0(parameter, "_xRE_COEFx"),"[g,i] = ",paste0(parameter, "_xRE_STDx"),"[i] * inprod(", lkj_module$cholesky_name, "[i,1:", n_par, "], ", paste0(parameter, "_xRE_Zx"), "[g,1:", n_par, "])\n",
        "   }\n",
        " }\n"
      ))
    }
  }else if(random_structure %in% c("cs", "hcs", "ar1", "car", "har")){
    block_prior <- .bt_random_prior_for_block(prior_random, random_term$block_name)
    centered_structure <- identical(parameterization$resolved, "centered")
    corr_module <- if(isTRUE(centered_structure)){
      .bt_JAGS_structured_corr_cholesky(
        node_prefix = parameter,
        prior_prefix = parameter_suffix,
        K = n_par,
        structure = random_structure,
        block_prior = block_prior,
        include_correlation = TRUE,
        require_rho = TRUE,
        distance_matrix = if(identical(random_structure, "car")){
          abs(outer(car_metadata$time_values, car_metadata$time_values, "-"))
        }else{
          NULL
        }
      )
    }else{
      .bt_JAGS_structured_corr_direct(
        node_prefix = parameter,
        prior_prefix = parameter_suffix,
        K = n_par,
        structure = random_structure,
        block_prior = block_prior,
        include_correlation = FALSE,
        require_rho = TRUE,
        distance_matrix = NULL
      )
    }
    random_syntax <- c(random_syntax, corr_module$syntax)
    for(corr_prior_name in names(corr_module$prior_list)){
      corr_module$prior_list[[corr_prior_name]] <- .bt_random_effect_set_prior_metadata(
        corr_module$prior_list[[corr_prior_name]],
        block = random_term$block_name,
        grouping = grouping_factor,
        type = "correlation",
        structure = random_structure,
        name = .bt_random_effect_public_name(random_term)
      )
    }
    new_prior_list <- c(new_prior_list, corr_module$prior_list)
    correlation_monitors <- corr_module$monitor
    if(isTRUE(centered_structure)){
      correlation_monitors <- corr_module$rho_name
    }
    add_parameters <- c(add_parameters, correlation_monitors)
    correlation_metadata <- corr_module$bridge
    if(identical(random_structure, "car")){
      correlation_metadata$time_variable <- car_metadata$time_variable
      correlation_metadata$time_values <- car_metadata$time_values
    }
    if(isTRUE(sampled_random_effect) && isTRUE(group_local)){
      random_syntax <- c(random_syntax, .bt_JAGS_structured_local_transform(
        parameter = parameter,
        layout = latent_layout,
        rho_name = corr_module$rho_name,
        sd_name = paste0(parameter, "_xRE_STDx"),
        row_indexed_external_sd = row_indexed_external_sd
      ))
    }else if(isTRUE(sampled_random_effect) && isTRUE(centered_structure)){
      correlation_expression <- function(row, column){
        paste0(corr_module$correlation_name, "[", row, ",", column, "]")
      }
      random_syntax <- c(random_syntax, .bt_JAGS_centered_correlated_random(
        parameter = parameter,
        K = n_par,
        n_groups = n_id,
        sd_name = paste0(parameter, "_xRE_STDx"),
        correlation_expression = correlation_expression,
        cholesky_name = corr_module$cholesky_name
      ))
    }else if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
      random_syntax <- c(random_syntax, paste0(
        .bt_JAGS_structured_dense_transform(
          parameter = parameter,
          structure = random_structure,
          K = n_par,
          n_groups = n_id,
          rho_name = corr_module$rho_name,
          sd_name = paste0(parameter, "_xRE_STDx"),
          row_indexed_external_sd = TRUE,
          car_time_values = if(identical(random_structure, "car")){
            car_metadata$time_values
          }else{
            NULL
          }
        )
      ))
    }else if(isTRUE(sampled_random_effect)){
      random_syntax <- c(random_syntax, .bt_JAGS_structured_dense_transform(
        parameter = parameter,
        structure = random_structure,
        K = n_par,
        n_groups = n_id,
        rho_name = corr_module$rho_name,
        sd_name = paste0(parameter, "_xRE_STDx"),
        row_indexed_external_sd = FALSE,
        car_time_values = if(identical(random_structure, "car")){
          car_metadata$time_values
        }else{
          NULL
        }
      ))
    }
  }

  # step 4
  if(isTRUE(sampled_random_effect) && !is.null(structured_layout) &&
     !isTRUE(row_indexed_external_sd)){
    random_syntax <- c(random_syntax, paste0(
      " for(i in 1:", nrow(model_matrix), "){\n",
      "   ", parameter, "[i] = ", parameter, "_xRE_COEFx[",
      parameter, "_xRE_MAPx[i],", parameter, "_xRE_COLx[i]]\n",
      " }\n"
    ))
  }else if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
    if(isTRUE(sd_binding$true_allocation)){
      allocation_target <- .bt_random_effect_allocation_target_metadata(
        sd_binding$allocations[[1L]],
        context = paste0(
          "Random-effect SD binding metadata for block '",
          random_term$block_name,
          "'"
        )
      )
      if(length(sd_binding$factors_by_column) > 0L &&
         !identical(allocation_target, "sd_component")){
        stop(
          "Random-effect SD binding metadata for block '",
          random_term$block_name,
          "' with row-indexed column factor chains require 'allocation$target' to be 'sd_component'.",
          call. = FALSE
        )
      }
    }
    if(isTRUE(sd_binding$true_allocation) &&
       length(sd_binding$allocations) > 0L &&
       identical(sd_binding$allocations[[1L]]$target, "sd_component")){
      .bt_check_random_sd_component_binding(
        binding = sd_binding,
        n_columns = n_par,
        context = paste0(
          "Random-effect SD binding metadata for block '",
          random_term$block_name,
          "'"
        )
      )
    }
    source_expression <- .bt_random_sd_binding_shared_source_expression(
      sd_binding,
      row_index = "i"
    )
    unit_expression <- if(!is.null(structured_layout)){
      paste0(
        parameter, "_xRE_UNIT_COEFx[", parameter, "_xRE_MAPx[i],",
        parameter, "_xRE_COLx[i]]"
      )
    }else{
      NULL
    }
    if(identical(sd_binding$application, "column")){
      if(length(sd_binding$factors_by_column) == 0L){
        stop(
          "Random-effect SD binding metadata for block '",
          random_term$block_name,
          "' are missing canonical 'binding$factors_by_column'.",
          call. = FALSE
        )
      }
      if(length(sd_binding$factors_by_column) != n_par){
        stop(
          "Random-effect SD binding metadata for block '",
          random_term$block_name,
          "' do not match the number of random-effect columns.",
          call. = FALSE
        )
      }
      if(!is.null(structured_layout)){
        column_scale_name <- paste0(parameter, "_xRE_ROW_COL_SCALEx")
        for(column in seq_len(n_par)){
          factor_expression <- .bt_random_sd_binding_factors_expression(
            sd_binding$factors_by_column[[column]]
          )
          random_syntax <- c(random_syntax, paste0(
            column_scale_name, "[", column, "] <- ", factor_expression, "\n"
          ))
        }
        row_contribution <- paste0(
          source_expression, " * ", column_scale_name, "[",
          parameter, "_xRE_COLx[i]] * ", unit_expression
        )
      }else{
        column_terms <- vapply(seq_len(n_par), function(column){
          factor_expression <- .bt_random_sd_binding_factors_expression(
            sd_binding$factors_by_column[[column]]
          )
          paste(
            c(
              if(!identical(factor_expression, "1")) factor_expression,
              paste0(parameter, "_xRE_UNIT_COEFx[", parameter,
                     "_xRE_MAPx[i],", column, "]"),
              paste0(parameter, "_xRE_DATAx[i,", column, "]")
            ),
            collapse = " * "
          )
        }, character(1))
        row_contribution <- paste0(
          source_expression, " * (", paste(column_terms, collapse = " + "), ")"
        )
      }
    }else{
      factor_expression <- .bt_random_sd_binding_factors_expression(sd_binding$factors)
      unit_contribution <- if(!is.null(structured_layout)){
        unit_expression
      }else{
        paste0(
          "inprod(",
          parameter, "_xRE_UNIT_COEFx[", parameter,
          "_xRE_MAPx[i], 1:", n_par, "], ",
          parameter, "_xRE_DATAx[i,1:", n_par, "]",
          ")"
        )
      }
      row_contribution <- paste0(
        source_expression, " * ", factor_expression, " * ", unit_contribution
      )
    }
    random_syntax <- c(random_syntax, paste0(
      " for(i in 1:",nrow(model_matrix),"){\n",
      "   ",parameter,"[i] = ", row_contribution, "\n",
      " }\n"
    ))
  }else if(isTRUE(sampled_random_effect)){
    random_syntax <- c(random_syntax, paste0(
      " for(i in 1:",nrow(model_matrix),"){\n",
      "   ",parameter,"[i] = inprod(", paste0(parameter, "_xRE_COEFx[", paste0(parameter, "_xRE_MAPx[i]"),", 1:",n_par,"]"), ", ", paste0(parameter, "_xRE_DATAx[i,1:", n_par,"]"),")\n",
      " }\n"
    ))
  }

  # create the JAGS data list
  if(isTRUE(sampled_random_effect) && !is.null(structured_layout)){
    JAGS_data[[paste0(parameter, "_xRE_MAPx")]] <- grouping_mapping
    JAGS_data[[paste0(parameter, "_xRE_COLx")]] <- structured_layout$row_column
  }else if(isTRUE(sampled_random_effect)){
    JAGS_data[[paste0(parameter, "_xRE_DATAx")]] <- model_matrix
    JAGS_data[[paste0(parameter, "_xRE_MAPx")]]  <- grouping_mapping
  }

  if(isTRUE(sampled_random_effect) && isTRUE(monitor_policy$latent)){
    add_parameters <- c(add_parameters, paste0(parameter, "_xRE_Zx"))
  }
  if(isTRUE(sampled_random_effect) &&
     isTRUE(monitor_policy$coefficients) && !isTRUE(row_indexed_external_sd)){
    add_parameters <- c(add_parameters, paste0(parameter, "_xRE_COEFx"))
  }

  random_term$model_matrix    <- model_matrix
  random_term$raw_column_names <- raw_column_names
  random_term$column_names     <- column_names
  random_term$contrasts        <- attr(model_matrix, "contrasts")
  random_term$contrast_matrices <- random_contrast_matrices
  random_term$contrast_owner   <- if(isTRUE(structure_owned_basis)){
    "structure"
  }else{
    "random_block"
  }
  random_term$xlevels          <- random_xlevels
  random_term$assign           <- attr(model_matrix, "assign")
  random_term$model_terms      <- model_terms
  random_term$model_terms_type <- model_terms_type
  random_term$group_levels     <- grouping_factor_levels
  random_term$group_map        <- grouping_mapping
  random_term$group_components <- grouping_metadata$components
  random_term$group_component_levels <- grouping_metadata$component_levels
  random_term$group_tuples     <- grouping_metadata$tuples
  random_term$group_labels     <- grouping_metadata$labels
  random_term$group_tuple_keys <- grouping_metadata$tuple_keys
  random_term$group_tuple_index <- grouping_metadata$tuple_index
  random_term$n_groups         <- n_id
  random_term$n_columns        <- n_par
  random_term$prior_terms      <- original_prior_names
  random_term$parameter_stem   <- parameter
  random_term$sd_parameter_names <- sd_parameter_names
  random_term$sd_leaves          <- sd_leaves
  random_jags_data_names <- names(JAGS_data)
  if(is.null(random_jags_data_names)){
    random_jags_data_names <- character()
  }
  random_term["jags_data_names"] <- list(random_jags_data_names)
  random_term$structure        <- random_structure
  random_term$homogeneous_sd   <- homogeneous_sd
  random_term$interface        <- "prior_random"
  random_term$sd_binding       <- sd_binding
  random_term$correlation      <- correlation_metadata
  random_term$group_covariance <- group_covariance
  random_term$car              <- car_metadata
  random_term$new_levels       <- block_prior$new_levels
  random_term$monitor          <- monitor_policy
  random_term$compile_mode     <- compile_mode
  random_term$parameterization_requested <- parameterization$requested
  random_term$parameterization_resolved  <- parameterization$resolved
  random_term$parameterization_reason    <- parameterization$reason
  random_term$parameterization_policy    <- parameterization$policy
  random_term$latent_layout <- latent_layout
  attr(random_term, "random_block") <- random_term$block_name
  attr(random_term, "compile_mode") <- compile_mode

  return(list(
    random_syntax  = random_syntax,
    formula_term   = if(isTRUE(sampled_random_effect)) paste0(parameter,"[i]") else character(),
    data           = JAGS_data,
    prior_list     = new_prior_list,
    random_scale_terms = random_scale_terms,
    add_parameters = unique(add_parameters),
    jags_modules   = unique(jags_modules),
    required_packages = unique(required_packages),
    random_effect  = random_term,
    formula        = formula
  ))
}
