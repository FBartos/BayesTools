.bt_random_effects_interface <- function(random_effects, prior_random = NULL){

  if(length(random_effects) == 0L){
    return("none")
  }
  if(!is.null(prior_random)){
    return("prior_random")
  }

  stop("Formula random effects require 'prior_random'.", call. = FALSE)
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
  grouping_values <- .bt_random_group_values(random_term, group_data)
  grouping_factor_levels <- levels(as.factor(grouping_values))
  grouping_mapping       <- as.numeric(factor(grouping_values, levels = grouping_factor_levels))

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
    if(is.factor(data[,predictor]) | is.character(data[,predictor])){
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
    prior_list = prior_list
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
  # step 1:
  if(isTRUE(sampled_random_effect)){
    random_syntax <- c(random_syntax, .add_JAGS_matrix(name = paste0(parameter, "_xRE_PRECx"), diag(1, n_par)))
    random_syntax <- c(random_syntax, paste0(
      " for(i in 1:",n_id,"){\n",
      "   ",paste0(parameter, "_xRE_Zx"),"[i,1:", n_par ,"] ~ dmnorm(rep(0, ", n_par,"), ", paste0(parameter, "_xRE_PRECx"), ")\n",
      " }\n"
    ))
  }

  # step 2
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
    if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
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
    if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
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
    corr_module <- .bt_JAGS_structured_corr_cholesky(
      node_prefix = parameter,
      prior_prefix = parameter_suffix,
      K = n_par,
      structure = random_structure,
      block_prior = block_prior,
      include_correlation = isTRUE(monitor_policy$correlation),
      require_rho = TRUE,
      distance_matrix = if(identical(random_structure, "car")) car_metadata$distance_matrix else NULL
    )
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
    add_parameters <- c(add_parameters, corr_module$monitor)
    correlation_metadata <- corr_module$bridge
    if(identical(random_structure, "car")){
      correlation_metadata$time_variable <- car_metadata$time_variable
      correlation_metadata$time_values <- car_metadata$time_values
    }
    if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
      random_syntax <- c(random_syntax, paste0(
        " for(g in 1:",n_id,"){\n",
        "   for(i in 1:",n_par,"){\n",
        "     ",paste0(parameter, "_xRE_UNIT_COEFx"),"[g,i] = inprod(", corr_module$cholesky_name, "[i,1:", n_par, "], ", paste0(parameter, "_xRE_Zx"), "[g,1:", n_par, "])\n",
        "   }\n",
        " }\n"
      ))
    }else if(isTRUE(sampled_random_effect)){
      random_syntax <- c(random_syntax, paste0(
        " for(g in 1:",n_id,"){\n",
        "   for(i in 1:",n_par,"){\n",
        "     ",paste0(parameter, "_xRE_COEFx"),"[g,i] = ",paste0(parameter, "_xRE_STDx"),"[i] * inprod(", corr_module$cholesky_name, "[i,1:", n_par, "], ", paste0(parameter, "_xRE_Zx"), "[g,1:", n_par, "])\n",
        "   }\n",
        " }\n"
      ))
    }
  }

  # step 4
  if(isTRUE(sampled_random_effect) && isTRUE(row_indexed_external_sd)){
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
      column_terms <- vapply(seq_len(n_par), function(column){
        factor_expression <- .bt_random_sd_binding_factors_expression(
          sd_binding$factors_by_column[[column]]
        )
        paste(
          c(
            if(!identical(factor_expression, "1")) factor_expression,
            paste0(parameter, "_xRE_UNIT_COEFx[", parameter, "_xRE_MAPx[i],", column, "]"),
            paste0(parameter, "_xRE_DATAx[i,", column, "]")
          ),
          collapse = " * "
        )
      }, character(1))
      row_contribution <- paste0(source_expression, " * (", paste(column_terms, collapse = " + "), ")")
    }else{
      factor_expression <- .bt_random_sd_binding_factors_expression(sd_binding$factors)
      unit_contribution <- paste0(
        "inprod(",
        paste0(parameter, "_xRE_UNIT_COEFx[", paste0(parameter, "_xRE_MAPx[i]"), ", 1:", n_par, "]"),
        ", ",
        paste0(parameter, "_xRE_DATAx[i,1:", n_par, "]"),
        ")"
      )
      row_contribution <- paste0(source_expression, " * ", factor_expression, " * ", unit_contribution)
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
  if(isTRUE(sampled_random_effect)){
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
  random_term$xlevels          <- random_xlevels
  random_term$assign           <- attr(model_matrix, "assign")
  random_term$model_terms      <- model_terms
  random_term$model_terms_type <- model_terms_type
  random_term$group_levels     <- grouping_factor_levels
  random_term$group_map        <- grouping_mapping
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
  random_term$car              <- car_metadata
  random_term$new_levels       <- block_prior$new_levels
  random_term$compile_mode     <- compile_mode
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

