.bt_JAGS_bridge_marginal_random_spec <- function(
    formula_design_list,
    formula_random_effects_marginalize_list,
    bridge_context){

  if(is.null(formula_random_effects_marginalize_list)){
    return(list())
  }
  check_list(
    formula_random_effects_marginalize_list,
    "formula_random_effects_marginalize_list",
    allow_NULL = FALSE
  )
  spec_names <- names(formula_random_effects_marginalize_list)
  if(is.null(spec_names) || anyNA(spec_names) || any(!nzchar(spec_names)) ||
     anyDuplicated(spec_names)){
    stop(
      "'formula_random_effects_marginalize_list' must be a uniquely named list.",
      call. = FALSE
    )
  }
  unknown_parameters <- setdiff(spec_names, names(formula_design_list))
  if(length(unknown_parameters) > 0L){
    stop(
      "'formula_random_effects_marginalize_list' contains unknown formula parameter(s): ",
      paste(unknown_parameters, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  out <- list()
  for(parameter in spec_names){
    request <- formula_random_effects_marginalize_list[[parameter]]
    if(is.character(request) || is.null(request)){
      blocks <- request
      row_blocks <- NULL
      factor_state <- FALSE
    }else if(is.list(request)){
      request_names <- names(request)
      if(is.null(request_names) || anyNA(request_names) ||
         any(!nzchar(request_names)) || anyDuplicated(request_names) ||
         !all(request_names %in% c("blocks", "row_blocks", "factor_state")) ||
         !"blocks" %in% request_names){
        stop(
          "Bridge-marginalized random-effect request for formula parameter '",
          parameter,
          "' must contain uniquely named 'blocks' and optional 'row_blocks' and 'factor_state' entries.",
          call. = FALSE
        )
      }
      blocks <- request$blocks
      row_blocks <- request$row_blocks
      factor_state <- if(is.null(request$factor_state)){
        FALSE
      }else{
        check_bool(
          request$factor_state,
          paste0(
            "formula_random_effects_marginalize_list[['",
            parameter,
            "']]$factor_state"
          )
        )
        request$factor_state
      }
    }else{
      stop(
        "Bridge-marginalized random-effect request for formula parameter '",
        parameter,
        "' must be a character vector or a named list.",
        call. = FALSE
      )
    }
    check_char(
      blocks,
      paste0("formula_random_effects_marginalize_list[['", parameter, "']]$blocks"),
      check_length = 0,
      allow_NULL = TRUE,
      allow_NA = FALSE
    )
    if(is.null(blocks) || length(blocks) == 0L){
      next
    }
    if(anyDuplicated(blocks)){
      stop(
        "Bridge-marginalized random-effect block names must be unique for formula parameter '",
        parameter,
        "'.",
        call. = FALSE
      )
    }

    design <- formula_design_list[[parameter]]
    if(!.bt_formula_design_has_any_random_effects(design)){
      stop(
        "Formula parameter '", parameter,
        "' has no random-effect blocks to marginalize during bridge sampling.",
        call. = FALSE
      )
    }
    if(!identical(design$random_effects_interface, "prior_random")){
      stop(
        "Bridge-only Gaussian random-effect marginalization requires the prior_random() interface for formula parameter '",
        parameter,
        "'.",
        call. = FALSE
      )
    }

    random_effects <- .bt_formula_design_random_effects(design)
    block_names <- vapply(random_effects, `[[`, character(1), "block_name")
    unknown_blocks <- setdiff(blocks, block_names)
    if(length(unknown_blocks) > 0L){
      stop(
        "Bridge-only random-effect marginalization contains unknown block(s) for formula parameter '",
        parameter,
        "': ",
        paste(unknown_blocks, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    selected <- random_effects[match(blocks, block_names)]
    compile_modes <- vapply(
      selected,
      .bt_random_effect_term_compile_mode,
      character(1)
    )
    if(any(compile_modes != "sampled")){
      invalid <- blocks[compile_modes != "sampled"]
      stop(
        "Bridge-only random-effect marginalization can select only fitted sampled blocks. Already marginalized block(s): ",
        paste(invalid, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    row_blocks <- .bt_JAGS_bridge_marginal_random_row_blocks(
      row_blocks = row_blocks,
      n_rows = nrow(selected[[1L]]$model_matrix),
      parameter = parameter
    )
    if(!is.null(row_blocks)){
      .bt_JAGS_bridge_validate_marginal_random_row_blocks(
        random_effects = selected,
        row_blocks = row_blocks,
        parameter = parameter
      )
    }
    if(isTRUE(factor_state) && is.null(row_blocks)){
      stop(
        "Bridge-marginalized 'factor_state' requires exact 'row_blocks' for formula parameter '",
        parameter,
        "'.",
        call. = FALSE
      )
    }

    out[[parameter]] <- list(
      blocks = blocks,
      row_blocks = row_blocks,
      factor_state = factor_state
    )
  }

  if(length(out) > 0L && !bridge_context %in% c("full", "marginal")){
    stop(
      "Bridge-only random-effect marginalization requires bridge_context = TRUE, 'full', or 'marginal' so the likelihood callback receives the exact covariance contribution.",
      call. = FALSE
    )
  }

  out
}

.bt_JAGS_bridge_marginal_random_row_blocks <- function(row_blocks, n_rows,
                                                        parameter){

  if(is.null(row_blocks)){
    return(NULL)
  }
  if(!is.list(row_blocks) || length(row_blocks) == 0L){
    stop(
      "Bridge-marginalized 'row_blocks' for formula parameter '",
      parameter,
      "' must be a non-empty list of row indices.",
      call. = FALSE
    )
  }

  out <- lapply(seq_along(row_blocks), function(block_i){
    index <- row_blocks[[block_i]]
    if(!is.numeric(index) || length(index) == 0L || anyNA(index) ||
       any(!is.finite(index)) || any(index != as.integer(index)) ||
       any(index < 1L) || any(index > n_rows) || anyDuplicated(index)){
      stop(
        "Bridge-marginalized row block ", block_i,
        " for formula parameter '", parameter,
        "' must contain unique integer indices between 1 and ", n_rows, ".",
        call. = FALSE
      )
    }
    as.integer(index)
  })
  if(!identical(sort(as.integer(unlist(out))), seq_len(n_rows))){
    stop(
      "Bridge-marginalized 'row_blocks' for formula parameter '",
      parameter,
      "' must partition every formula row exactly once.",
      call. = FALSE
    )
  }

  out
}

.bt_JAGS_bridge_validate_marginal_random_row_blocks <- function(
    random_effects, row_blocks, parameter){

  n_rows <- sum(lengths(row_blocks))
  membership <- integer(n_rows)
  for(block_i in seq_along(row_blocks)){
    membership[row_blocks[[block_i]]] <- block_i
  }

  for(random_term in random_effects){
    group_map <- as.integer(random_term$group_map)
    if(length(group_map) != n_rows){
      stop(
        "Selected random-effect blocks for formula parameter '", parameter,
        "' do not have a common row count.",
        call. = FALSE
      )
    }
    if(.bt_random_effect_has_known_group_covariance(random_term)){
      group_covariance <- .bt_random_effect_known_group_covariance(
        random_term,
        context = "Bridge-only random-effect row blocks"
      )
      cross_block <- outer(membership, membership, "!=")
      structural_covariance <-
        group_covariance$kernel[group_map, group_map, drop = FALSE] != 0
      separated <- any(cross_block & structural_covariance)
    }else{
      group_membership <- split(membership, group_map)
      separated <- any(vapply(
        group_membership,
        function(x) length(unique(x)) > 1L,
        logical(1)
      ))
    }
    if(separated){
      stop(
        "Bridge-marginalized 'row_blocks' separate a structurally nonzero covariance contribution from random-effect block '",
        random_term$block_name,
        "' for formula parameter '", parameter, "'.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_JAGS_bridge_marginal_random_design_list <- function(
    formula_design_list,
    marginal_random_spec){

  if(length(marginal_random_spec) == 0L){
    return(formula_design_list)
  }

  out <- formula_design_list
  for(parameter in names(marginal_random_spec)){
    design <- out[[parameter]]
    selected <- marginal_random_spec[[parameter]]$blocks
    random_effects <- .bt_formula_design_random_effects(design)
    for(term_i in seq_along(random_effects)){
      if(random_effects[[term_i]]$block_name %in% selected){
        random_effects[[term_i]]$compile_mode <- "marginalized"
        attr(random_effects[[term_i]], "compile_mode") <- "marginalized"
      }
    }
    out[[parameter]] <- .bt_formula_design_set_random_effects(
      design = design,
      random_effects = random_effects
    )
  }

  out
}

.bt_JAGS_bridge_marginal_random_latent_names <- function(
    formula_design_list,
    marginal_random_spec){

  out <- character()
  for(parameter in names(marginal_random_spec)){
    design <- formula_design_list[[parameter]]
    selected <- marginal_random_spec[[parameter]]$blocks
    random_effects <- .bt_formula_design_random_effects(design)
    for(random_term in random_effects){
      if(!random_term$block_name %in% selected){
        next
      }
      out <- c(out, as.vector(.bt_random_effect_latent_names(
        random_term = random_term,
        n_groups = random_term$n_groups,
        n_columns = random_term$n_columns
      )))
    }
  }

  unique(out)
}

.bt_JAGS_bridge_compile_marginal_random_evaluator <- function(
    formula_design_list,
    marginal_random_spec,
    formula_data_list,
    formula_prior_list,
    model_data,
    posterior_names = NULL){

  if(length(marginal_random_spec) == 0L){
    return(list(
      active = FALSE,
      covariance = function(samples, prior_parameters,
                            formula_prior_parameters, formula_parameters,
                            factor_covariance = TRUE,
                            factor_state = FALSE) list()
    ))
  }

  plans <- lapply(names(marginal_random_spec), function(parameter){
    design <- formula_design_list[[parameter]]
    selected <- marginal_random_spec[[parameter]]$blocks
    row_blocks <- marginal_random_spec[[parameter]]$row_blocks
    factor_state <- isTRUE(marginal_random_spec[[parameter]]$factor_state)
    random_effects <- .bt_formula_design_random_effects(design)
    block_names <- vapply(random_effects, `[[`, character(1), "block_name")
    random_effects <- random_effects[match(selected, block_names)]
    block_plans <- lapply(random_effects, function(random_term){
      row_indexed <- .bt_random_effect_has_row_indexed_external_sd(random_term)
      structure <- .bt_random_effect_structure(
        random_term,
        context = "Bridge-only random-effect marginal covariance"
      )
      group_covariance <- if(
        .bt_random_effect_has_known_group_covariance(random_term)
      ){
        .bt_random_effect_known_group_covariance(
          random_term,
          context = "Bridge-only random-effect marginal covariance"
        )$kernel
      }else{
        NULL
      }
      block_data <- .bt_random_effect_marginal_covariance_block_data(
        design = design,
        random_term = random_term,
        data = NULL
      )
      prediction_rows <- if(row_indexed){
        .bt_random_effect_prediction_fitted_rows(
          random_term = random_term,
          n_rows = nrow(block_data$model_matrix),
          data_supplied = FALSE,
          context = "Bridge-only random-effect marginal covariance"
        )
      }else{
        NULL
      }
      factor_plan <- list(
        type = if(row_indexed){
          "row_group"
        }else if(!is.null(group_covariance)){
          "known_group"
        }else{
          "group"
        },
        model_matrix = block_data$model_matrix,
        group_map = block_data$group_map,
        coefficient_structure = if(
          is.null(group_covariance) &&
          structure %in% c("ar1", "car", "har") &&
          ncol(block_data$model_matrix) > 1L
        ){
          "markov"
        }else if(structure %in% c("diag", "id") ||
                 ncol(block_data$model_matrix) == 1L){
          "diagonal"
        }else{
          "dense"
        }
      )
      if(!is.null(group_covariance)){
        factor_plan$group_covariance <- group_covariance
      }
      list(
        random_term = random_term,
        row_indexed = row_indexed,
        structure = structure,
        group_covariance = group_covariance,
        coefficient_cholesky_evaluator =
          .bt_random_effect_compile_cholesky_evaluator(
            random_term = random_term,
            n_columns = ncol(block_data$model_matrix),
            structure = structure,
            posterior_names = posterior_names
          ),
        sd_evaluator = .bt_JAGS_bridge_compile_random_sd_evaluator(
          random_term = random_term,
          prior_list = formula_prior_list[[parameter]],
          posterior_names = posterior_names
        ),
        model_matrix = block_data$model_matrix,
        group_map = block_data$group_map,
        source_data = .bt_JAGS_marglik_parameter_source_data(
          model_data = model_data,
          formula_data = if(!is.null(formula_data_list)){
            formula_data_list[[parameter]]
          }else{
            NULL
          },
          design = design
        ),
        prediction_rows = prediction_rows,
        posterior_names = posterior_names,
        factor_plan = factor_plan
      )
    })
    names(block_plans) <- selected
    factor_state_evaluators <- lapply(
      block_plans,
      .bt_JAGS_bridge_compile_marginal_random_block_factor_state_evaluator,
      prior_list = formula_prior_list[[parameter]]
    )
    row_names <- rownames(block_plans[[1L]]$model_matrix)
    if(is.null(row_names)){
      row_names <- as.character(seq_len(nrow(block_plans[[1L]]$model_matrix)))
    }
    list(
      parameter = parameter,
      blocks = block_plans,
      factor_state_evaluators = factor_state_evaluators,
      factor_plans = lapply(block_plans, `[[`, "factor_plan"),
      block_names = selected,
      structures = stats::setNames(vapply(
        block_plans,
        `[[`,
        character(1),
        "structure"
      ), selected),
      row_blocks = row_blocks,
      factor_state = factor_state,
      row_names = row_names,
      prior_list = formula_prior_list[[parameter]],
      contract_id = new.env(parent = emptyenv())
    )
  })
  names(plans) <- names(marginal_random_spec)
  forbidden_formula_parameters <- names(formula_design_list)[vapply(
    formula_design_list,
    .bt_formula_design_has_any_random_effects,
    logical(1)
  )]
  force(plans)
  force(forbidden_formula_parameters)

  list(
    active = TRUE,
    factor_states = function(posterior){
      posterior <- as.matrix(posterior)
      if(nrow(posterior) < 1L ||
         (ncol(posterior) > 0L && is.null(colnames(posterior)))){
        stop(
          "Random-effect factor-state samples must be a non-empty matrix with column names when coordinates are present.",
          call. = FALSE
        )
      }
      out <- lapply(plans, function(plan){
        .bt_JAGS_bridge_marginal_random_factor_states_batch(
          plan = plan,
          posterior = posterior
        )
      })
      if(any(vapply(out, is.null, logical(1)))){
        return(NULL)
      }
      out
    },
    factor_components = function(posterior){
      posterior <- as.matrix(posterior)
      if(nrow(posterior) < 1L ||
         (ncol(posterior) > 0L && is.null(colnames(posterior)))){
        stop(
          "Random-effect factor-component samples must be a non-empty matrix with column names when coordinates are present.",
          call. = FALSE
        )
      }
      out <- lapply(plans, function(plan){
        .bt_JAGS_bridge_marginal_random_factor_components_batch(
          plan = plan,
          posterior = posterior
        )
      })
      if(any(vapply(out, is.null, logical(1)))){
        return(NULL)
      }
      out
    },
    coefficient_scales = function(posterior, parameter, block){
      posterior <- as.matrix(posterior)
      plan <- plans[[parameter]]
      if(is.null(plan) || is.null(plan$blocks[[block]])){
        stop("Random-effect coefficient-scale block is unavailable.",
             call. = FALSE)
      }
      .bt_JAGS_bridge_marginal_random_block_coefficient_scales_batch(
        block_plan = plan$blocks[[block]],
        posterior = posterior,
        prior_list = plan$prior_list
      )
    },
    coefficient_cholesky = function(posterior, parameter, block){
      posterior <- as.matrix(posterior)
      plan <- plans[[parameter]]
      if(is.null(plan) || is.null(plan$blocks[[block]])){
        stop("Random-effect coefficient-correlation block is unavailable.",
             call. = FALSE)
      }
      .bt_JAGS_bridge_marginal_random_block_coefficient_cholesky_batch(
        block_plan = plan$blocks[[block]],
        posterior = posterior
      )
    },
    covariance = function(samples, prior_parameters,
                          formula_prior_parameters, formula_parameters,
                          factor_covariance = TRUE,
                          factor_state = FALSE){
      posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
      if(isTRUE(factor_state) && !isTRUE(factor_covariance)){
        compact <- lapply(plans, function(plan){
          value <- .bt_JAGS_bridge_marginal_random_factor_state_posterior(
            plan = plan,
            posterior = posterior
          )
          if(is.null(value)){
            return(NULL)
          }
          .bt_JAGS_bridge_marginal_random_parameter_value(plan, value)
        })
        if(!any(vapply(compact, is.null, logical(1)))){
          class(compact) <- c("BayesTools_bridge_marginal_random", "list")
          return(compact)
        }
      }
      source_parameters <- .bt_JAGS_bridge_context_source_parameters(
        samples = samples,
        prior_parameters = prior_parameters,
        formula_prior_parameters = formula_prior_parameters,
        formula_parameters = formula_parameters
      )
      source_parameters <- .bt_parameter_source_forbid_formula_parameters(
        source_parameters,
        forbidden_formula_parameters
      )

      out <- lapply(plans, function(plan){
        covariance <- .bt_JAGS_bridge_marginal_random_covariance(
          plan = plan,
          posterior = posterior,
          source_parameters = source_parameters,
          factor_covariance = factor_covariance,
          factor_state = factor_state
        )
        .bt_JAGS_bridge_marginal_random_parameter_value(
          plan = plan,
          covariance = covariance
        )
      })
      class(out) <- c("BayesTools_bridge_marginal_random", "list")
      out
    }
  )
}

.bt_JAGS_bridge_marginal_random_parameter_value <- function(plan,
                                                             covariance){

  value <- c(
    covariance,
    list(
      blocks = plan$block_names,
      structures = plan$structures,
      row_names = plan$row_names,
      dimension = length(plan$row_names)
    )
  )
  class(value) <- c(
    "BayesTools_bridge_marginal_random_parameter",
    "list"
  )
  value
}

.bt_JAGS_bridge_marginal_random_factor_state_posterior <- function(plan,
                                                                   posterior){

  if(!isTRUE(plan$factor_state) || is.null(plan$row_blocks)){
    return(NULL)
  }
  states <- lapply(
    plan$factor_state_evaluators,
    function(evaluator) evaluator(posterior)
  )
  if(any(vapply(states, is.null, logical(1)))){
    return(NULL)
  }
  names(states) <- names(plan$blocks)
  list(
    representation = "factor_state",
    contract_id = plan$contract_id,
    row_blocks = plan$row_blocks,
    factor_plans = plan$factor_plans,
    factor_states = states
  )
}

.bt_JAGS_bridge_compile_marginal_random_block_factor_state_evaluator <-
    function(block_plan, prior_list){

  if(isTRUE(block_plan$row_indexed)){
    return(function(posterior) NULL)
  }
  random_term <- block_plan$random_term
  n_columns   <- ncol(block_plan$model_matrix)
  sd_evaluator <- block_plan$sd_evaluator
  direct_sd_evaluators <- if(is.null(sd_evaluator)){
    sd_names <- random_term$sd_parameter_names
    if(is.null(sd_names) || length(sd_names) != n_columns || anyNA(sd_names)){
      NULL
    }else{
      lapply(
        sd_names,
        .bt_JAGS_bridge_compile_parameter_draw_evaluator,
        prior_list = prior_list,
        posterior_names = block_plan$posterior_names
      )
    }
  }else{
    NULL
  }
  include_markov <- identical(
    block_plan$factor_plan$coefficient_structure,
    "markov"
  )
  force(random_term)
  force(n_columns)
  force(sd_evaluator)
  force(direct_sd_evaluators)
  force(include_markov)

  function(posterior){
    column_scale <- if(!is.null(sd_evaluator)){
      sd_evaluator$posterior_values(posterior)
    }else if(!is.null(direct_sd_evaluators)){
      values <- lapply(
        direct_sd_evaluators,
        function(evaluator) evaluator(posterior)
      )
      if(any(vapply(values, is.null, logical(1)))){
        NULL
      }else{
        vapply(values, `[[`, numeric(1), 1L)
      }
    }else{
      NULL
    }
    if(is.null(column_scale)){
      return(NULL)
    }
    .bt_random_effect_marginal_covariance_validate_draw_matrix(
      draws = matrix(column_scale, nrow = 1L),
      n_draws = 1L,
      n_columns = n_columns,
      label = "SD",
      random_term = random_term,
      nonnegative = TRUE
    )
    coefficient <- .bt_JAGS_bridge_marginal_random_coefficient_geometry(
      random_term = random_term,
      posterior = posterior,
      column_scale = column_scale,
      covariance = FALSE,
      structure = block_plan$structure,
      cholesky_evaluator = block_plan$coefficient_cholesky_evaluator
    )
    .bt_JAGS_bridge_marginal_random_factor_value(
      coefficient = coefficient,
      include_markov = include_markov
    )
  }
}

.bt_JAGS_bridge_marginal_random_factor_components_batch <- function(
    plan, posterior){

  attr(
    posterior,
    "BayesTools_random_effect_dirichlet_draw_cache"
  ) <- new.env(parent = emptyenv())
  out <- lapply(
    plan$blocks,
    .bt_JAGS_bridge_marginal_random_block_factor_components_batch,
    posterior = posterior,
    prior_list = plan$prior_list
  )
  if(any(vapply(out, is.null, logical(1)))){
    return(NULL)
  }
  names(out) <- names(plan$blocks)
  out
}

.bt_JAGS_bridge_marginal_random_factor_states_batch <- function(plan,
                                                                 posterior){

  if(!isTRUE(plan$factor_state) || is.null(plan$row_blocks)){
    return(NULL)
  }
  attr(
    posterior,
    "BayesTools_random_effect_dirichlet_draw_cache"
  ) <- new.env(parent = emptyenv())

  block_states <- lapply(
    plan$blocks,
    .bt_JAGS_bridge_marginal_random_block_factor_states_batch,
    posterior = posterior,
    prior_list = plan$prior_list
  )
  if(any(vapply(block_states, is.null, logical(1)))){
    return(NULL)
  }

  states <- lapply(seq_len(nrow(posterior)), function(draw){
    out <- lapply(block_states, `[[`, draw)
    names(out) <- names(block_states)
    out
  })
  list(
    representation = "factor_state",
    contract_id = plan$contract_id,
    row_blocks = plan$row_blocks,
    factor_plans = plan$factor_plans,
    factor_states = states
  )
}

.bt_JAGS_bridge_marginal_random_block_factor_states_batch <- function(
    block_plan, posterior, prior_list){

  components <- .bt_JAGS_bridge_marginal_random_block_factor_components_batch(
    block_plan = block_plan,
    posterior = posterior,
    prior_list = prior_list
  )
  if(is.null(components)){
    return(NULL)
  }

  n_columns <- ncol(block_plan$model_matrix)
  out <- vector("list", nrow(posterior))
  for(draw in seq_len(nrow(posterior))){
    cholesky <- if(is.null(components$coefficient_cholesky)){
      NULL
    }else{
      matrix(
        components$coefficient_cholesky[draw, , ],
        nrow = n_columns,
        ncol = n_columns
      )
    }
    factor <- if(is.null(cholesky)){
      diag(
        components$coefficient_scale[draw, ],
        nrow = n_columns,
        ncol = n_columns
      )
    }else{
      cholesky * components$coefficient_scale[draw, ]
    }
    value <- list(coefficient_factor = factor)
    if(isTRUE(components$markov)){
      value$coefficient_scale <- components$coefficient_scale[draw, ]
      value$markov_transition <-
        cholesky[cbind(2:n_columns, seq_len(n_columns - 1L))] /
        diag(cholesky)[seq_len(n_columns - 1L)]
      value$markov_innovation_variance <- diag(cholesky)[2:n_columns]^2
    }
    if(!is.null(components$row_scale)){
      value$row_scale <- as.numeric(components$row_scale[draw, ])
    }
    out[[draw]] <- value
  }
  names(out) <- NULL
  out
}


.bt_JAGS_bridge_marginal_random_block_factor_components_batch <- function(
    block_plan, posterior, prior_list){

  random_term <- block_plan$random_term
  n_columns   <- ncol(block_plan$model_matrix)
  row_scale_draws <- NULL
  if(isTRUE(block_plan$row_indexed)){
    source <- .bt_random_effect_row_indexed_source(random_term)
    if(.bt_parameter_source_has_values(source$source)){
      return(NULL)
    }
    source_draws <- .bt_random_effect_row_indexed_source_draws(
      random_term = random_term,
      n_rows = nrow(block_plan$model_matrix),
      posterior = posterior,
      data = block_plan$source_data,
      prediction_rows = block_plan$prediction_rows,
      context = "Bridge-only random-effect marginal covariance"
    )
    .bt_random_effect_marginal_covariance_validate_draw_matrix(
      draws = source_draws,
      n_draws = nrow(posterior),
      n_columns = nrow(block_plan$model_matrix),
      label = "row-indexed SD source",
      random_term = random_term,
      nonnegative = TRUE,
      context = "Bridge-only random-effect marginal covariance"
    )
    column_allocation <-
      .bt_random_effect_row_indexed_column_allocation_draws(
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
        n_draws = nrow(posterior),
        n_columns = 1L,
        label = "row-indexed SD allocation",
        random_term = random_term,
        nonnegative = TRUE,
        context = "Bridge-only random-effect marginal covariance"
      )
      row_scale_draws <- source_draws * allocation
      sd_draws <- matrix(1, nrow = nrow(posterior), ncol = n_columns)
    }else{
      .bt_random_effect_marginal_covariance_validate_draw_matrix(
        draws = column_allocation,
        n_draws = nrow(posterior),
        n_columns = n_columns,
        label = "row-indexed column SD allocation",
        random_term = random_term,
        nonnegative = TRUE,
        context = "Bridge-only random-effect marginal covariance"
      )
      row_scale_draws <- source_draws
      sd_draws <- column_allocation
    }
  }else{
    sd_draws <- .bt_JAGS_bridge_marginal_random_block_coefficient_scales_batch(
      block_plan = block_plan,
      posterior = posterior,
      prior_list = prior_list
    )
  }
  if(is.null(sd_draws)){
    return(NULL)
  }

  structure <- block_plan$structure
  cholesky_draws <-
    .bt_JAGS_bridge_marginal_random_block_coefficient_cholesky_batch(
      block_plan = block_plan,
      posterior = posterior
    )
  if(!structure %in% c("diag", "id") && n_columns > 1L &&
     is.null(cholesky_draws)){
    return(NULL)
  }

  list(
    coefficient_scale = unname(sd_draws),
    coefficient_cholesky = if(is.null(cholesky_draws)){
      NULL
    }else{
      unname(cholesky_draws)
    },
    row_scale = if(is.null(row_scale_draws)){
      NULL
    }else{
      unname(row_scale_draws)
    },
    markov = structure %in% c("ar1", "car", "har") && n_columns > 1L
  )
}


.bt_JAGS_bridge_marginal_random_block_coefficient_scales_batch <- function(
    block_plan, posterior, prior_list){

  random_term <- block_plan$random_term
  n_columns   <- ncol(block_plan$model_matrix)
  sd_draws <- if(is.null(block_plan$sd_evaluator)){
    .bt_random_effect_sd_draws(
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior,
      prior_list = prior_list
    )
  }else{
    block_plan$sd_evaluator$posterior_draws(posterior)
  }
  if(is.null(sd_draws) && !is.null(block_plan$sd_evaluator)){
    sd_draws <- .bt_random_effect_sd_draws(
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior,
      prior_list = prior_list
    )
  }
  if(is.null(sd_draws)){
    return(NULL)
  }
  .bt_random_effect_marginal_covariance_validate_draw_matrix(
    draws = sd_draws,
    n_draws = nrow(posterior),
    n_columns = n_columns,
    label = "SD",
    random_term = random_term,
    nonnegative = TRUE
  )
  unname(sd_draws)
}


.bt_JAGS_bridge_marginal_random_block_coefficient_cholesky_batch <- function(
    block_plan, posterior){

  random_term <- block_plan$random_term
  n_columns   <- ncol(block_plan$model_matrix)
  if(block_plan$structure %in% c("diag", "id") || n_columns == 1L){
    return(NULL)
  }
  cholesky <- if(is.null(block_plan$coefficient_cholesky_evaluator)){
    .bt_random_effect_cholesky_draws(
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior
    )
  }else{
    block_plan$coefficient_cholesky_evaluator(posterior)
  }
  if(is.null(cholesky)){
    return(NULL)
  }
  .bt_random_effect_marginal_covariance_validate_correlation_cholesky(
    cholesky = cholesky,
    random_term = random_term,
    n_columns = n_columns,
    posterior = posterior
  )
  unname(cholesky)
}

.bt_JAGS_bridge_marginal_random_covariance <- function(
    plan, posterior, source_parameters, factor_covariance = TRUE,
    factor_state = FALSE){

  factor_covariance <- isTRUE(factor_covariance) || is.null(plan$row_blocks)

  if(isTRUE(factor_state) && isTRUE(plan$factor_state) &&
     !is.null(plan$row_blocks)){
    states <- lapply(
      plan$blocks,
      .bt_JAGS_bridge_marginal_random_factor_state,
      posterior = posterior,
      prior_list = plan$prior_list,
      source_parameters = source_parameters,
      factor_covariance = FALSE
    )
    return(list(
      representation = "factor_state",
      contract_id = plan$contract_id,
      row_blocks = plan$row_blocks,
      factor_plans = plan$factor_plans,
      factor_states = states
    ))
  }

  geometries <- lapply(
    plan$blocks,
    .bt_JAGS_bridge_marginal_random_geometry,
    posterior = posterior,
    prior_list = plan$prior_list,
    source_parameters = source_parameters,
    factor_covariance = factor_covariance
  )
  if(is.null(plan$row_blocks)){
    covariance <- NULL
    for(geometry in geometries){
      contribution <- .bt_JAGS_bridge_marginal_random_geometry_covariance(
        geometry = geometry,
        index = seq_len(length(geometry$group_map))
      )
      if(is.null(covariance)){
        covariance <- contribution
      }else{
        covariance <- covariance + contribution
      }
    }
    dimnames(covariance) <- list(plan$row_names, plan$row_names)
    return(list(
      representation = "dense",
      covariance = covariance
    ))
  }

  list(
    representation = "factor",
    row_blocks = plan$row_blocks,
    factors = geometries
  )
}

.bt_JAGS_bridge_marginal_random_geometry <- function(
    block_plan, posterior, prior_list, source_parameters,
    factor_covariance = TRUE){

  state <- .bt_JAGS_bridge_marginal_random_factor_state(
    block_plan = block_plan,
    posterior = posterior,
    prior_list = prior_list,
    source_parameters = source_parameters,
    factor_covariance = factor_covariance
  )
  c(block_plan$factor_plan, state)
}

.bt_JAGS_bridge_marginal_random_factor_state <- function(
    block_plan, posterior, prior_list, source_parameters,
    factor_covariance = TRUE){

  random_term <- block_plan$random_term
  model_matrix <- block_plan$model_matrix
  if(isTRUE(block_plan$row_indexed)){
    source_draws <- .bt_random_effect_row_indexed_source_draws(
      random_term = random_term,
      n_rows = nrow(model_matrix),
      posterior = posterior,
      data = block_plan$source_data,
      parameters = source_parameters,
      prediction_rows = block_plan$prediction_rows,
      context = "Bridge-only random-effect marginal covariance"
    )
    .bt_random_effect_marginal_covariance_validate_draw_matrix(
      draws = source_draws,
      n_draws = 1L,
      n_columns = nrow(model_matrix),
      label = "row-indexed SD source",
      random_term = random_term,
      nonnegative = TRUE,
      context = "Bridge-only random-effect marginal covariance"
    )
    column_allocation <-
      .bt_random_effect_row_indexed_column_allocation_draws(
        random_term = random_term,
        posterior = posterior,
        prior_list = prior_list,
        n_columns = ncol(model_matrix)
      )
    if(is.null(column_allocation)){
      allocation <- .bt_random_effect_row_indexed_allocation_draws(
        random_term = random_term,
        posterior = posterior,
        prior_list = prior_list
      )
      .bt_random_effect_marginal_covariance_validate_draw_matrix(
        draws = matrix(allocation, ncol = 1L),
        n_draws = 1L,
        n_columns = 1L,
        label = "row-indexed SD allocation",
        random_term = random_term,
        nonnegative = TRUE,
        context = "Bridge-only random-effect marginal covariance"
      )
      row_scale <- as.numeric(source_draws[1L, ]) * allocation[1L]
      column_scale <- rep(1, ncol(model_matrix))
    }else{
      .bt_random_effect_marginal_covariance_validate_draw_matrix(
        draws = column_allocation,
        n_draws = 1L,
        n_columns = ncol(model_matrix),
        label = "row-indexed column SD allocation",
        random_term = random_term,
        nonnegative = TRUE,
        context = "Bridge-only random-effect marginal covariance"
      )
      row_scale <- as.numeric(source_draws[1L, ])
      column_scale <- as.numeric(column_allocation[1L, ])
    }
    coefficient <- .bt_JAGS_bridge_marginal_random_coefficient_geometry(
      random_term = random_term,
      posterior = posterior,
      column_scale = column_scale,
      covariance = factor_covariance,
      structure = block_plan$structure,
      cholesky_evaluator = block_plan$coefficient_cholesky_evaluator
    )
    value <- .bt_JAGS_bridge_marginal_random_factor_value(
      coefficient = coefficient,
      row_scale = row_scale,
      include_markov = identical(
        block_plan$factor_plan$coefficient_structure,
        "markov"
      )
    )
    if(factor_covariance){
      value$coefficient_covariance <- coefficient$covariance
    }
    return(value)
  }

  n_columns <- ncol(model_matrix)
  sd_draws <- if(is.null(block_plan$sd_evaluator)) {
    .bt_random_effect_sd_draws(
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior,
      prior_list = prior_list
    )
  } else {
    sd_values <- block_plan$sd_evaluator$posterior_values(
      posterior,
      parameters = source_parameters
    )
    if(is.null(sd_values)){
      .bt_random_effect_sd_draws(
        random_term = random_term,
        n_columns = n_columns,
        posterior = posterior,
        prior_list = prior_list
      )
    }else{
      matrix(sd_values, nrow = 1L)
    }
  }
  if(is.null(sd_draws)){
    .bt_random_effect_marginal_covariance_missing_sd_stop(
      random_term = random_term,
      n_columns = n_columns
    )
  }
  .bt_random_effect_marginal_covariance_validate_draw_matrix(
    draws = sd_draws,
    n_draws = 1L,
    n_columns = n_columns,
    label = "SD",
    random_term = random_term,
    nonnegative = TRUE
  )

  coefficient <- .bt_JAGS_bridge_marginal_random_coefficient_geometry(
    random_term = random_term,
    posterior = posterior,
    column_scale = as.numeric(sd_draws[1L, ]),
    covariance = factor_covariance,
    structure = block_plan$structure,
    cholesky_evaluator = block_plan$coefficient_cholesky_evaluator
  )

  value <- .bt_JAGS_bridge_marginal_random_factor_value(
    coefficient = coefficient,
    include_markov = identical(
      block_plan$factor_plan$coefficient_structure,
      "markov"
    )
  )
  if(factor_covariance){
    value$coefficient_covariance <- coefficient$covariance
  }
  value
}

.bt_JAGS_bridge_marginal_random_factor_value <- function(
    coefficient, row_scale = NULL, include_markov = FALSE){

  value <- list(coefficient_factor = coefficient$factor)
  if(isTRUE(include_markov)){
    value$coefficient_scale <- coefficient$scale
    value$markov_transition <- coefficient$markov_transition
    value$markov_innovation_variance <-
      coefficient$markov_innovation_variance
  }
  if(!is.null(row_scale)){
    value$row_scale <- row_scale
  }
  value
}

.bt_JAGS_bridge_marginal_random_coefficient_geometry <- function(
    random_term, posterior, column_scale, covariance = TRUE,
    structure = NULL, cholesky_evaluator = NULL){

  n_columns <- length(column_scale)
  if(is.null(structure)){
    structure <- .bt_random_effect_structure(
      random_term,
      context = "Bridge-only random-effect marginal covariance"
    )
  }
  if(structure %in% c("diag", "id") || n_columns == 1L){
    factor <- diag(
      column_scale,
      nrow = n_columns,
      ncol = n_columns
    )
  }else{
    cholesky_draws <- if(is.null(cholesky_evaluator)){
      .bt_random_effect_cholesky_draws(
        random_term = random_term,
        n_columns = n_columns,
        posterior = posterior
      )
    }else{
      cholesky_evaluator(posterior)
    }
    if(is.null(cholesky_draws)){
      .bt_random_effect_marginal_covariance_missing_correlation_stop(
        random_term = random_term,
        n_columns = n_columns,
        posterior = posterior
      )
    }
    .bt_random_effect_marginal_covariance_validate_correlation_cholesky(
      cholesky = cholesky_draws,
      random_term = random_term,
      n_columns = n_columns,
      posterior = posterior
    )
    cholesky <- matrix(
      cholesky_draws[1L, , ],
      nrow = n_columns,
      ncol = n_columns
    )
    factor <- sweep(
      matrix(cholesky, nrow = n_columns, ncol = n_columns),
      MARGIN = 1L,
      STATS = column_scale,
      FUN = "*"
    )
  }
  if(any(!is.finite(factor))){
    stop(
      "Bridge-only random-effect coefficient factor for block '",
      random_term$block_name,
      "' must be finite.",
      call. = FALSE
    )
  }

  list(
    covariance = if(isTRUE(covariance)) tcrossprod(factor) else NULL,
    factor = factor,
    scale = if(structure %in% c("ar1", "car", "har") &&
               n_columns > 1L) column_scale else NULL,
    markov_transition = if(
      structure %in% c("ar1", "car", "har") && n_columns > 1L
    ){
      cholesky[cbind(2:n_columns, seq_len(n_columns - 1L))] /
        diag(cholesky)[seq_len(n_columns - 1L)]
    }else{
      NULL
    },
    markov_innovation_variance = if(
      structure %in% c("ar1", "car", "har") && n_columns > 1L
    ){
      diag(cholesky)[2:n_columns]^2
    }else{
      NULL
    }
  )
}

.bt_JAGS_bridge_marginal_random_geometry_covariance <- function(geometry,
                                                                 index){

  if(identical(geometry$type, "dense")){
    return(geometry$covariance[index, index, drop = FALSE])
  }

  Z <- geometry$model_matrix[index, , drop = FALSE]
  group_map <- geometry$group_map[index]
  if(identical(geometry$type, "row_group")){
    Z <- Z * geometry$row_scale[index]
  }
  coefficient_covariance <- geometry$coefficient_covariance
  if(is.null(coefficient_covariance)){
    coefficient_covariance <- tcrossprod(geometry$coefficient_factor)
  }
  if(identical(geometry$type, "known_group")){
    return(
      geometry$group_covariance[group_map, group_map, drop = FALSE] *
        tcrossprod(Z %*% coefficient_covariance, Z)
    )
  }
  if(!geometry$type %in% c("group", "row_group")){
    stop("Unknown bridge marginal random-effect geometry.", call. = FALSE)
  }

  same_group <- outer(group_map, group_map, "==")
  storage.mode(same_group) <- "double"
  same_group * tcrossprod(Z %*% coefficient_covariance, Z)
}
