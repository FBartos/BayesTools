#' Random-effect marginal covariance update plan
#'
#' @description
#' Resolves a selected public random-effect quantity to the exact scalar
#' dependence of its marginal covariance contribution. The result is compiled
#' exclusively from the fitted parameter map and formula random-effect
#' metadata; posterior draws are neither inspected nor used to infer an update
#' form. An unsupported plan can include a `conditional` entry containing
#' `required_active` and `plan`. The nested plan applies only to posterior rows
#' where every named inclusion indicator equals one. It does not declare an
#' unconditional update for the selected public quantity.
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
#'   affected random-effect blocks. For a plan in a source coordinate
#'   (`coefficient_input = "source"`), the blocks include every block whose SD
#'   binding uses that source, also through another variance allocation.
#'
#'   Fields shared by the supported families are `update` (the kind of
#'   update, below), `blocks`, `structure` and `component_index` (the
#'   covariance structure and the column of a per-column quantity),
#'   `coefficient_input`, `source_parameter`, `source_transform` and
#'   `source_transform_spec` (the fitted coordinate and its transform to the
#'   public quantity), `quantity_id`, `canonical_name`, `formula_parameter` and
#'   `dependencies`. Affine plans add `coefficient_transform` and, for the
#'   Dirichlet-weight updates of an allocation, `allocation` (`weight_name`,
#'   `index`, `n_targets`, and the allocation `scale` when the plan comes from
#'   the conditional route). `invariant_covariance` holds the
#'   covariance that does not depend on the selected coordinate when it is
#'   structurally available, and `NULL` otherwise. Factor plans of
#'   allocation-derived SDs add `coefficient_scale_transform`. An unsupported
#'   plan carries `reason` and `message` instead.
#'
#'   `update` is `"scale"` when the contribution of the affected blocks is
#'   scaled by one coefficient, `"column_scale"` when the SDs of some columns of
#'   a block change (`component_index`; the other column SDs of a direct
#'   per-column SD are held fixed, and the driven columns of an
#'   allocation-derived SD are listed in `driven_columns`), `"allocation"` when
#'   a Dirichlet weight of a variance allocation changes, and `"correlation"`
#'   for a correlation coordinate. A direct per-column SD with unequal column
#'   SDs reports `"column_scale"`.
#'
#'   The plan of an allocation-derived component SD or variance, and of an
#'   aggregate allocation SD or variance (`sd_common`, `var_common`, `sd_total`,
#'   `var_total`) that is a quantity-scale plan, has the field
#'   `held_fixed = "allocation_weights"`: whatever the covariance structure,
#'   setting the selected SD to `v` holds the allocation weights and every
#'   non-allocation coordinate (correlations included) fixed and moves the
#'   allocation scale `tau` of the SD `SD_k = tau * m_k(w)` to `v / m_k(w)`,
#'   where `m_k(w)` is the product of the square roots of the weights on the
#'   path from `tau` to the SD (`sqrt(w_j)` for a total-variance allocation,
#'   `sqrt(K_j w_j)` for a mean-variance one, over the nested allocations).
#'   This is the conditional of the variance prior declared on `(tau, w)`,
#'   which are independent a priori. For an SD value `v` the prior of the moved
#'   coordinate is the `tau` prior evaluated at `v / m_k(w)` with the Jacobian
#'   `1 / m_k(w)`; for a variance value `u` the SD is `sqrt(u)` and the Jacobian
#'   is `1 / m_k(w)^2` on the variance parameterization of the scale
#'   (`tau^2`, equivalently `1 / (2 m_k(w) sqrt(u))` from `tau`). The scale
#'   reaches every block whose SD binding uses it, also through nested
#'   allocations (a split root with one `<block>_components` allocation per
#'   block) or a shared [random_sd_source()]: `blocks` lists all of them, the
#'   selected block first. Every column SD that derives from `tau` is
#'   multiplied by the same factor `c = v / SD_k`, so a reached block whose
#'   columns are all driven by `tau` has its covariance contribution multiplied
#'   by `c^2` on every covariance structure: the plan is then the affine
#'   `update = "scale"` plan of the summed contribution of `blocks`. If `tau`
#'   drives only some columns of the selected block, which is the only block it
#'   reaches, the plan is the factor `"column_scale"` plan with
#'   `driven_columns` and the grid factor `column_factor` for them, when the
#'   selected column is among the driven ones (otherwise the plan is
#'   `unsupported` with the reason `"selected_column_not_driven"`); a scale
#'   that reaches several blocks and only part of their columns, or an
#'   aggregate whose scale drives only part of a block, has no exact scalar
#'   update and is `unsupported` with the reason `"partial_scale_reach"`. An
#'   aggregate of a nested allocation moves the root scale of its allocation
#'   chain as well, so its `blocks` are those of every sibling allocation of
#'   the root, not only its own.
#'   Plans of direct per-column SDs have no `held_fixed` field.
#'
#' @details
#' `family = "affine"` means that, conditional on every non-selected posterior
#' coordinate, the row covariance has the exact form `A + h(value) * B`, where
#' `value` is the coordinate declared by `coefficient_input` (`"source"` or
#' `"quantity"`) and `h` is recorded in `coefficient_transform`.
#' `family = "factor"` identifies an exact candidate-dependent covariance
#' factor, while `family = "markov"` identifies a compiled structured Markov
#' state. Allocation-derived component SDs are affine on every covariance
#' structure when the allocation scale drives all columns of every block it
#' reaches (`A + h(v) B` with `h` the square of the SD or the identity of the
#' variance, `B` the summed covariance of `blocks` per unit of `h`); they are
#' factor `column_scale` updates only when it drives part of the columns of the
#' one block it reaches. All of them hold the allocation weights fixed and scale
#' every driven column SD by the same factor (`held_fixed`, above), so the
#' `column_scale` plan of an allocation-derived SD differs from that of a direct
#' SD, which replaces one column scale only. Their
#' `coefficient_scale_transform` converts a selected variance to its square
#' root and leaves a selected SD unchanged. These
#' quantity-scale plans require the public SD to be on the fitted coefficient
#' scale; after formula scaling that changes it, the plan is unsupported.
#' Unsupported plans carry a structural reason. The accessor
#' never estimates affineness from evaluated covariance matrices.
#'
#' The invariant covariance of a block with several coefficient columns is a
#' dense rows x rows matrix of the fitted data. Its estimated peak memory is
#' checked against the `BayesTools.random_effects_memory_limit_bytes` option
#' (see [BayesTools]) before it is formed; above the ceiling the accessor stops
#' with an error of class `BayesTools_random_effects_memory_limit` (also
#' `BayesTools_random_effects_error`).
#'
#' The `conditional` entry of an unsupported plan is compiled only for a
#' variance allocation with exactly two single-column block targets, a scalar
#' total-variance source and inclusion gates on the blocks; an allocation with
#' more than two targets has no conditional plan and stays unsupported.
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
  if(identical(plan$family, "unsupported") && plan$reason %in% c(
    "gated_realized_allocation", "ambiguous_allocation_scale",
    "non_scalar_covariance_path"
  )){
    block_allocation <- .bt_random_effect_marginal_update_block_allocation(
      design = design,
      key = key
    )
    if(!is.null(block_allocation)){
      if(length(block_allocation$required_active) > 0L){
        plan$conditional <- block_allocation
      }else{
        plan <- block_allocation$plan
      }
    }
  }
  if(identical(plan$family, "affine") &&
     identical(plan$coefficient_input, "source")){
    plan$blocks <- .bt_random_effect_marginal_update_source_blocks(
      fit = fit,
      design = design,
      key = key,
      blocks = plan$blocks,
      source = if(is.null(plan$source_parameter)){
        key$source_parameter
      }else{
        plan$source_parameter
      }
    )
  }
  plan$invariant_covariance <-
    .bt_random_effect_marginal_update_invariant_covariance(
      design = design,
      plan = plan,
      key = key
    )
  plan$quantity_id       <- quantity$quantity_id
  plan$canonical_name    <- quantity$canonical_name
  plan$formula_parameter <- parameter
  if(is.null(plan$source_parameter)){
    plan$source_parameter <- key$source_parameter
  }
  if(is.null(plan$source_transform)){
    plan$source_transform <- key$source_transform
  }
  # `source_transform` names the transform; `source_transform_spec` is the
  # transform itself, so a consumer that needs to apply it to a draw does not
  # have to re-derive it from the catalog. NULL when the quantity has none.
  if(is.null(plan$source_transform_spec)){
    plan$source_transform_spec <- .bt_parameter_transform_from_quantity(
      fit, quantity)
  }
  plan$dependencies      <- key$dependencies
  class(plan) <- c(
    "BayesTools_random_effects_marginal_update_plan",
    "list"
  )
  if(!is.null(plan$conditional)){
    for(field in c("quantity_id", "canonical_name", "formula_parameter",
                   "dependencies")){
      plan$conditional$plan[[field]] <- plan[[field]]
    }
    class(plan$conditional$plan) <- class(plan)
  }
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
#' @param values finite candidate values on the scale declared by `update`.
#'   Source-coordinate factor updates replace `source_parameter` in the
#'   posterior. Quantity-scale factor updates accept non-negative values of
#'   the selected SD or variance; variance values are converted to SDs using
#'   the plan's `coefficient_scale_transform`.
#' @param posterior_samples optional posterior sample object. By default the
#'   posterior is obtained from `fit`.
#' @param prior_list optional formula prior list used to reconstruct random SD
#'   bindings. By default it is obtained from `fit`.
#'
#' @return A list of class
#'   `BayesTools_random_effects_marginal_update_grid`. Factor updates contain
#'   posterior coefficient scales, correlation Cholesky factors, and the
#'   selected candidate scales. The candidate scale replaces the scale of
#'   column `component_index` only for a plan that holds the sibling SDs fixed
#'   (a direct per-column SD). For a plan with `held_fixed =
#'   "allocation_weights"` (an allocation-derived component SD) the grid also
#'   contains `column_factor` (values by posterior draws), the factor
#'   `candidate_scale / coefficient_scale[, component_index]` by which the
#'   scale of every driven column (`driven_columns`, all columns when it is
#'   `NULL`) is multiplied; it is `NA` for the draws where the current
#'   selected SD is zero, which no allocation scale maps to a non-zero (or an
#'   unchanged zero) candidate. Markov updates contain posterior coefficient
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
  quantity_scale <- identical(update$coefficient_input, "quantity")
  if(quantity_scale && any(values < 0)){
    stop("'values' must be non-negative for a random-effect SD or variance update.",
         call. = FALSE)
  }
  source <- update$source_parameter
  if(!quantity_scale &&
     (!is.character(source) || length(source) != 1L ||
      is.na(source) || !nzchar(source))){
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
  if(!quantity_scale && !source %in% colnames(posterior)){
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
    if(quantity_scale){
      transform <- update$coefficient_scale_transform
      candidate_values <- if(identical(transform, list(type = "sqrt"))){
        sqrt(values)
      }else if(identical(transform, list(type = "identity"))){
        values
      }else{
        stop("Quantity-scale factor updates require a supported coefficient-scale transform.",
             call. = FALSE)
      }
      candidate_scale[] <- rep(candidate_values, times = nrow(posterior))
    }else{
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
    if(identical(update$held_fixed, "allocation_weights")){
      # the allocation scale moves with the weights fixed: every column scale
      # is multiplied by the factor of the selected one
      column_factor <- sweep(candidate_scale, 2L, scale[, component], "/")
      column_factor[!is.finite(column_factor)] <- NA_real_
      out$held_fixed <- update$held_fixed
      out$driven_columns <- update$driven_columns
      out$column_factor <- unname(column_factor)
    }
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
      # the rows x rows group kernel, design cross-product and their product
      .bt_random_effect_check_memory(
        estimate = .bt_random_effect_update_plan_memory_estimate(
          n_rows = nrow(block_data$model_matrix),
          n_blocks = length(random_terms)
        ),
        block_name = random_term$block_name,
        alternative = paste0(
          "The plan holds the dense rows x rows covariance of each block of the ",
          "fitted data: fit the model to fewer rows. Raise the option (or set ",
          "it to Inf) only after verifying the operation's memory budget."
        )
      )
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
    direct_source <- .bt_random_effect_marginal_update_single_source(key)
    coefficient_input <- if(direct_source) "source" else "quantity"
    coefficient_transform <- if(direct_source ||
                                    identical(evaluator, "allocation_sd")){
      list(type = "square")
    }else{
      list(type = "identity")
    }
    # The aggregate is the allocation scale times weights: the update moves the
    # root scale of the allocation chain with every weight fixed, which reaches
    # every block whose SD binding uses it, also those of sibling allocations.
    reach <- .bt_random_effect_marginal_update_scale_reach(
      design = .bt_random_effect_marginal_covariance_design(
        fit = fit,
        parameter = key$formula_parameter
      ),
      source = .bt_random_sd_binding_source_name(allocation$source),
      first = .bt_random_effect_marginal_update_allocation_blocks(
        fit,
        key,
        allocation
      )
    )
    if(identical(reach$family, "unsupported")){
      return(.bt_random_effect_marginal_update_unavailable(
        quantity = quantity,
        reason = reach$reason,
        message = reach$message
      ))
    }
    return(.bt_random_effect_marginal_update_affine(
      update = "scale",
      blocks = reach$blocks,
      coefficient_transform = coefficient_transform,
      coefficient_input = coefficient_input,
      held_fixed = if(!direct_source) "allocation_weights"
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
    if(identical(
      .bt_random_effect_allocation_target_metadata(allocation),
      "sd_component"
    )){
      structure <- .bt_random_effect_structure(
        random_term,
        context = "Random-effect marginal covariance update plan"
      )
      if(random_term$n_columns > 1L && !structure %in% c("id", "diag")){
        return(.bt_random_effect_marginal_update_unavailable(
          quantity = quantity,
          reason = "correlated_component_allocation",
          message = paste0(
            "Correlated component allocations require candidate-dependent ",
            "covariance evaluation."
          )
        ))
      }
    }
    # The catalog indexes sd_mult entries after the n_targets var_prop or
    # var_mult entries; the plan reports the Dirichlet weight coordinate.
    index <- key$index
    if(identical(quantity$quantity, "sd_mult") &&
       index > allocation$n_targets){
      index <- index - allocation$n_targets
    }
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
        index = index,
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
        coefficient_transform = parameter_transform(fit, selection),
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
        structure = structure,
        coefficient_input = "source"
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
        coefficient_transform = parameter_transform(fit, selection),
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
    if(!.bt_parameter_catalog_random_sd_is_direct(
      random_term = random_term,
      parameter = key$formula_parameter,
      formula_scale = attr(fit, "formula_scale", exact = TRUE)
    )){
      # Quantity-scale plans (affine and factor) place candidate values on the
      # fitted coefficient scale, which differs from the public SD after
      # formula scaling.
      return(.bt_random_effect_marginal_update_unavailable(
        quantity = quantity,
        reason = "scaled_component_sd",
        message = paste0(
          "The selected allocation-derived SD is reported on the original ",
          "predictor scale, not the fitted coefficient scale, so it has no ",
          "exact quantity-scale covariance update."
        )
      ))
    }
    # The update moves the allocation scale the SD derives from, with the
    # weights fixed: every block whose SD binding uses that scale moves with it.
    block <- random_term$block_name
    reach <- .bt_random_effect_marginal_update_scale_reach(
      design = .bt_random_effect_marginal_covariance_design(
        fit = fit,
        parameter = key$formula_parameter
      ),
      source = .bt_random_sd_binding_source_name(random_term$sd_binding$source),
      first = block,
      column = key$index
    )
    if(identical(reach$family, "unsupported")){
      return(.bt_random_effect_marginal_update_unavailable(
        quantity = quantity,
        reason = reach$reason,
        message = reach$message
      ))
    }
    if(identical(reach$family, "affine")){
      # the scale drives every column of every reached block: each block's
      # D R D scales by c^2, whatever the covariance structure
      return(.bt_random_effect_marginal_update_affine(
        update = "scale",
        blocks = reach$blocks,
        coefficient_transform = if(identical(key$evaluator, "sd_variance")){
          list(type = "identity")
        }else{
          list(type = "square")
        },
        coefficient_input = "quantity",
        component_index = key$index,
        structure = structure,
        held_fixed = "allocation_weights"
      ))
    }
    # the scale drives only part of the columns of the one block it reaches,
    # among them the selected one
    return(list(
      family = "factor",
      update = "column_scale",
      blocks = block,
      component_index = key$index,
      structure = structure,
      coefficient_input = "quantity",
      coefficient_scale_transform = if(identical(key$evaluator, "sd_variance")){
        list(type = "sqrt")
      }else{
        list(type = "identity")
      },
      held_fixed = "allocation_weights",
      driven_columns = reach$driven_columns
    ))
  }

  if(!.bt_random_effect_marginal_update_single_source(key)){
    return(.bt_random_effect_marginal_update_unavailable(
      quantity = quantity,
      reason = "composite_sd_source",
      message = paste0(
        "The selected random-effect SD or variance depends on several fitted ",
        "coordinates (for example, after formula scaling), so it has no exact ",
        "covariance update in one source coordinate."
      )
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
      structure = structure,
      coefficient_input = "source"
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


# The blocks that an allocation scale reaches: those whose SD binding uses the
# source (the root scale of nested allocations, a shared random_sd_source()),
# with the columns of each block that it drives.
.bt_random_effect_marginal_update_reach <- function(design, source){

  blocks <- character()
  driven <- list()
  for(random_term in .bt_formula_design_random_effects(design)){
    binding <- random_term$sd_binding
    column_sources <- if(length(binding$sources_by_column) > 0L){
      binding$sources_by_column
    }else{
      rep(list(binding$source), random_term$n_columns)
    }
    columns <- vapply(seq_len(random_term$n_columns), function(column){
      candidate <- if(column <= length(column_sources)) column_sources[[column]]
      !is.null(candidate) &&
        identical(.bt_random_sd_binding_source_name(candidate), source)
    }, logical(1))
    if(any(columns)){
      blocks <- c(blocks, random_term$block_name)
      driven[[random_term$block_name]] <- columns
    }
  }

  list(blocks = blocks, driven = driven)
}


# What moving an allocation scale does to the blocks it reaches. The scale is
# moved with every weight fixed, so every column it drives is multiplied by the
# same factor: when it drives all columns of all reached blocks, each block's
# D R D scales by c^2 (affine). The scale drives only part of the columns of one
# reached block (factor, which needs the selected 'column' to be among the driven
# ones), or of several blocks (unsupported). 'first' are the blocks of the
# selected quantity, listed first.
.bt_random_effect_marginal_update_scale_reach <- function(design, source, first,
                                                          column = NA_integer_){

  unsupported <- function(reason, message){
    list(family = "unsupported", reason = reason, message = message)
  }
  reach <- .bt_random_effect_marginal_update_reach(design, source)
  if(!all(first %in% reach$blocks)){
    return(unsupported(
      "ambiguous_allocation_scale",
      paste0(
        "The selected allocation quantity does not identify the allocation ",
        "scale it derives from."
      )
    ))
  }
  blocks <- c(first, setdiff(reach$blocks, first))
  if(all(vapply(reach$driven, all, logical(1)))){
    return(list(family = "affine", blocks = blocks))
  }
  if(length(blocks) == 1L && !is.na(column)){
    driven <- which(reach$driven[[blocks]])
    if(column %in% driven){
      return(list(family = "factor", blocks = blocks, driven_columns = driven))
    }
    return(unsupported(
      "selected_column_not_driven",
      paste0(
        "The allocation scale of the selected allocation-derived SD does not ",
        "drive the selected column, so it has no exact covariance update."
      )
    ))
  }

  unsupported(
    "partial_scale_reach",
    paste0(
      "The allocation scale of the selected allocation quantity drives only ",
      "part of the columns of a reached random-effect block, in several blocks ",
      "or without a selected column, so it has no exact scalar covariance update."
    )
  )
}


.bt_random_effect_marginal_update_single_source <- function(key){

  key$source_type %in% c("identity", "one_to_one_transform") &&
    is.character(key$source_parameter) &&
    length(key$source_parameter) == 1L &&
    !is.na(key$source_parameter) && nzchar(key$source_parameter)
}


.bt_random_effect_marginal_update_affine <- function(
    update, blocks, coefficient_transform, coefficient_input,
    component_index = NA_integer_, structure = "", allocation = NULL,
    held_fixed = NULL){

  out <- list(
    family = "affine",
    update = update,
    blocks = unname(blocks),
    coefficient_transform = coefficient_transform,
    coefficient_input = coefficient_input,
    component_index = component_index,
    structure = structure,
    allocation = allocation
  )
  if(!is.null(held_fixed)){
    out$held_fixed <- held_fixed
  }

  out
}


.bt_random_effect_marginal_update_block_allocation <- function(design, key){

  component <- key$evaluator %in% c("sd", "sd_variance") &&
    isTRUE(key$allocation_derived)
  aggregate <- key$evaluator %in% c("allocation_sd", "allocation_var")
  proportion <- identical(key$evaluator, "allocation") &&
    identical(key$source_transform, "var_prop")
  if(!component && !aggregate && !proportion){
    return(NULL)
  }

  random_terms <- .bt_formula_design_random_effects(design)
  if(length(random_terms) != 2L || any(vapply(random_terms, function(term){
    term$n_columns != 1L || !identical(term$sd_binding$application, "block") ||
      length(term$sd_binding$allocations) != 1L ||
      length(term$sd_binding$sources_by_column) > 0L ||
      length(term$sd_binding$factors_by_column) > 0L
  }, logical(1)))){
    return(NULL)
  }
  block_names <- vapply(random_terms, `[[`, character(1), "block_name")
  allocation <- random_terms[[1L]]$sd_binding$allocations[[1L]]
  if(!identical(.bt_random_effect_allocation_target_metadata(allocation),
                "block") ||
     !identical(allocation$n_targets, 2L) ||
     !setequal(unname(allocation$terms), block_names) ||
     length(allocation$inclusion) > 0L ||
     !identical(allocation$scale, "total_variance") ||
     !is.list(allocation$source) ||
     !identical(allocation$source$shape, "scalar")){
    return(NULL)
  }
  source <- .bt_random_sd_binding_source_name(allocation$source)
  weight <- allocation$weight_name
  if(!is.character(source) || length(source) != 1L || is.na(source) ||
     !nzchar(source) || !is.character(weight) || length(weight) != 1L ||
     is.na(weight) || !nzchar(weight)){
    return(NULL)
  }
  if((aggregate || proportion) &&
     !identical(key$allocation_label, allocation$label)){
    return(NULL)
  }

  term_metadata <- lapply(random_terms, function(term){

    binding <- term$sd_binding
    current <- binding$allocations[[1L]]
    if(!identical(current$label, allocation$label) ||
       !identical(current$weight_name, weight) ||
       !identical(current$n_targets, 2L) ||
       !identical(current$scale, allocation$scale) ||
       length(current$inclusion) > 0L ||
       !identical(binding$source$shape, "scalar") ||
       !identical(.bt_random_sd_binding_source_name(binding$source), source)){
      return(NULL)
    }
    factors <- .bt_random_effect_allocation_factors_metadata(binding)
    split <- vapply(factors, function(factor){
      identical(factor$weight_name, weight)
    }, logical(1))
    if(sum(split) != 1L){
      return(NULL)
    }
    factor <- factors[[which(split)]]
    if(!identical(factor$scale, allocation$scale) ||
       !identical(factor$n_targets, 2L) ||
       !is.numeric(factor$index) || length(factor$index) != 1L ||
       is.na(factor$index) || !factor$index %in% 1:2 ||
       !is.null(factor$inclusion_name)){
      return(NULL)
    }
    parents <- factors[!split]
    pure_gates <- vapply(parents, function(parent){
      is.null(parent$weight_name) && identical(parent$n_targets, 1L) &&
        identical(parent$scale, "total_variance") &&
        is.character(parent$inclusion_name) &&
        length(parent$inclusion_name) == 1L &&
        !is.na(parent$inclusion_name) && nzchar(parent$inclusion_name)
    }, logical(1))
    if(!all(pure_gates)){
      return(NULL)
    }
    list(
      index = as.integer(factor$index),
      required_active = sort(unique(vapply(
        parents, `[[`, character(1), "inclusion_name"
      )))
    )
  })
  if(any(vapply(term_metadata, is.null, logical(1))) ||
     !identical(sort(vapply(term_metadata, `[[`, integer(1), "index")), 1:2) ||
     !identical(term_metadata[[1L]]$required_active,
                term_metadata[[2L]]$required_active)){
    return(NULL)
  }
  required_active <- term_metadata[[1L]]$required_active
  if(length(required_active) == 0L && !component){
    return(NULL)
  }

  index <- if(component){
    block <- match(key$random_block, block_names)
    if(is.na(block)){
      return(NULL)
    }
    term_metadata[[block]]$index
  }else if(proportion){
    key$index
  }else{
    NA_integer_
  }
  if(proportion && (!is.numeric(index) || length(index) != 1L ||
                    is.na(index) || !index %in% 1:2)){
    return(NULL)
  }
  plan <- .bt_random_effect_marginal_update_affine(
    update = if(proportion) "allocation" else "scale",
    blocks = block_names,
    coefficient_transform = if(proportion ||
                               identical(key$evaluator, "sd_variance")){
      list(type = "identity")
    }else{
      list(type = "square")
    },
    coefficient_input = if(component) "quantity" else "source",
    component_index = if(component) key$index else NA_integer_,
    allocation = list(
      weight_name = weight,
      index = as.integer(index),
      n_targets = 2L,
      scale = allocation$scale
    ),
    held_fixed = if(component) "allocation_weights"
  )
  plan$source_parameter <- if(proportion) weight else source
  plan$source_transform <- key$source_transform
  list(required_active = required_active, plan = plan)
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


# Changing a source coordinate changes every block whose SD binding uses it,
# also through variance allocations: the reach of the source.
.bt_random_effect_marginal_update_source_blocks <- function(
    fit, design, key, blocks, source){

  if(!is.character(source) || length(source) != 1L || is.na(source) ||
     !nzchar(source)){
    return(blocks)
  }

  unique(c(
    blocks,
    .bt_random_effect_marginal_update_reach(design, source)$blocks
  ))
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
