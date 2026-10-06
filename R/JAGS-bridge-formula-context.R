.bt_JAGS_bridge_formula_input_supplied <- function(formula_list,
                                                   formula_data_list,
                                                   formula_prior_list,
                                                   formula_scale_list,
                                                   formula_random_prior_list,
                                                   formula_random_effects_compile_list = NULL){

  !is.null(formula_list) ||
    !is.null(formula_data_list) ||
    !is.null(formula_prior_list) ||
    !is.null(formula_scale_list) ||
    !is.null(formula_random_prior_list) ||
    !is.null(formula_random_effects_compile_list)
}

.bt_JAGS_bridge_formula_context <- function(fit,
                                            formula_list,
                                            formula_data_list,
                                            formula_prior_list,
                                            formula_scale_list,
                                            formula_random_prior_list,
                                            formula_random_effects_compile_list = NULL){

  formula_input_supplied <- .bt_JAGS_bridge_formula_input_supplied(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_scale_list = formula_scale_list,
    formula_random_prior_list = formula_random_prior_list,
    formula_random_effects_compile_list = formula_random_effects_compile_list
  )
  fitted_formula_design <- .bt_JAGS_bridge_formula_design_list(
    attr(fit, "formula_design")
  )
  has_fitted_formula_design <- length(fitted_formula_design) > 0L

  if(!has_fitted_formula_design){
    if(formula_input_supplied){
      .bt_stop_refit_required(
        "JAGS_bridgesampling() cannot reconstruct formula parameters because ",
        "the fitted formula-design metadata are missing. Refit the model with ",
        "this version of BayesTools; supplied formula inputs cannot replace ",
        "fitted replay metadata."
      )
    }
    return(.bt_JAGS_bridge_empty_formula_context())
  }

  fitted_context <- .bt_JAGS_bridge_formula_context_from_design(
    fitted_formula_design
  )
  if(!formula_input_supplied){
    return(fitted_context)
  }

  rebuildability <- tryCatch(.bt_JAGS_bridge_formula_inputs_rebuildability(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_scale_list = formula_scale_list,
    formula_random_prior_list = formula_random_prior_list,
    formula_random_effects_compile_list = formula_random_effects_compile_list
  ), error = function(e){
    list(
      rebuildable = FALSE,
      detail = conditionMessage(e)
    )
  })
  if(!isTRUE(rebuildability$rebuildable)){
    .bt_JAGS_bridge_stop_formula_mismatch(rebuildability$detail)
  }

  rebuilt_context <- tryCatch(.bt_JAGS_bridge_rebuild_formula_context(
    fit = fit,
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_scale_list = formula_scale_list,
    formula_random_prior_list = formula_random_prior_list,
    formula_random_effects_compile_list = formula_random_effects_compile_list
  ), error = function(e){
    e
  })
  if(inherits(rebuilt_context, "error")){
    .bt_JAGS_bridge_stop_formula_mismatch(
      paste0("formula inputs could not be rebuilt: ", conditionMessage(rebuilt_context))
    )
  }
  mismatches <- .bt_JAGS_bridge_formula_design_mismatches(
    fitted_formula_design = fitted_formula_design,
    rebuilt_formula_design = rebuilt_context$formula_design_list
  )
  if(length(mismatches) > 0L){
    .bt_JAGS_bridge_stop_formula_mismatch(mismatches)
  }

  fitted_context
}

.bt_JAGS_bridge_formula_inputs_rebuildability <- function(formula_list,
                                                          formula_data_list,
                                                          formula_prior_list,
                                                          formula_scale_list,
                                                          formula_random_prior_list,
                                                          formula_random_effects_compile_list = NULL){

  if(is.null(formula_list) || is.null(formula_data_list) || is.null(formula_prior_list)){
    return(list(
      rebuildable = FALSE,
      detail = "formula-related inputs are incomplete"
    ))
  }

  check_list(formula_list, "formula_list", allow_NULL = FALSE)
  check_list(formula_data_list, "formula_data_list", allow_NULL = FALSE)
  check_list(formula_prior_list, "formula_prior_list", allow_NULL = FALSE)
  check_list(formula_scale_list, "formula_scale_list", allow_NULL = TRUE)
  check_list(formula_random_prior_list, "formula_random_prior_list", allow_NULL = TRUE)
  check_list(formula_random_effects_compile_list, "formula_random_effects_compile_list", allow_NULL = TRUE)
  .bt_validate_jags_formula_lists(
    formula_list = formula_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    formula_random_prior_list = formula_random_prior_list,
    formula_random_effects_compile_list = formula_random_effects_compile_list,
    formula_scale_list = formula_scale_list
  )
  if(!is.null(formula_random_prior_list)){
    for(parameter in names(formula_random_prior_list)){
      .bt_check_prior_random(formula_random_prior_list[[parameter]])
    }
  }
  if(!is.null(formula_random_effects_compile_list)){
    for(parameter in names(formula_random_effects_compile_list)){
      .bt_check_random_effects_compile(formula_random_effects_compile_list[[parameter]])
    }
  }

  for(parameter in names(formula_list)){
    if(.has_random_effects(formula_list[[parameter]]) &&
       (is.null(formula_random_prior_list) || is.null(formula_random_prior_list[[parameter]]))){
      detail <- paste0(
        "formula random effects for parameter '",
        parameter,
        "' cannot be rebuilt without a matching 'formula_random_prior_list' entry"
      )
      return(list(rebuildable = FALSE, detail = detail))
    }
  }

  list(rebuildable = TRUE, detail = NULL)
}

.bt_JAGS_bridge_rebuild_formula_context <- function(fit,
                                                    formula_list,
                                                    formula_data_list,
                                                    formula_prior_list,
                                                    formula_scale_list,
                                                    formula_random_prior_list,
                                                    formula_random_effects_compile_list = NULL){

  if(is.null(formula_scale_list)){
    formula_scale_list <- .JAGS_formula_scale_list_from_fit(
      fit,
      names(formula_list)
    )
  }

  formula_output <- list()
  for(parameter in names(formula_list)){
    formula_output[[parameter]] <- JAGS_formula(
      formula       = formula_list[[parameter]],
      parameter     = parameter,
      data          = formula_data_list[[parameter]],
      prior_list    = formula_prior_list[[parameter]],
      formula_scale = if(!is.null(formula_scale_list)) formula_scale_list[[parameter]] else NULL,
      prior_random  = if(!is.null(formula_random_prior_list)) formula_random_prior_list[[parameter]] else NULL,
      random_effects_compile = if(!is.null(formula_random_effects_compile_list)) formula_random_effects_compile_list[[parameter]] else NULL
    )
    priors <- attr(fit, "prior_list", exact = TRUE)
    rebuilt_priors <- formula_output[[parameter]]$prior_list
    priors <- c(rebuilt_priors, priors[setdiff(names(priors), names(rebuilt_priors))])
    fitted_design <- attr(fit, "formula_design", exact = TRUE)[[parameter]]
    fitted_spec <- attr(fitted_design$formula_scale, "unscale_design", exact = TRUE)
    completed_scale <- .bt_formula_scale_finalize(formula_output[[parameter]]$formula_scale,
      formula_output[[parameter]]$formula_design, priors,
      model_data = as.list(fitted_spec$state_constants), owner_scope = "fit")
    formula_output[[parameter]]$formula_scale <- completed_scale
    formula_output[[parameter]]$formula_design$formula_scale <- completed_scale
  }
  formula_data_output <- lapply(names(formula_output), function(parameter){
    .bt_JAGS_bridge_formula_data_with_raw(
      generated_data = formula_output[[parameter]][["data"]],
      raw_data = formula_data_list[[parameter]]
    )
  })
  names(formula_data_output) <- names(formula_output)

  list(
    formula_design_list = lapply(formula_output, function(output) output[["formula_design"]]),
    formula_list        = lapply(formula_output, function(output) output[["formula"]]),
    formula_data_list   = formula_data_output,
    formula_prior_list  = lapply(formula_output, function(output) output[["prior_list"]])
  )
}

.bt_JAGS_bridge_formula_data_with_raw <- function(generated_data, raw_data){

  out <- list()
  out <- .bt_JAGS_marglik_merge_source_data(out, raw_data)
  out <- .bt_JAGS_marglik_merge_source_data(out, generated_data)

  out
}

.bt_JAGS_bridge_empty_formula_context <- function(){

  list(
    formula_design_list = NULL,
    formula_list        = NULL,
    formula_data_list   = NULL,
    formula_prior_list  = list()
  )
}

.bt_JAGS_bridge_formula_context_from_design <- function(formula_design_list,
                                                       formula_data_list = NULL){

  for(design in formula_design_list){
    .bt_validate_formula_design_replay_schema(
      design,
      context = "JAGS_bridgesampling()"
    )
    .bt_formula_scale_finalized_check(design$formula_scale, design$parameter, require_owner = TRUE)
  }

  list(
    formula_design_list = formula_design_list,
    formula_list        = .bt_JAGS_bridge_formula_list_from_design(formula_design_list),
    formula_data_list   = formula_data_list,
    formula_prior_list  = .bt_JAGS_bridge_formula_prior_list_from_design(formula_design_list)
  )
}

.bt_JAGS_bridge_prior_list_from_fit <- function(fit, formula_design_list){

  fit_prior_list <- attr(fit, "prior_list")
  .bt_JAGS_bridge_non_formula_prior_list(
    prior_list = fit_prior_list,
    formula_design_list = formula_design_list,
    check_duplicates = FALSE
  )
}

# Formula priors come from the fitted formula design. A supplied 'prior_list'
# entry of a formula prior is dropped when it specifies the same prior as the
# fitted one (see .bt_JAGS_bridge_prior_semantics()) and rejected otherwise.
.bt_JAGS_bridge_non_formula_prior_list <- function(prior_list,
                                                   formula_design_list,
                                                   check_duplicates = FALSE){

  if(is.null(prior_list)){
    return(list())
  }
  if(!is.list(prior_list)){
    return(prior_list)
  }

  formula_prior_names <- .bt_JAGS_bridge_formula_prior_names(formula_design_list)
  if(length(formula_prior_names) == 0L){
    return(prior_list)
  }

  overlap <- intersect(names(prior_list), formula_prior_names)
  if(length(overlap) > 0L){
    if(isTRUE(check_duplicates)){
      fitted_priors <- do.call(c, unname(
        .bt_JAGS_bridge_formula_prior_list_from_design(formula_design_list)
      ))
      different <- overlap[!vapply(overlap, function(name){
        identical(
          .bt_JAGS_bridge_prior_semantics(prior_list[[name]]),
          .bt_JAGS_bridge_prior_semantics(fitted_priors[[name]])
        )
      }, logical(1))]
      if(length(different) > 0L){
        .bt_stop_refit_required(
          "JAGS_bridgesampling() received formula priors in 'prior_list' that ",
          "differ from the fitted formula priors: ",
          paste(utils::head(different, 8L), collapse = ", "),
          if(length(different) > 8L) ", ..." else "",
          ". Formula priors are taken from the fitted formula design; remove ",
          "them from 'prior_list' or refit the model with the intended priors."
        )
      }
    }
    prior_list <- prior_list[setdiff(names(prior_list), overlap)]
  }

  prior_list
}

# The part of a prior that defines its distribution: the class (prior type and
# contrast), the fields (distribution, parameters, truncation, and the fields of
# ordered priors), and the attributes 'multiply_by', 'components',
# 'prior_weights' (mixing weights), and 'inclusion_prior', with nested priors
# reduced the same way. Bookkeeping attributes that JAGS_formula() adds (formula
# parameter, levels and level names, factor design, interaction terms, ...) and
# the model prior weight are ignored.
.bt_JAGS_bridge_prior_semantics <- function(prior){

  if(!is.prior(prior)){
    return(prior)
  }

  fields <- unclass(prior)
  attributes(fields) <- list(names = names(fields))
  fields <- lapply(fields, .bt_JAGS_bridge_prior_semantics)
  fields[["prior_weights"]] <- NULL

  semantic_attributes <- c("multiply_by", "components", "prior_weights", "inclusion_prior")
  list(
    class      = class(prior),
    fields     = fields,
    attributes = lapply(
      stats::setNames(nm = semantic_attributes),
      function(name) .bt_JAGS_bridge_prior_semantics(attr(prior, name, exact = TRUE))
    )
  )
}

.bt_JAGS_bridge_formula_prior_names <- function(formula_design_list){

  formula_prior_list <- .bt_JAGS_bridge_formula_prior_list_from_design(
    formula_design_list
  )
  unique(unlist(lapply(formula_prior_list, names), use.names = FALSE))
}

.bt_JAGS_bridge_formula_list_from_design <- function(formula_design_list){

  if(length(formula_design_list) == 0L){
    return(NULL)
  }

  formula_list <- lapply(formula_design_list, function(design){
    if(inherits(design, "BayesTools_formula_design")){
      return(design$formula)
    }
    NULL
  })
  formula_list[!vapply(formula_list, is.null, logical(1))]
}

.bt_JAGS_bridge_formula_prior_list_from_design <- function(formula_design_list){

  if(length(formula_design_list) == 0L){
    return(list())
  }

  prior_list <- lapply(formula_design_list, function(design){
    if(inherits(design, "BayesTools_formula_design") &&
       is.list(design$prior_list)){
      return(design$prior_list)
    }
    list()
  })
  prior_list[!vapply(prior_list, function(x) length(x) == 0L, logical(1))]
}

.bt_JAGS_bridge_stop_formula_mismatch <- function(mismatches){

  if(length(mismatches) == 0L){
    mismatches <- "unknown formula-design mismatch"
  }

  stop(
    "JAGS_bridgesampling() supplied formula-related inputs do not fully match ",
    "the fitted formula design. First mismatch: ",
    mismatches[[1L]],
    if(length(mismatches) > 1L) paste0(" Additional mismatches: ", length(mismatches) - 1L, ".") else "",
    call. = FALSE
  )
}

.bt_JAGS_bridge_formula_design_mismatches <- function(fitted_formula_design,
                                                      rebuilt_formula_design){

  fitted_formula_design <- .bt_JAGS_bridge_formula_design_list(fitted_formula_design)
  rebuilt_formula_design <- .bt_JAGS_bridge_formula_design_list(rebuilt_formula_design)
  mismatches <- character()

  if(!setequal(names(fitted_formula_design), names(rebuilt_formula_design))){
    mismatches <- c(
      mismatches,
      paste0(
        "formula parameter names differ; fitted: ",
        paste(names(fitted_formula_design), collapse = ", "),
        "; supplied: ",
        paste(names(rebuilt_formula_design), collapse = ", ")
      )
    )
  }

  for(parameter in intersect(names(fitted_formula_design), names(rebuilt_formula_design))){
    fitted <- fitted_formula_design[[parameter]]
    rebuilt <- rebuilt_formula_design[[parameter]]
    mismatches <- c(
      mismatches,
      .bt_JAGS_bridge_formula_design_parameter_mismatches(
        parameter = parameter,
        fitted = fitted,
        rebuilt = rebuilt
      )
    )
  }

  mismatches
}

.bt_JAGS_bridge_formula_design_parameter_mismatches <- function(parameter,
                                                               fitted,
                                                               rebuilt){

  mismatches <- character()
  if(!inherits(fitted, "BayesTools_formula_design") ||
     !inherits(rebuilt, "BayesTools_formula_design")){
    return(paste0("formula design metadata for parameter '", parameter, "' are incomplete"))
  }

  if(!.bt_JAGS_bridge_formulas_equal(fitted$formula, rebuilt$formula)){
    mismatches <- c(mismatches, paste0("formula differs for parameter '", parameter, "'"))
  }
  if(!identical(isTRUE(fitted$log_intercept), isTRUE(rebuilt$log_intercept))){
    mismatches <- c(mismatches, paste0("log-intercept metadata differ for parameter '", parameter, "'"))
  }
  if(!identical(dim(fitted$model_matrix), dim(rebuilt$model_matrix)) ||
     !identical(colnames(fitted$model_matrix), colnames(rebuilt$model_matrix))){
    mismatches <- c(mismatches, paste0("fixed-effect model matrix shape or columns differ for parameter '", parameter, "'"))
  }else if(!identical(
    unname(fitted$model_matrix),
    unname(rebuilt$model_matrix)
  )){
    mismatches <- c(mismatches, paste0("fixed-effect model matrix values differ for parameter '", parameter, "'"))
  }
  if(!identical(fitted$assign, rebuilt$assign)){
    mismatches <- c(mismatches, paste0("fixed-effect model term assignments differ for parameter '", parameter, "'"))
  }
  if(!identical(fitted$contrasts, rebuilt$contrasts) ||
     !identical(
       fitted$contrast_matrices,
       rebuilt$contrast_matrices
     ) ||
     !identical(fitted$xlevels, rebuilt$xlevels)){
    mismatches <- c(mismatches, paste0("factor contrasts or levels differ for parameter '", parameter, "'"))
  }
  if(!identical(fitted$source_data, rebuilt$source_data)){
    mismatches <- c(mismatches, paste0(
      "original formula source data differ for parameter '", parameter, "'"
    ))
  }
  if(!identical(fitted$transformed_terms, rebuilt$transformed_terms)){
    mismatches <- c(mismatches, paste0(
      "literal expression() terms differ for parameter '", parameter, "'"
    ))
  }
  if(!identical(names(fitted$prior_list), names(rebuilt$prior_list))){
    mismatches <- c(mismatches, paste0("formula prior names differ for parameter '", parameter, "'"))
  }else{
    for(prior_name in names(fitted$prior_list)){
      if(!.bt_JAGS_bridge_formula_prior_metadata_equal(
        fitted$prior_list[[prior_name]],
        rebuilt$prior_list[[prior_name]]
      )){
        mismatches <- c(mismatches, paste0("formula prior metadata differ for parameter '", parameter, "', prior '", prior_name, "'"))
        break
      }
    }
  }

  random_mismatch <- tryCatch({
    .bt_JAGS_bridge_validate_formula_random_design(
      parameter = parameter,
      fitted = fitted,
      rebuilt = rebuilt
    )
    NULL
  }, error = function(e) conditionMessage(e))
  if(!is.null(random_mismatch)){
    mismatches <- c(mismatches, random_mismatch)
  }
  if(!identical(fitted$formula_scale, rebuilt$formula_scale)){
    mismatches <- c(mismatches, paste0("formula scaling metadata differ for parameter '", parameter, "'"))
  }

  mismatches
}

.bt_JAGS_bridge_formula_prior_metadata_equal <- function(x, y){

  semantic_attributes <- c(
    "class", "names", "multiply_by", "levels", "level_names",
    "interaction", "interaction_terms", "term_components", "factor_terms",
    "factor_contrasts", "factor_design", "factor_cell_names",
    "ordered_metadata", "random_factor", "random_grouping_factor",
    "random_allocation", "random_allocation_terms",
    "random_allocation_parent", "random_allocation_inclusion", "K",
    "components", "prior_weights", "model_prior_weights",
    "inclusion_prior", "component"
  )
  prior_core <- function(prior){
    prior_attributes <- attributes(prior)
    attributes(prior) <- prior_attributes[
      intersect(names(prior_attributes), semantic_attributes)
    ]
    prior
  }

  identical(prior_core(x), prior_core(y))
}

.bt_JAGS_bridge_formulas_equal <- function(x, y){

  isTRUE(all.equal(x, y, check.environment = FALSE))
}
