# Reconstruct ordered-factor term data from the design stored by JAGS_formula().
# stats::model.matrix() expands the first no-intercept factor to level indicators,
# even when a full-rank ordered contrast was supplied.
.marginal_posterior_term_data <- function(model_matrix, terms_indexes, term_index,
                                          data, prior_info, term_name){

  term_data <- model_matrix[, terms_indexes == term_index, drop = FALSE]
  if(!isTRUE(prior_info[["ordered"]])){
    return(term_data)
  }

  factor_terms <- prior_info[["factor_terms"]]
  factor_design <- prior_info[["factor_design"]]
  level_names <- prior_info[["level_names"]]

  if(is.null(factor_terms) || length(factor_terms) == 0L ||
     is.null(factor_design) || is.null(level_names) ||
     any(!factor_terms %in% names(data))){
    stop("Ordered factor metadata for '", term_name, "' are incomplete.", call. = FALSE)
  }

  if(!is.list(level_names)){
    if(length(factor_terms) != 1L){
      stop("Ordered factor metadata for '", term_name, "' are incomplete.", call. = FALSE)
    }
    level_names <- stats::setNames(list(level_names), factor_terms)
  }else{
    if(is.null(names(level_names))){
      if(length(level_names) != length(factor_terms)){
        stop("Ordered factor metadata for '", term_name, "' are incomplete.", call. = FALSE)
      }
      names(level_names) <- factor_terms
    }
    level_names <- level_names[factor_terms]
  }

  if(any(vapply(level_names, is.null, logical(1)))){
    stop("Ordered factor metadata for '", term_name, "' are incomplete.", call. = FALSE)
  }

  factor_design <- as.matrix(factor_design)
  cell_grid <- .factor_cell_grid(level_names)
  if(nrow(factor_design) != nrow(cell_grid) ||
     ncol(factor_design) != prior_info[["levels"]]){
    stop(
      "Ordered factor metadata for '", term_name,
      "' do not match its coefficient shape.",
      call. = FALSE
    )
  }

  term_data <- matrix(0, nrow = nrow(data), ncol = ncol(factor_design))
  for(row_i in seq_len(nrow(data))){
    factor_values <- data[row_i, factor_terms, drop = FALSE]
    if(anyNA(factor_values)){
      next
    }

    cell_match <- rep(TRUE, nrow(cell_grid))
    for(factor_term in factor_terms){
      cell_match <- cell_match &
        as.character(cell_grid[[factor_term]]) == as.character(factor_values[[factor_term]])
    }
    cell_index <- which(cell_match)
    if(length(cell_index) != 1L){
      stop(
        "Factor values for ordered term '", term_name,
        "' do not match its stored levels.",
        call. = FALSE
      )
    }
    term_data[row_i, ] <- factor_design[cell_index, ]
  }

  continuous_terms <- setdiff(prior_info[["term_components"]], factor_terms)
  for(continuous_term in continuous_terms){
    if(!continuous_term %in% names(data) || !is.numeric(data[[continuous_term]])){
      stop(
        "Continuous component '", continuous_term, "' of ordered term '",
        term_name, "' cannot be evaluated from 'at'.",
        call. = FALSE
      )
    }
    term_data <- term_data * matrix(
      data[[continuous_term]],
      nrow = nrow(term_data),
      ncol = ncol(term_data)
    )
  }

  term_data
}

#' @title Model-average marginal posterior distributions
#'
#' @description Creates marginal model-averages posterior distributions for a given
#' parameter based on model-averaged posterior samples and parameter name
#' (and formula with at specification).
#'
#' @param samples model-averaged posterior samples created by \code{mix_posteriors()}.
#' Mixed factor posteriors that lack the factor metadata of this version of
#' BayesTools (such as those created by BayesTools 0.3.0) must be recreated
#' from models fitted with this version: they stop with an error of class
#' \code{BayesTools_refit_required} (see [JAGS_validate_fit_contract()]).
#' @param parameter parameter of interest
#' @param formula model formula (needs to be specified if \code{parameter} was part of a formula)
#' @param at named list with predictor levels of the formula for which marginalization
#' should be performed. If a predictor level is missing, \code{0} is used for continuous
#' predictors, the baseline factor level is used for factors with \code{contrast = "treatment"} prior
#' distributions, and the parameter is completely omitted for factors with
#' \code{contrast = "meandif"}, \code{contrast = "orthonormal"},
#' \code{contrast = "independent"}, and ordered-factor levels. A predictor
#' without its own main-effect term (e.g., \code{x} in \code{~ g + g:x}) is
#' handled by the type, levels, and contrast recorded for it in the
#' interaction terms that contain it.
#' @param prior_samples whether marginal prior distributions should be generated
#' @param use_formula whether the parameter should be evaluated as a part of supplied formula
#' @param n_samples controls the numerical grid used for model-averaged
#' prior densities. For a simple parameter, an explicitly attached
#' \code{prior_density} is authoritative and is propagated instead of being
#' reconstructed from \code{prior_list}.
#' @param transformation,transformation_arguments optional transformation of
#' the marginal posterior, applied with [posterior_transform()] to the
#' untransformed marginal posterior and its metadata.
#' @inheritParams density.prior
#'
#' @details When the mixed posterior samples carry deterministic
#' \code{posterior_density}, \code{posterior_ordinate}, or \code{support}
#' metadata ([posterior_metadata()]), \code{marginal_posterior()} propagates
#' matching metadata to the returned marginal posterior. Matching uses the
#' parameter name, list/level names such as \code{theta[A]}, and conditional
#' metadata when present. Exact support metadata is propagated even when
#' \code{prior_samples = FALSE}; requesting prior samples adds prior-density
#' metadata but does not replace already attached posterior support.
#' A \code{transformation} is applied by [posterior_transform()] to the
#' marginal posterior with all its metadata (supports, atoms, prior densities,
#' stored posterior densities and ordinates, and component supports); the
#' mixed posterior samples themselves must be untransformed (samples
#' transformed with [posterior_transform()] whose prior is their
#' \code{prior_list} are refused). If support metadata is
#' absent, support is inferred from prior metadata only when the posterior
#' samples are on the raw, unconditioned prior scale; otherwise the deterministic
#' prior-density context is used so formula-scale transformations and
#' conditional model restrictions are respected. Marginal posteriors with
#' prior samples also record each draw's mixture component and each
#' component's exact support (\code{components} metadata): the model
#' of \code{mix_posteriors()} ensembles, or the combination of the mixture and
#' spike-and-slab component indicators of a single fit
#' (\code{as_mixed_posteriors()}). \code{Savage_Dickey_BF()} uses them when the
#' components' supports differ. With \code{log(intercept)} formula scaling, the
#' unscaled intercept (\code{transform_scaled} samples) is the exp of a linear
#' combination of the fitted coefficients with the log of the fitted
#' intercept: its simple marginal posterior (\code{use_formula = FALSE}) takes
#' its prior density and support from the exp of that combination, and the
#' linear predictors of a formula marginal posterior are linear in the same
#' coefficients.
#'
#' @return \code{marginal_posterior} returns a named list of mixed marginal posterior
#' distributions (either vectors or matrices).
#'
#' @export
marginal_posterior <- function(samples, parameter, formula = NULL, at = NULL, prior_samples = FALSE, use_formula = TRUE,
                               transformation = NULL, transformation_arguments = NULL, transformation_settings = FALSE,
                               n_samples = 10000, ...){

  check_list(samples, "samples")
  if(!inherits(samples, "mixed_posteriors"))
    stop("'samples' must be a be an object generated by 'mix_posteriors' function.")
  check_char(parameter, "parameter", allow_values = names(samples))
  if(is.numeric(samples[[parameter]]) && !inherits(samples[[parameter]], "mixed_posteriors")){
    .bt_draws_stop_plain(paste0(
      "The posterior samples of '", parameter,
      "' must be created by 'mix_posteriors' or 'as_mixed_posteriors', not plain numeric draws"
    ))
  }
  if(!is.null(formula) && !is.language(formula))
    stop("'formula' must be a formula")
  if(!is.null(at) && !is.list(at))
    stop("'at' must be a list")
  check_bool(prior_samples, "prior_samples")
  check_bool(use_formula, "use_formula")
  .check_transformation_input(transformation, transformation_arguments, transformation_settings)
  .marginal_posterior_check_untransformed(
    samples,
    if(use_formula && inherits(samples[[parameter]], "mixed_posteriors.formula")) names(samples) else parameter
  )


  # deal formula vs non-formula marginal posterior
  if(use_formula && inherits(samples[[parameter]], "mixed_posteriors.formula")){

      # remove the specified response (would crash the model.frame if not included)
      formula_log_intercept <- attr(formula, "log(intercept)", exact = TRUE)
      formula <- .remove_response(formula)
      formula_parameter <- .bt_meta_get(samples[[parameter]], "formula_parameter")
      log_intercept <- .marginal_posterior_log_intercept(
        samples               = samples,
        formula_log_intercept = formula_log_intercept,
        formula_parameter     = formula_parameter
      )

      ### extract the terms information from the formula
      formula_terms          <- stats::terms(formula)
      has_intercept          <- attr(formula_terms, "intercept") == 1
      predictors             <- as.character(attr(formula_terms, "variables"))[-1]
      model_terms            <- c(if(has_intercept) "intercept", attr(formula_terms, "term.labels"))

      JAGS_model_terms <- JAGS_parameter_names(parameters = model_terms, formula_parameter = formula_parameter)
      # an intercept-only formula has no predictors
      JAGS_predictors  <- if(length(predictors) > 0L){
        JAGS_parameter_names(parameters = predictors, formula_parameter = formula_parameter)
      }else{
        character()
      }


      ### obtain posterior samples and check that all are present
      if(!all(JAGS_model_terms %in% names(samples)))
        stop(paste0("The posterior samples for the ", paste0("'", JAGS_model_terms[!JAGS_model_terms %in% names(samples)], "'", collapse = ", ")," term is missing in the samples."))


      ### obtain prior list and information
      prior_list <- lapply(names(samples), function(model_term) attr(samples[[model_term]], "prior_list"))
      names(prior_list) <- names(samples)

      # get parameter information
      priors_info <- lapply(names(prior_list), function(model_term) list(
        term              = model_term,
        intercept         = model_term == "intercept",
        factor            = inherits(samples[[model_term]], "mixed_posteriors.factor"),
        levels            = attr(samples[[model_term]], "levels"),
        level_names       = attr(samples[[model_term]], "level_names"),
        interaction       = attr(samples[[model_term]], "interaction", exact = TRUE),
        interaction_terms = attr(samples[[model_term]], "interaction_terms"),
        term_components   = attr(samples[[model_term]], "term_components"),
        factor_terms      = attr(samples[[model_term]], "factor_terms"),
        factor_contrasts  = attr(samples[[model_term]], "factor_contrasts"),
        factor_design     = attr(samples[[model_term]], "factor_design"),
        treatment         = attr(samples[[model_term]], "treatment"),
        independent       = attr(samples[[model_term]], "independent"),
        orthonormal       = attr(samples[[model_term]], "orthonormal"),
        meandif           = attr(samples[[model_term]], "meandif"),
        ordered           = attr(samples[[model_term]], "ordered", exact = TRUE)
        ))
      names(priors_info) <- names(prior_list)
      model_terms_type <- sapply(JAGS_model_terms, function(model_term){
        if(priors_info[[model_term]][["factor"]]){
          return("factor")
        }else{
          return("continuous")
        }
      })
      # a predictor without its own main-effect term (e.g., `x` in
      # `~ g + g:x`) takes its type, fitted levels, and contrast from the
      # interaction terms that contain it
      predictors_info <- lapply(seq_along(predictors), function(i){
        if(JAGS_predictors[i] %in% JAGS_model_terms){
          return(priors_info[[JAGS_predictors[i]]])
        }
        .marginal_posterior_interaction_predictor_info(
          predictor   = predictors[i],
          parameter   = JAGS_predictors[i],
          priors_info = priors_info[JAGS_model_terms]
        )
      })
      names(predictors_info) <- JAGS_predictors
      predictors_type <- vapply(predictors_info, function(predictor_info){
        if(isTRUE(predictor_info[["factor"]])) "factor" else "continuous"
      }, character(1))


      ### prepare at specification
      # in case of an interaction, all levels need to be set
      # (the manipulated predictors are the term's components)
      if(!is.null(priors_info[[parameter]][["interaction"]]) && priors_info[[parameter]][["interaction"]]){
        manipulated_predictors <- priors_info[[parameter]][["interaction_terms"]]
        at_manipulated <- JAGS_parameter_names(manipulated_predictors, formula_parameter = formula_parameter)
      }else{
        at_manipulated <- parameter
        manipulated_predictors <- .bt_label_parts_coefficient(parameter, formula_parameter)$components
      }
      intercept_term <- identical(manipulated_predictors, "intercept")

      if(!all(names(at) %in% predictors))
        stop(paste0("The following values passed via the 'at' argument do not correspond to the specified model: ", paste0("'", names(at)[!names(at) %in% predictors], "'", collapse = ", ")))
      if(any(manipulated_predictors %in% names(at)))
        stop("Values of the parameter of interested cannot be specified via the 'at' argument.")

      # fill in with default values if needed
      for(i in seq_along(predictors)){
        if(JAGS_predictors[i] %in% at_manipulated){
          # specify levels for the parameter of interest
          if(predictors_type[[JAGS_predictors[i]]] == "continuous"){
            at[[predictors[i]]] <- c(-1, 0, 1)
          }else{
            at[[predictors[i]]] <- predictors_info[[JAGS_predictors[i]]][["level_names"]]
          }
        }else if(is.null(at[[predictors[i]]])){
          # specify levels for the remaining parameters
          if(predictors_type[[JAGS_predictors[i]]] == "continuous"){
            # fill in zeroes for unspecified continuous predictors
            at[[predictors[i]]] <- 0
          }else if(predictors_info[[JAGS_predictors[i]]][["treatment"]]){
            # fill in the default category for unspecified treatment factors
            at[[predictors[i]]] <- predictors_info[[JAGS_predictors[i]]][["level_names"]][1]
          }else{
            # fill in NA for any other factor type
            at[[predictors[i]]] <- NA
          }
        }
      }

      # transform to a data.frame (one row without predictors)
      data <- if(length(at) > 0L){
        as.data.frame(expand.grid(at))
      }else{
        data.frame(row.names = 1L)
      }

      # check the specified data
      if(any(predictors_type == "factor")){

        # check the proper data input for each factor variable
        for(i in seq_along(predictors_type)[predictors_type == "factor"]){

          if(is.factor(data[,predictors[i]])){
            if(all(levels(data[,predictors[i]]) %in% predictors_info[[JAGS_predictors[i]]][["level_names"]])){
              # either the formatting is correct, or the supplied levels are a subset of the original levels
              # reformat to check ordering and etc...
              data[,predictors[i]] <- factor(data[,predictors[i]], levels = predictors_info[[JAGS_predictors[i]]][["level_names"]])
            }else{
              # there are some additional levels
              stop(paste0("Levels specified in the '", predictors[i], "' factor variable do not match the levels used for model specification."))
            }
          }else if(all(stats::na.omit(unique(data[,predictors[i]])) %in% predictors_info[[JAGS_predictors[i]]][["level_names"]])){
            # the variable was not passed as a factor but the values matches the factor levels
            data[,predictors[i]] <- factor(data[,predictors[i]], levels = predictors_info[[JAGS_predictors[i]]][["level_names"]])
          }else{
            # there are some additional mismatching values
            stop(paste0("Levels specified in the '", predictors[i], "' factor variable do not match the levels used for model specification."))
          }

          # set the contrast
          factor_flag <- function(field){
            .marginal_posterior_factor_flag(predictors_info[[JAGS_predictors[i]]], field, JAGS_predictors[i])
          }
          if(factor_flag("orthonormal")){
            stats::contrasts(data[,predictors[i]]) <- "contr.orthonormal"
          }else if(factor_flag("meandif")){
            stats::contrasts(data[,predictors[i]]) <- "contr.meandif"
          }else if(factor_flag("independent")){
            stats::contrasts(data[,predictors[i]]) <- "contr.independent"
          }else if(factor_flag("ordered")){
            factor_contrasts <- unlist(
              predictors_info[[JAGS_predictors[i]]][["factor_contrasts"]],
              use.names = FALSE
            )
            ordered_contrasts <- factor_contrasts[
              .prior_ordered_is_contrast_name(factor_contrasts)
            ]
            if(length(ordered_contrasts) != 1L){
              stop(
                "Ordered contrast metadata for '", predictors[i],
                "' are incomplete.",
                call. = FALSE
              )
            }
            stats::contrasts(data[,predictors[i]]) <- ordered_contrasts
          }else if(factor_flag("treatment")){
            stats::contrasts(data[,predictors[i]]) <- "contr.treatment"
            if(anyNA(data[,predictors[i]]))
              stop("Unspecified levels in the '", predictors[i], "' factor (NAs not allowed for 'treatment' factors).")
          }
        }
      }
      if(any(predictors_type == "continuous")){

        # check the proper data input for each continuous variable
        for(i in seq_along(predictors_type)[predictors_type == "continuous"]){
          if(anyNA(data[,predictors[i]]))
            stop("Unspecified levels in the '", predictors[i], "' variable (NAs not allowed for continuous variables).")
          if(!is.numeric(data[,predictors[i]]))
            stop("Nonnumeric values in the '", predictors[i], "' continuous variable.")
        }
      }


      ### get the design matrix
      model_frame  <- stats::model.frame(formula, data = data, na.action = NULL)
      model_matrix <- .bt_model_matrix(model_frame, data = model_frame, formula = formula)

      # replaces NAs by zero to omit the corresponding coefficients
      model_matrix[is.na(model_matrix)] <- 0


      ### prepare posterior samples and information
      for(i in seq_along(samples)){
        # de-name factor levels
        if(priors_info[[i]][["factor"]]){
          if(priors_info[[i]][["levels"]] == 1){
            colnames(samples[[i]]) <- priors_info[[i]][["term"]]
          }else{
            colnames(samples[[i]]) <- paste0(priors_info[[i]][["term"]], "[", 1:priors_info[[i]][["levels"]], "]")
          }
        }
      }
      posterior_samples_matrix <- do.call(cbind, samples)


      # the model of each draw of a model-averaged ensemble (the terms' draws
      # must be aligned by model and draw); single fits have no model index
      model_component <- NULL
      if(!inherits(samples, "as_mixed_posteriors")){
        term_metadata <- function(field){
          do.call(cbind, lapply(model_terms, function(x){
            .bt_meta_get(samples[[JAGS_parameter_names(x, formula_parameter = formula_parameter)]], field)
          }))
        }
        model_component <- term_metadata("component")
        draw_index <- term_metadata("draw_index")
        if(is.null(model_component) || is.null(draw_index) ||
           !all(model_component[,1] == model_component) || !all(draw_index[,1] == draw_index))
          stop("the posterior samples are not alligned across models/draws")
        model_component <- model_component[,1]
      }


      ### evaluate the design matrix on the samples -> output[data, posterior]
      # The linear predictor is the fixed part of the registered
      # 'linear_predictor' node: the intercept (log(intercept) when declared,
      # never scaled) and the terms with their per-model multipliers.
      terms <- list()
      if(has_intercept){

        terms_indexes    <- attr(model_matrix, "assign") + 1
        terms_indexes[1] <- 0

        intercept_values <- posterior_samples_matrix[,JAGS_parameter_names("intercept", formula_parameter = formula_parameter)]
        if(log_intercept && any(!(intercept_values > 0))){
          stop(
            "The formula for '", formula_parameter, "' uses log(intercept), ",
            "but some intercept samples are not positive.",
            call. = FALSE
          )
        }
        terms[[1L]] <- .bt_dnode_linear_predictor_term(
          parameter  = formula_parameter,
          model_term = "intercept",
          type       = "intercept",
          columns    = 1L,
          prior      = NULL,
          log        = log_intercept
        )
        terms[[1L]]$coefficient_names <- JAGS_parameter_names("intercept", formula_parameter = formula_parameter)

      }else{

        terms_indexes <- attr(model_matrix, "assign")

      }

      # the remaining terms (omitting the intercept indexed as 0)
      for(i in unique(terms_indexes[terms_indexes > 0])){
        term <- .bt_dnode_linear_predictor_term(
          parameter  = formula_parameter,
          model_term = model_terms[i],
          type       = if(model_terms_type[i] == "factor") "factor" else "continuous",
          columns    = which(terms_indexes == i),
          prior      = NULL
        )
        term$term_index <- i
        term$coefficient_names <- paste0(
          JAGS_model_terms[i],
          if(model_terms_type[i] == "factor" && priors_info[[JAGS_model_terms[i]]][["levels"]] > 1) paste0("[", 1:priors_info[[JAGS_model_terms[i]]][["levels"]], "]")
        )
        terms[[length(terms) + 1L]] <- term
      }

      marginal_posterior_samples <- .bt_dnode_linear_predictor_fixed(
        terms        = terms,
        model_matrix = model_matrix,
        n_draws      = nrow(posterior_samples_matrix),
        values_of    = function(term){
          posterior_samples_matrix[, term$coefficient_names, drop = FALSE]
        },
        multiplier_of = function(term){
          # the multipliers of the term's prior in each model's draws
          .marginal_posterior_term_multiplier(
            JAGS_model_terms[term$term_index],
            prior_list  = prior_list,
            posterior   = posterior_samples_matrix,
            model_component = model_component,
            simple_list = inherits(samples, "as_mixed_posteriors")
          )
        },
        data_of = function(term){
          .marginal_posterior_term_data(
            model_matrix  = model_matrix,
            terms_indexes = terms_indexes,
            term_index    = term$term_index,
            data          = data,
            prior_info    = priors_info[[JAGS_model_terms[term$term_index]]],
            term_name     = JAGS_model_terms[term$term_index]
          )
        }
      )


      ### split the output into lists based on specification
      # create indexing and names for the manipulated predictors
      if(intercept_term){

        class(marginal_posterior_samples)             <- c(class(marginal_posterior_samples), "marginal_posterior.simple")
        attr(marginal_posterior_samples, "parameter") <- parameter
        attr(marginal_posterior_samples, "level")     <- "intercept"
        attr(marginal_posterior_samples, "data")      <- data
        marginal_posterior_samples <- .bt_meta_set(
          marginal_posterior_samples,
          "quantities",
          .marginal_posterior_level_quantities(
            formula_parameter = formula_parameter,
            components        = "intercept",
            level_frame       = NULL
          )
        )

        marginal_posterior_samples <- list("intercept" = marginal_posterior_samples)

        attr(marginal_posterior_samples, "data")        <- data
        attr(marginal_posterior_samples, "level_at")    <- NULL
        attr(marginal_posterior_samples, "level_names") <- "intercept"
        attr(marginal_posterior_samples, "parameter")   <- parameter

      }else{

        at_index_output        <- at[manipulated_predictors]
        at_index_output.names  <- at_index_output

        # rename continuous predictors levels
        for(i in seq_along(at_manipulated)){
          if(predictors_type[[at_manipulated[i]]] == "continuous"){
            at_index_output.names[[manipulated_predictors[i]]] <- paste0(at_index_output.names[[manipulated_predictors[i]]], "SD")
          }
        }

        at_index_output_frame       <- expand.grid(at_index_output)
        at_index_output.names_frame <- expand.grid(at_index_output.names)
        # the level names are the rendered level labels of the marginal means
        level_quantities <- .marginal_posterior_level_quantities(
          formula_parameter = formula_parameter,
          components        = manipulated_predictors,
          level_frame       = at_index_output.names_frame
        )
        level_names <- level_quantities$column

        # split the output samples
        data_split <- lapply(1:nrow(at_index_output_frame), function(i){
          apply(do.call(cbind, lapply(colnames(at_index_output_frame), function(pred){
            data[, pred] == at_index_output_frame[i, pred]
          })), 1, all)
        })
        marginal_posterior_samples <- lapply(seq_along(data_split), function(lvl){
          temp_marginal_posterior_samples <- marginal_posterior_samples[data_split[[lvl]],]
          temp_data                       <- data[data_split[[lvl]],]
          class(temp_marginal_posterior_samples)             <- c(class(temp_marginal_posterior_samples), "marginal_posterior.simple")
          attr(temp_marginal_posterior_samples, "parameter") <- parameter
          attr(temp_marginal_posterior_samples, "level")     <- level_names[lvl]
          attr(temp_marginal_posterior_samples, "data")      <- temp_data
          temp_marginal_posterior_samples <- .bt_meta_set(
            temp_marginal_posterior_samples,
            "quantities",
            level_quantities[lvl, , drop = FALSE]
          )
          return(temp_marginal_posterior_samples)
        })
        names(marginal_posterior_samples) <- level_names
        class(marginal_posterior_samples) <- c(class(marginal_posterior_samples), "marginal_posterior.factor")

        attr(marginal_posterior_samples, "data")        <- data
        attr(marginal_posterior_samples, "level_at")    <- at_index_output_frame
        attr(marginal_posterior_samples, "level_names") <- level_names
        attr(marginal_posterior_samples, "parameter")   <- parameter

      }


      prior_density_context <- .marginal_posterior_prior_density_context(
        samples          = samples,
        prior_list       = prior_list,
        column_names     = colnames(posterior_samples_matrix),
        n_samples        = n_samples,
        condition_source = samples[[parameter]]
      )

      if(!is.null(prior_density_context)){

      linear_weights <- matrix(
        0,
        nrow = nrow(data),
        ncol = length(prior_density_context$column_names),
        dimnames = list(NULL, prior_density_context$column_names)
      )

      if(has_intercept){
        terms_indexes    <- attr(model_matrix, "assign") + 1
        terms_indexes[1] <- 0
        intercept_name <- JAGS_parameter_names("intercept", formula_parameter = formula_parameter)
        if(intercept_name %in% colnames(linear_weights)){
          linear_weights[, intercept_name] <- 1
        }
      }else{
        terms_indexes <- attr(model_matrix, "assign")
      }

      for(i in unique(terms_indexes[terms_indexes > 0])){
        temp_data <- .marginal_posterior_term_data(
          model_matrix  = model_matrix,
          terms_indexes = terms_indexes,
          term_index    = i,
          data          = data,
          prior_info    = priors_info[[JAGS_model_terms[i]]],
          term_name     = JAGS_model_terms[i]
        )
        temp_all_columns <- paste0(
          JAGS_model_terms[i],
          if(model_terms_type[i] == "factor" && priors_info[[JAGS_model_terms[i]]][["levels"]] > 1) paste0("[", 1:priors_info[[JAGS_model_terms[i]]][["levels"]], "]")
        )
        temp_columns_keep <- temp_all_columns %in% colnames(linear_weights)
        temp_columns <- temp_all_columns[temp_columns_keep]
        temp_data <- temp_data[, temp_columns_keep, drop = FALSE]
        if(length(temp_columns) == 0)
          next

        linear_weights[, temp_columns] <- linear_weights[, temp_columns, drop = FALSE] + temp_data[, seq_along(temp_columns), drop = FALSE]
      }

      log_source_transforms <- NULL
      if(has_intercept && log_intercept){
        log_source_transforms <- stats::setNames("log", intercept_name)
      }
      # draws of model-averaged ensembles are aligned across terms by model
      ensemble_model_component <- model_component

      if(intercept_term){

        prior_weights <- linear_weights
        marginal_posterior_samples[["intercept"]] <- .marginal_posterior_formula_level_metadata(
          marginal                 = marginal_posterior_samples[["intercept"]],
          samples                  = samples,
          prior_list               = prior_list,
          prior_density_context    = prior_density_context,
          weights                  = prior_weights,
          column_name              = "intercept",
          source_transforms        = log_source_transforms,
          model_component          = ensemble_model_component
        )

      }else{

        level_prior_weights <- vector("list", length(level_names))
        names(level_prior_weights) <- level_names
        for(lvl in seq_along(level_names)){
          prior_weights <- linear_weights[data_split[[lvl]], , drop = FALSE]
          level_prior_weights[[level_names[lvl]]] <- prior_weights
          marginal_posterior_samples[[level_names[lvl]]] <- .marginal_posterior_formula_level_metadata(
            marginal                 = marginal_posterior_samples[[level_names[lvl]]],
            samples                  = samples,
            prior_list               = prior_list,
            prior_density_context    = prior_density_context,
            weights                  = prior_weights,
            column_name              = level_names[lvl],
            source_transforms        = log_source_transforms,
            model_component          = ensemble_model_component
          )
        }
      }

      # add priors
      if(prior_samples){

        if(sum(grepl(":", model_terms, fixed = TRUE)) > 5){
          warning(
            "Deterministic marginal prior densities with more than five interaction terms can be slow.",
            call. = FALSE
          )
        }

        if(intercept_term){

          prior_weights <- linear_weights
          prior_density <- .prior_density_from_context_rows(
            prior_density_context,
            prior_weights,
            source_transforms = log_source_transforms
          )
          marginal_posterior_samples[["intercept"]] <- .bt_meta_update(
            marginal_posterior_samples[["intercept"]],
            linear_weights = prior_weights,
            prior_density  = prior_density,
            prior_context  = prior_density_context
          )

        }else{

          for(lvl in seq_along(level_names)){
            prior_weights <- level_prior_weights[[level_names[lvl]]]
            prior_density <- .prior_density_from_context_rows(
              prior_density_context,
              prior_weights,
              source_transforms = log_source_transforms
            )
            marginal_posterior_samples[[level_names[lvl]]] <- .bt_meta_update(
              marginal_posterior_samples[[level_names[lvl]]],
              linear_weights = prior_weights,
              prior_density  = prior_density,
              prior_context  = prior_density_context
            )
          }
        }

        marginal_posterior_samples <- .bt_meta_set(marginal_posterior_samples, "prior_context", prior_density_context)
      }

      }

      marginal_posterior_samples <- .bt_meta_set(marginal_posterior_samples, "formula_parameter", formula_parameter)
      class(marginal_posterior_samples) <- c(class(marginal_posterior_samples), "marginal_posterior.formula")

  }else{

    if(!is.null(formula))
      stop("'formula' is supposed to be NULL when dealing with simple posteriors")
    if(!is.null(at))
      stop("'at' is supposed to be NULL when dealing with simple posteriors")


    ### obtain prior list and information
    prior_list <- lapply(names(samples), function(model_term) attr(samples[[model_term]], "prior_list"))
    names(prior_list) <- names(samples)
    prior_density_context <- NULL
    factor_weights <- NULL
    can_use_raw_prior_support <- .marginal_posterior_can_use_raw_prior_support(
      samples,
      parameter
    )


    ### extract the corresponding samples
    if(inherits(samples[[parameter]], "mixed_posteriors.factor")){

      parameter_samples <- samples[[parameter]]
      posterior_density_sources <- .posterior_density_sources(samples, parameter_samples)
      posterior_ordinate_sources <- .posterior_ordinate_sources(samples, parameter_samples)
      posterior_density_conditional <- .bt_meta_condition(parameter_samples, "conditional")
      posterior_density_conditional_rule <- .bt_meta_condition(parameter_samples, "conditional_rule")
      posterior_density_condition_key <- .bt_meta_condition(parameter_samples, "condition_key")

      # transform factor levels
      marginal_posterior_samples <- transform_factor_samples(samples)
      marginal_posterior_samples <- transform_treatment_samples(marginal_posterior_samples)[[parameter]]
      marginal_factor_atoms <- .posterior_atoms_get(marginal_posterior_samples)
      marginal_posterior_samples <- .bt_meta_set(marginal_posterior_samples, "posterior_density", NULL)
      marginal_posterior_samples <- .bt_meta_set(marginal_posterior_samples, "posterior_ordinate", NULL)
      marginal_factor_metadata <- marginal_posterior_samples
      factor_weights <- .prior_factor_level_weight_matrix(
        sample_metadata = marginal_factor_metadata,
        parameter       = parameter,
        samples         = samples
      )
      marginal_factor_support <- .bt_meta_get(marginal_factor_metadata, "support")
      marginal_factor_support <- .marginal_posterior_support_for_context(
        marginal_factor_support,
        can_use_raw_prior_support
      )

      level_names <- attr(marginal_posterior_samples, "level_names")
      if(is.null(level_names) || is.list(level_names)){
        level_names <- .factor_cell_labels(.factor_level_list(marginal_posterior_samples))
      }
      if(is.null(marginal_factor_support)){
        prior_density_context <- .marginal_posterior_prior_density_context(
          samples    = samples,
          prior_list = prior_list,
          n_samples  = n_samples,
          condition_source = parameter_samples,
          raw_coefficients = TRUE
        )
      }

      # the fitted coefficient combination of each level
      level_quantities <- .bt_draws_quantities(marginal_posterior_samples)

      # create output object
      marginal_posterior_samples <- lapply(seq_along(level_names), function(lvl_i){
        temp_marginal_posterior_samples <- marginal_posterior_samples[,lvl_i]
        class(temp_marginal_posterior_samples) <- c(class(temp_marginal_posterior_samples), "marginal_posterior.factor")
        attr(temp_marginal_posterior_samples, "parameter")  <- parameter
        attr(temp_marginal_posterior_samples, "level_name") <- level_names[lvl_i]
        if(!is.null(level_quantities)){
          temp_quantities <- level_quantities[lvl_i, , drop = FALSE]
          temp_quantities$column <- level_names[lvl_i]
          temp_marginal_posterior_samples <- .bt_meta_set(
            temp_marginal_posterior_samples,
            "quantities",
            temp_quantities
          )
        }
        temp_support <- NULL
        if(!is.null(marginal_factor_support)){
          temp_support <- .posterior_support_get(
            marginal_factor_metadata,
            colnames(marginal_posterior_samples)[lvl_i]
          )
          if(is.null(temp_support)){
            temp_support <- .posterior_support_get(
              marginal_factor_metadata,
              level_names[lvl_i]
            )
          }
        }
        if(is.null(temp_support) && !is.null(factor_weights) &&
           !is.null(prior_density_context) &&
           lvl_i <= nrow(factor_weights)){
          weights <- rep(0, length(prior_density_context$column_names))
          names(weights) <- prior_density_context$column_names
          weights[colnames(factor_weights)] <- factor_weights[lvl_i, ]
          temp_support <- .posterior_support_from_prior_context_weights(
            prior_density_context,
            weights
          )
        }
        temp_marginal_posterior_samples <- .posterior_support_set(
          temp_marginal_posterior_samples,
          temp_support
        )
        zero_design <- lvl_i <= nrow(factor_weights) &&
          all(factor_weights[lvl_i, ] == 0)
        if(zero_design || !is.null(marginal_factor_atoms)){
          temp_atoms <- if(zero_design){
            .posterior_atoms_new(
              locations = matrix(0, nrow = 1L, ncol = 1L),
              mass = 1,
              column_names = colnames(marginal_posterior_samples)[lvl_i],
              source = "zero_factor_design",
              declared = TRUE
            )
          }else{
            .posterior_atoms_for_column(
              marginal_factor_atoms,
              lvl_i
            )
          }
          temp_marginal_posterior_samples <- .posterior_atoms_set(
            temp_marginal_posterior_samples,
            temp_atoms
          )
        }
        posterior_density <- .posterior_density_from_sources(
          sources          = posterior_density_sources,
          aliases          = .posterior_density_aliases(
            level_names[lvl_i],
            colnames(marginal_posterior_samples)[lvl_i],
            paste0(parameter, "[", level_names[lvl_i], "]")
          ),
          conditional      = posterior_density_conditional,
          conditional_rule = posterior_density_conditional_rule,
          condition_key    = posterior_density_condition_key
        )
        if(!is.null(posterior_density)){
          temp_marginal_posterior_samples <- .bt_meta_set(temp_marginal_posterior_samples, "posterior_density", posterior_density)
        }
        posterior_ordinate <- .posterior_ordinate_from_sources(
          sources          = posterior_ordinate_sources,
          aliases          = .posterior_density_aliases(
            level_names[lvl_i],
            colnames(marginal_posterior_samples)[lvl_i],
            paste0(parameter, "[", level_names[lvl_i], "]")
          ),
          conditional      = posterior_density_conditional,
          conditional_rule = posterior_density_conditional_rule,
          condition_key    = posterior_density_condition_key
        )
        if(!is.null(posterior_ordinate)){
          temp_marginal_posterior_samples <- .bt_meta_set(temp_marginal_posterior_samples, "posterior_ordinate", posterior_ordinate)
        }
        return(temp_marginal_posterior_samples)
      })
      names(marginal_posterior_samples) <- level_names
      attr(marginal_posterior_samples, "level_names") <- level_names
      class(marginal_posterior_samples) <- c(class(marginal_posterior_samples), "marginal_posterior.factor")

    }else if(inherits(samples[[parameter]], "mixed_posteriors.simple")){

      marginal_posterior_samples <- samples[[parameter]]
      marginal_support <- .marginal_posterior_support_for_context(
        .posterior_support_get(marginal_posterior_samples),
        can_use_raw_prior_support
      )
      if(is.null(marginal_support) && can_use_raw_prior_support){
        marginal_support <- .posterior_support_from_prior_list(prior_list[[parameter]])
      }
      if(is.null(marginal_support)){
        prior_density_context <- .marginal_posterior_prior_density_context(
          samples    = samples,
          prior_list = prior_list,
          n_samples  = n_samples,
          condition_source = samples[[parameter]],
          raw_coefficients = TRUE
        )
        if(!is.null(prior_density_context)){
          weights <- rep(0, length(prior_density_context$column_names))
          names(weights) <- prior_density_context$column_names
          if(parameter %in% names(weights)){
            weights[[parameter]] <- 1
            target <- .marginal_posterior_simple_target(prior_density_context, parameter)
            marginal_support <- .posterior_support_from_prior_context_weights(
              prior_density_context,
              weights,
              output_transformation = target$output_transformation,
              source_transforms     = target$source_transforms
            )
          }
        }
      }

      marginal_posterior_samples <- .posterior_density_attach(
        samples          = marginal_posterior_samples,
        sources          = .posterior_density_sources(samples[[parameter]]),
        parameter        = parameter,
        conditional      = .bt_meta_condition(samples[[parameter]], "conditional"),
        conditional_rule = .bt_meta_condition(samples[[parameter]], "conditional_rule"),
        condition_key    = .bt_meta_condition(samples[[parameter]], "condition_key"),
        allow_unlabeled  = TRUE
      )
      marginal_posterior_samples <- .posterior_ordinate_attach(
        samples          = marginal_posterior_samples,
        sources          = .posterior_ordinate_sources(samples[[parameter]]),
        parameter        = parameter,
        conditional      = .bt_meta_condition(samples[[parameter]], "conditional"),
        conditional_rule = .bt_meta_condition(samples[[parameter]], "conditional_rule"),
        condition_key    = .bt_meta_condition(samples[[parameter]], "condition_key"),
        allow_unlabeled  = TRUE
      )
      marginal_posterior_samples <- .posterior_density_attach(
        samples          = marginal_posterior_samples,
        sources          = .posterior_density_sources(samples),
        parameter        = parameter,
        conditional      = .bt_meta_condition(samples[[parameter]], "conditional"),
        conditional_rule = .bt_meta_condition(samples[[parameter]], "conditional_rule"),
        condition_key    = .bt_meta_condition(samples[[parameter]], "condition_key"),
        allow_unlabeled  = FALSE
      )
      marginal_posterior_samples <- .posterior_ordinate_attach(
        samples          = marginal_posterior_samples,
        sources          = .posterior_ordinate_sources(samples),
        parameter        = parameter,
        conditional      = .bt_meta_condition(samples[[parameter]], "conditional"),
        conditional_rule = .bt_meta_condition(samples[[parameter]], "conditional_rule"),
        condition_key    = .bt_meta_condition(samples[[parameter]], "condition_key"),
        allow_unlabeled  = FALSE
      )

      marginal_posterior_samples <- .posterior_support_set(
        marginal_posterior_samples,
        marginal_support
      )
      class(marginal_posterior_samples) <- c(class(marginal_posterior_samples), "marginal_posterior.simple")

    }else{

      stop(
        "'marginal_posterior()' is not supported for ",
        .marginal_posterior_samples_type(samples[[parameter]]),
        " posterior samples ('", parameter, "'). Marginal posterior ",
        "distributions are available for simple, factor, and formula parameters.",
        call. = FALSE
      )
    }


    # add prior densities
    if(prior_samples){

      if(is.null(prior_density_context)){
        prior_density_context <- .marginal_posterior_prior_density_context(
          samples    = samples,
          prior_list = prior_list,
          n_samples  = n_samples,
          condition_source = samples[[parameter]],
          raw_coefficients = TRUE
        )
      }

      if(inherits(samples[[parameter]], "mixed_posteriors.factor")){

        if(is.null(factor_weights)){
          factor_weights <- .prior_factor_level_weight_matrix(
            sample_metadata = marginal_factor_metadata,
            parameter       = parameter,
            samples         = samples
          )
        }

        for(lvl_i in seq_along(level_names)){
          weights <- rep(0, length(prior_density_context$column_names))
          names(weights) <- prior_density_context$column_names
          weights[colnames(factor_weights)] <- factor_weights[lvl_i, ]

          prior_density <- .prior_density_from_context(
            prior_density_context,
            weights
          )
          marginal_posterior_samples[[level_names[lvl_i]]] <- .bt_meta_update(
            marginal_posterior_samples[[level_names[lvl_i]]],
            linear_weights = weights,
            prior_density  = prior_density,
            prior_context  = prior_density_context,
            components     = .marginal_posterior_components(
              context         = prior_density_context,
              weights         = weights,
              model_component = .bt_draws_model_component(samples[[parameter]]),
              n_values        = length(marginal_posterior_samples[[level_names[lvl_i]]]),
              samples         = samples
            )
          )
        }

      }else if(inherits(samples[[parameter]], "mixed_posteriors.simple")){

        weights <- rep(0, length(prior_density_context$column_names))
        names(weights) <- prior_density_context$column_names
        weights[[parameter]] <- 1
        target <- .marginal_posterior_simple_target(prior_density_context, parameter)

        prior_density <- .prior_density_from_context(
          prior_density_context,
          weights,
          source_transforms     = target$source_transforms,
          output_transformation = target$output_transformation
        )
        marginal_posterior_samples <- .bt_meta_update(
          marginal_posterior_samples,
          linear_weights = weights,
          prior_density  = prior_density,
          prior_context  = prior_density_context,
          components     = .marginal_posterior_components(
            context               = prior_density_context,
            weights               = weights,
            model_component       = .bt_draws_model_component(samples[[parameter]]),
            n_values              = length(marginal_posterior_samples),
            samples               = samples,
            source_transforms     = target$source_transforms,
            output_transformation = target$output_transformation
          )
        )
      }

      if(inherits(samples[[parameter]], "mixed_posteriors.factor")){
        # the list of levels (the draws of each level carry it already)
        marginal_posterior_samples <- .bt_meta_set(marginal_posterior_samples, "prior_context", prior_density_context)
      }
    }
  }

  if(inherits(samples[[parameter]], "mixed_posteriors.simple")){
    stored_prior_density <- .bt_meta_get(samples[[parameter]], "prior_density")
    if(!is.null(stored_prior_density)){
      marginal_posterior_samples <- .bt_meta_set(marginal_posterior_samples, "prior_density", stored_prior_density)
      stored_prior_context <- .bt_meta_get(samples[[parameter]], "prior_context")
      if(!is.null(stored_prior_context)){
        marginal_posterior_samples <- .bt_meta_set(marginal_posterior_samples, "prior_context", stored_prior_context)
      }
    }
  }
  marginal_posterior_samples <- .marginal_posterior_attach_precomputed_metadata(
    marginal              = marginal_posterior_samples,
    samples               = samples,
    parameter             = parameter,
    condition_source      = samples[[parameter]]
  )

  marginal_posterior_samples <- .marginal_posterior_set_condition_attributes(
    marginal_posterior_samples,
    samples,
    condition_source = samples[[parameter]]
  )
  class(marginal_posterior_samples) <- c(class(marginal_posterior_samples), "marginal_posterior")

  # the transformed marginal posterior: the draws with their metadata
  if(!is.null(transformation)){
    marginal_posterior_samples <- .bt_posterior_transform(
      marginal_posterior_samples,
      .bt_posterior_transformation(transformation, transformation_arguments)
    )
  }
  return(marginal_posterior_samples)
}

# Mixed posteriors transformed by posterior_transform() carry the prior of
# their untransformed values ('prior_list'), from which marginal_posterior()
# would build their prior densities and linear predictors; only draws whose
# prior density is their stored metadata (prior_none(), e.g.
# parameter_mixed_posterior()) remain valid inputs.
.marginal_posterior_check_untransformed <- function(samples, parameters){

  for(parameter in parameters){
    prior <- attr(samples[[parameter]], "prior_list", exact = TRUE)
    if(length(.bt_draws_output_transformations(samples[[parameter]])) > 0L &&
       !(is.prior(prior) && is.prior.none(prior))){
      stop(
        "The posterior samples of '", parameter, "' were transformed with ",
        "'posterior_transform()', but 'marginal_posterior()' builds their prior ",
        "densities from the untransformed prior distributions. Pass the ",
        "transformation to 'marginal_posterior(transformation = )' instead.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

# The column table of estimated marginal means: each level is a prediction at
# the stated values of the manipulated predictors (factor levels, or
# continuous values in standard deviations), not a combination of fitted
# coefficients, so it declares no fitted coordinates.
.marginal_posterior_level_quantities <- function(formula_parameter, components,
                                                 level_frame){

  if(is.null(level_frame)){
    parts <- list(.bt_label_parts(
      components        = components,
      formula_parameter = formula_parameter,
      marginal          = TRUE
    ))
  }else{
    parts <- lapply(seq_len(nrow(level_frame)), function(level){
      .bt_label_parts(
        components        = components,
        formula_parameter = formula_parameter,
        levels            = stats::setNames(
          vapply(level_frame, function(column) as.character(column[[level]]), character(1)),
          names(level_frame)
        ),
        marginal          = TRUE
      )
    })
  }
  .bt_draws_quantity_table(
    columns      = .bt_label(parts, style = "plot"),
    quantity_ids = rep("", length(parts)),
    dependencies = rep(list(character()), length(parts)),
    weights      = rep(list(numeric()), length(parts)),
    label_parts  = parts
  )
}

.marginal_posterior_prior_density_context <- function(samples, prior_list,
                                                       column_names = NULL,
                                                       n_samples = 10000,
                                                       condition_source = NULL,
                                                       raw_coefficients = FALSE){

  condition_metadata <- .marginal_posterior_condition_metadata(
    samples,
    condition_source = condition_source
  )
  prior_density_context <- .bt_meta_get(samples, "prior_context")
  if(!is.null(prior_density_context) &&
     .marginal_posterior_context_matches_condition(prior_density_context, condition_metadata)){
    return(.marginal_posterior_canonical_context(
      prior_density_context,
      raw_coefficients = raw_coefficients
    ))
  }

  prior_list <- .marginal_posterior_canonical_prior_list(
    prior_list,
    raw_coefficients = raw_coefficients
  )

  if(is.null(column_names)){
    column_names <- unique(unlist(lapply(names(prior_list), function(parameter_name){
      parameter_prior <- prior_list[[parameter_name]]
      if(is.prior(parameter_prior)){
        .prior_linear_prior_columns(parameter_name, parameter_prior)
      }else{
        .prior_linear_prior_columns(parameter_name, parameter_prior[[1]])
      }
    }), use.names = FALSE))
  }

  # build failures propagate: the context carries the support, atoms and
  # mixture components of the marginal posterior
  .prior_density_build_context(
    prior_list       = prior_list,
    column_names     = column_names,
    n_grid           = max(16L, n_samples),
    conditional      = condition_metadata[["conditional"]],
    conditional_rule = condition_metadata[["conditional_rule"]],
    condition_event  = condition_metadata[["condition_event"]]
  )
}

# Whether the formula's linear predictor uses log(intercept). Persisted fit
# metadata (formula design, formula-scale information) is authoritative; an
# explicit "log(intercept)" attribute on 'formula' must agree with it.
.marginal_posterior_log_intercept <- function(samples, formula_log_intercept,
                                              formula_parameter){

  declared <- unlist(lapply(samples, function(parameter_samples){
    if(!identical(.bt_meta_get(parameter_samples, "formula_parameter"), formula_parameter)){
      return(NULL)
    }
    .bt_meta_get(parameter_samples, "log_intercept")
  }), use.names = FALSE)

  formula_scale <- .bt_meta_get(samples, "formula_scale")
  if(is.list(formula_scale) && length(formula_parameter) == 1L &&
     !is.null(formula_scale[[formula_parameter]])){
    scale_log_intercept <- attr(formula_scale[[formula_parameter]], "log_intercept", exact = TRUE)
    if(!is.null(scale_log_intercept)){
      declared <- c(declared, isTRUE(scale_log_intercept))
    }
  }

  declared <- unique(declared)
  if(anyNA(declared) || length(declared) > 1L){
    stop(
      "The mixed models disagree on whether the formula for '",
      formula_parameter, "' uses log(intercept).",
      call. = FALSE
    )
  }

  if(!is.null(formula_log_intercept)){
    formula_log_intercept <- isTRUE(formula_log_intercept)
    if(length(declared) == 1L && !identical(declared, formula_log_intercept)){
      stop(
        "The 'log(intercept)' attribute of 'formula' does not match the fitted ",
        "formula for '", formula_parameter, "'.",
        call. = FALSE
      )
    }
    return(formula_log_intercept)
  }

  length(declared) == 1L && isTRUE(declared)
}

# Support, posterior-atom and component metadata for one formula level. Build
# failures propagate. A log(intercept) term enters the linear predictor through
# the log of the intercept (a log source), in which the level is linear in the
# fitted coefficients, also for unscaled (transform_scaled) coefficients: its
# support, components and atoms are derived through that log source. Its
# linear weights are not those of a linear combination of the coefficients
# themselves ('joint_prior_transformation').
.marginal_posterior_formula_level_metadata <- function(marginal, samples, prior_list,
                                                       prior_density_context, weights,
                                                       column_name,
                                                       source_transforms = NULL,
                                                       model_component = NULL){

  log_columns <- intersect(names(source_transforms)[source_transforms == "log"], colnames(weights))
  if(length(log_columns) > 0L && any(weights[, log_columns] != 0)){
    marginal <- .bt_meta_set(marginal, "joint_prior_transformation", "log_intercept")
  }
  support <- .posterior_support_from_prior_context_weights(
    prior_density_context,
    weights,
    source_transforms = source_transforms
  )
  components <- .marginal_posterior_components(
    context           = prior_density_context,
    weights           = weights,
    model_component   = model_component,
    n_values          = length(marginal),
    samples           = samples,
    source_transforms = source_transforms
  )
  marginal <- .posterior_support_set(marginal, support)
  marginal <- .posterior_components_set(marginal, components)

  atoms <- .posterior_atoms_formula(
    samples,
    prior_list,
    weights,
    column_name       = column_name,
    source_transforms = source_transforms
  )
  if(!is.null(atoms)){
    marginal <- .posterior_atoms_set(marginal, atoms)
  }
  fitted_weights <- weights
  if(isTRUE(.bt_meta_get(samples,"transform_scaled"))){
    context <- .prior_density_context(prior_list,colnames(weights),formula_scale=.bt_meta_get(samples,"formula_scale"))
    fitted_weights <- matrix(0,nrow(weights),ncol(weights),dimnames=dimnames(weights))
    for(i in seq_len(nrow(weights))){
      standardized <- .prior_density_context_standardized_weights(context,weights[i,],source_transforms)
      fitted_weights[i,names(standardized)] <- standardized
    }
  }
  ordered <- .bt_ordered_formula_projections(samples,fitted_weights,source_transforms)
  if(!is.null(ordered)){
    values <- as.vector(do.call(rbind,lapply(ordered,`[[`,"values")))
    atom <- as.vector(do.call(rbind,lapply(ordered,`[[`,"atom")))
    exact <- as.vector(do.call(rbind,lapply(ordered,`[[`,"exact")))
    if(length(values)!=length(marginal)) .bt_ordered_stop("Ordered formula sources do not align with marginal draw rows. Recreate marginal posteriors from the source fits.")
    replace <- !is.na(atom) | exact
    marginal <- .bt_draws_transform_values(marginal,function(old){
      old[replace] <- values[replace]
      old
    })
  }

  marginal
}

# Information on a formula predictor that has no main-effect term and enters
# the formula only through interactions (e.g., `x` and `h` in
# `~ g + g:x + g:h`), in the shape of the term information of
# `marginal_posterior()`: a factor component of a factor term takes the
# fitted levels and contrast stored for it on that term, and any other
# component is continuous.
.marginal_posterior_interaction_predictor_info <- function(predictor, parameter,
                                                           priors_info){

  containing <- Filter(function(term_info){
    predictor %in% c(term_info[["term_components"]], term_info[["interaction_terms"]])
  }, priors_info)
  if(length(containing) == 0L){
    stop(
      "The formula predictor '", predictor, "' is not part of any term of the ",
      "mixed posterior samples.",
      call. = FALSE
    )
  }
  factor_info <- Filter(function(term_info){
    predictor %in% term_info[["factor_terms"]]
  }, containing)
  if(length(factor_info) == 0L){
    return(list(term = parameter, factor = FALSE))
  }

  level_names <- factor_info[[1L]][["level_names"]]
  if(is.list(level_names)){
    level_names <- level_names[[predictor]]
  }
  # A term that codes the predictor by level indicators records the
  # independent coding for it; the fitted contrast of the predictor is the
  # one recorded by a term that codes it by that contrast.
  contrasts <- unlist(lapply(factor_info, function(term_info){
    factor_contrasts <- term_info[["factor_contrasts"]]
    if(predictor %in% names(factor_contrasts)){
      unname(factor_contrasts[[predictor]])
    }
  }), use.names = FALSE)
  contrast <- if(length(contrasts) > 0L){
    coded <- contrasts[contrasts != "contr.independent"]
    if(length(coded) > 0L) coded[[1L]] else contrasts[[1L]]
  }
  if(is.null(level_names) || is.null(contrast)){
    .bt_stop_refit_required(
      "The mixed posterior samples lack the fitted levels or contrast of the ",
      "factor predictor '", predictor, "'. Recreate them with mix_posteriors() ",
      "or as_mixed_posteriors() from models fitted with the current version."
    )
  }

  list(
    term             = parameter,
    factor           = TRUE,
    level_names      = level_names,
    factor_contrasts = stats::setNames(contrast, predictor),
    treatment        = identical(contrast, "contr.treatment"),
    independent      = identical(contrast, "contr.independent"),
    orthonormal      = identical(contrast, "contr.orthonormal"),
    meandif          = identical(contrast, "contr.meandif"),
    ordered          = .prior_ordered_is_contrast_name(contrast)
  )
}

# Contrast flag of a mixed factor posterior. Objects created before ordered
# factors were supported (BayesTools 0.3.0) lack part of this metadata; their
# remedy is models fitted with this version (class BayesTools_refit_required).
.marginal_posterior_factor_flag <- function(prior_info, field, parameter){

  value <- prior_info[[field]]
  if(!is.logical(value) || length(value) != 1L || is.na(value)){
    .bt_stop_refit_required(
      "The mixed posterior samples of '", parameter, "' lack the factor ",
      "metadata recorded by the current BayesTools version (missing: '", field,
      "'); they were likely created by BayesTools 0.3.0. Recreate them with ",
      "mix_posteriors() or as_mixed_posteriors() from models fitted with the ",
      "current version."
    )
  }

  value
}

.marginal_posterior_samples_type <- function(parameter_samples){

  if(inherits(parameter_samples, "mixed_posteriors.weightfunction")){
    return("weightfunction")
  }
  if(inherits(parameter_samples, "mixed_posteriors.bias")){
    return("publication-bias")
  }
  if(inherits(parameter_samples, "mixed_posteriors.phacking")){
    return("p-hacking")
  }
  if(inherits(parameter_samples, "mixed_posteriors.vector")){
    return("vector (e.g., Dirichlet simplex)")
  }

  "these"
}

# Per-draw component index and per-component exact supports of a mixture
# marginal (Savage-Dickey mixes per-component posterior ordinates when the
# components' supports differ). The component of a draw is its model in a
# model-mixture ensemble (mix_posteriors()), or, for a single fit
# (as_mixed_posteriors()), the combination of the component indices of the
# mixture and spike-and-slab priors entering the quantity. A marginal evaluated
# at several design rows stores each draw's rows consecutively, so the index
# repeats per row. Build failures propagate; marginals whose quantity involves
# no mixture have no components (the pooled posterior ordinate is used).
.marginal_posterior_components <- function(context, weights, model_component, n_values,
                                           samples = NULL, source_transforms = NULL,
                                           output_transformation = NULL){

  n_rows <- if(is.null(dim(weights))) 1L else nrow(weights)
  if(n_rows < 1L){
    return(NULL)
  }

  if(inherits(context, "prior_density_model_mixture_context")){
    if(is.null(model_component) || length(model_component) * n_rows != n_values){
      return(NULL)
    }
    return(.posterior_components_new(
      index    = rep(model_component, each = n_rows),
      supports = .posterior_support_model_components(
        context, weights,
        output_transformation = output_transformation,
        source_transforms     = source_transforms
      ),
      keys     = matrix(
        seq_along(context$model_weights),
        ncol = 1L,
        dimnames = list(NULL, ".model")
      )
    ))
  }

  if(!inherits(samples, "as_mixed_posteriors") ||
     !(inherits(context, "prior_density_context") ||
       inherits(context, "prior_density_conditional_context"))){
    return(NULL)
  }

  .marginal_posterior_indicator_components(
    context               = context,
    weights               = weights,
    samples               = samples,
    n_values              = n_values,
    source_transforms     = source_transforms,
    output_transformation = output_transformation
  )
}

.marginal_posterior_indicator_components <- function(context, weights, samples,
                                                     n_values, source_transforms = NULL,
                                                     output_transformation = NULL){

  parameters <- .posterior_components_mixture_parameters(context, weights, source_transforms)
  if(length(parameters) == 0L){
    return(NULL)
  }

  n_rows <- if(is.null(dim(weights))) 1L else nrow(weights)
  indicators <- lapply(parameters, function(parameter){
    component <- .bt_meta_get(samples[[parameter]], "component")
    if(is.null(component) ||
       !.bt_meta_get(samples[[parameter]], "component_source") %in% c("mixture", "spike_and_slab") ||
       length(component) * n_rows != n_values){
      stop("Mixture component indices are unavailable.", call. = FALSE)
    }
    as.numeric(component)
  })
  indicators <- do.call(cbind, indicators)
  colnames(indicators) <- parameters

  draw_keys <- do.call(paste, c(as.data.frame(indicators), sep = "\r"))
  unique_keys <- unique(draw_keys)
  index <- match(draw_keys, unique_keys)
  keys <- indicators[match(unique_keys, draw_keys), , drop = FALSE]
  rownames(keys) <- NULL

  .posterior_components_new(
    index    = rep(index, each = n_rows),
    supports = .posterior_components_supports(
      context               = context,
      keys                  = keys,
      weights               = weights,
      output_transformation = output_transformation,
      source_transforms     = source_transforms
    ),
    keys     = keys
  )
}

# The prior-density target of the simple marginal posterior of 'parameter' in
# 'context': the coefficient itself, or for the unscaled intercept of a
# log-intercept formula scaling, which is not linear in the fitted
# coefficients, the exp of its log, which is (a log source with an exp output
# transformation).
.marginal_posterior_simple_target <- function(context, parameter){

  target <- list(source_transforms = NULL, output_transformation = NULL)
  for(transform in context$transforms){
    if(isTRUE(transform$log_intercept) && identical(transform$intercept, parameter)){
      target$source_transforms     <- stats::setNames("log", parameter)
      target$output_transformation <- "exp"
    }
  }
  target
}

# Prior lists used by marginal prior-density contexts. Structural zeros are
# expressed on every coefficient column, and monitored coefficient columns
# (raw_coefficients) are the raw JAGS nodes: a formula prior's 'multiply_by'
# scales only the linear predictor, never the coefficient itself.
.marginal_posterior_canonical_context <- function(context, raw_coefficients = FALSE){

  if(is.null(context)){
    return(context)
  }

  if(!is.null(context$prior_list)){
    context$prior_list <- .marginal_posterior_canonical_prior_list(
      context$prior_list,
      raw_coefficients = raw_coefficients
    )
  }
  if(!is.null(context$prior_lists)){
    context$prior_lists <- lapply(
      context$prior_lists,
      .marginal_posterior_canonical_prior_list,
      raw_coefficients = raw_coefficients
    )
  }

  context
}

.marginal_posterior_canonical_prior_list <- function(prior_list, raw_coefficients = FALSE){

  if(!is.list(prior_list)){
    return(prior_list)
  }

  for(parameter in names(prior_list)){
    entry <- prior_list[[parameter]]
    if(is.null(entry)){
      next
    }
    if(is.prior(entry)){
      prior_list[[parameter]] <- .marginal_posterior_structural_zero_prior(entry)
    }else if(is.list(entry)){
      K <- .marginal_posterior_model_list_dimension(entry)
      for(i in seq_along(entry)){
        if(!is.prior(entry[[i]])){
          next
        }
        location <- entry[[i]]$parameters[["location"]]
        if(K > 1L && is.prior.point(entry[[i]]) && !is.prior.vector(entry[[i]]) &&
           is.numeric(location) && length(location) == 1L){
          # a scalar point prior (e.g. the spike(0) filled in for a model that
          # omits the term) fixes every coefficient column of the term
          entry[[i]] <- .marginal_posterior_zero_vector_prior(
            entry[[i]], K,
            location = entry[[i]]$parameters[["location"]]
          )
        }else{
          entry[[i]] <- .marginal_posterior_structural_zero_prior(entry[[i]])
        }
      }
      prior_list[[parameter]] <- entry
    }
  }

  if(isTRUE(raw_coefficients)){
    prior_list <- .marginal_posterior_strip_multiply_by(prior_list)
  }

  prior_list
}

# A structurally zero coefficient vector (e.g., an ordered prior with a point
# total at zero) as a point at zero on every coefficient column.
.marginal_posterior_structural_zero_prior <- function(prior){

  if(!.posterior_atoms_is_ordered_zero_total(prior)){
    return(prior)
  }

  .marginal_posterior_zero_vector_prior(prior, .prior_linear_prior_dimension(prior))
}

# Number of coefficient columns of a model-averaged term, from the models whose
# prior is not a scalar point.
.marginal_posterior_model_list_dimension <- function(priors){

  for(prior in priors){
    if(!is.prior(prior) || (is.prior.point(prior) && !is.prior.vector(prior))){
      next
    }
    K <- .prior_linear_prior_dimension(prior)
    if(length(K) == 1L && !is.na(K)){
      return(as.integer(K))
    }
  }

  1L
}

.marginal_posterior_zero_vector_prior <- function(prior, K, location = 0){

  zero_prior <- prior("mpoint", list(location = location, K = K))
  model_weight <- .prior_model_weight(prior)
  if(!is.null(model_weight)){
    zero_prior <- .set_prior_model_weight(zero_prior, model_weight)
  }
  for(attribute in c("parameter", "multiply_by")){
    attr(zero_prior, attribute) <- attr(prior, attribute, exact = TRUE)
  }

  zero_prior
}

.marginal_posterior_strip_multiply_by <- function(x){

  if(!is.list(x)){
    return(x)
  }

  attr(x, "multiply_by") <- NULL
  x_is_prior <- is.prior(x)
  for(i in seq_along(x)){
    # Recurse into prior components (mixtures, ordered totals) and into
    # containers of priors, never into prior parameter lists.
    if(is.list(x[[i]]) && (is.prior(x[[i]]) || !x_is_prior)){
      x[[i]] <- .marginal_posterior_strip_multiply_by(x[[i]])
    }
  }

  x
}

.marginal_posterior_context_condition_metadata <- function(context){

  .marginal_posterior_normalize_condition_metadata(
    conditional      = context[["conditional"]],
    conditional_rule = context[["conditional_rule"]],
    condition_key    = context[["condition_key"]],
    condition_event  = context[["condition_event"]]
  )
}

.marginal_posterior_normalize_condition_metadata <- function(conditional = NULL,
                                                             conditional_rule = NULL,
                                                             condition_key = NULL,
                                                             condition_event = NULL){

  if(is.null(conditional) && !is.null(condition_event)){
    conditional <- condition_event[["conditional"]]
  }
  if(is.null(conditional_rule) && !is.null(condition_event)){
    conditional_rule <- condition_event[["conditional_rule"]]
  }
  if(is.null(condition_key) && !is.null(condition_event)){
    condition_key <- condition_event[["condition_key"]]
  }
  if(is.null(conditional_rule) && length(conditional) > 0L){
    conditional_rule <- "AND"
  }
  if(is.null(condition_key) && length(conditional) > 0L){
    condition_key <- .condition_event_key(conditional, conditional_rule)
  }

  list(
    conditional      = conditional,
    conditional_rule = conditional_rule,
    condition_key    = condition_key,
    condition_event  = condition_event
  )
}

.marginal_posterior_context_matches_condition <- function(context,
                                                          condition_metadata){

  context_metadata <- .marginal_posterior_context_condition_metadata(context)

  requested_conditional <- .posterior_density_normalize_condition(
    condition_metadata[["conditional"]]
  )
  context_conditional <- .posterior_density_normalize_condition(
    context_metadata[["conditional"]]
  )

  if(length(requested_conditional) == 0L &&
     length(context_conditional) == 0L){
    return(TRUE)
  }
  if(length(requested_conditional) == 0L ||
     length(context_conditional) == 0L){
    return(FALSE)
  }

  requested_key <- condition_metadata[["condition_key"]]
  if(is.null(requested_key)){
    requested_key <- .condition_event_key(
      requested_conditional,
      condition_metadata[["conditional_rule"]]
    )
  }
  context_key <- context_metadata[["condition_key"]]
  if(is.null(context_key)){
    context_key <- .condition_event_key(
      context_conditional,
      context_metadata[["conditional_rule"]]
    )
  }

  identical(as.character(context_key), as.character(requested_key))
}

.marginal_posterior_condition_metadata <- function(samples, condition_source = NULL){

  sources <- list(samples, condition_source)
  sources <- sources[!vapply(sources, is.null, logical(1))]
  if(length(sources) == 0L){
    return(list(
      conditional      = NULL,
      conditional_rule = NULL,
      condition_key    = NULL,
      condition_event  = NULL
    ))
  }

  source_metadata <- lapply(sources, function(source){
    # one read (and check) of the conditioning metadata of each source
    condition <- .bt_meta_get(source, "condition")
    condition_event <- condition[["resolved_condition_event"]]
    if(is.null(condition_event)){
      condition_event <- condition[["condition_event"]]
    }

    .marginal_posterior_normalize_condition_metadata(
      conditional      = condition[["conditional"]],
      conditional_rule = condition[["conditional_rule"]],
      condition_key    = condition[["condition_key"]],
      condition_event  = condition_event
    )
  })

  conditioned <- vapply(
    source_metadata,
    function(metadata) length(metadata[["conditional"]]) > 0L,
    logical(1)
  )
  if(any(conditioned)){
    return(source_metadata[[which(conditioned)[1L]]])
  }

  present <- vapply(
    source_metadata,
    function(metadata){
      !is.null(metadata[["conditional"]]) ||
        !is.null(metadata[["conditional_rule"]]) ||
        !is.null(metadata[["condition_key"]]) ||
        !is.null(metadata[["condition_event"]])
    },
    logical(1)
  )
  if(any(present)){
    return(source_metadata[[which(present)[1L]]])
  }

  source_metadata[[1L]]
}

.marginal_posterior_can_use_raw_prior_support <- function(samples, parameter){

  if(isTRUE(.bt_meta_get(samples, "transform_scaled")) ||
     isTRUE(.bt_meta_get(samples[[parameter]], "transform_scaled"))){
    return(FALSE)
  }

  condition_metadata <- .marginal_posterior_condition_metadata(
    samples,
    condition_source = samples[[parameter]]
  )

  length(condition_metadata[["conditional"]]) == 0L
}

.marginal_posterior_support_for_context <- function(support,
                                                    can_use_raw_prior_support){

  if(isTRUE(can_use_raw_prior_support)){
    return(support)
  }
  if(.posterior_support_is_raw_prior(support)){
    return(NULL)
  }

  support
}

.marginal_posterior_attach_precomputed_metadata <- function(marginal, samples,
                                                            parameter,
                                                            condition_source = NULL){

  if(!is.list(marginal)){
    return(marginal)
  }

  density_sources <- .posterior_density_sources(samples)
  ordinate_sources <- .posterior_ordinate_sources(samples)
  if(length(density_sources) == 0L && length(ordinate_sources) == 0L){
    return(marginal)
  }

  condition_metadata <- .marginal_posterior_condition_metadata(
    samples,
    condition_source = condition_source
  )
  sample_names <- names(marginal)

  for(i in seq_along(marginal)){
    sample_name <- if(!is.null(sample_names)) sample_names[[i]] else NULL
    child_level <- attr(marginal[[i]], "level_name", exact = TRUE)
    if(is.null(child_level)){
      child_level <- attr(marginal[[i]], "level", exact = TRUE)
    }
    aliases <- .posterior_density_aliases(
      sample_name,
      child_level,
      if(!is.null(sample_name)) paste0(parameter, "[", sample_name, "]"),
      if(!is.null(child_level)) paste0(parameter, "[", child_level, "]"),
      attr(marginal[[i]], "factor_cell_names", exact = TRUE)
    )

    if(length(density_sources) > 0L){
      density <- .posterior_density_from_sources(
        sources          = density_sources,
        aliases          = aliases,
        conditional      = condition_metadata[["conditional"]],
        conditional_rule = condition_metadata[["conditional_rule"]],
        condition_key    = condition_metadata[["condition_key"]],
        allow_unlabeled  = FALSE
      )
      if(!is.null(density) &&
         is.null(.bt_meta_get(marginal[[i]], "posterior_density"))){
        marginal[[i]] <- .bt_meta_set(marginal[[i]], "posterior_density", density)
      }
    }

    if(length(ordinate_sources) > 0L){
      ordinate <- .posterior_ordinate_from_sources(
        sources          = ordinate_sources,
        aliases          = aliases,
        conditional      = condition_metadata[["conditional"]],
        conditional_rule = condition_metadata[["conditional_rule"]],
        condition_key    = condition_metadata[["condition_key"]],
        allow_unlabeled  = FALSE
      )
      if(!is.null(ordinate) &&
         is.null(.bt_meta_get(marginal[[i]], "posterior_ordinate"))){
        marginal[[i]] <- .bt_meta_set(marginal[[i]], "posterior_ordinate", ordinate)
      }
    }
  }

  marginal
}

.marginal_posterior_set_condition_attributes <- function(marginal, samples,
                                                         condition_source = NULL){

  condition_metadata <- .marginal_posterior_condition_metadata(
    samples,
    condition_source = condition_source
  )
  conditional <- condition_metadata[["conditional"]]
  conditional_rule <- condition_metadata[["conditional_rule"]]
  condition_key <- condition_metadata[["condition_key"]]
  condition_event <- condition_metadata[["condition_event"]]

  if(is.null(conditional) && is.null(conditional_rule) &&
     is.null(condition_key) && is.null(condition_event)){
    return(marginal)
  }
  if(is.null(conditional_rule)){
    conditional_rule <- "AND"
  }
  if(is.null(condition_key)){
    condition_key <- .condition_event_key(conditional, conditional_rule)
  }

  set_attrs <- function(x){
    condition <- .bt_meta_get(x, "condition")
    if(is.null(condition)){
      condition <- list()
    }
    condition[["conditional"]]      <- conditional
    condition[["conditional_rule"]] <- conditional_rule
    condition[["condition_key"]]    <- condition_key
    condition[["averaged"]]         <- .condition_is_averaged(conditional)
    if(!is.null(condition_event)){
      condition[["condition_event"]]          <- condition_event
      condition[["resolved_condition_event"]] <- condition_event
    }
    .bt_meta_set(x, "condition", if(length(condition) > 0L) condition)
  }

  if(is.list(marginal)){
    marginal_attributes <- attributes(marginal)
    for(i in seq_along(marginal)){
      marginal[[i]] <- set_attrs(marginal[[i]])
    }
    attributes(marginal) <- marginal_attributes
  }
  marginal <- set_attrs(marginal)

  marginal
}

.marginal_posterior_parameter_samples <- function(samples, parameter){

  parameter_samples <- samples[[parameter]]

  if(is.list(parameter_samples)){
    out <- lapply(parameter_samples, as.numeric)
  }else{
    out <- list(as.numeric(parameter_samples))
    names(out) <- parameter
  }

  out
}

.marginal_posterior_parameter_prior_densities <- function(samples, parameter){

  parameter_samples <- samples[[parameter]]

  if(is.list(parameter_samples)){
    out <- lapply(parameter_samples, .bt_meta_get, field = "prior_density")
  }else{
    out <- list(.bt_meta_get(parameter_samples, "prior_density"))
    names(out) <- parameter
  }

  out
}

.marginal_posterior_parameter_posterior_densities <- function(samples, parameter){

  parameter_samples <- samples[[parameter]]

  if(is.list(parameter_samples)){
    out <- .posterior_density_child_attributes(parameter_samples)
  }else{
    out <- list(.bt_meta_get(parameter_samples, "posterior_density"))
    names(out) <- parameter
  }

  out
}

# The per-draw multipliers of a formula term (the 'multiply_by' of the term's
# prior in the model of each draw), or NULL when no draw scales the term.
.marginal_posterior_term_multiplier <- function(term, prior_list, posterior, model_component, simple_list = FALSE){

  multiplier <- .get_combined_parameter_scaling_factor_matrix(
    term,
    prior_list  = prior_list,
    posterior   = posterior,
    model_component = model_component,
    nrow        = 1L,
    simple_list = simple_list
  )[1L, ]
  if(isTRUE(all(multiplier == 1))){
    return(NULL)
  }

  multiplier
}

.get_combined_parameter_scaling_factor_matrix <- function(term, prior_list, posterior, model_component, nrow, simple_list = FALSE){

  if(simple_list){
    temp_multiply_by <- .get_parameter_scaling_factor_matrix(term, prior_list, posterior, nrow = nrow, ncol = nrow(posterior))
  }else{
    temp_multiply_by <- do.call(cbind, lapply(unique(model_component), function(m){
      temp_prior_list <- lapply(prior_list, function(parameter_priors) parameter_priors[[m]])
      temp_posterior  <- posterior[model_component == m,,drop=FALSE]
      return(.get_parameter_scaling_factor_matrix(term, temp_prior_list, temp_posterior, nrow = nrow, ncol = sum(model_component == m)))
    }))
  }

  return(temp_multiply_by)
}
