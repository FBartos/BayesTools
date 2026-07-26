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
#' @param samples model-averaged posterior samples created by \code{mix_posteriors()}
#' @param parameter parameter of interest
#' @param formula model formula (needs to be specified if \code{parameter} was part of a formula)
#' @param at named list with predictor levels of the formula for which marginalization
#' should be performed. If a predictor level is missing, \code{0} is used for continuous
#' predictors, the baseline factor level is used for factors with \code{contrast = "treatment"} prior
#' distributions, and the parameter is completely omitted for factors with
#' \code{contrast = "meandif"}, \code{contrast = "orthonormal"},
#' \code{contrast = "independent"}, and ordered-factor levels.
#' @param prior_samples whether marginal prior distributions should be generated
#' @param use_formula whether the parameter should be evaluated as a part of supplied formula
#' @param n_samples controls the numerical grid used for model-averaged
#' prior densities
#' @inheritParams density.prior
#'
#' @details When the mixed posterior samples carry deterministic
#' \code{posterior_density}, \code{posterior_ordinate}, or
#' \code{posterior_support} metadata, \code{marginal_posterior()} propagates
#' matching metadata to the returned marginal posterior. Matching uses the
#' parameter name, list/level names such as \code{theta[A]}, and conditional
#' metadata when present. Exact support metadata is propagated even when
#' \code{prior_samples = FALSE}; requesting prior samples adds prior-density
#' metadata but does not replace already attached posterior support.
#' Transformations drop stored posterior density and ordinate metadata because
#' those estimates are no longer on the returned scale; exact support metadata
#' is transformed when the transformation is supported. If support metadata is
#' absent, support is inferred from prior metadata only when the posterior
#' samples are on the raw, unconditioned prior scale; otherwise the deterministic
#' prior-density context is used so formula-scale transformations and
#' conditional model restrictions are respected.
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
  if(!is.null(formula) && !is.language(formula))
    stop("'formula' must be a formula")
  if(!is.null(at) && !is.list(at))
    stop("'at' must be a list")
  check_bool(prior_samples, "prior_samples")
  check_bool(use_formula, "use_formula")
  .check_transformation_input(transformation, transformation_arguments, transformation_settings)


  # deal formula vs non-formula marginal posterior
  if(use_formula && inherits(samples[[parameter]], "mixed_posteriors.formula")){

      # remove the specified response (would crash the model.frame if not included)
      formula <- .remove_response(formula)
      formula_parameter <- attr(samples[[parameter]], "formula_parameter")

      ### extract the terms information from the formula
      formula_terms          <- stats::terms(formula)
      has_intercept          <- attr(formula_terms, "intercept") == 1
      predictors             <- as.character(attr(formula_terms, "variables"))[-1]
      model_terms            <- c(if(has_intercept) "intercept", attr(formula_terms, "term.labels"))

      JAGS_model_terms <- JAGS_parameter_names(parameters = model_terms, formula_parameter = formula_parameter)
      JAGS_predictors  <- JAGS_parameter_names(parameters = predictors, formula_parameter = formula_parameter)


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
        interaction       = attr(samples[[model_term]], "interaction"),
        interaction_terms = attr(samples[[model_term]], "interaction_terms"),
        term_components   = attr(samples[[model_term]], "term_components"),
        factor_terms      = attr(samples[[model_term]], "factor_terms"),
        factor_contrasts  = attr(samples[[model_term]], "factor_contrasts"),
        factor_design     = attr(samples[[model_term]], "factor_design"),
        treatment         = attr(samples[[model_term]], "treatment"),
        independent       = attr(samples[[model_term]], "independent"),
        orthonormal       = attr(samples[[model_term]], "orthonormal"),
        meandif           = attr(samples[[model_term]], "meandif"),
        ordered           = attr(samples[[model_term]], "ordered")
        ))
      names(priors_info) <- names(prior_list)
      model_terms_type <- sapply(JAGS_model_terms, function(model_term){
        if(priors_info[[model_term]][["factor"]]){
          return("factor")
        }else{
          return("continuous")
        }
      })
      predictors_type  <- model_terms_type[JAGS_parameter_names(parameters = predictors, formula_parameter = formula_parameter)]


      ### prepare at specification
      # in case of an interaction, all levels need to be set
      if(!is.null(priors_info[[parameter]][["interaction"]]) && priors_info[[parameter]][["interaction"]]){
        at_manipulated <- JAGS_parameter_names(priors_info[[parameter]][["interaction_terms"]], formula_parameter = formula_parameter)
      }else{
        at_manipulated <- parameter
      }

      if(!all(names(at) %in% predictors))
        stop(paste0("The following values passed via the 'at' argument do not correspond to the specified model: ", paste0("'", names(at)[!names(at) %in% predictors], "'", collapse = ", ")))
      if(any(format_parameter_names(at_manipulated, formula_parameters = formula_parameter, formula_prefix = FALSE) %in% names(at)))
        stop("Values of the parameter of interested cannot be specified via the 'at' argument.")

      # fill in with default values if needed
      for(i in seq_along(predictors)){
        if(JAGS_predictors[i] %in% at_manipulated){
          # specify levels for the parameter of interest
          if(model_terms_type[[JAGS_predictors[i]]] == "continuous"){
            at[[predictors[i]]] <- c(-1, 0, 1)
          }else{
            at[[predictors[i]]] <- priors_info[[JAGS_predictors[i]]][["level_names"]]
          }
        }else if(is.null(at[[predictors[i]]])){
          # specify levels for the remaining parameters
          if(model_terms_type[[JAGS_predictors[i]]] == "continuous"){
            # fill in zeroes for unspecified continuous predictors
            at[[predictors[i]]] <- 0
          }else if(priors_info[[JAGS_predictors[i]]][["treatment"]]){
            # fill in the default category for unspecified treatment factors
            at[[predictors[i]]] <- priors_info[[JAGS_predictors[i]]][["level_names"]][1]
          }else{
            # fill in NA for any other factor type
            at[[predictors[i]]] <- NA
          }
        }
      }

      # transform to a data.frame
      data <- as.data.frame(expand.grid(at))

      # check the specified data
      if(any(predictors_type == "factor")){

        # check the proper data input for each factor variable
        for(i in seq_along(predictors_type)[predictors_type == "factor"]){

          if(is.factor(data[,predictors[i]])){
            if(all(levels(data[,predictors[i]]) %in% priors_info[[JAGS_predictors[i]]][["level_names"]])){
              # either the formatting is correct, or the supplied levels are a subset of the original levels
              # reformat to check ordering and etc...
              data[,predictors[i]] <- factor(data[,predictors[i]], levels = priors_info[[JAGS_predictors[i]]][["level_names"]])
            }else{
              # there are some additional levels
              stop(paste0("Levels specified in the '", predictors[i], "' factor variable do not match the levels used for model specification."))
            }
          }else if(all(stats::na.omit(unique(data[,predictors[i]])) %in% priors_info[[JAGS_predictors[i]]][["level_names"]])){
            # the variable was not passed as a factor but the values matches the factor levels
            data[,predictors[i]] <- factor(data[,predictors[i]], levels = priors_info[[JAGS_predictors[i]]][["level_names"]])
          }else{
            # there are some additional mismatching values
            stop(paste0("Levels specified in the '", predictors[i], "' factor variable do not match the levels used for model specification."))
          }

          # set the contrast
          if(priors_info[[JAGS_predictors[i]]][["orthonormal"]]){
            stats::contrasts(data[,predictors[i]]) <- "contr.orthonormal"
          }else if(priors_info[[JAGS_predictors[i]]][["meandif"]]){
            stats::contrasts(data[,predictors[i]]) <- "contr.meandif"
          }else if(priors_info[[JAGS_predictors[i]]][["independent"]]){
            stats::contrasts(data[,predictors[i]]) <- "contr.independent"
          }else if(priors_info[[JAGS_predictors[i]]][["ordered"]]){
            factor_contrasts <- unlist(
              priors_info[[JAGS_predictors[i]]][["factor_contrasts"]],
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
          }else if(priors_info[[JAGS_predictors[i]]][["treatment"]]){
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


      # obtain samples information
      models_ind <- do.call(cbind, lapply(c(if(has_intercept) "intercept", model_terms), function(x) attr(samples[[JAGS_parameter_names(x, formula_parameter = formula_parameter)]], "models_ind")))
      sample_ind <- do.call(cbind, lapply(c(if(has_intercept) "intercept", model_terms), function(x) attr(samples[[JAGS_parameter_names(x, formula_parameter = formula_parameter)]], "sample_ind")))
      if(!inherits(samples, "as_mixed_posteriors") && (!all(models_ind[,1] == models_ind) || !all(sample_ind[,1] == sample_ind)))
        stop("the posterior samples are not alligned across models/draws")
      models_ind <- models_ind[,1]


      ### evaluate the design matrix on the samples -> output[data, posterior]
      if(has_intercept){

        terms_indexes    <- attr(model_matrix, "assign") + 1
        terms_indexes[1] <- 0

        # get model/sample indices and check for scaling factors
        temp_multiply_by <- .get_combined_parameter_scaling_factor_matrix(
          JAGS_parameter_names("intercept", formula_parameter = formula_parameter),
          prior_list  = prior_list,
          posterior   = posterior_samples_matrix,
          models_ind  = models_ind,
          nrow        = nrow(data),
          simple_list = inherits(samples, "as_mixed_posteriors")
        )


        marginal_posterior_samples <- temp_multiply_by * matrix(posterior_samples_matrix[,JAGS_parameter_names("intercept", formula_parameter = formula_parameter)],
                                                       nrow = nrow(data), ncol = nrow(posterior_samples_matrix), byrow = TRUE)

      }else{

        terms_indexes <- attr(model_matrix, "assign")
        marginal_posterior_samples <- matrix(0, nrow = nrow(data), ncol = nrow(posterior_samples_matrix))

      }

      # add remaining terms (omitting the intercept indexed as 0)
      for(i in unique(terms_indexes[terms_indexes > 0])){

        # subset the model matrix
        temp_data <- .marginal_posterior_term_data(
          model_matrix  = model_matrix,
          terms_indexes = terms_indexes,
          term_index    = i,
          data          = data,
          prior_info    = priors_info[[JAGS_model_terms[i]]],
          term_name     = JAGS_model_terms[i]
        )

        temp_posterior <- posterior_samples_matrix[,paste0(
          JAGS_model_terms[i],
          if(model_terms_type[i] == "factor" && priors_info[[JAGS_model_terms[i]]][["levels"]] > 1) paste0("[", 1:priors_info[[JAGS_model_terms[i]]][["levels"]], "]"))
          ,drop = FALSE]

        # check for scaling factors
        temp_multiply_by <- .get_combined_parameter_scaling_factor_matrix(
          JAGS_model_terms[i],
          prior_list  = prior_list,
          posterior   = posterior_samples_matrix,
          models_ind  = models_ind,
          nrow        = nrow(data),
          simple_list = inherits(samples, "as_mixed_posteriors")
        )

        marginal_posterior_samples <- marginal_posterior_samples + temp_multiply_by * (temp_data %*% t(temp_posterior))

      }


      # apply transformations
      if(!is.null(transformation)){
        marginal_posterior_samples <- .density.prior_transformation_x(marginal_posterior_samples, transformation, transformation_arguments)
      }


      ### split the output into lists based on specification
      # create indexing and names for the manipulated predictors
      if(length(at_manipulated) == 1 && format_parameter_names(at_manipulated, formula_parameters = formula_parameter, formula_prefix = FALSE) == "intercept"){

        class(marginal_posterior_samples)             <- c(class(marginal_posterior_samples), "marginal_posterior.simple")
        attr(marginal_posterior_samples, "parameter") <- parameter
        attr(marginal_posterior_samples, "level")     <- "intercept"
        attr(marginal_posterior_samples, "data")      <- data

        marginal_posterior_samples <- list("intercept" = marginal_posterior_samples)

        attr(marginal_posterior_samples, "data")        <- data
        attr(marginal_posterior_samples, "level_at")    <- NULL
        attr(marginal_posterior_samples, "level_names") <- "intercept"
        attr(marginal_posterior_samples, "parameter")   <- parameter

      }else{

        manipulated_predictors <- format_parameter_names(at_manipulated, formula_parameters = formula_parameter, formula_prefix = FALSE)
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
        level_names                 <- apply(at_index_output.names_frame, 1, paste0, collapse = ", ")

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
        allow_failure    = !prior_samples,
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

      if(length(at_manipulated) == 1 && format_parameter_names(at_manipulated, formula_parameters = formula_parameter, formula_prefix = FALSE) == "intercept"){

        prior_weights <- linear_weights
        marginal_posterior_samples[["intercept"]] <- .posterior_support_set(
          marginal_posterior_samples[["intercept"]],
          .posterior_support_from_prior_context_weights(
            prior_density_context,
            prior_weights,
            output_transformation           = transformation,
            output_transformation_arguments = transformation_arguments
          )
        )
        intercept_atoms <- .posterior_atoms_formula(
          samples,
          prior_list,
          prior_weights,
          transformation = transformation,
          transformation_arguments = transformation_arguments,
          column_name = "intercept"
        )
        if(!is.null(intercept_atoms)){
          marginal_posterior_samples[["intercept"]] <- .posterior_atoms_set(
            marginal_posterior_samples[["intercept"]],
            intercept_atoms
          )
        }

      }else{

        level_prior_weights <- vector("list", length(level_names))
        names(level_prior_weights) <- level_names
        for(lvl in seq_along(level_names)){
          prior_weights <- linear_weights[data_split[[lvl]], , drop = FALSE]
          level_prior_weights[[level_names[lvl]]] <- prior_weights
          marginal_posterior_samples[[level_names[lvl]]] <- .posterior_support_set(
            marginal_posterior_samples[[level_names[lvl]]],
            .posterior_support_from_prior_context_weights(
              prior_density_context,
              prior_weights,
              output_transformation           = transformation,
              output_transformation_arguments = transformation_arguments
            )
          )
          level_atoms <- .posterior_atoms_formula(
            samples,
            prior_list,
            prior_weights,
            transformation = transformation,
            transformation_arguments = transformation_arguments,
            column_name = level_names[lvl]
          )
          if(!is.null(level_atoms)){
            marginal_posterior_samples[[level_names[lvl]]] <-
              .posterior_atoms_set(
                marginal_posterior_samples[[level_names[lvl]]],
                level_atoms
              )
          }
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

        if(length(at_manipulated) == 1 && format_parameter_names(at_manipulated, formula_parameters = formula_parameter, formula_prefix = FALSE) == "intercept"){

          prior_weights <- linear_weights
          prior_density <- .prior_density_from_context_rows(
            prior_density_context,
            prior_weights,
            output_transformation           = transformation,
            output_transformation_arguments = transformation_arguments
          )
          attr(marginal_posterior_samples[["intercept"]], "linear_weights") <- prior_weights
          attr(marginal_posterior_samples[["intercept"]], "prior_density") <- prior_density
          attr(marginal_posterior_samples[["intercept"]], "prior_density_context") <- prior_density_context

        }else{

          for(lvl in seq_along(level_names)){
            prior_weights <- level_prior_weights[[level_names[lvl]]]
            prior_density <- .prior_density_from_context_rows(
              prior_density_context,
              prior_weights,
              output_transformation           = transformation,
              output_transformation_arguments = transformation_arguments
            )
            attr(marginal_posterior_samples[[level_names[lvl]]], "linear_weights") <- prior_weights
            attr(marginal_posterior_samples[[level_names[lvl]]], "prior_density") <- prior_density
            attr(marginal_posterior_samples[[level_names[lvl]]], "prior_density_context") <- prior_density_context
          }
        }

        attr(marginal_posterior_samples, "prior_density_context") <- prior_density_context
      }

      }

      attr(marginal_posterior_samples, "formula_parameter") <- formula_parameter
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
      posterior_density_conditional <- attr(parameter_samples, "conditional", exact = TRUE)
      posterior_density_conditional_rule <- attr(parameter_samples, "conditional_rule", exact = TRUE)
      posterior_density_condition_key <- attr(parameter_samples, "condition_key", exact = TRUE)

      # transform factor levels
      marginal_posterior_samples <- transform_factor_samples(samples)
      marginal_posterior_samples <- transform_treatment_samples(marginal_posterior_samples)[[parameter]]
      marginal_factor_atoms <- .posterior_atoms_get(marginal_posterior_samples)
      attr(marginal_posterior_samples, "posterior_density") <- NULL
      attr(marginal_posterior_samples, "posterior_ordinate") <- NULL
      marginal_factor_metadata <- marginal_posterior_samples
      marginal_factor_support <- attr(marginal_factor_metadata, "posterior_support", exact = TRUE)
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
          allow_failure = TRUE,
          condition_source = parameter_samples
        )
        factor_weights <- .prior_factor_level_weight_matrix(
          sample_metadata = marginal_factor_metadata,
          parameter       = parameter,
          samples         = samples
        )
      }

      # apply transformations
      if(!is.null(transformation)){
        marginal_posterior_samples <- .density.prior_transformation_x(marginal_posterior_samples, transformation, transformation_arguments)
      }

      # create output object
      marginal_posterior_samples <- lapply(seq_along(level_names), function(lvl_i){
        temp_marginal_posterior_samples <- marginal_posterior_samples[,lvl_i]
        class(temp_marginal_posterior_samples) <- c(class(temp_marginal_posterior_samples), "marginal_posterior.factor")
        attr(temp_marginal_posterior_samples, "parameter")  <- parameter
        attr(temp_marginal_posterior_samples, "level_name") <- level_names[lvl_i]
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
            weights,
            output_transformation           = transformation,
            output_transformation_arguments = transformation_arguments
          )
        }else if(!is.null(temp_support) && !is.null(transformation)){
          temp_support <- .posterior_support_transform(
            temp_support,
            transformation,
            transformation_arguments
          )
        }
        temp_marginal_posterior_samples <- .posterior_support_set(
          temp_marginal_posterior_samples,
          temp_support
        )
        if(!is.null(marginal_factor_atoms)){
          temp_atoms <- .posterior_atoms_for_column(
            marginal_factor_atoms,
            lvl_i
          )
          if(!is.null(transformation)){
            temp_atoms <- .posterior_atoms_transform(
              temp_atoms,
              transformation,
              transformation_arguments
            )
          }
          temp_marginal_posterior_samples <- .posterior_atoms_set(
            temp_marginal_posterior_samples,
            temp_atoms
          )
        }
        if(is.null(transformation)){
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
            attr(temp_marginal_posterior_samples, "posterior_density") <- posterior_density
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
            attr(temp_marginal_posterior_samples, "posterior_ordinate") <- posterior_ordinate
          }
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
          allow_failure = TRUE,
          condition_source = samples[[parameter]]
        )
        if(!is.null(prior_density_context)){
          weights <- rep(0, length(prior_density_context$column_names))
          names(weights) <- prior_density_context$column_names
          if(parameter %in% names(weights)){
            weights[[parameter]] <- 1
            marginal_support <- .posterior_support_from_prior_context_weights(
              prior_density_context,
              weights
            )
          }
        }
      }

      # apply transformations
      if(!is.null(transformation)){
        marginal_atoms <- .posterior_atoms_get(marginal_posterior_samples)
        marginal_posterior_samples <- .density.prior_transformation_x(marginal_posterior_samples, transformation, transformation_arguments)
        attr(marginal_posterior_samples, "posterior_density") <- NULL
        attr(marginal_posterior_samples, "posterior_ordinate") <- NULL
        if(!is.null(marginal_atoms)){
          marginal_posterior_samples <- .posterior_atoms_set(
            marginal_posterior_samples,
            .posterior_atoms_transform(
              marginal_atoms,
              transformation,
              transformation_arguments
            )
          )
        }
        marginal_support <- .posterior_support_transform(
          marginal_support,
          transformation,
          transformation_arguments
        )
      }else{
        marginal_posterior_samples <- .posterior_density_attach(
          samples          = marginal_posterior_samples,
          sources          = .posterior_density_sources(samples[[parameter]]),
          parameter        = parameter,
          conditional      = attr(samples[[parameter]], "conditional", exact = TRUE),
          conditional_rule = attr(samples[[parameter]], "conditional_rule", exact = TRUE),
          condition_key    = attr(samples[[parameter]], "condition_key", exact = TRUE),
          allow_unlabeled  = TRUE
        )
        marginal_posterior_samples <- .posterior_ordinate_attach(
          samples          = marginal_posterior_samples,
          sources          = .posterior_ordinate_sources(samples[[parameter]]),
          parameter        = parameter,
          conditional      = attr(samples[[parameter]], "conditional", exact = TRUE),
          conditional_rule = attr(samples[[parameter]], "conditional_rule", exact = TRUE),
          condition_key    = attr(samples[[parameter]], "condition_key", exact = TRUE),
          allow_unlabeled  = TRUE
        )
        marginal_posterior_samples <- .posterior_density_attach(
          samples          = marginal_posterior_samples,
          sources          = .posterior_density_sources(samples),
          parameter        = parameter,
          conditional      = attr(samples[[parameter]], "conditional", exact = TRUE),
          conditional_rule = attr(samples[[parameter]], "conditional_rule", exact = TRUE),
          condition_key    = attr(samples[[parameter]], "condition_key", exact = TRUE),
          allow_unlabeled  = FALSE
        )
        marginal_posterior_samples <- .posterior_ordinate_attach(
          samples          = marginal_posterior_samples,
          sources          = .posterior_ordinate_sources(samples),
          parameter        = parameter,
          conditional      = attr(samples[[parameter]], "conditional", exact = TRUE),
          conditional_rule = attr(samples[[parameter]], "conditional_rule", exact = TRUE),
          condition_key    = attr(samples[[parameter]], "condition_key", exact = TRUE),
          allow_unlabeled  = FALSE
        )
      }

      marginal_posterior_samples <- .posterior_support_set(
        marginal_posterior_samples,
        marginal_support
      )
      class(marginal_posterior_samples) <- c(class(marginal_posterior_samples), "marginal_posterior.simple")

    }


    # add prior densities
    if(prior_samples){

      if(is.null(prior_density_context)){
        prior_density_context <- .marginal_posterior_prior_density_context(
          samples    = samples,
          prior_list = prior_list,
          n_samples  = n_samples,
          condition_source = samples[[parameter]]
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
            weights,
            output_transformation           = transformation,
            output_transformation_arguments = transformation_arguments
          )
          attr(marginal_posterior_samples[[level_names[lvl_i]]], "linear_weights") <- weights
          attr(marginal_posterior_samples[[level_names[lvl_i]]], "prior_density") <- prior_density
          attr(marginal_posterior_samples[[level_names[lvl_i]]], "prior_density_context") <- prior_density_context
        }

      }else if(inherits(samples[[parameter]], "mixed_posteriors.simple")){

        weights <- rep(0, length(prior_density_context$column_names))
        names(weights) <- prior_density_context$column_names
        weights[[parameter]] <- 1

        prior_density <- .prior_density_from_context(
          prior_density_context,
          weights,
          output_transformation           = transformation,
          output_transformation_arguments = transformation_arguments
        )
        attr(marginal_posterior_samples, "linear_weights") <- weights
        attr(marginal_posterior_samples, "prior_density") <- prior_density
        attr(marginal_posterior_samples, "prior_density_context") <- prior_density_context
      }

      attr(marginal_posterior_samples, "prior_density_context") <- prior_density_context
    }
  }

  if(is.null(transformation)){
    marginal_posterior_samples <- .marginal_posterior_attach_precomputed_metadata(
      marginal              = marginal_posterior_samples,
      samples               = samples,
      parameter             = parameter,
      condition_source      = samples[[parameter]]
    )
  }

  marginal_posterior_samples <- .marginal_posterior_set_condition_attributes(
    marginal_posterior_samples,
    samples,
    condition_source = samples[[parameter]]
  )
  class(marginal_posterior_samples) <- c(class(marginal_posterior_samples), "marginal_posterior")
  return(marginal_posterior_samples)
}

.marginal_posterior_prior_density_context <- function(samples, prior_list,
                                                       column_names = NULL,
                                                       n_samples = 10000,
                                                       allow_failure = FALSE,
                                                       condition_source = NULL){

  condition_metadata <- .marginal_posterior_condition_metadata(
    samples,
    condition_source = condition_source
  )
  prior_density_context <- attr(samples, "prior_density_context")
  if(!is.null(prior_density_context) &&
     .marginal_posterior_context_matches_condition(prior_density_context, condition_metadata)){
    return(prior_density_context)
  }

  tryCatch(
    {
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

      .prior_density_build_context(
        prior_list       = prior_list,
        column_names     = column_names,
        n_grid           = max(16L, n_samples),
        conditional      = condition_metadata[["conditional"]],
        conditional_rule = condition_metadata[["conditional_rule"]],
        condition_event  = condition_metadata[["condition_event"]]
      )
    },
    error = function(e){
      if(isTRUE(allow_failure)){
        return(NULL)
      }
      stop(conditionMessage(e), call. = FALSE)
    }
  )
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
    condition_event <- attr(source, "resolved_condition_event", exact = TRUE)
    if(is.null(condition_event)){
      condition_event <- attr(source, "condition_event", exact = TRUE)
    }

    .marginal_posterior_normalize_condition_metadata(
      conditional      = attr(source, "conditional", exact = TRUE),
      conditional_rule = attr(source, "conditional_rule", exact = TRUE),
      condition_key    = attr(source, "condition_key", exact = TRUE),
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

  if(isTRUE(attr(samples, "transform_scaled", exact = TRUE)) ||
     isTRUE(attr(samples[[parameter]], "transform_scaled", exact = TRUE))){
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
         is.null(attr(marginal[[i]], "posterior_density", exact = TRUE))){
        attr(marginal[[i]], "posterior_density") <- density
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
         is.null(attr(marginal[[i]], "posterior_ordinate", exact = TRUE))){
        attr(marginal[[i]], "posterior_ordinate") <- ordinate
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
    attr(x, "conditional")      <- conditional
    attr(x, "conditional_rule") <- conditional_rule
    attr(x, "condition_key")    <- condition_key
    if(!is.null(condition_event)){
      attr(x, "condition_event")          <- condition_event
      attr(x, "resolved_condition_event") <- condition_event
    }
    x
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
    out <- lapply(parameter_samples, attr, which = "prior_density")
  }else{
    out <- list(attr(parameter_samples, "prior_density"))
    names(out) <- parameter
  }

  out
}

.marginal_posterior_parameter_posterior_densities <- function(samples, parameter){

  parameter_samples <- samples[[parameter]]

  if(is.list(parameter_samples)){
    out <- .posterior_density_child_attributes(parameter_samples)
  }else{
    out <- list(attr(parameter_samples, "posterior_density"))
    names(out) <- parameter
  }

  out
}

.get_combined_parameter_scaling_factor_matrix <- function(term, prior_list, posterior, models_ind, nrow, simple_list = FALSE){

  if(simple_list){
    temp_multiply_by <- .get_parameter_scaling_factor_matrix(term, prior_list, posterior, nrow = nrow, ncol = nrow(posterior))
  }else{
    model_samples <- table(models_ind)

    temp_multiply_by <- do.call(cbind, lapply(unique(models_ind), function(m){
      temp_prior_list <- lapply(prior_list, function(parameter_priors) parameter_priors[[m]])
      temp_posterior  <- posterior[models_ind == m,,drop=FALSE]
      return(.get_parameter_scaling_factor_matrix(term, temp_prior_list, temp_posterior, nrow = nrow, ncol = sum(models_ind == m)))
    }))
  }

  return(temp_multiply_by)
}
