#' @title Create BayesTools model tables
#'
#' @description Creates model summary based on a model objects or
#' provides estimates table for a runjags fit.
#'
#' @param model model object containing a list of \code{priors}
#' and \code{inference} object, The \code{inference} must be a
#' named list with information about the model: model number
#' \code{m_number}, marginal likelihood \code{marglik}, prior and
#' posterior probability \code{prior_prob} and \code{post_prob},
#' and model inclusion Bayes factor \code{inclusion_BF}
#' @param fit runjags model fit
#' @param conditional summarizes estimates conditional on being included
#' in the model for spike and slab priors. Defaults to \code{FALSE}.
#' @param transformations named list of transformations to be applied
#' to specific parameters
#' @param model_description named list with additional description
#' to be added to the table
#' @param remove_inclusion whether estimates of the inclusion probabilities
#' should be excluded from the summary table. Defaults to \code{FALSE}.
#' @param remove_parameters parameters to be removed from the summary.
#' Can be \code{NULL} (default, no removal), a character vector of parameter
#' names to remove, or \code{TRUE} to remove all parameters that are not
#' part of any formula. For formula random effects, character filters also
#' accept semantic aliases \code{"random"}, \code{"random_sd"},
#' \code{"random_rho"}, \code{"random_correlation"},
#' \code{"random_variance_fraction"}, \code{"random_variance_ratio"},
#' \code{"random_allocation"}, and \code{"random_sd_multiplier"}.
#' @param remove_formulas character vector of formula names whose parameters
#' should be removed from the summary. Defaults to \code{NULL}.
#' @param keep_parameters character vector of parameter names to keep.
#' All other parameters will be removed unless they belong to formulas
#' specified in \code{keep_formulas}. The random-effect aliases listed for
#' \code{remove_parameters} can also be used here, for example
#' \code{keep_parameters = c("random_sd", "random_rho")}.
#' @param keep_formulas character vector of formula names whose parameters
#' should be kept. All other parameters will be removed unless they are
#' specified in \code{keep_parameters}. Defaults to \code{NULL}.
#' @param random_effects_summary random-effect reporting mode for JAGS estimates
#' tables. \code{"standard"} replaces raw random-effect implementation
#' parameters with semantic SD, rho/correlation, true variance-fraction, and
#' mean-variance SD-component variance-ratio summaries. \code{"full"} also
#' includes heterogeneous SD multipliers.
#' \code{"raw"} keeps the historical raw monitored parameters, and
#' \code{"none"} removes random-effect parameters from the table. When used
#' together with \code{transform_scaled = TRUE}, SD and correlation summaries
#' for formula random effects are computed after applying the original-scale
#' formula transformation.
#' @param random_effects_metadata whether to add random-effect metadata columns
#' to JAGS estimates tables. When \code{TRUE}, the table includes the
#' user-facing random-effect name, grouping label, and covariance structure
#' type for random-effect rows. Defaults to \code{FALSE}.
#' @param remove_random_effects,keep_random_effects character vectors of
#' random-effect names/blocks used to remove or keep random-effect rows.
#' These match explicit \code{name = } values, generated block names, and
#' grouping labels. Non-random parameters are left unaffected; combine with
#' \code{keep_parameters = "random"} to return only the selected random-effect
#' rows.
#' @param remove_random_structures,keep_random_structures character vectors of
#' random-effect covariance structures (for example \code{"diag"}, \code{"us"},
#' \code{"ar1"}, or \code{"hcs"}) used to remove or keep random-effect rows.
#' Non-random parameters are left unaffected.
#' @param return_samples whether to return the transoformed and formated samples
#' instead of the table. Defaults to \code{FALSE}.
#' @param remove_diagnostics whether to exclude MCMC diagnostics (MCMC error,
#' ESS, R-hat) from the output table. Defaults to \code{FALSE}. Setting to
#' \code{TRUE} will exclude diagnostics columns regardless of the
#' \code{conditional} setting.
#' @param diagnostic_columns MCMC diagnostic columns to display in JAGS
#' estimates tables. Can be \code{"all"}, \code{"none"}, \code{TRUE},
#' \code{FALSE}, or a character vector containing any subset of
#' \code{"MCMC_error"}, \code{"MCMC_SD_error"}, \code{"ESS"}, and
#' \code{"R_hat"}. Defaults to the
#' \code{BayesTools.JAGS_estimates_diagnostic_columns} option, or all
#' diagnostics unless \code{remove_diagnostics = TRUE}.
#' @param BF_diagnostics whether to add MCMC diagnostics for Bayes factors
#' computed from model indicator frequencies. The Bayes factor error is
#' reported as a relative Monte Carlo standard error percentage. Defaults to
#' \code{FALSE}.
#' @param BF_diagnostic_columns MCMC diagnostic columns to display in JAGS
#' inclusion Bayes factor tables. Can be \code{"all"}, \code{"none"},
#' \code{TRUE}, \code{FALSE}, or a character vector containing any subset of
#' \code{"ESS"}, \code{"MCMC_error"}, and \code{"BF_error_percent"}.
#' Defaults to the \code{BayesTools.JAGS_BF_diagnostic_columns} option, or all
#' diagnostics when \code{BF_diagnostics = TRUE} and none otherwise.
#' @inheritParams BayesTools_ensemble_tables
#'
#'
#' @return \code{model_summary_table} returns a table with
#' overview of the fitted model, \code{runjags_estimates_table} returns
#' a table with MCMC estimates, and \code{runjags_estimates_empty_table}
#' returns an empty estimates table. All of the tables are objects of
#' class 'BayesTools_table'.
#'
#' @details For product-space JAGS inclusion Bayes factors, posterior
#' inclusion probabilities of exactly 0 or 1 cannot produce a finite
#' point estimate of the model odds ratio. In that case, inclusion BFs
#' marked with \code{"<"} or \code{">"} are finite-sample bounds: posterior
#' inclusion probabilities of 0 or 1 were replaced by \code{1/S} or
#' \code{(S - 1)/S}, where \code{S} is the number of posterior samples.
#' This is a reporting convention, not an unbiased finite Bayes-factor
#' estimate. If the prior inclusion probability is exactly 0 or 1, the
#' inclusion Bayes factor is undefined and reported as \code{NA}, because
#' the corresponding inclusion/exclusion comparison was not tested.
#'
#' @export JAGS_summary_table
#' @export JAGS_estimates_table
#' @export JAGS_inference_table
#' @export model_summary_table
#' @export runjags_estimates_table
#' @export runjags_inference_table
#' @export model_summary_empty_table
#' @export JAGS_inference_empty_table
#' @export JAGS_estimates_empty_table
#' @export runjags_estimates_empty_table
#' @export runjags_inference_empty_table
#' @export stan_estimates_table
#' @name BayesTools_model_tables
#'
#' @seealso [BayesTools_ensemble_tables]
NULL

#' @rdname BayesTools_model_tables
model_summary_table <- function(model, model_description = NULL, title = NULL, footnotes = NULL, warnings = NULL,
                                remove_spike_0 = TRUE, short_name = FALSE, formula_prefix = TRUE, remove_parameters = NULL){

  # check input
  check_list(model, "model", check_names = "inference", allow_other = TRUE, all_objects = TRUE)
  prior_list <- attr(model[["fit"]], "prior_list")
  check_list(prior_list, "model:priors")
  if(!all(sapply(prior_list, is.prior)))
    stop("'model:priors' must be a list of priors.")
  model_inference <- model[["inference"]]
  check_list(model_inference, "model:inference", check_names = c("m_number", "marglik", "prior_prob", "post_prob", "inclusion_BF"), allow_other = TRUE, all_objects = TRUE)
  check_int(model_inference[["m_number"]],      "model_inference:model_number")
  check_real(model_inference[["marglik"]],      "model_inference:marglik")
  check_real(model_inference[["prior_prob"]],   "model_inference:prior_prob",   lower = 0, upper = 1)
  check_real(model_inference[["post_prob"]],    "model_inference:post_prob",   lower = 0, upper = 1)
  check_real(model_inference[["inclusion_BF"]], "model_inference:inclusion_BF", lower = 0)
  check_list(model_description, "model_description", allow_NULL = TRUE)
  check_bool(remove_spike_0, "remove_spike_0", allow_NA = FALSE)
  check_bool(short_name, "short_name", allow_NA = FALSE)
  check_char(title, "title", allow_NULL = TRUE)
  check_char(footnotes, "footnotes", check_length = 0, allow_NULL = TRUE)
  check_char(warnings, "warnings", check_length = 0, allow_NULL = TRUE)
  check_bool(formula_prefix, "formula_prefix")
  check_char(remove_parameters, "remove_parameters", allow_NULL = TRUE, check_length = 0)

  # prepare the columns
  summary_names  <- c(
    "Model",
    if(!is.null(model_description)) names(model_description),
    "Prior prob.",
    "log(marglik)",
    "Post. prob.",
    "Inclusion BF")
  summary_values <- c(
    model_inference[["m_number"]],
    if(!is.null(model_description)) unlist(model_description),
    .format_column(model_inference[["prior_prob"]],   "probability"),
    .format_column(model_inference[["marglik"]],      "marglik"),
    .format_column(model_inference[["post_prob"]],    "probability"),
    .format_column(model_inference[["inclusion_BF"]], "BF"))

  summary_priors  <- "Parameter prior distributions"
  for(i in seq_along(prior_list)){
    # get the prior name
    if(is.prior.none(prior_list[[i]])){
      next
    }else if(remove_spike_0 && is.prior.point(prior_list[[i]]) && prior_list[[i]][["parameters"]][["location"]] == 0 || (names(prior_list)[[i]] %in% remove_parameters)){
      next
    }else if(is.prior.weightfunction(prior_list[[i]]) | is.prior.PET(prior_list[[i]]) | is.prior.PEESE(prior_list[[i]]) |
             is_prior_phacking(prior_list[[i]]) | is_prior_bias(prior_list[[i]])){
      temp_prior <- print(prior_list[[i]], silent = TRUE, short_name = short_name)
    }else if(is.prior.simple(prior_list[[i]]) | is.prior.vector(prior_list[[i]]) | is.prior.factor(prior_list[[i]]) | is.prior.spike_and_slab(prior_list[[i]]) | is.prior.mixture(prior_list[[i]])){
      temp_prior <- paste0(names(prior_list)[i], " ~ " , print(prior_list[[i]], silent = TRUE, short_name = short_name, inline = TRUE))
    }else if(is.prior.point(prior_list[[i]])){
      temp_prior <- paste0(names(prior_list)[i], " = " , print(prior_list[[i]], silent = TRUE, short_name = short_name))
    }
    # change the formula formatting
    if(!is.null(attr(prior_list[[i]], "parameter"))){
      temp_prior <- gsub(
        paste0(attr(prior_list[[i]], "parameter"), "_"),
        if(formula_prefix) paste0("(", attr(prior_list[[i]], "parameter"), ") ") else "",
        temp_prior)
      temp_prior <- gsub("__xXx__", ":", temp_prior)
    }
    summary_priors <- c(summary_priors, temp_prior)
  }

  if(length(summary_names) > length(summary_priors)){
    summary_priors <- c(summary_priors, rep("", length(summary_names) - length(summary_priors)))
  }else if(length(summary_names) < length(summary_priors)){
    summary_names  <- c(summary_names,  rep("", length(summary_priors) - length(summary_names)))
    summary_values <- c(summary_values, rep("", length(summary_priors) - length(summary_values)))
  }
  summary_names <- paste0(summary_names, "  ")

  summary_table <- data.frame(cbind(
    summary_names,
    summary_values,
    rep("           ", length(summary_names)),
    summary_priors
  ))
  names(summary_table) <- NULL

  # prepare output
  class(summary_table)             <- c("BayesTools_table", class(summary_table))
  attr(summary_table, "type")      <- c("string_left", "string", "string", "prior")
  attr(summary_table, "rownames")  <- FALSE
  attr(summary_table, "as.matrix") <- TRUE
  attr(summary_table, "title")     <- title
  attr(summary_table, "footnotes") <- footnotes
  attr(summary_table, "warnings")  <- warnings

  return(summary_table)
}

#' @rdname BayesTools_model_tables
runjags_estimates_table  <- function(fit, transformations = NULL, title = NULL, footnotes = NULL, warnings = NULL, conditional = FALSE,
                                     probs = c(0.025, 0.5, 0.975), remove_spike_0 = TRUE, transform_factors = FALSE, transform_orthonormal = FALSE,
                                     formula_prefix = TRUE, remove_inclusion = FALSE, remove_parameters = NULL, remove_formulas = NULL,
                                     keep_parameters = NULL, keep_formulas = NULL, return_samples = FALSE, transform_scaled = FALSE,
                                     random_effects_summary = c("standard", "full", "raw", "none"),
                                     random_effects_metadata = FALSE,
                                     remove_random_effects = NULL, keep_random_effects = NULL,
                                     remove_random_structures = NULL, keep_random_structures = NULL,
                                     remove_diagnostics = FALSE,
                                     diagnostic_columns = getOption("BayesTools.JAGS_estimates_diagnostic_columns", if(remove_diagnostics) "none" else "all")){

  .check_runjags()
  # most of the code is shared with .diagnostics_plot_data function (keep them in sync on update)

  # check fits
  if(!inherits(fit, "runjags"))
    stop("'fit' must be a runjags fit")
  if(!inherits(fit, "BayesTools_fit"))
    stop("'fit' must be a BayesTools fit")
  parameter_registry <- JAGS_parameter_registry(fit)
  prior_list <- attr(fit, "prior_list")
  check_list(prior_list, "prior_list")
  if(!all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_list(transformations, "transformations", allow_NULL = TRUE)
  if(!is.null(transformations) && any(!sapply(transformations, function(trans)is.function(trans[["fun"]]))))
    stop("'transformations' must be list of functions in the 'fun' element.")
  check_char(title, "title", allow_NULL = TRUE)
  check_char(footnotes, "footnotes", check_length = 0, allow_NULL = TRUE)
  check_char(warnings, "warnings", check_length = 0, allow_NULL = TRUE)
  check_real(probs, "probs", lower = 0, upper = 1, check_length = 0)
  check_bool(remove_spike_0, "remove_spike_0", allow_NA = FALSE)
  check_bool(conditional, "conditional", allow_NA = FALSE)
  check_bool(transform_factors, "transform_factors")
  check_bool(transform_orthonormal, "transform_orthonormal")
  check_bool(formula_prefix, "formula_prefix")
  check_bool(transform_scaled, "transform_scaled")
  random_effects_summary <- match.arg(random_effects_summary)
  check_bool(random_effects_metadata, "random_effects_metadata")
  check_bool(remove_diagnostics, "remove_diagnostics")
  diagnostic_columns <- .normalize_diagnostic_columns(diagnostic_columns, .JAGS_estimates_diagnostic_columns(), "diagnostic_columns")
  if(remove_diagnostics){
    diagnostic_columns <- character()
  }
  summary_diagnostic_columns <- if(conditional) character() else diagnostic_columns
  if(!is.null(remove_parameters) && !is.logical(remove_parameters))
    check_char(remove_parameters, "remove_parameters", allow_NULL = TRUE, check_length = 0)
  if(is.logical(remove_parameters))
    check_bool(remove_parameters, "remove_parameters")
  check_char(remove_formulas, "remove_formulas", allow_NULL = TRUE, check_length = 0)
  check_char(keep_parameters, "keep_parameters", allow_NULL = TRUE, check_length = 0)
  check_char(keep_formulas, "keep_formulas", allow_NULL = TRUE, check_length = 0)
  check_char(remove_random_effects, "remove_random_effects", allow_NULL = TRUE, check_length = 0)
  check_char(keep_random_effects, "keep_random_effects", allow_NULL = TRUE, check_length = 0)
  check_char(remove_random_structures, "remove_random_structures", allow_NULL = TRUE, check_length = 0)
  check_char(keep_random_structures, "keep_random_structures", allow_NULL = TRUE, check_length = 0)

  # depreciate
  transform_factors <- .depreciate.transform_orthonormal(transform_orthonormal, transform_factors)

  # get model samples
  model_samples <- .extract_posterior_samples(fit, as_list = FALSE)

  formula_scale <- NULL
  if(transform_scaled){
    formula_scale <- attr(fit, "formula_scale")
  }

  random_summary <- .bt_random_effect_summary_samples(
    model_samples = model_samples,
    prior_list = prior_list,
    formula_design = attr(fit, "formula_design"),
    parameter_registry = parameter_registry,
    mode = random_effects_summary,
    formula_scale = if(transform_scaled) formula_scale else NULL
  )
  model_samples <- random_summary$model_samples
  prior_list <- random_summary$prior_list

  # Transform scaled coefficients after deriving random-effect summaries. This
  # lets summaries reconstruct point-prior/allocation SDs before covariance
  # transformations, while fixed effects still retain all interaction columns.
  if(transform_scaled && !is.null(formula_scale) && length(formula_scale) > 0){
    model_samples <- transform_scale_samples(model_samples, formula_scale)
  }

  model_samples <- .materialize_missing_point_prior_samples(model_samples, prior_list)
  if(remove_inclusion){
    random_inclusion <- vapply(
      prior_list,
      function(prior) identical(attr(prior, "random_summary", exact = TRUE), "inclusion"),
      logical(1)
    )
    if(any(random_inclusion)){
      random_inclusion_names <- names(prior_list)[random_inclusion]
      model_samples <- model_samples[, !colnames(model_samples) %in% random_inclusion_names, drop = FALSE]
      prior_list <- prior_list[!random_inclusion]
    }
  }

  ### remove un-wanted estimates (or support values) - spike and slab priors already dealt with later (also remove the item from prior list)
  # compute filtered parameters using the helper function
  remove_params_vec <- .filter_parameters(
    prior_list        = prior_list,
    remove_parameters = remove_parameters,
    remove_formulas   = remove_formulas,
    keep_parameters   = keep_parameters,
    keep_formulas     = keep_formulas,
    remove_random_effects = remove_random_effects,
    keep_random_effects = keep_random_effects,
    remove_random_structures = remove_random_structures,
    keep_random_structures = keep_random_structures,
    remove_spike_0    = remove_spike_0
  )

  cleaned       <- .remove_auxiliary_parameters(model_samples, prior_list, remove_params_vec)
  model_samples <- cleaned$model_samples
  prior_list    <- cleaned$prior_list
  model_samples <- .bt_JAGS_estimates_filter_raw_random_columns(
    model_samples = model_samples,
    prior_list = prior_list,
    parameter_registry = parameter_registry,
    remove_parameters = remove_parameters,
    remove_formulas = remove_formulas,
    keep_parameters = keep_parameters,
    keep_formulas = keep_formulas,
    remove_random_effects = remove_random_effects,
    keep_random_effects = keep_random_effects,
    remove_random_structures = remove_random_structures,
    keep_random_structures = keep_random_structures
  )

  # simplify mixture and spike and slab priors to simple priors
  # the samples and summary can be dealt with as any other prior (i.e., transformations later)
  for(par in names(prior_list)){
    if(is.prior.spike_and_slab(prior_list[[par]])){

      # process spike and slab using helper function
      processed     <- .process_spike_and_slab(model_samples, prior_list, par, conditional, remove_inclusion, warnings)
      model_samples <- processed$model_samples
      prior_list    <- processed$prior_list
      warnings      <- processed$warnings

    }else if(is.prior.mixture(prior_list[[par]])){

      # check for publication bias component
      is_bias_mixture <- inherits(prior_list[[par]], "prior.bias_mixture")
      is_PET          <- sapply(prior_list[[par]], is.prior.PET)
      is_PEESE        <- sapply(prior_list[[par]], is.prior.PEESE)
      if(is_bias_mixture){
        branch_info   <- .selection_prior_branch_info(prior_list[[par]])
        has_selection <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
        has_phacking  <- vapply(branch_info, function(x) !is.null(x$phacking), logical(1))
      }else{
        branch_info   <- vector("list", length(prior_list[[par]]))
        has_selection <- rep(FALSE, length(prior_list[[par]]))
        has_phacking  <- rep(FALSE, length(prior_list[[par]]))
      }

      # distinguish between null/alternative and component type notations
      components <- attr(prior_list[[par]], "components")

      if(any(is_PET | is_PEESE | has_selection | has_phacking)){

        # change the samples between conditional/averaged based on the preferences
        if(conditional){

          if(any(is_PET)){
            # compute the number of conditional samples
            n_conditional_samples <- sum(model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% which(is_PET))

            # replace null samples with NAs (important for later transformations)
            model_samples[!model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% which(is_PET), "PET"] <- NA

            # add warnings about conditional summary
            warnings <- c(warnings, .runjags_conditional_warning("PET", n_conditional_samples))
          }

          if(any(is_PEESE)){
            # compute the number of conditional samples
            n_conditional_samples <- sum(model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% which(is_PEESE))

            # replace null samples with NAs (important for later transformations)
            model_samples[!model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% which(is_PEESE), "PEESE"] <- NA

            # add warnings about conditional summary
            warnings <- c(warnings, .runjags_conditional_warning("PEESE", n_conditional_samples))
          }

          if(any(has_selection)){
            # compute the number of conditional samples
            n_conditional_samples <- sum(model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% which(has_selection))

            # replace null samples with NAs (important for later transformations)
            model_samples[!model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% which(has_selection), grepl("omega", colnames(model_samples))] <- NA

            # add warnings about conditional summary
            warnings <- c(warnings, .runjags_conditional_warning("omega", n_conditional_samples))
          }

          if(any(has_phacking)){
            # compute the number of conditional samples
            n_conditional_samples <- sum(model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% which(has_phacking))

            # replace non-p-hacking samples with NAs
            phacking_report_columns <- .selection_phacking_report_parameters(lapply(branch_info[has_phacking], function(x) x$phacking))
            phacking_columns <- intersect(c(phacking_report_columns, "phack_kind"), colnames(model_samples))
            if(length(phacking_columns) > 0){
              model_samples[!model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% which(has_phacking), phacking_columns] <- NA
            }

            # add warnings about conditional summary
            warning_columns <- intersect(phacking_report_columns, phacking_columns)
            if(length(warning_columns) > 0){
              warnings <- c(warnings, .runjags_conditional_warning(warning_columns, n_conditional_samples))
            }
          }

        }

        # re-format the weightfunctions
        if(any(has_selection) || any(has_phacking)){

          # rename
          omega_cuts      <- if(any(has_selection)){
            selection_priors <- lapply(branch_info[has_selection], function(x) x$selection)
            weightfunctions_mapping(selection_priors, cuts_only = TRUE, one_sided = TRUE)
          }else{
            c(0, 1)
          }
          omega_names_old <- paste0("omega[", 1:(length(omega_cuts)-1),"]")
          omega_names     <- sapply(1:(length(omega_cuts)-1), function(i)paste0("omega[",omega_cuts[i],",",omega_cuts[i+1],"]"))
          colnames(model_samples)[which(colnames(model_samples) %in% omega_names_old)] <- omega_names

          # remove if requested
          if("omega" %in% remove_parameters){
            model_samples <- model_samples[,!colnames(model_samples) %in% omega_names,drop=FALSE]
            if(any(has_selection)){
              prior_list[[par]][has_selection] <- NULL
            }
          }
        }

        # add the simpler priors to the prior list
        if(any(is_PET)){
          prior_list[["PET"]] <- prior_list[[par]][is_PET][1]
        }
        if(any(is_PEESE)){
          prior_list[["PEESE"]] <- prior_list[[par]][is_PEESE][1]
        }
        if(any(has_selection)){
          prior_list[["omega"]] <- branch_info[has_selection][[1]]$selection
        }
        if(any(has_phacking)){
          phacking_priors <- lapply(branch_info[has_phacking], function(x) x$phacking)
          drop_phacking <- .selection_phacking_unreported_parameters(phacking_priors)
          model_samples <- model_samples[, !colnames(model_samples) %in% drop_phacking, drop = FALSE]

          phacking_report_columns <- .selection_phacking_report_parameters(phacking_priors)
          for(phacking_report_column in phacking_report_columns){
            prior_list[[phacking_report_column]] <- phacking_priors[[1]]$alpha
          }
        }

      }else{

        # prepare parameter names
        if(inherits(prior_list[[par]], "prior.factor_mixture")){
          par_names <- .JAGS_prior_factor_names(par, prior_list[[par]])
        }else{
          par_names <- par
        }

        # change the samples between conditional/averaged based on the preferences
        if(conditional){

          if(all(components %in% c("null", "alternative"))){

            # select the corresponding indicators
            this_component_indicator <- which(components == "alternative")

            # compute the number of conditional samples
            n_conditional_samples <- sum(model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% this_component_indicator)

            # replace null samples with NAs (important for later transformations)
            model_samples[!model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% this_component_indicator, par_names] <- NA

            # add warnings about conditional summary
            warnings <- c(warnings, .runjags_conditional_warning(par_names, n_conditional_samples))

          }else{

            # remove the join samples and replace with individual conditional samples
            temp_position    <- min(which(colnames(model_samples) %in% par))
            temp_all_samples <- model_samples[, colnames(model_samples) %in% par,drop=FALSE]
            temp_new_samples <- list()
            model_samples    <- model_samples[,!colnames(model_samples) %in% par,drop=FALSE]

             # component-by-component replacement
            for(component in unique(components[components != "null"])){

              # create component specific samples
              temp_par_names <- paste0(par_names, "[", component, "]")
              temp_new_samples[[component]]           <- temp_all_samples
              colnames(temp_new_samples[[component]]) <- temp_par_names

              # select the corresponding indicators
              this_component_indicator <- which(components == component)

              # compute the number of conditional samples
              n_conditional_samples <- sum(model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% this_component_indicator)

              # replace null samples with NAs (important for later transformations)
              temp_new_samples[[component]][!model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% this_component_indicator,] <- NA

              # add warnings about conditional summary
              warnings <- c(warnings, .runjags_conditional_warning(temp_par_names, n_conditional_samples))

              # forward transformations to the conditional estimates
              if(!is.null(transformations[[par]])){
                transformations[[temp_par_names]] <- transformations[[par]]
                attr(prior_list[[par]][which(components == component)][1], "parameter") <- attr(prior_list[[par]], "parameter")
                prior_list[[temp_par_names]]      <- prior_list[[par]][which(components == component)][1]
              }
            }

            # place the transformed samples back
            model_samples <- cbind(
              if(temp_position > 1) model_samples[,1:(temp_position-1),drop=FALSE],
              do.call(cbind, temp_new_samples),
              if(temp_position <= ncol(model_samples)) model_samples[,temp_position:ncol(model_samples),drop=FALSE]
            )

            # remove the original parameter transformations
            if(!is.null(transformations[[par]])){
              transformations[[par]] <- NULL
              prior_list[[par]]      <- NULL
            }
          }
        }
      }

      # remove/rename the inclusions probabilities
      if(remove_inclusion){
        model_samples   <- model_samples[,colnames(model_samples) != paste0(par, "_indicator"),drop=FALSE]
      }else{
        if(all(components %in% c("null", "alternative"))){
          # replace and rename in the samples
          model_samples[,colnames(model_samples) == paste0(par, "_indicator")] <- ifelse(
            model_samples[,colnames(model_samples) == paste0(par, "_indicator")] %in% which(components == "alternative"), 1, 0)
          colnames(model_samples)[colnames(model_samples) == paste0(par, "_indicator")] <- paste0(par, " (inclusion)")
        }else{
          # extract
          temp_position <- min(which(colnames(model_samples) %in% paste0(par, "_indicator")))
          temp_samples  <- model_samples[,colnames(model_samples) == paste0(par, "_indicator")]
          model_samples <- model_samples[,colnames(model_samples) != paste0(par, "_indicator"),drop=FALSE]

          # compute component specific indicators
          temp_new_samples <- lapply(unique(components), function(component) ifelse(temp_samples %in% which(components == component), 1, 0))
          temp_new_samples <- do.call(cbind, temp_new_samples)
          colnames(temp_new_samples) <- paste0(par, " (inclusion: ", unique(components),")")

          # place the transformed samples back
          model_samples <- cbind(
            if(temp_position > 1) model_samples[,1:(temp_position-1),drop=FALSE],
            temp_new_samples,
            if(temp_position <= ncol(model_samples)) model_samples[,temp_position:ncol(model_samples),drop=FALSE]
          )
        }
      }
    }
  }

  # remove transformations for removed variables
  if(!is.null(transformations)){
    transformations <-  transformations[names(transformations) %in% names(prior_list)]
  }

  # apply transformations (not orthornormal if they are to be returned transformed to diffs)
  model_samples <- .apply_parameter_transformations(model_samples, transformations, prior_list, transform_factors)

  # transform orthonormal factors to differences from mean
  model_samples <- .transform_factor_contrasts(model_samples, prior_list, transform_factors, transformations)

  # rename factor levels
  model_samples <- .rename_factor_levels(model_samples, prior_list)

  # store parameter names before removing formula attachments
  parameter_names <- colnames(model_samples)

  # rename formula parameters
  if(any(!sapply(lapply(prior_list, attr, which = "parameter"), is.null))){
    raw_parameter_names <- colnames(model_samples)
    colnames(model_samples) <- format_parameter_names(
      parameters         = colnames(model_samples),
      formula_parameters = unique(unlist(lapply(prior_list, attr, which = "parameter"))),
      formula_random     = unique(unlist(lapply(prior_list, attr, which = "random_factor"))),
      formula_prefix     = formula_prefix,
      formula_scale      = if(transform_scaled) formula_scale else NULL)
    colnames(model_samples) <- .bt_random_effect_summary_display_names(
      names = colnames(model_samples),
      raw_names = raw_parameter_names,
      prior_list = prior_list,
      formula_prefix = formula_prefix,
      parameter_registry = parameter_registry
    )
  }

  # return samples if requested
  if(return_samples){
    attr(model_samples, "prior_list") <- prior_list
    return(model_samples)
  }

  # compute the summary
  if(ncol(model_samples) == 0){
    empty_table <- runjags_estimates_empty_table(
      probs              = probs,
      title              = title,
      footnotes          = footnotes,
      warnings           = warnings,
      remove_diagnostics = remove_diagnostics,
      diagnostic_columns = summary_diagnostic_columns
    )
    if(random_effects_metadata){
      empty_table <- .bt_random_effect_summary_add_metadata_columns(
        table = empty_table,
        parameter_names = character(),
        prior_list = prior_list,
        parameter_registry = parameter_registry
      )
    }
    return(empty_table)
  }else{
    runjags_summary <- .runjags_summary_fast(
      model_samples       = model_samples,
      n_samples           = fit$sample,
      n_chains            = length(fit$mcmc),
      conditional         = conditional,
      probs               = probs,
      remove_diagnostics  = remove_diagnostics,
      diagnostic_columns  = diagnostic_columns
    )
  }

  # prepare output
  n_estimate_cols <- 2 + length(probs)  # Mean, SD, quantiles
  class(runjags_summary)              <- c("BayesTools_table", "BayesTools_runjags_summary", class(runjags_summary))
  attr(runjags_summary, "type")       <- c(rep("estimate", n_estimate_cols), summary_diagnostic_columns)
  attr(runjags_summary, "parameters") <- parameter_names
  attr(runjags_summary, "rownames")   <- TRUE
  attr(runjags_summary, "title")      <- title
  attr(runjags_summary, "footnotes")  <- footnotes
  attr(runjags_summary, "warnings")   <- warnings
  if(random_effects_metadata){
    runjags_summary <- .bt_random_effect_summary_add_metadata_columns(
      table = runjags_summary,
      parameter_names = parameter_names,
      prior_list = prior_list,
      parameter_registry = parameter_registry
    )
  }

  return(runjags_summary)
}

.bt_JAGS_estimates_filter_raw_random_columns <- function(model_samples,
                                                         prior_list,
                                                         parameter_registry = NULL,
                                                         formula_design = NULL,
                                                         remove_parameters = NULL,
                                                         remove_formulas = NULL,
                                                         keep_parameters = NULL,
                                                         keep_formulas = NULL,
                                                         remove_random_effects = NULL,
                                                         keep_random_effects = NULL,
                                                         remove_random_structures = NULL,
                                                         keep_random_structures = NULL){

  if(is.null(parameter_registry)){
    parameter_registry <- .bt_build_parameter_registry(
      columns = colnames(model_samples),
      prior_list = prior_list,
      formula_design = formula_design
    )
  }
  model_samples <- .bt_random_effect_summary_filter_raw_columns(
    model_samples = model_samples,
    parameter_registry = parameter_registry,
    remove_random_effects = remove_random_effects,
    keep_random_effects = keep_random_effects,
    remove_random_structures = remove_random_structures,
    keep_random_structures = keep_random_structures
  )

  column_names <- colnames(model_samples)
  if(length(column_names) == 0L){
    return(model_samples)
  }
  .bt_validate_parameter_registry(parameter_registry)
  registry_rows <- match(column_names, parameter_registry$canonical_name)
  registered <- !is.na(registry_rows)
  row_roles <- rep("", length(column_names))
  row_roles[registered] <- parameter_registry$role[registry_rows[registered]]
  raw_random <- registered &
    (startsWith(row_roles, "random_") | row_roles == "allocation")
  if(!any(raw_random)){
    return(model_samples)
  }

  remove_aliases <- .bt_JAGS_estimates_random_aliases(remove_parameters)
  keep_aliases <- .bt_JAGS_estimates_random_aliases(keep_parameters)
  keep_active <- !is.null(keep_parameters) || !is.null(keep_formulas)
  if(length(remove_aliases) == 0L && length(remove_formulas) == 0L && !keep_active){
    return(model_samples)
  }

  random_prior_columns <- .bt_JAGS_estimates_random_prior_columns(
    column_names = column_names,
    prior_list = prior_list
  )
  remove_columns <- rep(FALSE, length(column_names))
  random_blocks <- unique(
    parameter_registry$random_block[registry_rows[raw_random]]
  )

  for(random_block in random_blocks){
      term_columns <- raw_random &
        parameter_registry$random_block[registry_rows] == random_block
      if(!any(term_columns)){
        next
      }
      formula_parameter <- unique(
        parameter_registry$formula_parameter[registry_rows[term_columns]]
      )
      formula_parameter <- formula_parameter[nzchar(formula_parameter)]
      remove_formula <- any(formula_parameter %in% remove_formulas)
      keep_formula <- any(formula_parameter %in% keep_formulas)

      if(remove_formula || .bt_JAGS_estimates_random_alias_has_all(remove_aliases)){
        remove_columns <- remove_columns | term_columns
      }else if(.bt_JAGS_estimates_random_alias_has_correlation(remove_aliases)){
        remove_columns <- remove_columns |
          (term_columns & row_roles == "random_correlation")
      }

      if(keep_active){
        keep_columns <- random_prior_columns
        if(keep_formula || .bt_JAGS_estimates_random_alias_has_all(keep_aliases)){
          keep_columns <- keep_columns | term_columns
        }
        if(.bt_JAGS_estimates_random_alias_has_correlation(keep_aliases)){
          keep_columns <- keep_columns |
            (term_columns & row_roles == "random_correlation")
        }
        if(!is.null(keep_random_effects) || !is.null(keep_random_structures)){
          term_matches <- TRUE
          if(!is.null(keep_random_effects)){
            term_registry_rows <- registry_rows[term_columns]
            term_matches <- term_matches &&
              (random_block %in% keep_random_effects ||
                 any(parameter_registry$random_name[
                   term_registry_rows
                 ] %in% keep_random_effects) ||
                 any(parameter_registry$random_grouping[
                   term_registry_rows
                 ] %in% keep_random_effects))
          }
          if(!is.null(keep_random_structures)){
            term_matches <- term_matches && any(
              parameter_registry$random_structure[
                registry_rows[term_columns]
              ] %in% keep_random_structures
            )
          }
          if(term_matches){
            keep_columns <- keep_columns | term_columns
          }
        }
        remove_columns <- remove_columns | (term_columns & !keep_columns)
      }
  }

  model_samples[, !remove_columns, drop = FALSE]
}

.bt_JAGS_estimates_random_prior_columns <- function(column_names, prior_list){

  if(length(prior_list) == 0L){
    return(rep(FALSE, length(column_names)))
  }
  random_flags <- .bt_random_effect_prior_flags(prior_list)
  random_prior_names <- random_flags$name[random_flags$any]

  .bt_random_effect_summary_parameter_columns(
    column_names = column_names,
    parameter_names = random_prior_names
  )
}

.bt_JAGS_estimates_random_aliases <- function(parameters){

  if(is.null(parameters) || !is.character(parameters)){
    return(character())
  }

  intersect(
    parameters,
    c(
      "random", "random_effects", "random_sd", "random_rho",
      "random_cor", "random_correlation", "random_variance_fraction",
      "random_variance_ratio", "random_allocation", "random_sd_multiplier"
    )
  )
}

.bt_JAGS_estimates_random_alias_has_all <- function(aliases){

  any(aliases %in% c("random", "random_effects"))
}

.bt_JAGS_estimates_random_alias_has_correlation <- function(aliases){

  any(aliases %in% c("random_rho", "random_cor", "random_correlation"))
}

.bt_JAGS_estimates_raw_random_correlation_columns <- function(column_names){

  grepl("_xRE_CORx", column_names, fixed = TRUE) |
    grepl("_rho(_z|_logit)?(\\[|$)", column_names)
}

#' @rdname BayesTools_model_tables
runjags_inference_table  <- function(fit, title = NULL, footnotes = NULL, warnings = NULL, formula_prefix = TRUE,
                                     logBF = FALSE, BF01 = FALSE, BF_diagnostics = FALSE,
                                     BF_diagnostic_columns = getOption("BayesTools.JAGS_BF_diagnostic_columns", if(BF_diagnostics) "all" else "none")){

  # check fits
  if(!inherits(fit, "runjags"))
    stop("'fit' must be a runjags fit")
  if(!inherits(fit, "BayesTools_fit"))
    stop("'fit' must be a BayesTools fit")
  prior_list <- attr(fit, "prior_list")
  check_list(prior_list, "prior_list")
  if(!all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_char(title, "title", allow_NULL = TRUE)
  check_char(footnotes, "footnotes", check_length = 0, allow_NULL = TRUE)
  check_char(warnings, "warnings", check_length = 0, allow_NULL = TRUE)
  check_bool(formula_prefix, "formula_prefix")
  check_bool(logBF, "logBF", allow_NA = FALSE)
  check_bool(BF01,  "BF01",  allow_NA = FALSE)
  check_bool(BF_diagnostics, "BF_diagnostics")
  BF_diagnostic_columns <- .normalize_diagnostic_columns(BF_diagnostic_columns, .JAGS_BF_diagnostic_columns(), "BF_diagnostic_columns")
  BF_diagnostics        <- length(BF_diagnostic_columns) > 0
  BF_error_diagnostics  <- "BF_error_percent" %in% BF_diagnostic_columns

  # return empty table if none of the priors is spike and slab
  if(!any(sapply(prior_list, function(p) is.prior.spike_and_slab(p) | is.prior.mixture(p)))){
    runjags_summary <- runjags_inference_empty_table(
      title          = title,
      footnotes      = footnotes,
      warnings       = warnings,
      logBF          = logBF,
      BF01           = BF01,
      BF_diagnostics = BF_diagnostics,
      BF_diagnostic_columns = BF_diagnostic_columns
    )
    return(runjags_summary)
  }

  # extract samples
  model_samples   <- .extract_posterior_samples(fit, as_list = FALSE)
  if(BF_diagnostics){
    model_samples_list <- .extract_posterior_samples(fit, as_list = TRUE)
  }else{
    model_samples_list <- NULL
  }
  runjags_summary <- data.frame(matrix(nrow = 0, ncol = 4 + length(BF_diagnostic_columns)))
  colnames(runjags_summary) <- c("Parameter", "prior_prob", "post_prob", "inclusion_BF", BF_diagnostic_columns)
  BF_bound_operators <- character()

  for(par in names(prior_list)){
    if(is.prior.spike_and_slab(prior_list[[par]])){

      temp_prior_prob <- mean(.get_spike_and_slab_inclusion(prior_list[[par]]))
      temp_post_prob  <- mean(model_samples[,paste0(par, "_indicator")])
      temp_BF         <- inclusion_BF(
        prior_probs = c(null = 1 - temp_prior_prob, alternative = temp_prior_prob),
        post_probs  = c(null = 1 - temp_post_prob,  alternative = temp_post_prob),
        is_null     = c(TRUE, FALSE)
      )
      temp_BF_reporting <- .indicator_BF_reporting_value(temp_BF, temp_post_prob, temp_prior_prob, nrow(model_samples))
      temp_row        <- data.frame(
        Parameter    = par,
        prior_prob   = temp_prior_prob,
        post_prob    = temp_post_prob,
        inclusion_BF = temp_BF_reporting$value
      )
      if(BF_diagnostics){
        temp_indicator_list <- .runjags_indicator_list(model_samples_list, paste0(par, "_indicator"), 1)
        temp_diagnostics    <- .indicator_BF_diagnostics(temp_indicator_list, temp_prior_prob, temp_BF)
        temp_row            <- cbind(temp_row, .indicator_BF_diagnostic_row(temp_diagnostics)[, BF_diagnostic_columns, drop = FALSE])
        if(BF_error_diagnostics){
          warnings <- c(warnings, .indicator_BF_warnings(par, temp_diagnostics))
        }
      }

      runjags_summary <- rbind(runjags_summary, temp_row)
      BF_bound_operators <- c(BF_bound_operators, temp_BF_reporting$operator)
    }else if(is.prior.mixture(prior_list[[par]])){

      # extract the components and prior probabilities
      components      <- attr(prior_list[[par]], "components")
      temp_prior_prob <- attr(prior_list[[par]], "prior_weights")
      temp_prior_prob <- sapply(unique(components), function(component) sum(temp_prior_prob[which(components == component)])) / sum(temp_prior_prob)
      temp_post_prob  <- sapply(unique(components), function(component) mean(model_samples[,paste0(par, "_indicator")] %in% which(components == component)))

      # if only null and alternative are specified, removed the null component
      if(all(components %in% c("null", "alternative"))){

        if(all(components == "null")){
          temp_prior_prob <- c(temp_prior_prob, "alternative" = 0)
          temp_post_prob  <- c(temp_post_prob,  "alternative" = 0)
        }
        if(all(components == "alternative")){
          temp_prior_prob <- c("null" = 0, temp_prior_prob)
          temp_post_prob  <- c("null" = 0,  temp_post_prob)
        }

        temp_BF  <- inclusion_BF(prior_probs = temp_prior_prob, post_probs = temp_post_prob, is_null = names(temp_post_prob) != "alternative")
        temp_BF_reporting <- .indicator_BF_reporting_value(temp_BF, temp_post_prob[["alternative"]], temp_prior_prob[["alternative"]], nrow(model_samples))
        temp_row <- data.frame(
          Parameter    = par,
          prior_prob   = temp_prior_prob[["alternative"]],
          post_prob    = temp_post_prob[["alternative"]],
          inclusion_BF = temp_BF_reporting$value
        )
        if(BF_diagnostics){
          temp_indicator_list <- .runjags_indicator_list(model_samples_list, paste0(par, "_indicator"), which(components == "alternative"))
          temp_diagnostics    <- .indicator_BF_diagnostics(temp_indicator_list, temp_prior_prob[["alternative"]], temp_BF)
          temp_row            <- cbind(temp_row, .indicator_BF_diagnostic_row(temp_diagnostics)[, BF_diagnostic_columns, drop = FALSE])
          if(BF_error_diagnostics){
            warnings <- c(warnings, .indicator_BF_warnings(par, temp_diagnostics))
          }
        }

        runjags_summary <- rbind(runjags_summary, temp_row)
        BF_bound_operators <- c(BF_bound_operators, temp_BF_reporting$operator)

      }else{

        # compute summary for each component
        for(component in unique(components)){
          temp_parameter <- paste0(par, " [", component, "]")
          temp_BF        <- inclusion_BF(prior_probs = temp_prior_prob, post_probs = temp_post_prob, is_null = names(temp_post_prob) != component)
          temp_BF_reporting <- .indicator_BF_reporting_value(temp_BF, temp_post_prob[[component]], temp_prior_prob[[component]], nrow(model_samples))
          temp_row       <- data.frame(
            Parameter    = paste0(par, " [", component, "]"),
            prior_prob   = temp_prior_prob[[component]],
            post_prob    = temp_post_prob[[component]],
            inclusion_BF = temp_BF_reporting$value
          )
          if(BF_diagnostics){
            temp_indicator_list <- .runjags_indicator_list(model_samples_list, paste0(par, "_indicator"), which(components == component))
            temp_diagnostics    <- .indicator_BF_diagnostics(temp_indicator_list, temp_prior_prob[[component]], temp_BF)
            temp_row            <- cbind(temp_row, .indicator_BF_diagnostic_row(temp_diagnostics)[, BF_diagnostic_columns, drop = FALSE])
            if(BF_error_diagnostics){
              warnings <- c(warnings, .indicator_BF_warnings(temp_parameter, temp_diagnostics))
            }
          }

          runjags_summary <- rbind(runjags_summary, temp_row)
          BF_bound_operators <- c(BF_bound_operators, temp_BF_reporting$operator)
        }

      }
    }
  }

  # store parameter names before removing formula attachments
  parameter_names           <- runjags_summary$Parameter
  rownames(runjags_summary) <- parameter_names
  runjags_summary           <- runjags_summary[,-1]

  # format BF and BF MC error on the requested scale
  temp_inclusion_BF <- runjags_summary[,"inclusion_BF"]
  attr(temp_inclusion_BF, "bound_operator") <- BF_bound_operators
  runjags_summary[["inclusion_BF"]] <- format_BF(temp_inclusion_BF, logBF = logBF, BF01 = BF01, inclusion = TRUE)
  if("MCMC_error" %in% BF_diagnostic_columns){
    attr(runjags_summary[["MCMC_error"]], "name") <- "error(Post. prob.)"
  }
  if("BF_error_percent" %in% BF_diagnostic_columns){
    attr(runjags_summary[["BF_error_percent"]], "name") <- .BF_error_column_name(BF01)
  }

  # rename formula parameters
  if(any(!sapply(lapply(prior_list, attr, which = "parameter"), is.null))){
    rownames(runjags_summary) <- format_parameter_names(
      parameters         = rownames(runjags_summary),
      formula_parameters = unique(unlist(lapply(prior_list, attr, which = "parameter"))),
      formula_random     = unique(unlist(lapply(prior_list, attr, which = "random_factor"))),
      formula_prefix     = formula_prefix,
      formula_scale      = NULL)
  }

  class(runjags_summary)               <- c("BayesTools_table", "BayesTools_runjags_inference", class(runjags_summary))
  attr(runjags_summary, "type")        <- c("prior_prob", "post_prob", "inclusion_BF", .JAGS_BF_diagnostic_column_types(BF_diagnostic_columns))
  attr(runjags_summary, "parameters")  <- parameter_names
  attr(runjags_summary, "rownames")    <- TRUE
  attr(runjags_summary, "title")       <- title
  attr(runjags_summary, "footnotes")   <- footnotes
  attr(runjags_summary, "warnings")    <- warnings

  return(runjags_summary)
}

#' @rdname BayesTools_model_tables
JAGS_estimates_table <- runjags_estimates_table

#' @rdname BayesTools_model_tables
JAGS_inference_table <- runjags_inference_table

#' @rdname BayesTools_model_tables
JAGS_summary_table   <- model_summary_table

#' @rdname BayesTools_model_tables
model_summary_empty_table <- function(model_description = NULL, title = NULL, footnotes = NULL, warnings = NULL){

  check_list(model_description, "model_description", allow_NULL = TRUE)

  summary_names  <- c(
    "Model",
    if(!is.null(model_description)) names(model_description),
    "Prior prob.",
    "log(marglik)",
    "Post. prob.",
    "Inclusion BF")


  summary_names <- paste0(summary_names, "  ")

  empty_table <- data.frame(cbind(
    summary_names,
    rep("",            length(summary_names)),
    rep("           ", length(summary_names)),
    c("Parameter prior distributions", rep("", length(summary_names) - 1))
  ))
  names(empty_table) <- NULL

  # prepare output
  class(empty_table)             <- c("BayesTools_table", class(empty_table))
  attr(empty_table, "type")      <- c("string_left", "string", "string", "prior")
  attr(empty_table, "rownames")  <- FALSE
  attr(empty_table, "as.matrix") <- TRUE
  attr(empty_table, "title")     <- title
  attr(empty_table, "footnotes") <- footnotes
  attr(empty_table, "warnings")  <- warnings

  return(empty_table)
}

#' @rdname BayesTools_model_tables
runjags_estimates_empty_table <- function(probs = c(0.025, 0.5, 0.975), title = NULL, footnotes = NULL, warnings = NULL,
                                          remove_diagnostics = FALSE,
                                          diagnostic_columns = getOption("BayesTools.JAGS_estimates_diagnostic_columns", if(remove_diagnostics) "none" else "all")){

  check_bool(remove_diagnostics, "remove_diagnostics")
  diagnostic_columns <- .normalize_diagnostic_columns(diagnostic_columns, .JAGS_estimates_diagnostic_columns(), "diagnostic_columns")
  if(remove_diagnostics){
    diagnostic_columns <- character()
  }
  n_estimate_cols <- 2 + length(probs)  # Mean, SD, quantiles
  empty_table <- data.frame(matrix(nrow = 0, ncol = n_estimate_cols + length(diagnostic_columns)), check.names = FALSE)
  colnames(empty_table) <- c("Mean", "SD", as.character(probs), diagnostic_columns)

  class(empty_table)             <- c("BayesTools_table", "BayesTools_runjags_summary", class(empty_table))
  attr(empty_table, "type")      <- c(rep("estimate", n_estimate_cols), diagnostic_columns)
  attr(empty_table, "rownames")  <- FALSE
  attr(empty_table, "title")     <- title
  attr(empty_table, "footnotes") <- footnotes
  attr(empty_table, "warnings")  <- warnings

  return(empty_table)
}

#' @rdname BayesTools_model_tables
runjags_inference_empty_table <- function(title = NULL, footnotes = NULL, warnings = NULL,
                                          logBF = FALSE, BF01 = FALSE, BF_diagnostics = FALSE,
                                          BF_diagnostic_columns = getOption("BayesTools.JAGS_BF_diagnostic_columns", if(BF_diagnostics) "all" else "none")){

  check_char(title, "title", allow_NULL = TRUE)
  check_char(footnotes, "footnotes", check_length = 0, allow_NULL = TRUE)
  check_char(warnings, "warnings", check_length = 0, allow_NULL = TRUE)
  check_bool(logBF, "logBF", allow_NA = FALSE)
  check_bool(BF01,  "BF01",  allow_NA = FALSE)
  check_bool(BF_diagnostics, "BF_diagnostics")
  BF_diagnostic_columns <- .normalize_diagnostic_columns(BF_diagnostic_columns, .JAGS_BF_diagnostic_columns(), "BF_diagnostic_columns")
  BF_diagnostics        <- length(BF_diagnostic_columns) > 0

  empty_table <- data.frame(matrix(nrow = 0, ncol = 3 + length(BF_diagnostic_columns)))
  colnames(empty_table) <- c("prior_prob", "post_prob", "inclusion_BF", BF_diagnostic_columns)
  attr(empty_table[["inclusion_BF"]], "name") <- .BF_column_name(logBF = logBF, BF01 = BF01, inclusion = TRUE)
  if("MCMC_error" %in% BF_diagnostic_columns){
    attr(empty_table[["MCMC_error"]], "name") <- "error(Post. prob.)"
  }
  if("BF_error_percent" %in% BF_diagnostic_columns){
    attr(empty_table[["BF_error_percent"]], "name") <- .BF_error_column_name(BF01)
  }

  class(empty_table)             <- c("BayesTools_table", "BayesTools_runjags_inference", class(empty_table))
  attr(empty_table, "type")      <- c("prior_prob", "post_prob", "inclusion_BF", .JAGS_BF_diagnostic_column_types(BF_diagnostic_columns))
  attr(empty_table, "rownames")  <- FALSE
  attr(empty_table, "title")     <- title
  attr(empty_table, "footnotes") <- footnotes
  attr(empty_table, "warnings")  <- warnings

  return(empty_table)
}

#' @rdname BayesTools_model_tables
JAGS_estimates_empty_table <- runjags_estimates_empty_table

#' @rdname BayesTools_model_tables
JAGS_inference_empty_table <- runjags_inference_empty_table

#' @rdname BayesTools_model_tables
stan_estimates_table  <- function(fit, transformations = NULL, title = NULL, footnotes = NULL, warnings = NULL){

  # this is a simplification of the runjags_estimates_table function for stan
  .check_rstan()

  # check fits
  if(!inherits(fit, "stanfit"))
    stop("'fit' must be a rstan fit")
  prior_list <- attr(fit, "prior_list")
  check_list(prior_list, "prior_list")
  if(!all(sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")
  check_list(transformations, "transformations", allow_NULL = TRUE)
  if(!is.null(transformations) && any(!sapply(transformations, function(trans)is.function(trans[["fun"]]))))
    stop("'transformations' must be list of functions in the 'fun' element.")
  check_char(title, "title", allow_NULL = TRUE)
  check_char(footnotes, "footnotes", check_length = 0, allow_NULL = TRUE)
  check_char(warnings, "warnings", check_length = 0, allow_NULL = TRUE)

  # obtain model information
  stan_summary  <- data.frame(rstan::summary(fit)$summary)
  model_samples <- .extract_stan(fit, drop = FALSE)

  # remove un-wanted columns
  stan_summary <- stan_summary[,!colnames(stan_summary) %in% c("X25.", "X75."), drop = FALSE]

  # remove un-wanted rows
  stan_summary <- stan_summary[-nrow(stan_summary),, drop = FALSE]

  # rename columns to match runjags output
  colnames(stan_summary) <- c("Mean", "MCerr", "SD", "Lower95", "Median", "Upper95", "SSeff", "psrf")

  # add MC.ofSD estimates
  stan_summary[, "MC.ofSD"] <- stan_summary[, "MCerr"] / stan_summary[, "SD"]


  # apply transformations
  if(!is.null(transformations)){
    for(par in names(transformations)){
      model_samples[,par] <- do.call(transformations[[par]][["fun"]], c(list(model_samples[,par]), transformations[[par]][["arg"]]))
      transformed_summary <- .stan_transformed_summary(model_samples[,par], stan_summary[par, "SSeff"])
      stan_summary[par, names(transformed_summary)] <- transformed_summary
      stan_summary[par, "MC.ofSD"] <- stan_summary[par, "MCerr"] / stan_summary[par, "SD"]
    }
  }


  # rename the rest
  colnames(stan_summary)[colnames(stan_summary) == "Lower95"] <- "0.025"
  colnames(stan_summary)[colnames(stan_summary) == "Median"]  <- "0.5"
  colnames(stan_summary)[colnames(stan_summary) == "Upper95"] <- "0.975"
  colnames(stan_summary)[colnames(stan_summary) == "MCerr"]   <- "MCMC_error"
  colnames(stan_summary)[colnames(stan_summary) == "MC.ofSD"] <- "MCMC_SD_error"
  colnames(stan_summary)[colnames(stan_summary) == "SSeff"]   <- "ESS"
  colnames(stan_summary)[colnames(stan_summary) == "psrf"]    <- "R_hat"

  # reorder the columns
  stan_summary <- stan_summary[,c("Mean", "SD", "0.025", "0.5", "0.975", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"), drop = FALSE]

  # store parameter names
  parameter_names <- rownames(stan_summary)

  # prepare output
  class(stan_summary)              <- c("BayesTools_table", "BayesTools_stan_summary", class(stan_summary))
  attr(stan_summary, "type")       <- c(rep("estimate", 5), "MCMC_error", "MCMC_SD_error", "ESS", "R_hat")
  attr(stan_summary, "parameters") <- parameter_names
  attr(stan_summary, "rownames")   <- TRUE
  attr(stan_summary, "title")      <- title
  attr(stan_summary, "footnotes")  <- footnotes
  attr(stan_summary, "warnings")   <- warnings

  return(stan_summary)
}

.stan_transformed_summary <- function(samples, SSeff){

  samples <- as.numeric(samples)
  qs      <- stats::quantile(samples, probs = c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)
  sd      <- stats::sd(samples, na.rm = TRUE)
  MCerr   <- if(is.finite(SSeff) && SSeff > 0) sd / sqrt(SSeff) else NA_real_

  c(
    "Mean"    = mean(samples, na.rm = TRUE),
    "MCerr"   = MCerr,
    "SD"      = sd,
    "Lower95" = qs[1],
    "Median"  = qs[2],
    "Upper95" = qs[3]
  )
}
