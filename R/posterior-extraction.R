#' @title Helper functions for extracting and formatting posterior distributions
#'
#' @description Internal helper functions to extract posterior samples from JAGS
#' fits and reformat them for further processing (summary tables, diagnostics, plots).
#' These functions consolidate common logic that was duplicated across
#' \code{runjags_estimates_table}, \code{.diagnostics_plot_data}, and plotting functions.
#'
#' @name posterior_extraction_helpers
#' @keywords internal
NULL


#' @rdname posterior_extraction_helpers
#' @param fit a JAGS model fit object
#' @param as_list whether to return samples as mcmc.list (TRUE) or merged matrix (FALSE)
#' @return matrix or mcmc.list of posterior samples
.extract_posterior_samples <- function(fit, as_list = FALSE) {
  
  if (as_list) {
    # Use generic function to allow S3 method dispatch (runjags has its own as.mcmc.list method)
    model_samples <- coda::as.mcmc.list(fit)
  } else {
    # Use generic function to allow S3 method dispatch (runjags has its own as.mcmc method)
    model_samples <- suppressWarnings(coda::as.mcmc(fit))
  }
  
  return(model_samples)
}


#' @rdname posterior_extraction_helpers
#' @param model_samples matrix of posterior samples
#' @param prior_list list of prior objects
#' @param remove_parameters character vector of parameter names to remove
#' @return list with cleaned model_samples and updated prior_list
.remove_auxiliary_parameters <- function(model_samples, prior_list, remove_parameters = NULL) {
  
  for (i in rev(seq_along(prior_list))) {
    
    # invgamma support parameter
    if (is.prior.simple(prior_list[[i]]) && prior_list[[i]][["distribution"]] == "invgamma") {
      model_samples <- model_samples[, colnames(model_samples) != paste0("inv_", names(prior_list)[i]), drop = FALSE]
    }
    
    # weightfunction parameters
    if (is.prior.weightfunction(prior_list[[i]])) {
      # remove etas
      if (prior_list[[i]][["distribution"]] %in% c("one.sided", "two.sided")) {
        model_samples <- model_samples[, !grepl("eta", colnames(model_samples)), drop = FALSE]
      }
      
      # rename the omegas
      omega_cuts      <- weightfunctions_mapping(prior_list[i], cuts_only = TRUE)
      omega_names_old <- paste0("omega[", 1:(length(omega_cuts) - 1), "]")
      omega_names     <- sapply(1:(length(omega_cuts) - 1), function(j) paste0("omega[", omega_cuts[j], ",", omega_cuts[j + 1], "]"))
      
      # change the order of omegas
      model_samples[, which(colnames(model_samples) %in% omega_names_old)] <- model_samples[, rev(which(colnames(model_samples) %in% omega_names_old)), drop = FALSE]
      colnames(model_samples)[which(colnames(model_samples) %in% omega_names_old)] <- omega_names
      
      # remove omegas if requested
      if ("omega" %in% remove_parameters) {
        model_samples <- model_samples[, !colnames(model_samples) %in% omega_names, drop = FALSE]
        prior_list[i] <- NULL
      }
      
    } else if (names(prior_list)[[i]] %in% remove_parameters) {
      # remove parameters to be excluded (note: spike_0 removal is handled by caller)
      if (is.prior.factor(prior_list[[i]])) {
        model_samples <- model_samples[, !colnames(model_samples) %in% .JAGS_prior_factor_names(names(prior_list)[i], prior_list[[i]]), drop = FALSE]
      } else {
        model_samples <- model_samples[, colnames(model_samples) != names(prior_list)[i], drop = FALSE]
      }
      prior_list[i] <- NULL
    }
  }
  
  return(list(model_samples = model_samples, prior_list = prior_list))
}


#' @rdname posterior_extraction_helpers
#' @param par parameter name
#' @param conditional whether to compute conditional summary
#' @param remove_inclusion whether to remove inclusion indicators
#' @param warnings character vector for collecting warnings
#' @return list with updated model_samples, prior_list, and warnings
.process_spike_and_slab <- function(model_samples, prior_list, par, conditional = FALSE, remove_inclusion = FALSE, warnings = NULL) {
  
  # prepare parameter names
  if (is.prior.factor(.get_spike_and_slab_variable(prior_list[[par]]))) {
    if (.get_prior_factor_levels(.get_spike_and_slab_variable(prior_list[[par]])) == 1) {
      par_names <- par
    } else {
      par_names <- paste0(par, "[", 1:.get_prior_factor_levels(.get_spike_and_slab_variable(prior_list[[par]])), "]")
    }
  } else {
    par_names <- par
  }
  
  # change the samples between conditional/averaged based on the preferences
  if (conditional) {
    # compute the number of conditional samples
    n_conditional_samples <- sum(model_samples[, colnames(model_samples) == paste0(par, "_indicator")] == 1)
    
    # replace null samples with NAs (important for later transformations)
    model_samples[model_samples[, colnames(model_samples) == paste0(par, "_indicator")] != 1, par_names] <- NA
    
    # add warnings about conditional summary
    warnings <- c(warnings, .runjags_conditional_warning(par_names, n_conditional_samples))
  }
  
  # remove the inclusion
  model_samples <- model_samples[, colnames(model_samples) != paste0(par, "_inclusion"), drop = FALSE]
  
  # remove the latent variable
  model_samples <- model_samples[, !colnames(model_samples) %in% gsub(par, paste0(par, "_variable"), par_names), drop = FALSE]
  
  # remove/rename the inclusions probabilities
  if (remove_inclusion) {
    model_samples <- model_samples[, colnames(model_samples) != paste0(par, "_indicator"), drop = FALSE]
  } else {
    colnames(model_samples)[colnames(model_samples) == paste0(par, "_indicator")] <- paste0(par, " (inclusion)")
  }
  
  # modify the parameter list (forward the parameter attribute)
  variable_component <- .get_spike_and_slab_variable(prior_list[[par]])
  attr(variable_component, "parameter") <- attr(prior_list[[par]], "parameter")
  prior_list[[par]] <- variable_component
  
  return(list(model_samples = model_samples, prior_list = prior_list, warnings = warnings))
}


#' @rdname posterior_extraction_helpers
#' @param transformations list of transformations to apply
#' @param transform_factors whether orthonormal/meandif will be transformed later
#' @return updated model_samples matrix
.apply_parameter_transformations <- function(model_samples, transformations, prior_list, transform_factors = FALSE) {
  
  if (is.null(transformations)) {
    return(model_samples)
  }
  
  for (par in names(transformations)) {
    if (!is.prior.factor(prior_list[[par]])) {
      # non-factor priors
      model_samples[, par] <- do.call(transformations[[par]][["fun"]], c(list(model_samples[, par]), transformations[[par]][["arg"]]))
    } else if ((!transform_factors && (is.prior.orthonormal(prior_list[[par]]) | is.prior.meandif(prior_list[[par]]))) || is.prior.treatment(prior_list[[par]])) {
      # treatment priors, or orthonormal/meandif that won't be transformed to differences
      par_names <- .JAGS_prior_factor_names(par, prior_list[[par]])
      
      for (i in seq_along(par_names)) {
        model_samples[, par_names[i]] <- do.call(transformations[[par]][["fun"]], c(list(model_samples[, par_names[i]]), transformations[[par]][["arg"]]))
      }
    }
  }
  
  return(model_samples)
}


#' @rdname posterior_extraction_helpers
#' @param transform_factors whether to transform orthonormal/meandif to differences
#' @return updated model_samples matrix
.transform_factor_contrasts <- function(model_samples, prior_list, transform_factors = FALSE, transformations = NULL) {
  
  if (!transform_factors || !any(sapply(prior_list, function(x) is.prior.orthonormal(x) | is.prior.meandif(x)))) {
    return(model_samples)
  }
  
  message("The transformation was applied to the differences from the mean. Note that non-linear transformations do not map from the orthonormal/meandif contrasts to the differences from the mean.")
  
  for (par in names(prior_list)[sapply(prior_list, function(x) is.prior.orthonormal(x) | is.prior.meandif(x))]) {
    
    par_names <- .JAGS_prior_factor_names(par, prior_list[[par]])
    
    temp_position <- min(which(colnames(model_samples) %in% par_names))
    temp_samples  <- model_samples[, colnames(model_samples) %in% par_names, drop = FALSE]
    model_samples <- model_samples[, !colnames(model_samples) %in% par_names, drop = FALSE]
    
    if (is.prior.orthonormal(prior_list[[par]])) {
      transformed_samples <- temp_samples %*% t(contr.orthonormal(1:(.get_prior_factor_levels(prior_list[[par]]) + 1)))
    } else if (is.prior.meandif(prior_list[[par]])) {
      transformed_samples <- temp_samples %*% t(contr.meandif(1:(.get_prior_factor_levels(prior_list[[par]]) + 1)))
    }
    
    # apply transformation if specified
    if (!is.null(transformations[par])) {
      for (i in 1:ncol(transformed_samples)) {
        transformed_samples[, i] <- do.call(transformations[[par]][["fun"]], c(list(transformed_samples[, i]), transformations[[par]][["arg"]]))
      }
    }
    
    if (.is_prior_interaction(prior_list[[par]])) {
      if (length(.get_prior_factor_level_names(prior_list[[par]])) == 1) {
        transformed_names <- paste0(par, " [dif: ", .get_prior_factor_level_names(prior_list[[par]])[[1]], "]")
      } else {
        stop("orthonormal/meandif de-transformation for interaction of multiple factors is not implemented.")
      }
    } else {
      transformed_names <- paste0(par, " [dif: ", .get_prior_factor_level_names(prior_list[[par]]), "]")
    }
    colnames(transformed_samples) <- transformed_names
    
    # place the transformed samples back
    model_samples <- cbind(
      if (temp_position > 1) model_samples[, 1:(temp_position - 1), drop = FALSE],
      transformed_samples,
      if (temp_position <= ncol(model_samples)) model_samples[, temp_position:ncol(model_samples), drop = FALSE]
    )
  }
  
  return(model_samples)
}


#' @rdname posterior_extraction_helpers
#' @return updated model_samples matrix with renamed columns
.rename_factor_levels <- function(model_samples, prior_list) {
  
  # rename treatment factor levels
  if (any(sapply(prior_list, is.prior.treatment))) {
    for (par in names(prior_list)[sapply(prior_list, is.prior.treatment)]) {
      if (!.is_prior_interaction(prior_list[[par]])) {
        if (.get_prior_factor_levels(prior_list[[par]]) == 1) {
          colnames(model_samples)[colnames(model_samples) == par] <-
            paste0(par, "[", .get_prior_factor_level_names(prior_list[[par]])[-1], "]")
        } else {
          colnames(model_samples)[colnames(model_samples) %in% paste0(par, "[", 1:.get_prior_factor_levels(prior_list[[par]]), "]")] <-
            paste0(par, "[", .get_prior_factor_level_names(prior_list[[par]])[-1], "]")
        }
      } else if (length(attr(prior_list[[par]], "levels")) == 1) {
        colnames(model_samples)[colnames(model_samples) %in% paste0(par, "[", 1:.get_prior_factor_levels(prior_list[[par]]), "]")] <-
          paste0(par, "[", .get_prior_factor_level_names(prior_list[[par]])[[1]][-1], "]")
      }
    }
  }
  
  # rename independent factor levels
  if (any(sapply(prior_list, is.prior.independent))) {
    for (par in names(prior_list)[sapply(prior_list, is.prior.independent)]) {
      if (!.is_prior_interaction(prior_list[[par]])) {
        if (.get_prior_factor_levels(prior_list[[par]]) == 1) {
          colnames(model_samples)[colnames(model_samples) == par] <-
            paste0(par, "[", .get_prior_factor_level_names(prior_list[[par]]), "]")
        } else {
          colnames(model_samples)[colnames(model_samples) %in% paste0(par, "[", 1:.get_prior_factor_levels(prior_list[[par]]), "]")] <-
            paste0(par, "[", .get_prior_factor_level_names(prior_list[[par]]), "]")
        }
      } else if (length(attr(prior_list[[par]], "levels")) == 1) {
        colnames(model_samples)[colnames(model_samples) %in% paste0(par, "[", 1:.get_prior_factor_levels(prior_list[[par]]), "]")] <-
          paste0(par, "[", .get_prior_factor_level_names(prior_list[[par]])[[1]], "]")
      }
    }
  }
  
  return(model_samples)
}
