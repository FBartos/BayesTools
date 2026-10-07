# Plots use declared posterior atoms; point masses are never inferred from the
# prior list and the per-draw model indicators.
.plot_data_stop_unknown_atoms <- function(){

  stop(
    "Posterior atom status is unknown. Plotting posterior samples requires an ",
    "explicit atom/no-atom declaration; attach posterior_atom_attribute() ",
    "metadata or use a BayesTools posterior producer that records it. ",
    "Posteriors created by BayesTools 0.3.0 do not record it: recompute them ",
    "with the current version (refitting models fitted with 0.3.0).",
    call. = FALSE
  )
}
.bias_prior_list_for_condition <- function(prior_list, condition_event){

  .model_probability_plot_check(prior_list)
  prior_list_fallback <- .weightfunction_expand_bias_mixture_priors(prior_list)

  if(is.null(condition_event) ||
     length(.posterior_density_normalize_condition(condition_event[["conditional"]])) == 0L){
    return(prior_list_fallback)
  }

  bias_family <- condition_event[["families"]][["bias"]]
  if(is.null(bias_family) || length(bias_family[["labels"]]) == 0L){
    return(prior_list_fallback)
  }

  values <- .condition_event_bias_label_values(prior_list, bias_family[["labels"]])
  if(is.null(dim(values))){
    keep <- as.logical(values)
  }else{
    rule_fun <- .condition_event_rule_function(condition_event[["conditional_rule"]])
    keep <- apply(values, 1L, rule_fun)
  }
  keep[is.na(keep)] <- FALSE
  if(length(keep) != length(prior_list_fallback) || !any(keep)){
    return(prior_list_fallback)
  }

  prior_weights <- vapply(prior_list_fallback, .prior_model_weight, numeric(1))
  prior_weights[!is.finite(prior_weights) | prior_weights < 0] <- 0
  if(sum(prior_weights[keep]) <= 0){
    return(prior_list_fallback)
  }
  prior_weights[!keep] <- 0
  prior_weights <- prior_weights / sum(prior_weights)

  for(i in seq_along(prior_list_fallback)){
    prior_list_fallback[[i]] <- .set_prior_model_weight(prior_list_fallback[[i]], prior_weights[i])
  }

  prior_list_fallback
}
.bias_samples_condition_event <- function(samples){

  condition_event <- .bt_meta_condition(samples[["bias"]], "resolved_condition_event")
  if(is.null(condition_event)){
    condition_event <- .bt_meta_condition(samples[["bias"]], "condition_event")
  }
  if(is.null(condition_event)){
    condition_event <- .bt_meta_condition(samples, "resolved_condition_event")
  }
  if(is.null(condition_event)){
    condition_event <- .bt_meta_condition(samples, "condition_event")
  }

  condition_event
}
.bias_samples_prior_list <- function(samples){

  # bias branches with the prior weights implied by the samples' condition
  prior_list <- .bias_prior_list_for_condition(
    attr(samples[["bias"]], "prior_list"),
    .bias_samples_condition_event(samples)
  )

  # with the model prior list, the weights are P(branch | event) of the joint
  # condition event, which may combine bias labels with other parameters
  branch_weights <- .bias_samples_branch_weights(samples)
  if(!is.null(branch_weights) && length(branch_weights) == length(prior_list)){
    for(i in seq_along(prior_list)){
      prior_list[[i]] <- .set_prior_model_weight(prior_list[[i]], branch_weights[i])
    }
  }

  attr(prior_list, "omega_context") <- attr(samples[["bias"]], "omega_context")

  prior_list
}
.bias_samples_branch_weights <- function(samples){

  condition_event <- .bias_samples_condition_event(samples)
  if(is.null(condition_event) ||
     length(.posterior_density_normalize_condition(condition_event[["conditional"]])) == 0L){
    return(NULL)
  }

  model_prior_list <- attr(samples, "prior_list", exact = TRUE)
  if(!is.list(model_prior_list) || is.prior(model_prior_list) ||
     !is.prior.mixture(model_prior_list[["bias"]])){
    return(NULL)
  }

  # tag the branches to recover them from the model options
  bias_prior <- model_prior_list[["bias"]]
  n_branches <- length(bias_prior)
  for(k in seq_len(n_branches)){
    branch <- bias_prior[[k]]
    attr(branch, "bias_branch_index") <- k
    bias_prior[[k]] <- branch
  }
  model_prior_list[["bias"]] <- bias_prior

  models <- .condition_event_model_options(model_prior_list, condition_event)
  if(is.null(models) || length(models[["prior_lists"]]) == 0L){
    return(NULL)
  }

  prior_weights <- attr(bias_prior, "prior_weights")
  if(is.null(prior_weights)){
    prior_weights <- vapply(seq_len(n_branches), function(k) .prior_model_weight(bias_prior[[k]]), numeric(1))
  }
  prior_weights <- prior_weights / sum(prior_weights)

  weights <- numeric(n_branches)
  for(i in seq_along(models[["prior_lists"]])){
    branch <- attr(models[["prior_lists"]][[i]][["bias"]], "bias_branch_index", exact = TRUE)
    if(is.null(branch)){
      # the bias prior is not part of this option: all branches in proportion
      weights <- weights + models[["weights"]][i] * prior_weights
    }else{
      weights[branch] <- weights[branch] + models[["weights"]][i]
    }
  }

  weights / sum(weights)
}
.petpeese_samples_prior_pairs <- function(samples, mu_name){

  # Joint prior of the effect and the bias branch under the samples'
  # condition. A condition may involve mu, the bias event, other parameters,
  # or combine them with OR, so mu and bias are paired within each model
  # option of the condition event rather than conditioned separately.
  prior_list <- attr(samples, "prior_list", exact = TRUE)
  if(is.list(prior_list) && is.prior(prior_list[[mu_name]]) && is.prior(prior_list[["bias"]])){
    .model_probability_plot_check(prior_list[[mu_name]])
    .model_probability_plot_check(prior_list[["bias"]])
  }
  if(!is.list(prior_list) || is.prior(prior_list) ||
     !is.prior(prior_list[[mu_name]]) || !is.prior(prior_list[["bias"]])){
    return(NULL)
  }

  condition_event <- .bias_samples_condition_event(samples)
  models <- NULL
  if(!is.null(condition_event) &&
     length(.posterior_density_normalize_condition(condition_event[["conditional"]])) > 0L){
    models <- .condition_event_model_options(prior_list, condition_event)
  }
  if(is.null(models)){
    models <- list(prior_lists = list(prior_list), weights = 1)
  }
  if(length(models[["prior_lists"]]) == 0L){
    stop("The condition of the samples has zero prior probability.", call. = FALSE)
  }
  if(!is.null(models$log_weights)) .model_probability_measure_stop(list(
    probabilities=models$weights, logs=models$log_weights, declaration=models$model_probability_declaration))

  mu_priors   <- list()
  bias_priors <- list()
  weights     <- numeric()
  for(i in seq_along(models[["prior_lists"]])){
    model_prior_list <- models[["prior_lists"]][[i]]
    bias_prior       <- model_prior_list[["bias"]]
    if(is.prior.mixture(bias_prior)){
      branch_weights <- attr(bias_prior, "prior_weights")
      if(is.null(branch_weights)){
        branch_weights <- vapply(bias_prior, .prior_model_weight, numeric(1))
      }
      branches <- lapply(seq_along(bias_prior), function(k) bias_prior[[k]])
    }else{
      branch_weights <- 1
      branches       <- list(bias_prior)
    }
    branch_weights <- branch_weights / sum(branch_weights)

    for(k in seq_along(branches)){
      mu_priors[[length(mu_priors) + 1L]]     <- model_prior_list[[mu_name]]
      bias_priors[[length(bias_priors) + 1L]] <- .petpeese_bias_branch(branches[[k]])
      weights <- c(weights, models[["weights"]][i] * branch_weights[k])
    }
  }

  # merge repeated (mu, bias) pairs
  keep <- rep(TRUE, length(weights))
  for(i in seq_along(weights)){
    if(!keep[i]){
      next
    }
    for(j in seq_along(weights)[-seq_len(i)]){
      if(keep[j] && identical(mu_priors[[i]], mu_priors[[j]]) &&
         identical(bias_priors[[i]], bias_priors[[j]])){
        weights[i] <- weights[i] + weights[j]
        keep[j]    <- FALSE
      }
    }
  }
  keep <- keep & weights > 0

  bias_priors <- lapply(which(keep), function(i){
    .set_prior_model_weight(bias_priors[[i]], weights[i])
  })

  list(
    prior_list    = bias_priors,
    prior_list_mu = mu_priors[keep]
  )
}
.petpeese_bias_branch <- function(branch){

  # branches without PET or PEESE terms (no bias, weightfunctions, ...)
  # imply PET = PEESE = 0
  if(!(is.prior.PET(branch) || is.prior.PEESE(branch))){
    return(prior("point", parameters = list(location = 0)))
  }

  # weights of the plotted pairs are set explicitly
  attr(branch, "model_prior_weights") <- NULL
  branch
}
.simplify_as_mixed_posterior_bias <- function(samples, parameter) {

  ### replace all remaining priors by null prior
  prior_list <- .bias_samples_prior_list(samples)

  if (parameter == "PET") {
    prior_ind <- which(sapply(prior_list, \(x) !is.prior.PET(x)))
  } else if (parameter == "PEESE") {
    prior_ind <- which(sapply(prior_list, \(x) !is.prior.PEESE(x)))
  } else if (parameter == "omega") {
    prior_ind <- which(sapply(prior_list, \(x) !.weightfunction_prior_has_selection(x)))
  }
  if (length(prior_ind) > 0) {
    for (i in prior_ind) {
      parent_prior <- prior_list[[i]]
      prior_list[[i]] <- .prior_density_copy_parent_attributes(
        if(parameter == "omega") prior_none() else prior("point", parameters = list(0)), parent_prior)
      if(is.null(attr(parent_prior, "model_probability_declaration", exact = TRUE))){
        prior_list[[i]] <- .set_prior_model_weight(prior_list[[i]], .prior_model_weight(parent_prior))
      }
    }
  }

  ### create new samples
  original <- samples[["bias"]]
  selected <- grepl(parameter, colnames(original))
  new_samples <- .bt_draws_plain(original)[, selected, drop = FALSE]
  metadata <- .bt_meta_current_container(original)
  columns <- colnames(new_samples)
  if(!is.null(metadata)){
    if(!is.null(metadata$atoms) && any(selected)){
      atoms <- .posterior_atoms_from_attribute(metadata$atoms)
      metadata$atoms <- .posterior_atoms_new(atoms$locations[, selected, drop = FALSE], atoms$mass,
        column_names = columns, source = atoms$source,
        component_probabilities = atoms$component_probabilities,
        component_log_probabilities = atoms$component_log_probabilities,
        model_probability_declaration = atoms$model_probability_declaration,
        marginals = if(!is.null(atoms$marginals)) atoms$marginals[selected])
    }else metadata$atoms <- NULL
    for(field in c("support", "prior_densities")){
      if(.posterior_metadata_is_container(metadata[[field]])) metadata[[field]] <- metadata[[field]][columns]
    }
    for(field in c("quantities", "original_scale_quantities", "level_quantities")){
      if(!is.null(metadata[[field]]) && nrow(metadata[[field]]) == ncol(original)){
        metadata[[field]] <- metadata[[field]][selected, , drop = FALSE]
      }
    }
    if(!is.null(metadata$measure_unavailable)){
      metadata$measure_unavailable <- metadata$measure_unavailable[
        metadata$measure_unavailable$column %in% c(columns, parameter), , drop = FALSE]
      if(nrow(metadata$measure_unavailable) == 0L) metadata$measure_unavailable <- NULL
    }
  }
  if(parameter %in% c("PET", "PEESE") && ncol(new_samples) == 0L){
    indicator <- .bt_draws_component(samples[["bias"]])
    atoms <- .posterior_atoms_get(samples[["bias"]])
    probabilities <- if(is.null(atoms)) NULL else atoms$component_probabilities
    active <- unique(c(which(probabilities > 0), indicator))
    if(length(indicator) != nrow(new_samples) ||
       any(!indicator %in% seq_along(prior_list)) ||
       any(!active %in% seq_along(prior_list))){
      stop("Bias posterior component metadata are unavailable for scalar plotting.", call. = FALSE)
    }
    if(!all(vapply(prior_list[active], is.prior.point, logical(1)))){
      stop("Posterior samples for '", parameter,
           "' are unavailable because an active bias branch is not a point prior.",
           call. = FALSE)
    }
    locations <- vapply(prior_list[active], function(prior){
      prior$parameters[["location"]]
    }, numeric(1))
    new_samples <- matrix(
      locations[match(indicator, active)], ncol = 1L,
      dimnames = list(NULL, parameter)
    )
  }

  ### store attribute
  std_attrs  <- c("dim", "dimnames", "names", "prior_list", "mcpar", .bt_meta_attribute)
  all_attrs  <- attributes(samples[["bias"]])
  to_restore <- setdiff(names(all_attrs), std_attrs)

  ### re-assign attributes
  for (a in to_restore) {
    attr(new_samples, a) <- all_attrs[[a]]
  }
  if(!is.null(metadata)) new_samples <- .bt_meta_write(new_samples, metadata, .bt_meta_fingerprint(new_samples))
  new_samples <- .bt_meta_refresh(new_samples)

  # remove `mixed_posteriors.bias` class
  class(new_samples) <- class(new_samples)[!class(new_samples) %in% "mixed_posteriors.bias"]
  if(parameter == "omega") class(new_samples) <- unique(c(class(new_samples), "mixed_posteriors.weightfunction"))

  ### assign prior list and model indicator
  attr(prior_list, "omega_context") <- attr(samples[["bias"]], "omega_context")
  attr(new_samples, "prior_list") <- prior_list
  if(parameter %in% c("PET", "PEESE") && ncol(new_samples) == 1L){
    atoms <- .posterior_atoms_get(samples[["bias"]])
    if(!is.null(atoms)){
      probabilities <- atoms$component_probabilities
      scalar_atoms <- if(is.null(probabilities)){
        .posterior_atoms_from_components(
          prior_list, .bt_draws_component(new_samples),
          n_columns = 1L, column_names = parameter
        )
      }else{
        pair <- .bt_meta_get(original, "model_probabilities")$posterior
        if(is.null(pair) && !is.null(atoms$model_probability_declaration)) pair <- list(
          probabilities = atoms$component_probabilities, logs = atoms$component_log_probabilities,
          declaration = atoms$model_probability_declaration)
        .posterior_atoms_from_priors(
          prior_list, probabilities, n_columns = 1L,
          column_names = parameter, source = atoms$source, posterior_pair = pair
        )
      }
      new_samples <- .posterior_atoms_set(new_samples, scalar_atoms)
    }
  }

  ### remove the old samples & store new samples
  samples[["bias"]]    <- NULL
  samples[[parameter]] <- new_samples

  return(samples)
}
