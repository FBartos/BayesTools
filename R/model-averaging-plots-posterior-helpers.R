.simplify_spike_samples           <- function(samples, prior_list){

  # Check if we're dealing with spike_and_slab or mixture (which are single priors) vs list of priors
  is_spike_and_slab <- is.prior.spike_and_slab(prior_list)
  is_mixture <- is.prior.mixture(prior_list)

  # If we have a spike_and_slab or mixture prior, we need to iterate over their components
  # Otherwise, we have a list of individual priors
  if(is_spike_and_slab || is_mixture) {
    # For spike_and_slab and mixture, iterate over the components
    components_to_iterate <- prior_list
    component_indices <- seq_along(prior_list)
  } else {
    # For lists of priors, iterate over the list
    components_to_iterate <- prior_list
    component_indices <- seq_along(prior_list)
  }

  # aggregate for each spike
  priors_point_map <- data.frame(do.call(rbind, lapply(component_indices, function(i) {
    current_component <- components_to_iterate[[i]]
    if(is.prior.point(current_component)){
      if(is_spike_and_slab) {
        # For spike_and_slab: dbern() generates 0 (null) and 1 (alternative)
        # We need to determine which component this is
        component_name <- attr(current_component, "component")
        model_index <- if(component_name == "null") 0 else 1
      } else {
        # For mixture or list of priors: dcat() generates 1, 2, 3... so index i maps to JAGS index i
        model_index <- i
      }
      c("location" = current_component$parameters[["location"]], "frequency" = sum(attr(samples, "models_ind") == model_index))
    }
  })))


  # return the input with fewer than 2 inputs
  if(nrow(priors_point_map) < 2){
    spike_probability = data.frame(cbind(
      "location"    = priors_point_map[, "location"],
      "probability" = priors_point_map[, "frequency"] / if(!is.matrix(samples)) length(samples) else nrow(samples) ))
    spike_probability <- spike_probability[priors_point_map[, "frequency"] != 0, ]
    return(spike_probability)
  }

  # find unique spikes
  unique_map <- cbind("location" = unique(priors_point_map[, "location"]), "frequency" = 0)

  # collect them
  for(i in 1:nrow(unique_map)){
    unique_map[i, "frequency"] <- sum(priors_point_map[sapply(priors_point_map[, "location"], function(l) isTRUE(all.equal(l, unname(unique_map[i, "location"])))), "frequency"])
  }

  spike_probability = data.frame(cbind(
    "location"    = unique_map[, "location"],
    "probability" = unique_map[, "frequency"] / if(!is.matrix(samples)) length(samples) else nrow(samples) ))
  spike_probability <- spike_probability[unique_map[, "frequency"] != 0, ]

  return(spike_probability)
}
.bias_prior_list_for_condition <- function(prior_list, condition_event){

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
.simplify_as_mixed_posterior_bias <- function(samples, parameter) {

  ### replace all remaining priors by null prior
  prior_list <- attr(samples[["bias"]], "prior_list")
  condition_event <- attr(samples[["bias"]], "resolved_condition_event", exact = TRUE)
  if(is.null(condition_event)){
    condition_event <- attr(samples[["bias"]], "condition_event", exact = TRUE)
  }
  if(is.null(condition_event)){
    condition_event <- attr(samples, "resolved_condition_event", exact = TRUE)
  }
  if(is.null(condition_event)){
    condition_event <- attr(samples, "condition_event", exact = TRUE)
  }
  prior_list <- .bias_prior_list_for_condition(prior_list, condition_event)

  if (parameter == "PET") {
    prior_ind <- which(sapply(prior_list, \(x) !is.prior.PET(x)))
  } else if (parameter == "PEESE") {
    prior_ind <- which(sapply(prior_list, \(x) !is.prior.PEESE(x)))
  } else if (parameter == "omega") {
    prior_ind <- which(sapply(prior_list, \(x) !.weightfunction_prior_has_selection(x)))
  }
  if (length(prior_ind) > 0) {
    for (i in prior_ind) {
      temp_weight     <- prior_list[[i]][["prior_weights"]]
      prior_list[[i]] <- if (parameter == "omega") prior_none() else prior("point", parameters = list(0))
      prior_list[[i]][["prior_weights"]] <- temp_weight
    }
  }

  ### create new samples
  new_samples <- samples[["bias"]][, grepl(parameter, colnames(samples[["bias"]])),drop=FALSE]

  ### store attribute
  std_attrs  <- c("dim", "dimnames", "names", "prior_list", "mcpar")
  all_attrs  <- attributes(samples[["bias"]])
  to_restore <- setdiff(names(all_attrs), std_attrs)

  ### re-assign attributes
  for (a in to_restore) {
    attr(new_samples, a) <- all_attrs[[a]]
  }

  # remove `mixed_posteriors.bias` class
  class(new_samples) <- class(new_samples)[!class(new_samples) %in% "mixed_posteriors.bias"]

  ### assign prior list and model indicator
  attr(prior_list, "omega_context") <- attr(samples[["bias"]], "omega_context")
  attr(new_samples, "prior_list") <- prior_list

  ### remove the old samples & store new samples
  samples[["bias"]]    <- NULL
  samples[[parameter]] <- new_samples

  return(samples)
}
