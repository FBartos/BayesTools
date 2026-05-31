.condition_special_labels <- function(){

  c("PET", "PEESE", "PETPEESE", "omega", "phacking", "alpha", "pi_null")
}

.condition_normalize_labels <- function(conditional){

  conditional <- unlist(conditional, use.names = FALSE)
  conditional <- as.character(conditional)
  conditional <- conditional[!is.na(conditional) & nzchar(conditional)]

  unique(conditional)
}

.condition_labels_key <- function(conditional){

  conditional <- sort(unique(.condition_key_labels(conditional)))
  paste0(c(length(conditional), conditional), collapse = "\r")
}

.condition_event_key <- function(conditional, conditional_rule = "AND"){

  if(is.null(conditional_rule)){
    conditional_rule <- "AND"
  }
  conditional <- sort(unique(.condition_key_labels(conditional)))
  if(length(conditional) == 0L){
    return("<averaged>")
  }

  paste0(c(conditional_rule, length(conditional), conditional), collapse = "\r")
}

.condition_key_labels <- function(conditional){

  conditional <- .condition_normalize_labels(conditional)
  vapply(conditional, .condition_key_label, character(1), USE.NAMES = FALSE)
}

.condition_key_label <- function(label){

  if(label %in% c("alpha", "pi_null")){
    return("phacking")
  }

  label
}

.condition_event <- function(prior_list, conditional = NULL,
                             conditional_rule = "AND"){

  conditional <- .condition_normalize_labels(conditional)

  out <- list(
    conditional      = conditional,
    conditional_rule = conditional_rule,
    families         = .condition_event_families(prior_list, conditional),
    condition_key    = .condition_event_key(conditional, conditional_rule)
  )
  class(out) <- "BayesTools_condition_event"

  return(out)
}

.condition_event_rule_function <- function(conditional_rule){

  if(identical(conditional_rule, "AND")){
    return(all)
  }
  if(identical(conditional_rule, "OR")){
    return(any)
  }

  stop("'conditional_rule' must be either 'AND' or 'OR'.", call. = FALSE)
}

.condition_event_set_attributes <- function(x, condition_event, effective = FALSE){

  attr(x, "conditional")              <- condition_event[["conditional"]]
  attr(x, "conditional_rule")         <- condition_event[["conditional_rule"]]
  attr(x, "condition_key")            <- condition_event[["condition_key"]]
  attr(x, "condition_event")          <- condition_event
  attr(x, "resolved_condition_event") <- condition_event

  if(effective){
    attr(x, "effective_conditional")      <- condition_event[["conditional"]]
    attr(x, "effective_conditional_rule") <- condition_event[["conditional_rule"]]
  }

  x
}

.condition_event_families <- function(prior_list, conditional){

  families <- list()
  for(label in conditional){
    family <- .condition_event_label_family(prior_list, label)
    if(is.null(families[[family[["name"]]]])){
      families[[family[["name"]]]] <- family
    }
    families[[family[["name"]]]][["labels"]] <- c(
      families[[family[["name"]]]][["labels"]],
      label
    )
  }

  families
}

.condition_event_label_family <- function(prior_list, label){

  if(label %in% .condition_special_labels()){
    if(.condition_event_has_bias_prior(prior_list)){
      return(list(name = "bias", type = "bias", labels = character()))
    }
    if(.condition_event_special_is_tautology(prior_list, label)){
      return(list(name = paste0("<tautology>:", label), type = "tautology", labels = character()))
    }
  }

  if(label %in% names(prior_list)){
    if(label == "bias" && .condition_event_has_bias_prior(prior_list)){
      return(list(name = "bias", type = "bias", labels = character()))
    }
    return(list(name = label, type = "parameter", labels = character()))
  }

  list(name = paste0("<unknown>:", label), type = "unknown", labels = character())
}

.condition_event_has_bias_prior <- function(prior_list){

  !is.null(prior_list[["bias"]]) &&
    (is.prior.mixture(prior_list[["bias"]]) ||
       is.prior.PET(prior_list[["bias"]]) ||
       is.prior.PEESE(prior_list[["bias"]]) ||
       is.prior.weightfunction(prior_list[["bias"]]) ||
       is_prior_phacking(prior_list[["bias"]]) ||
       is_prior_bias(prior_list[["bias"]]) ||
       is.prior.none(prior_list[["bias"]]))
}

.condition_event_special_is_tautology <- function(prior_list, label){

  if(label == "omega"){
    return(any(vapply(prior_list, function(x){
      is.prior.weightfunction(x) || (is_prior_bias(x) && !is.null(x$selection))
    }, logical(1))))
  }

  if(label %in% c("phacking", "alpha", "pi_null")){
    return(any(vapply(prior_list, function(x){
      is_prior_phacking(x) || (is_prior_bias(x) && !is.null(x$phacking))
    }, logical(1))))
  }

  FALSE
}

.condition_event_label_posterior_mask <- function(prior_list, model_samples,
                                                  label){

  if(label %in% .condition_special_labels() &&
     .condition_event_has_bias_prior(prior_list)){
    prior <- prior_list[["bias"]]
    values <- .condition_event_bias_label_values(prior, label)
    if(is.prior.mixture(prior)){
      return(model_samples[, "bias_indicator"] %in% which(values))
    }
    return(rep(isTRUE(values[[1]]), nrow(model_samples)))
  }

  if(.condition_event_special_is_tautology(prior_list, label)){
    return(rep(TRUE, nrow(model_samples)))
  }

  if(!label %in% names(prior_list)){
    warning(
      sprintf(
        "The parameter '%s' is not a conditional parameter. All samples are assumed to come from the conditional posterior distribution.",
        label
      ),
      call. = FALSE,
      immediate. = TRUE
    )
    return(rep(TRUE, nrow(model_samples)))
  }

  prior <- prior_list[[label]]
  if(is.prior.spike_and_slab(prior)){
    return(model_samples[, paste0(label, "_indicator")] == 1)
  }

  if(is.prior.mixture(prior)){
    components <- attr(prior, "components")
    if(!all(components %in% c("null", "alternative"))){
      stop(
        "conditional mixture posterior distributions are available only for 'null' and 'alternative' components",
        call. = FALSE
      )
    }
    return(model_samples[, paste0(label, "_indicator")] %in%
             which(components == "alternative"))
  }

  warning(
    sprintf(
      "The parameter '%s' is not a conditional parameter. All samples are assumed to come from the conditional posterior distribution.",
      label
    ),
    call. = FALSE,
    immediate. = TRUE
  )
  rep(TRUE, nrow(model_samples))
}

.condition_event_posterior_mask <- function(event, prior_list, model_samples){

  conditional <- event[["conditional"]]
  if(length(conditional) == 0L){
    return(rep(TRUE, nrow(model_samples)))
  }

  conditioning_samples <- do.call(cbind, lapply(conditional, function(label){
    .condition_event_label_posterior_mask(prior_list, model_samples, label)
  }))
  rule_fun <- .condition_event_rule_function(event[["conditional_rule"]])

  apply(conditioning_samples, 1, rule_fun)
}

.condition_event_bias_label_values <- function(prior, labels){

  labels <- .condition_normalize_labels(labels)
  branches <- if(is.prior.mixture(prior)) prior else list(prior)
  branch_info <- lapply(branches, .selection_branch_info)

  is_PET        <- vapply(branches, is.prior.PET, logical(1))
  is_PEESE      <- vapply(branches, is.prior.PEESE, logical(1))
  has_selection <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
  has_phacking  <- vapply(branch_info, function(x) !is.null(x$phacking), logical(1))

  values <- matrix(FALSE, nrow = length(branches), ncol = length(labels))
  colnames(values) <- labels

  for(label in labels){
    values[, label] <- switch(
      label,
      "bias" = .condition_event_bias_alternative(prior, branches),
      "PET" = is_PET,
      "PEESE" = is_PEESE,
      "PETPEESE" = is_PET | is_PEESE,
      "omega" = has_selection,
      "phacking" = has_phacking,
      "alpha" = has_phacking,
      "pi_null" = has_phacking,
      rep(TRUE, length(branches))
    )
  }

  if(length(labels) == 1L){
    return(values[, 1])
  }

  values
}

.condition_event_bias_alternative <- function(prior, branches){

  if(is.prior.mixture(prior)){
    components <- attr(prior, "components")
    if(!is.null(components) && length(components) == length(branches)){
      if(!all(components %in% c("null", "alternative"))){
        stop(
          "conditional mixture posterior distributions are available only for 'null' and 'alternative' components",
          call. = FALSE
        )
      }
      return(components == "alternative")
    }
  }

  !vapply(branches, is.prior.none, logical(1))
}

.condition_event_family_options <- function(prior_list, family){

  labels <- unique(family[["labels"]])

  if(family[["type"]] %in% c("unknown", "tautology")){
    values <- rep(TRUE, length(labels))
    names(values) <- labels
    return(list(list(
      family      = family[["name"]],
      prior       = NULL,
      probability = 1,
      values      = values,
      replace     = FALSE
    )))
  }

  if(family[["type"]] == "bias"){
    return(.condition_event_bias_options(prior_list[["bias"]], labels))
  }

  prior <- prior_list[[family[["name"]]]]
  options <- .prior_density_condition_component(prior)
  lapply(options, function(option){
    values <- rep(option[["alternative"]], length(labels))
    names(values) <- labels
    list(
      family      = family[["name"]],
      prior       = option[["prior"]],
      probability = option[["probability"]],
      values      = values,
      replace     = TRUE
    )
  })
}

.condition_event_bias_options <- function(prior, labels){

  branches <- if(is.prior.mixture(prior)) prior else list(prior)
  probabilities <- if(is.prior.mixture(prior)){
    prior_weights <- attr(prior, "prior_weights")
    prior_weights / sum(prior_weights)
  }else{
    1
  }
  values <- .condition_event_bias_label_values(prior, labels)
  if(is.null(dim(values))){
    values <- matrix(values, ncol = 1, dimnames = list(NULL, labels))
  }

  lapply(seq_along(branches), function(i){
    list(
      family      = "bias",
      prior       = branches[[i]],
      probability = probabilities[i],
      values      = values[i, ],
      replace     = TRUE
    )
  })
}

.condition_event_model_options <- function(prior_list, event){

  if(length(event[["conditional"]]) == 0L){
    return(NULL)
  }

  families <- event[["families"]]
  if(length(families) == 0L){
    return(NULL)
  }

  options <- lapply(families, function(family){
    .condition_event_family_options(prior_list, family)
  })
  names(options) <- names(families)

  option_grid <- expand.grid(lapply(options, seq_along))
  keep <- logical(nrow(option_grid))
  model_weights <- numeric(nrow(option_grid))
  rule_fun <- .condition_event_rule_function(event[["conditional_rule"]])

  for(i in seq_len(nrow(option_grid))){
    values <- logical(length(event[["conditional"]]))
    names(values) <- event[["conditional"]]
    probabilities <- numeric(length(families))

    for(j in seq_along(families)){
      option <- options[[j]][[option_grid[i, j]]]
      values[names(option[["values"]])] <- option[["values"]]
      probabilities[j] <- option[["probability"]]
    }

    keep[i] <- rule_fun(values)
    model_weights[i] <- prod(probabilities)
  }

  event_probability <- sum(model_weights[keep & model_weights > 0])

  option_grid <- option_grid[keep & model_weights > 0, , drop = FALSE]
  model_weights <- model_weights[keep & model_weights > 0]

  if(nrow(option_grid) == 0L){
    return(list(
      prior_lists       = list(),
      weights           = numeric(),
      event_probability = 0
    ))
  }

  prior_lists <- lapply(seq_len(nrow(option_grid)), function(i){
    model_prior_list <- prior_list
    for(j in seq_along(families)){
      option <- options[[j]][[option_grid[i, j]]]
      if(isTRUE(option[["replace"]])){
        family_name <- option[["family"]]
        model_prior_list[[family_name]] <- .prior_density_copy_parent_attributes(
          component = option[["prior"]],
          parent    = prior_list[[family_name]]
        )
      }
    }
    model_prior_list
  })

  list(
    prior_lists       = prior_lists,
    weights           = model_weights / sum(model_weights),
    event_probability = event_probability
  )
}

.condition_event_active_labels <- function(prior_list, weights, conditional){

  conditional <- .condition_normalize_labels(conditional)
  if(length(conditional) == 0L){
    return(character())
  }

  active <- .prior_linear_active_parameters(prior_list, weights)
  if(length(active) == 0L){
    return(character())
  }

  keep <- vapply(conditional, function(label){
    if(label %in% .condition_special_labels()){
      return("bias" %in% active || label %in% active)
    }
    label %in% active
  }, logical(1))

  conditional[keep]
}
